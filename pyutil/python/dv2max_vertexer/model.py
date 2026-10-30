"""DV2Net: residual sparse U-Net (spconv or the in-repo torch backend) + SPVCNN point branch + FiLM event
context + DeepVtx-2 heads.  No attention, no transformer block anywhere (tests/test_no_attention.py).

Backbone (per level l = 0..L-1, voxel size voxel * 2**l):
  stem       BN(inputs) -> SubM 3^3 (Fv -> C0)
  encoder    `blocks` pre-activation residual blocks (BN, LeakyReLU, SubM 3^3, BN, LeakyReLU, SubM 3^3, + skip);
             at the `large_kernel.levels` coarsest levels the first conv of each block is the 13^3
             spatial-partition large kernel (+ parallel SubM 3^3); then BN, LeakyReLU, k2 s2 down conv
  context    mean + max pooling of the bottleneck per event, concatenated with event-level features -> MLP ->
             ctx; FiLM (gamma, beta from ctx, zero-initialised) on the bottleneck and every decoder level
  decoder    inverse conv (paired with the down conv), concat skip, linear fuse, `blocks` residual blocks, FiLM
  points     per-point MLP of the point features; at the bottleneck and every decoder level the level's voxel
             features are gathered to the points (devoxelise) and the points' mean is scattered back (SPVCNN)
Precision: sparse convs run in `cdtype` (bf16/fp16 under AMP on CUDA); BatchNorm, the residual stream and all
geometry (coordinates, offsets, log-scales) are fp32.

Heads (level 0 unless noted; `heads` config switches):
  heat        pass 1: P(within 3 cm of the target) (penalty-reduced focal on a Gaussian); pass 2: P(< 1 cm)
  any         vote weight, any visible neutrino
  offset      Laplace mean (3) + log-scale (3), fp32, anchor = voxel centroid
  cluster     per (event, clustering_pr cluster): scatter mean+max of voxel features -> MLP -> P(host)
  event       ctx -> MLP: pass 1 P(visible nu in event); pass 2 P(target inside the crop)
  point       per imaging point: vertexness P(< 1 cm) + ordinal distance-to-vertex classes (shell edges)
  shell_cone  (C2 arm) per voxel: ordinal shell classes of the distance to the vertex + von Mises-Fisher
              direction to the vertex (cone)
  listwise    (C3 arm) per cluster and per voxel scores for a hierarchical listwise softmax
  deep        per decoder level l >= 1: heat logit (deep supervision)
"""
import torch
from torch import nn
from .scatter import scatter_max, scatter_mean

from . import sparse
from .inputs import N_EVENT_FEATURES


def act(x):
    return nn.functional.leaky_relu(x, 0.1)


def mlp(cin, hidden, cout, zero_last=False):
    m = nn.Sequential(nn.Linear(cin, hidden), nn.LeakyReLU(0.1), nn.Linear(hidden, cout))
    if zero_last:
        nn.init.zeros_(m[-1].weight)
        nn.init.zeros_(m[-1].bias)
    return m


class ResBlock(nn.Module):
    def __init__(self, backend, C, key, lk_size=None):
        super().__init__()
        self.bn1, self.bn2 = nn.BatchNorm1d(C), nn.BatchNorm1d(C)
        self.c1 = sparse.LargeKernelPartition3d(backend, C, key, lk_size) if lk_size else sparse.subm(backend, C, C,
                                                                                                        key)
        self.c2 = sparse.subm(backend, C, C, key)

    def forward(self, x, cdt):
        h = act(self.bn1(x.feats.float())).to(cdt)
        h = self.c1(x.replace(h))
        h2 = act(self.bn2(h.feats.float())).to(cdt)
        h2 = self.c2(h.replace(h2))
        return x.replace(x.feats.float() + h2.feats.float())


class Down(nn.Module):
    def __init__(self, backend, cin, cout, key):
        super().__init__()
        self.bn = nn.BatchNorm1d(cin)
        self.conv = sparse.down(backend, cin, cout, key)

    def forward(self, x, cdt):
        y = self.conv(x.replace(act(self.bn(x.feats.float())).to(cdt)))
        return y.replace(y.feats.float())


class Up(nn.Module):
    def __init__(self, backend, cin, cout, key):
        super().__init__()
        self.bn = nn.BatchNorm1d(cin)
        self.conv = sparse.up(backend, cin, cout, key)

    def forward(self, x, cdt):
        y = self.conv(x.replace(act(self.bn(x.feats.float())).to(cdt)))
        return y.replace(y.feats.float())


class FiLM(nn.Module):
    def __init__(self, ctx_dim, C):
        super().__init__()
        self.lin = nn.Linear(ctx_dim, 2 * C)
        nn.init.zeros_(self.lin.weight)
        nn.init.zeros_(self.lin.bias)

    def forward(self, h, ctx, batch):
        g, b = self.lin(ctx).float().chunk(2, dim=1)
        return h * (1 + g[batch]) + b[batch]


class DV2Net(nn.Module):
    def __init__(self, pcfg, n_vox_feat, n_pt_feat, backend):
        super().__init__()
        self.pcfg = pcfg
        self.backend = backend
        L = pcfg["levels"]
        C = list(pcfg["channels"])
        nb = pcfg["blocks"]
        lk = pcfg["large_kernel"]
        lk_from = L - lk["levels"] if lk["enabled"] else L
        self.L, self.C = L, C
        self.cdtype = torch.float32
        self.return_features = False          # SSL pretraining reads the level-0 decoder features
        # input standardisation: eps 0.1 so that rare, near-constant input channels (variance ~1e-7 measured on
        # an overfit run) cannot be amplified ~300x at test time; normalised inputs are also clamped to +-20
        self.in_bn = nn.BatchNorm1d(n_vox_feat, eps=0.1)
        self.stem = sparse.subm(backend, n_vox_feat, C[0], "L0")
        self.enc = nn.ModuleList()
        self.down = nn.ModuleList()
        for l in range(L):
            self.enc.append(nn.ModuleList([ResBlock(backend, C[l], f"L{l}", lk["size"] if (l >= lk_from and b == 0)
                                                    else None) for b in range(nb)]))
            if l < L - 1:
                self.down.append(Down(backend, C[l], C[l + 1], f"D{l}"))
        self.up = nn.ModuleList([Up(backend, C[l + 1], C[l], f"D{l}") for l in range(L - 1)])
        self.fuse = nn.ModuleList([nn.Linear(2 * C[l], C[l]) for l in range(L - 1)])
        self.dec = nn.ModuleList([nn.ModuleList([ResBlock(backend, C[l], f"L{l}") for _ in range(nb)])
                                  for l in range(L - 1)])
        cd = pcfg["ctx_dim"]
        self.ctx = mlp(2 * C[-1] + N_EVENT_FEATURES, cd, cd)
        self.use_film = pcfg["film"]
        if self.use_film:
            self.film = nn.ModuleList([FiLM(cd, C[l]) for l in range(L)])
        pb = pcfg["point_branch"]
        self.use_points = pb["enabled"]
        D = pb["dim"]
        if self.use_points:
            self.pt_bn = nn.BatchNorm1d(n_pt_feat, eps=0.1)
            self.pt_in = mlp(n_pt_feat, D, D)
            self.v2p = nn.ModuleList([nn.Linear(C[l], D) for l in range(L)])
            self.p2v = nn.ModuleList([nn.Linear(D, C[l]) for l in range(L)])
            for m in self.p2v:
                nn.init.zeros_(m.weight)
                nn.init.zeros_(m.bias)
        h = pcfg["heads"]
        tg = pcfg["target"]
        C0 = C[0]
        self.head_heat = mlp(C0, C0, 1)
        self.head_any = mlp(C0, C0, 1)
        self.heads_cfg = h
        self.offset_scale = float(tg["offset_radius"]) / 2.0
        if h["offset"]:
            self.head_off = mlp(C0, C0, 6, zero_last=True)
        if h["cluster"] or h["listwise"]:
            self.head_cluster = mlp(2 * C0, C0, 2)            # [P(host) logit, C3 listwise cluster score]
        if h["event"]:
            self.head_event = mlp(cd, cd, 1)
        K = len(tg["shell_edges"])
        if h["point"] and self.use_points:
            self.head_point = mlp(D, D, 1 + K)
        if h["shell_cone"]:
            self.head_shell = mlp(C0, C0, K)
            self.head_cone = mlp(C0, C0, 4)
        if h["listwise"]:
            self.head_lw = mlp(C0, C0, 1)
        if h["deep_supervision"]:
            self.head_deep = nn.ModuleList([nn.Linear(C[l], 1) for l in range(1, L - 1)])
        # CenterNet prior: sparse positives -> start the vertexness logits at P = 0.01
        prior = -4.595
        for m in [self.head_heat[-1], self.head_any[-1]] + (list(self.head_deep) if h["deep_supervision"] else []):
            nn.init.constant_(m.bias, prior)
        if h["point"] and self.use_points:
            with torch.no_grad():
                self.head_point[-1].bias[0] = prior

    # ------------------------------------------------------------------ helpers
    @staticmethod
    def _pool_event(h, batch, B):
        mean = scatter_mean(h, batch, dim=0, dim_size=B)
        mx = scatter_max(h, batch, dim=0, dim_size=B)[0]
        return torch.cat([mean, mx], 1)

    def _level_rows(self, st, grid0, l):
        """(N0,) row of each level-0 voxel in the level-l tensor."""
        q = grid0.clone()
        q[:, 1:] = torch.div(q[:, 1:], 2 ** l, rounding_mode="floor")
        rows = sparse.lookup_rows(q, st)
        if bool((rows < 0).any()):
            raise RuntimeError(f"level {l}: {int((rows < 0).sum())} voxels without a parent")
        return rows

    # ------------------------------------------------------------------ forward
    def forward(self, b):
        cdt = self.cdtype
        grid0 = b["vox_coords"].long()                         # (N0, 4) batch, grid
        B = int(b["batch_size"])
        x = self.in_bn(b["vox_feat"].float()).clamp(-20.0, 20.0)
        st = sparse.STensor.create(x.to(cdt), grid0, b["spatial_shape"], B, self.backend)
        st = self.stem(st)
        st = st.replace(st.feats.float())
        skips = []
        for l in range(self.L):
            for blk in self.enc[l]:
                st = blk(st, cdt)
            skips.append(st)
            if l < self.L - 1:
                st = self.down[l](st, cdt)
        # levels: rows of each level-0 voxel at every level (for devoxelisation and level coordinates)
        v2l = [torch.arange(len(grid0), device=grid0.device)] + [self._level_rows(skips[l], grid0, l)
                                                                   for l in range(1, self.L)]
        bott = skips[-1]
        bb = bott.batch
        ev = b["ev_feat"].float()          # bounded logs / flags; no BatchNorm (it breaks on a batch of one event)
        ctx = self.ctx(torch.cat([self._pool_event(bott.feats.float(), bb, B), ev], 1)).float()
        pt = None
        if self.use_points:
            pt = self.pt_in(self.pt_bn(b["pt_feat"].float()).clamp(-20.0, 20.0)).float()
            p2v0 = b["pt2vox"].long()

        def fuse_level(h, l, st_l, pt):
            if self.use_film:
                h = self.film[l](h, ctx, st_l.batch)
            if pt is not None:
                rows = v2l[l][p2v0]
                pt = pt + self.v2p[l](h[rows])
                h = h + self.p2v[l](scatter_mean(pt, rows, dim=0, dim_size=len(h)))
            return h, pt

        h, pt = fuse_level(bott.feats.float(), self.L - 1, bott, pt)
        st = bott.replace(h)
        deep = []
        for l in range(self.L - 2, -1, -1):
            upst = self.up[l](st, cdt)
            skip = skips[l]
            h = self.fuse[l](torch.cat([upst.feats.float(), skip.feats.float()], 1))
            st = skip.replace(h)
            for blk in self.dec[l]:
                st = blk(st, cdt)
            h, pt = fuse_level(st.feats.float(), l, st, pt)
            st = st.replace(h)
            if l >= 1 and self.heads_cfg["deep_supervision"]:
                deep.append((l, st, self.head_deep[l - 1](h).squeeze(1)))
        h0 = st.feats.float()
        # level-0 tensor rows are in input order for both backends (submanifold keeps indices; checked in tests)
        out = {"coord": b["vox_coord"].float(), "batch": grid0[:, 0], "ctx": ctx}
        out["heat"] = self.head_heat(h0).squeeze(1).float()
        out["any"] = self.head_any(h0).squeeze(1).float()
        hc = self.heads_cfg
        if hc["offset"]:
            with torch.autocast(device_type=h0.device.type, enabled=False):
                raw = self.head_off(h0.float())
            out["offset"] = raw[:, :3] * self.offset_scale
            out["log_scale"] = raw[:, 3:].clamp(-6.0, 4.0)
        if hc["cluster"] or hc["listwise"]:
            cl = b["vox_cluster"].long()
            key = grid0[:, 0] * (int(cl.max().item()) + 2 if len(cl) else 1) + (cl + 1)
            uk, inv = torch.unique(key, return_inverse=True)
            pooled = torch.cat([scatter_mean(h0, inv, dim=0, dim_size=len(uk)),
                                scatter_max(h0, inv, dim=0, dim_size=len(uk))[0]], 1)
            s = self.head_cluster(pooled).float()
            cb = torch.zeros(len(uk), dtype=torch.long, device=h0.device).scatter_(0, inv, grid0[:, 0])
            cid = torch.zeros(len(uk), dtype=torch.long, device=h0.device).scatter_(0, inv, cl)
            out.update(cluster_logit=s[:, 0], cluster_lw=s[:, 1], cluster_batch=cb, cluster_id=cid,
                       vox_group=inv)
        if hc["event"]:
            out["event_logit"] = self.head_event(ctx).squeeze(1).float()
        if hc["point"] and self.use_points:
            ph = self.head_point(pt).float()
            out["point_logit"] = ph[:, 0]
            out["point_shell"] = ph[:, 1:]
            out["pt_xyz"] = b["pt_xyz"].float()
            out["pt_batch"] = grid0[b["pt2vox"].long(), 0]
        if hc["shell_cone"]:
            out["shell"] = self.head_shell(h0).float()
            out["cone"] = self.head_cone(h0).float()
        if hc["listwise"]:
            out["lw"] = self.head_lw(h0).squeeze(1).float()
        if deep:
            dl = []
            for l, st_l, logit in deep:
                rows = v2l[l]
                n_l = len(st_l.feats)
                coord_l = scatter_mean(b["vox_coord"].float(), rows, dim=0, dim_size=n_l)
                dl.append({"level": l, "coord": coord_l, "batch": st_l.batch, "logit": logit.float()})
            out["deep"] = dl
        if self.return_features:
            out["feat0"] = h0
        return out

    def set_compute_dtype(self, dtype):
        self.cdtype = dtype
