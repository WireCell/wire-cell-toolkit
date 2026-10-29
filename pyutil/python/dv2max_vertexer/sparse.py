"""Sparse 3D convolutions with two interchangeable backends.  No attention anywhere.
(dv2max_vertexer: vendored from dv2max/sparse.py with the spconv adapters removed; only the torch backend
below is used, with the parameter names / shapes of spconv-trained state_dicts.)

Backends (config `backend`, 'auto' = spconv if importable):
  'spconv'  spconv 2.x (spconv-cu121).  On CUDA the module's default implicit-GEMM algorithm is used; on CPU the
            'Native' algorithm (the only one spconv runs without a GPU; measured to work on the login node).
  'torch'   in-repo fallback in plain PyTorch (hash lookup with sort + searchsorted, per-offset matmul and
            index_add).  Same parameter names, shapes and conventions as spconv, so a state_dict trained with
            one backend loads into the other (tests/test_sparse.py checks equal outputs on CPU).
Conventions (measured on spconv 2.3.8, see tests): weights are [Cout, kx, ky, kz, Cin] for every algorithm;
  submanifold 3^3:   out[p]      = sum_o W[o + 1] x[p + o]           (o in {-1,0,1}^3)
  down k=2 s=2:      out[q]      = sum_o W[o] x[2q + o]               (o in {0,1}^3)
  inverse (up):      out[2q + o] = W[o] x[q]                          (fine sites of the matching down conv)
Coordinates are (batch, x, y, z) integer grid indices; spatial_shape is padded to a multiple of 2**(levels-1)
so every stride-2 level keeps all voxels (spconv drops the last odd index otherwise).

LargeKernelPartition3d: LargeKernel3D-style spatial-partition convolution (Chen et al. 2023) over a size^3
window: offsets are grouped into 3 bins per axis ({<=-2}, {-1,0,1}, {>=2}), 27 shared weight matrices, summed
over the neighbours of each group; a position-specific 3^3 submanifold conv runs in parallel.  It is used only
at the coarsest levels (few voxels), implemented with explicit neighbour pairs in torch for both backends.
"""
import math

import torch
from torch import nn

# dv2max_vertexer: torch backend only (spconv is not installed in the in-chain environment)
spconv = None
ConvAlgo = None
HAVE_SPCONV = False

OFFSETS3 = torch.tensor([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)], dtype=torch.long)
OFFSETS2 = torch.tensor([[i, j, k] for i in (0, 1) for j in (0, 1) for k in (0, 1)], dtype=torch.long)


# ---------------------------------------------------------------------------- tensor wrapper
class STensor:
    """Backend-agnostic sparse tensor.  feats (N,C); coords (N,4) int (b,x,y,z); spatial_shape [3]; batch_size.
    `cache` is shared by every tensor derived at the same resolution (neighbour maps, level lookups)."""

    def __init__(self, feats, coords, spatial_shape, batch_size, backend, sp=None, cache=None):
        self.backend = backend
        self._feats = feats
        self.coords = coords
        self.spatial_shape = [int(s) for s in spatial_shape]
        self.batch_size = int(batch_size)
        self.sp = sp
        self.cache = {} if cache is None else cache

    @classmethod
    def create(cls, feats, coords, spatial_shape, batch_size, backend):
        coords = coords.to(feats.device)
        if backend == "spconv":
            sp = spconv.SparseConvTensor(feats, coords.int().contiguous(), list(spatial_shape), int(batch_size))
            return cls(feats, sp.indices, spatial_shape, batch_size, backend, sp=sp)
        return cls(feats, coords.long().contiguous(), spatial_shape, batch_size, backend)

    @property
    def feats(self):
        return self.sp.features if self.sp is not None else self._feats

    def replace(self, feats):
        if self.sp is not None:
            sp = self.sp.replace_feature(feats)
            return STensor(feats, self.coords, self.spatial_shape, self.batch_size, self.backend, sp, self.cache)
        return STensor(feats, self.coords, self.spatial_shape, self.batch_size, self.backend, None, self.cache)

    @property
    def batch(self):
        return self.coords[:, 0].long()

    def __len__(self):
        return len(self.coords)


def grid_keys(coords, spatial_shape):
    """int64 keys of (b,x,y,z) rows; -1 where any spatial index is outside [0, shape)."""
    c = coords.long()
    s0, s1, s2 = (int(v) for v in spatial_shape)
    ok = (c[:, 1] >= 0) & (c[:, 1] < s0) & (c[:, 2] >= 0) & (c[:, 2] < s1) & (c[:, 3] >= 0) & (c[:, 3] < s2)
    k = ((c[:, 0] * s0 + c[:, 1]) * s1 + c[:, 2]) * s2 + c[:, 3]
    return torch.where(ok, k, torch.full_like(k, -1))


class KeyIndex:
    """Sorted key table for exact row lookup of grid coordinates."""

    def __init__(self, coords, spatial_shape):
        self.shape = spatial_shape
        keys = grid_keys(coords, spatial_shape)
        self.sorted, self.perm = torch.sort(keys)

    def lookup(self, query_coords):
        """Row index of each query coordinate (-1 if absent / out of range)."""
        q = grid_keys(query_coords, self.shape)
        if len(self.sorted) == 0:
            return torch.full_like(q, -1)
        pos = torch.searchsorted(self.sorted, q).clamp(max=len(self.sorted) - 1)
        hit = (self.sorted[pos] == q) & (q >= 0)
        return torch.where(hit, self.perm[pos], torch.full_like(q, -1))


def lookup_rows(query_coords, st):
    """Rows of `st` holding the given (b,x,y,z) coordinates (-1 if absent).  The key table is cached on st."""
    ki = st.cache.get("_keyindex")
    if ki is None or ki[0] is not st.coords:
        ki = (st.coords, KeyIndex(st.coords, st.spatial_shape))
        st.cache["_keyindex"] = ki
    return ki[1].lookup(query_coords.to(st.coords.device))


# ---------------------------------------------------------------------------- torch backend ops
def _kaiming(weight, fan_in):
    nn.init.normal_(weight, 0.0, math.sqrt(2.0 / max(1, fan_in)))


class TorchSubMConv3d(nn.Module):
    """Submanifold k^3 convolution (k odd), spconv conventions, plain torch."""

    def __init__(self, cin, cout, k=3, key=None, bias=False):
        super().__init__()
        self.cin, self.cout, self.k, self.key = cin, cout, k, key
        self.weight = nn.Parameter(torch.empty(cout, k, k, k, cin))
        _kaiming(self.weight, cin * k ** 3)
        self.bias = nn.Parameter(torch.zeros(cout)) if bias else None
        r = k // 2
        self.register_buffer("offsets", torch.tensor([[i, j, l] for i in range(-r, r + 1) for j in range(-r, r + 1)
                                                      for l in range(-r, r + 1)], dtype=torch.long), persistent=False)

    def neighbours(self, st):
        ck = ("subm_nbr", self.k, self.key)
        if ck not in st.cache:
            ki = KeyIndex(st.coords, st.spatial_shape)
            c = st.coords.long()
            nb = []
            for o in self.offsets.to(c.device):
                q = c.clone()
                q[:, 1:] += o
                nb.append(ki.lookup(q))
            st.cache[ck] = torch.stack(nb, 0)          # (K, N)
        return st.cache[ck]

    def forward(self, st):
        x = st.feats
        nb = self.neighbours(st)
        W = self.weight.reshape(self.cout, -1, self.cin)          # (Cout, K, Cin), K in (i, j, l) order
        out = x.new_zeros(len(x), self.cout, dtype=torch.float32)  # fp32 accumulation (bf16 / fp16 inputs)
        for kk in range(nb.shape[0]):
            j = nb[kk]
            m = j >= 0
            if not bool(m.any()):
                continue
            rows = torch.nonzero(m).squeeze(1)
            out.index_add_(0, rows, (x[j[rows]] @ W[:, kk, :].t().to(x.dtype)).to(out.dtype))
        if self.bias is not None:
            out = out + self.bias.float()
        return st.replace(out.to(x.dtype))


class TorchDownConv3d(nn.Module):
    """k=2, s=2 sparse convolution (output voxels = occupied parents)."""

    def __init__(self, cin, cout, key=None, bias=False):
        super().__init__()
        self.cin, self.cout, self.key = cin, cout, key
        self.weight = nn.Parameter(torch.empty(cout, 2, 2, 2, cin))
        _kaiming(self.weight, cin * 8)
        self.bias = nn.Parameter(torch.zeros(cout)) if bias else None

    def forward(self, st):
        c = st.coords.long()
        parent = c.clone()
        parent[:, 1:] = torch.div(c[:, 1:], 2, rounding_mode="floor")
        off = (c[:, 1] - 2 * parent[:, 1]) * 4 + (c[:, 2] - 2 * parent[:, 2]) * 2 + (c[:, 3] - 2 * parent[:, 3])
        shape = [(s + 1) // 2 for s in st.spatial_shape]
        keys = grid_keys(parent, shape)
        ukeys, inv = torch.unique(keys, sorted=True, return_inverse=True)
        M = len(ukeys)
        pc = torch.zeros(M, 4, dtype=torch.long, device=c.device)
        pc[inv] = parent
        x = st.feats
        W = self.weight.reshape(self.cout, 8, self.cin)
        out = x.new_zeros(M, self.cout, dtype=torch.float32)
        for o in range(8):
            rows = torch.nonzero(off == o).squeeze(1)
            if len(rows):
                out.index_add_(0, inv[rows], (x[rows] @ W[:, o, :].t().to(x.dtype)).to(out.dtype))
        if self.bias is not None:
            out = out + self.bias.float()
        st.cache[("down", self.key)] = (inv, off)
        return STensor(out.to(x.dtype), pc, shape, st.batch_size, "torch", None, {"_fine": st})


class TorchInverseConv3d(nn.Module):
    """Inverse of the TorchDownConv3d with the same key: back to the fine voxels of that conv's input."""

    def __init__(self, cin, cout, key=None, bias=False):
        super().__init__()
        self.cin, self.cout, self.key = cin, cout, key
        self.weight = nn.Parameter(torch.empty(cout, 2, 2, 2, cin))
        _kaiming(self.weight, cin)
        self.bias = nn.Parameter(torch.zeros(cout)) if bias else None

    def forward(self, st):
        fine = st.cache["_fine"]
        inv, off = fine.cache[("down", self.key)]
        x = st.feats
        W = self.weight.reshape(self.cout, 8, self.cin)
        out = x.new_zeros(len(inv), self.cout)
        for o in range(8):
            rows = torch.nonzero(off == o).squeeze(1)
            if len(rows):
                out.index_copy_(0, rows, (x[inv[rows]] @ W[:, o, :].t().to(x.dtype)).to(out.dtype))
        if self.bias is not None:
            out = out + self.bias.to(out.dtype)
        return fine.replace(out)


# ---------------------------------------------------------------------------- torch-backend wrappers (same names)
class _TorchWrap(nn.Module):
    def __init__(self, conv):
        super().__init__()
        self.conv = conv

    def forward(self, st):
        return self.conv(st)


def subm(backend, cin, cout, key, k=3, bias=False):
    if backend != "torch":
        raise ValueError(f"dv2max_vertexer has only the torch backend, got {backend!r}")
    return _TorchWrap(TorchSubMConv3d(cin, cout, k, key, bias))


def down(backend, cin, cout, key, bias=False):
    if backend != "torch":
        raise ValueError(f"dv2max_vertexer has only the torch backend, got {backend!r}")
    return _TorchWrap(TorchDownConv3d(cin, cout, key, bias))


def up(backend, cin, cout, key, bias=False):
    if backend != "torch":
        raise ValueError(f"dv2max_vertexer has only the torch backend, got {backend!r}")
    return _TorchWrap(TorchInverseConv3d(cin, cout, key, bias))


# ---------------------------------------------------------------------------- large kernel (both backends)
def partition_index(offset, size):
    """Group index 0..26 of integer offsets (…,3) in a size^3 window: per axis {<=-2: 0, -1..1: 1, >=2: 2}."""
    r = size // 2
    if torch.any(offset.abs() > r):
        raise ValueError("offset outside the kernel window")
    g = torch.where(offset <= -2, 0, torch.where(offset >= 2, 2, 1))
    return g[..., 0] * 9 + g[..., 1] * 3 + g[..., 2]


def window_pairs(coords, radius, chunk=4096):
    """All (i, j) row pairs of the same batch whose Chebyshev grid distance is <= radius (i == j included).
    Returns (i, j, offset (P,3) = c_j - c_i)."""
    c = coords.long()
    b = c[:, 0]
    I, J = [], []
    for bb in torch.unique(b):
        rows = torch.nonzero(b == bb).squeeze(1)
        xyz = c[rows, 1:]
        for a in range(0, len(rows), chunk):
            q = xyz[a:a + chunk]
            d = (q[:, None, :] - xyz[None, :, :]).abs().amax(-1)
            ii, jj = torch.nonzero(d <= radius, as_tuple=True)
            I.append(rows[a + ii])
            J.append(rows[jj])
    if not I:
        e = torch.zeros(0, dtype=torch.long, device=c.device)
        return e, e, torch.zeros(0, 3, dtype=torch.long, device=c.device)
    i, j = torch.cat(I), torch.cat(J)
    return i, j, c[j, 1:] - c[i, 1:]


class LargeKernelPartition3d(nn.Module):
    """Spatial-partition large-kernel conv (27 shared weights over a size^3 window) + parallel 3^3 submanifold
    conv (position-specific centre).  Output = LK(x) + SubM3(x)."""

    def __init__(self, backend, channels, key, size=13):
        super().__init__()
        self.size = size
        self.key = key
        self.weight = nn.Parameter(torch.empty(27, channels, channels))
        nn.init.normal_(self.weight, 0.0, math.sqrt(1.0 / (27 * channels)))
        self.small = subm(backend, channels, channels, key)

    def pairs(self, st):
        ck = ("lk_pairs", self.size)
        if ck not in st.cache:
            i, j, off = window_pairs(st.coords, self.size // 2)
            st.cache[ck] = (i, j, partition_index(off, self.size))
        return st.cache[ck]

    def forward(self, st):
        x = st.feats
        i, j, g = self.pairs(st)
        N, C = x.shape
        A = x.new_zeros(N * 27, C, dtype=torch.float32)          # fp32 sums over up to 5^3 neighbours per group
        A.index_add_(0, i * 27 + g, x[j].float())
        lk = A.view(N, 27 * C).to(x.dtype) @ self.weight.reshape(27 * C, -1).to(x.dtype)
        sm = self.small(st)
        return sm.replace(sm.feats + lk.to(sm.feats.dtype))
