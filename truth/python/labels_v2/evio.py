"""labels_v2.evio -- one campaign event's inputs for the v2 labeler (pr_ml doc 54 sec 1).

Loads, read only:
  deposits   sim-depos.tar.{bz2,zst,gz,} : gen-0 drifted rows (t ns, q e, x y z mm, sigma_L, sigma_T) and, through the gen-1
             prior, the true position P3 and G4 time T0 of each (as d50d_real.load_depos; one prior level)
  truth      truth_<k>.npz / truth.npz (truth_table.py): truth_pf, sed_ids (tid, orig), sed_energy, truth_nu, gen_labels
  pctree     mlpr_v1/pr-pctree.tar.gz: the /live tree's "3d" and "truth" clouds in train_targets.read_event's order
  settings   mlpr_v1/glabels/<event>.json "inputs" (wires basename, depo_toff_us, reco_toff_us, speed_mmus, tie_order) and
             mlpr_v1/truth.json (time_offset_us, speed_mmus, extra.real_assoc smear widths, extra.extra_smear_long)
  frames     sp-frames.tar* or sim-frames-anode*.tar*: tickinfo and chanmask_bad; else chanmask-apa*.json (static), else none
  wires      <wires dir>/<basename of inputs.wires> (validate_2d.load_wires)
Both event layouts (plain; event.tar.zst + repack.txt) are read: packed members are streamed once into a scratch dir and
checked against repack.txt (a copy of pr_ml/s0/event_files.py's reader, kept local so the toolkit package has no wcp import).
"""
import glob, hashlib, io, json, os, re, shutil, subprocess, sys, tarfile, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TP = os.path.dirname(HERE)                     # truth/python
if TP not in sys.path:
    sys.path.insert(0, TP)
import truth_table as TT          # noqa: E402
import validate_2d as V2          # noqa: E402  (load_wires, assign_faces)
from xres_check import load_faces, assign_depo_faces   # noqa: E402

ZSTD = os.environ.get("ZSTD") or shutil.which("zstd") or "/nfs/data/1/xqian/toolkit-dev/local/bin/zstd"
WIRES_DIR_DEFAULT = os.environ.get("LABELS_V2_WIRES_DIR", "/home/xqian/toolkit-dev/wire-cell-data")
US = 1000.0


# ----------------------------------------------------------------------------------------------- event layouts
class EventFiles:
    """one event dir, plain or packed (event.tar.zst + repack.txt); path(rel) gives a readable file"""
    ARCH, REC = "event.tar.zst", "repack.txt"

    def __init__(self, d, tmpdir=None):
        self.d = d.rstrip("/")
        self.packed = os.path.isfile(os.path.join(self.d, self.ARCH))
        self.rec = self._read_rec() if self.packed else {}
        self.tmpdir = tmpdir or os.environ.get("TMPDIR", "/home/xqian/tmp")
        self._tmp = None; self._got = {}

    def _read_rec(self):
        rec = {}
        for l in open(os.path.join(self.d, self.REC)):
            if l.startswith("#") or l.startswith(("verified ", "complete ")) or not l.strip():
                continue
            k, s, n, r = l.rstrip("\n").split(" ", 3)
            rec[r] = (k, s, int(n))
        return rec

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()

    def close(self):
        if self._tmp:
            shutil.rmtree(self._tmp, ignore_errors=True); self._tmp = None
        self._got = {}

    def has(self, rel):
        p = os.path.join(self.d, rel)
        return rel in self.rec or (os.path.isfile(p) and not (self.packed and os.path.islink(p)))

    def find(self, prefix):
        if self.packed:
            return sorted(r for r in self.rec if r.startswith(prefix))
        return sorted(os.path.relpath(p, self.d) for p in glob.glob(os.path.join(self.d, glob.escape(prefix) + "*"))
                      if os.path.isfile(p) and not os.path.islink(p))

    def fetch(self, rels):
        if not self.packed:
            return
        want = {}
        for r in rels:
            if r in self._got:
                continue
            if r not in self.rec:
                raise FileNotFoundError(f"{r} not in {self.d}/{self.REC}")
            want[r + ".gunzipped" if self.rec[r][0] == "gz" else r] = r
        if not want:
            return
        if self._tmp is None:
            self._tmp = tempfile.mkdtemp(prefix="lv2.", dir=self.tmpdir)
        z = subprocess.Popen([ZSTD, "-q", "-dc", "--long=31", os.path.join(self.d, self.ARCH)], stdout=subprocess.PIPE)
        try:
            with tarfile.open(fileobj=z.stdout, mode="r|") as t:
                for m in t:
                    rel = want.get(m.name)
                    if rel is None or not m.isfile():
                        continue
                    dst = os.path.join(self._tmp, m.name); os.makedirs(os.path.dirname(dst), exist_ok=True)
                    h = hashlib.sha256(); src = t.extractfile(m)
                    with open(dst, "wb") as o:
                        for b in iter(lambda: src.read(1 << 20), b""):
                            h.update(b); o.write(b)
                    if h.hexdigest() != self.rec[rel][1]:
                        raise RuntimeError(f"{self.d}: member {m.name} differs from {self.REC}")
                    self._got[rel] = dst
        finally:
            z.stdout.close(); rc = z.wait()
        if rc not in (0, -13):
            raise RuntimeError(f"zstd -d failed ({rc}) on {self.d}")
        miss = [r for r in want.values() if r not in self._got]
        if miss:
            raise FileNotFoundError(f"{self.d}: {miss[0]} listed in {self.REC} but not in {self.ARCH}")

    def path(self, rel):
        if not self.packed:
            p = os.path.join(self.d, rel)
            if not os.path.isfile(p):
                raise FileNotFoundError(p)
            return p
        self.fetch([rel]); return self._got[rel]

    def stream_sha(self, rel):
        if self.packed:
            return self.rec[rel][1]
        return sha_file(os.path.join(self.d, rel), rel.endswith(".tar.gz"))


def sha_file(p, gz=False):
    h = hashlib.sha256()
    if gz:
        with subprocess.Popen(["gzip", "-dc", p], stdout=subprocess.PIPE) as g:
            for b in iter(lambda: g.stdout.read(1 << 20), b""):
                h.update(b)
    else:
        with open(p, "rb") as f:
            for b in iter(lambda: f.read(1 << 20), b""):
                h.update(b)
    return h.hexdigest()


def resolve_truth(inp_truth, root, det, ev, files=None, stagea_root=None):
    """the event's truth table in either layout or under stageA (both namings); -> (path, identity)"""
    if files is not None and files.has("truth.npz"):
        return files.path("truth.npz"), "truth.npz"
    for r in [root] + ([stagea_root] if stagea_root else []):
        p = re.sub(r"^/prml/", r.rstrip("/") + "/", inp_truth or "")
        cands = [p] if p else []
        m = re.match(r"^(.*)_e(\d+)$", ev)
        if m:
            cands.append(f"{r}/{det}/stageA/{m.group(1)}/truth_{m.group(2)}.npz")
        cands.append(f"{r}/{det}/stageA/{ev}/truth.npz")
        for c in cands:
            if os.path.isfile(c) and not (os.path.islink(c) and not os.path.exists(c)):
                return c, "stageA/" + os.path.basename(os.path.dirname(c)) + "/" + os.path.basename(c)
    raise FileNotFoundError(f"no truth table for {det} {ev} (asked {inp_truth})")


# ----------------------------------------------------------------------------------------------- tars
def open_tar(fn):
    """tarfile of .tar / .tar.gz / .tar.bz2 / .tar.zst (zst streamed through zstd into memory)"""
    if fn.endswith(".zst"):
        b = subprocess.run([ZSTD, "-dcq", fn], check=True, capture_output=True).stdout
        return tarfile.open(fileobj=io.BytesIO(b))
    return tarfile.open(fn)


def load_depos(fn):
    """-> D (n, 7) gen-0 rows, I (n, 4), P3 true position (mm), T0 G4 time (ns); one prior level (d50d_real.load_depos)"""
    with open_tar(fn) as t:
        mem = {m.name: m for m in t.getmembers()}
        D = np.load(io.BytesIO(t.extractfile(mem["depo_data_0.npy"]).read()), allow_pickle=False).astype(np.float64)
        I = np.load(io.BytesIO(t.extractfile(mem["depo_info_0.npy"]).read()), allow_pickle=False).astype(np.int64)
    gen = I[:, 2]
    act = np.flatnonzero(gen == 0); pri = np.flatnonzero(gen == 1)
    ch = I[pri, 3]; ok = (ch >= 0) & (ch < len(D))
    prior = np.full(len(D), -1, np.int64); prior[ch[ok]] = pri[ok]
    pa = prior[act]; have = pa >= 0
    P3 = np.where(have[:, None], D[np.maximum(pa, 0), 2:5], D[act, 2:5])
    T0 = np.where(have, D[np.maximum(pa, 0), 0], D[act, 0])
    return D[act], I[act], P3, T0


NTICKS_DEFAULT = {"sbnd": 3427, "pdhd": 6000, "pdvd": 6400, "fdhd": 6000, "fdvd": 8500}   # the g25 frames (doc 54 sec 2); used when no frames are kept


def load_frames_info(fns):
    """tickinfo (3,) and the tick count of the SP charge frame (dnnsp else gauss), and the bad channels {channel: [(t0, t1), ...]}"""
    ti, bad, nfr, nticks = None, {}, 0, None
    for fn in fns:
        with open_tar(fn) as t:
            mem = {m.name: m for m in t.getmembers()}
            tag = "dnnsp" if any(n.startswith("frame_dnnsp") for n in mem) else "gauss"
            for n in sorted(mem):
                if re.match(r"tickinfo_" + tag + r"\d*_\d+\.npy", n):
                    ti = np.load(io.BytesIO(t.extractfile(mem[n]).read()), allow_pickle=False); nfr += 1
                    h = t.extractfile(mem[n.replace("tickinfo_", "frame_")]); v = np.lib.format.read_magic(h)      # the header only: (nchannels, nticks)
                    sh = (np.lib.format.read_array_header_1_0(h) if v == (1, 0) else np.lib.format.read_array_header_2_0(h))[0]
                    nticks = int(sh[1]) if nticks is None else max(nticks, int(sh[1]))
            for n in mem:
                if n.startswith("chanmask_bad"):
                    b = np.load(io.BytesIO(t.extractfile(mem[n]).read()), allow_pickle=False)
                    for c, t0, t1 in np.asarray(b).reshape(-1, 3).tolist():
                        bad.setdefault(int(c), []).append((int(t0), int(t1)))
    return ti, bad, nfr, nticks


def load_static_bad(files):
    """chanmask-apa*.json (SBND): {"chanmask_bad_<id>": [[ch, t0, t1], ...]} -> {channel: [(t0, t1)]}"""
    bad = {}
    for rel in files.find("chanmask-apa"):
        if not rel.endswith(".json"):
            continue
        d = json.load(open(files.path(rel)))
        for k, rows in d.items():
            if not k.startswith("chanmask_bad"):
                continue
            for c, t0, t1 in rows:
                bad.setdefault(int(c), []).append((int(t0), int(t1)))
    return bad


# ----------------------------------------------------------------------------------------------- the pctree points
def load_points(pctree):
    """the PR points in train_targets.read_event's order (nodes ascending, points in node order): X (reco frame, mm),
    cluster index per point, and the truth cloud columns tid, orig_tid, pf_node, cls, instance_valid, d_traj (defaults
    -3 / -3 / -1 / -1 / 0 / 0 where a node has no truth cloud)"""
    from steiner_truth import load_tree, per_node
    tree = load_tree(pctree)
    P3, TR = per_node(tree, "3d"), per_node(tree, "truth")
    co = tree["cluster_of"]
    cols = {k: [] for k in ("x", "y", "z", "cl", "node", "tid", "orig_tid", "pf_node", "cls", "instance_valid", "d_traj")}
    for nd in sorted(P3):
        p3 = P3[nd]; m = len(p3["x"]); tr = TR.get(nd)
        g = (lambda k, dflt: np.asarray(tr[k]) if tr is not None and k in tr else np.full(m, dflt))
        for k, v in (("x", p3["x"]), ("y", p3["y"]), ("z", p3["z"]), ("cl", np.full(m, int(co[nd]))), ("node", np.full(m, nd)),
                     ("tid", g("tid", -3)), ("orig_tid", g("orig_tid", -3)), ("pf_node", g("pf_node", -1)), ("cls", g("cls", -1)),
                     ("instance_valid", g("instance_valid", 0)), ("d_traj", g("d_traj", 0.0))):
            cols[k].append(np.asarray(v))
    P = {k: (np.concatenate(v) if v else np.zeros(0)) for k, v in cols.items()}
    out = dict(X=np.column_stack([P.pop("x"), P.pop("y"), P.pop("z")]).astype(np.float64) if len(P["cl"]) else np.zeros((0, 3)))
    for k in ("cl", "node", "tid", "orig_tid", "pf_node", "cls", "instance_valid"):
        out[k] = P[k].astype(np.int64)
    out["d_traj"] = P["d_traj"].astype(np.float64)
    return out


# ----------------------------------------------------------------------------------------------- the event
class Event:
    """everything the labeler reads for one event, loaded once"""

    @classmethod
    def synthetic(cls, det, wires, D, I, P3, T0, truth, points, inputs, tj, tickinfo=(0.0, 500.0, 0.0), bad=None):
        """an event from arrays (the selftests): D (n, 7) drifted rows, I (n, 4) [id, pdg, gen 0, child], P3 true positions,
        T0 G4 times, truth = truth_table.read_truth_npz-like dict (truth_pf, sed_ids, sed_energy, gen_labels, truth_nu),
        points = dict(X, tid, orig_tid, ...), inputs = dict(depo_toff_us, reco_toff_us, speed_mmus), tj = truth.json dict"""
        self = cls.__new__(cls)
        self.d = f"/synthetic/{det}/synthetic"; self.ev = "synthetic"; self.det = det; self.root = "/synthetic"
        self.files = None; self.G = {"inputs": inputs}; self.inputs = inputs; self.tj = tj; self.wires = wires
        self.depo_file = self.truth_file = self.pctree_file = None; self.truth_id = "synthetic"
        self.tickinfo = np.asarray(tickinfo, np.float64); self.bad = bad or {}; self.bad_source = "synthetic"
        self.nticks = NTICKS_DEFAULT.get(det, 6000); self.nticks_source = "default"
        self.D, self.I, self.P3, self.T0 = np.asarray(D, np.float64), np.asarray(I, np.int64), np.asarray(P3, np.float64), np.asarray(T0, np.float64)
        self.truth = truth
        self._finish()
        n = len(points["X"])
        self.pts = dict(X=np.asarray(points["X"], np.float64), cl=np.asarray(points.get("cl", np.zeros(n)), np.int64), node=np.asarray(points.get("node", np.zeros(n)), np.int64),
                        tid=np.asarray(points["tid"], np.int64), orig_tid=np.asarray(points["orig_tid"], np.int64), pf_node=np.asarray(points.get("pf_node", np.full(n, -1)), np.int64),
                        cls=np.asarray(points.get("cls", np.full(n, -1)), np.int64), instance_valid=np.asarray(points.get("instance_valid", np.ones(n)), np.int64),
                        d_traj=np.asarray(points.get("d_traj", np.zeros(n)), np.float64))
        self.meta_inputs = {"synthetic": "synthetic"}
        return self

    def __init__(self, d, root=None, wires_dir=None, stagea_root=None, tmpdir=None):
        self.d = d.rstrip("/")
        self.ev = os.path.basename(self.d); self.det = os.path.basename(os.path.dirname(self.d))
        self.root = root or os.path.dirname(os.path.dirname(self.d))
        self.wires_dir = wires_dir or WIRES_DIR_DEFAULT
        self.files = EventFiles(self.d, tmpdir)
        F = self.files
        gj = f"mlpr_v1/glabels/{self.ev}.json"
        need = [r for r in (gj, "mlpr_v1/truth.json", "mlpr_v1/pr-pctree.tar.gz", "truth.npz", *F.find("sim-depos.tar"),
                            *F.find("sp-frames.tar"), *F.find("sim-frames-anode"), *F.find("chanmask-apa")) if F.has(r)]
        F.fetch(need)
        self.G = json.load(open(F.path(gj)))
        if "error" in self.G:
            raise RuntimeError("glabels: " + str(self.G["error"])[:80])
        self.inputs = self.G["inputs"]
        self.tj = json.load(open(F.path("mlpr_v1/truth.json")))
        self.wires = os.path.join(self.wires_dir, os.path.basename(self.inputs["wires"]))
        self.depo_file = F.path(F.find("sim-depos.tar")[0])
        self.truth_file, self.truth_id = resolve_truth(self.inputs.get("truth"), self.root, self.det, self.ev, F, stagea_root)
        self.pctree_file = F.path("mlpr_v1/pr-pctree.tar.gz")
        # frames: tickinfo + bad channels
        fr = F.find("sp-frames.tar") or F.find("sim-frames-anode")
        self.tickinfo, self.bad, nfr, self.nticks = load_frames_info([F.path(r) for r in fr]) if fr else (None, {}, 0, None)
        self.nticks_source = "frames" if self.nticks else "default"
        if not self.nticks:
            self.nticks = NTICKS_DEFAULT.get(self.det, 6000)
        self.bad_source = "frames:chanmask_bad" if (nfr and self.bad) else ("frames:no_mask" if nfr else "none")
        if not nfr:
            sb = load_static_bad(F)
            if sb:
                self.bad, self.bad_source = sb, "static:chanmask-apa.json"
        if self.tickinfo is None:
            self.tickinfo = np.array([0.0, 500.0, 0.0])
        # deposits + truth
        self.D, self.I, self.P3, self.T0 = load_depos(self.depo_file)
        self.truth = TT.read_truth_npz(self.truth_file)
        self._finish()
        # pctree points
        self.pts = load_points(self.pctree_file)
        self.meta_inputs = {os.path.basename(self.depo_file): F.stream_sha(F.find("sim-depos.tar")[0]),
                            self.truth_id: sha_file(self.truth_file), "pr-pctree.tar.gz": F.stream_sha("mlpr_v1/pr-pctree.tar.gz")}

    def _finish(self):
        """deposit ids, charge, energy, geometry, faces and the reco-frame placement"""
        pf, c = self.truth["truth_pf"], TT.pf_col
        sed = self.truth["sed_ids"]
        ids = self.I[:, 0]; inr = (ids >= 0) & (ids < len(sed)); j = np.clip(ids, 0, max(len(sed) - 1, 0))
        if len(sed):
            self.tid = np.abs(np.where(inr, sed[j, 0], ids)).astype(np.int64)
            self.otid = np.abs(np.where(inr, sed[j, 1] if sed.shape[1] > 1 else sed[j, 0], ids)).astype(np.int64)
        else:
            self.tid = np.abs(ids).astype(np.int64); self.otid = self.tid.copy()
        se = self.truth.get("sed_energy", np.zeros(0))
        self.energy = np.where(inr & (ids < len(se)), se[np.clip(ids, 0, max(len(se) - 1, 0))], 0.0).astype(np.float64) if len(se) else np.zeros(len(ids))
        self.q = np.abs(self.D[:, 1])
        self.n = len(self.D)
        # geometry and the reco-frame placement (steiner_truth.build_trajectories' convention)
        self.geo = V2.load_wires(self.wires)
        self.faces = load_faces(self.wires)
        self.da, self.df = V2.assign_faces(self.D, self.geo)
        self.placed = self.da >= 0
        # inside the active volume: within its face's wire extent in y and z with no margin (assign_faces allows 50 mm)
        self.inside = np.zeros(self.n, bool)
        for an, faces in self.geo.items():
            for fi, fg in faces.items():
                g0 = fg[0]; m = (self.da == an) & (self.df == fi)
                self.inside[m] = (self.D[m, 3] >= g0["ymin"]) & (self.D[m, 3] <= g0["ymax"]) & (self.D[m, 4] >= g0["zmin"]) & (self.D[m, 4] <= g0["zmax"])
        self.Xreco = self._place()

    def _place(self):
        """reco-frame position of every deposit: x from the arrival time on its own face (as build_trajectories), y z drifted"""
        keys, dface = assign_depo_faces(self.D, self.faces)
        x = np.full(self.n, np.nan)
        toff = float(self.inputs["depo_toff_us"]); roff = float(self.inputs["reco_toff_us"]); speed = float(self.inputs["speed_mmus"])
        for ki, k in enumerate(keys):
            sel = dface == ki
            if not sel.any():
                continue
            xo = self.faces[k]["xorig"]
            dirx = float(np.sign(np.median(self.D[sel, 2] - xo))) or 1.0
            x[sel] = xo + dirx * (self.D[sel, 0] / US + toff + roff) * speed
        return np.column_stack([x, self.D[:, 3], self.D[:, 4]])

    def close(self):
        if self.files is not None:
            self.files.close()

    # ---- settings for the projection
    def projection_settings(self):
        tj = self.tj; ra = tj["extra"]["real_assoc"]
        return dict(toff_ns=float(tj["time_offset_us"]) * 1000.0, speed=float(tj["speed_mmus"]),
                    ext_tr=(float(ra["extra_smear_tran_induction"]), float(ra["extra_smear_tran_induction"]), float(ra["extra_smear_tran_collection"])),
                    ext_long_us=float(tj["extra"]["extra_smear_long"]) / 1000.0,
                    tick=float(self.tickinfo[1]), tbin=float(self.tickinfo[2]), nticks=int(self.nticks))
