#!/usr/bin/env python3
"""labels_v2.py -- the truth labeler v2 for one campaign event (pr_ml docs 50a / 50b secs 2.1-2.8; built in doc 54).

  labels_v2.py <event dir> <out dir> [--root R] [--wires-dir D] [--stagea-root S] [--slim F] [--det DET] [--tmpdir T]
      -> <out dir>/<event>.labels_v2.npz; prints "<event> ok|empty|fail <why>"
  labels_v2.py --selftest                                 (see selftest_v2.py)

Per point (pt_*, in the pctree / slim / term order): cls7 (0 ghost, 1 track, 2 stub, 3 EM segment, 4 EM blip, 5 Michel,
6 delta, -1 other), cls_obj (0 track, 1 stub, 2 delta, 3 michel, 4 shower, 5 LE blip cluster, 6 ghost, 7 other), inst
(compact instance id, -1 none), piece_id, piece_type (0 none, 1 segment, 2 blip), d_own_seg (mm; -1 none), overlap bits
(1 not free, 2 shower-ringed, 4 folded delta / Michel on its track, 8 absorbed EM object, 16 swallowed stub), stretch
(0 none, 1 traversed, 2 hidden end, 3 swallowed, 4 buried), kp_type / kp_dist / kp_index, ghost_dist, match_dist,
other_reason bits (1 ghost, 2 shared, 4 no truth, 8 far from own trajectory, 16 no deposit, 32 / 64 fallback match, 128 the
point's instance holds < 2,400 e),
and the v1 columns lab_tid / lab_orig_tid / lab_cls / lab_pf_node / lab_d_traj for the shift accounting.
Per instance (inst_*): root_row, g4_id, kind (0 track kind, 1 delta, 2 Michel, 3 EM object, 4 other), cls_obj, pdg,
parent_inst, parent_g4, edge, interaction, origin (index into inst_origin_names), flags (1 orphan, 2 buried, 4 hidden start, 8 hidden end, 16 vertex stub,
32 absorbed, 64 swallowed, 128 resolved, 256 hidden attach, 512 primary, 2048 visible, 4096 unresolved long track), n_dep, nseg, npiece, q, e_dep_mev,
ke_mev, vis_frac, ratio, chord_mm, run_mm, wait_ns, main_ext_mm (the delta / Michel main free piece), max_ext_mm (an EM instance's largest piece, >= 2-view extent).  Key points: kp_x/y/z, kp_kind, kp_inst.  meta: json.
"""
import ctypes, hashlib, json, os, sys, time
if __name__ == "__main__":
    try:
        ctypes.CDLL(None).prctl(41, 1, 0, 0, 0)     # PR_SET_THP_DISABLE: THP page faults are ~200x slow on wcgpu1 (pr_ml doc 50e)
    except Exception:
        pass
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from labels_v2 import VERSION, evio, rules, points, keypoints, writer   # noqa: E402


def peak_rss_mb():
    try:
        for l in open("/proc/self/status"):
            if l.startswith("VmHWM"):
                return int(l.split()[1]) // 1024
    except Exception:
        pass
    return -1


def label(d, outdir, root=None, wires_dir=None, stagea_root=None, slim=None, det=None, tmpdir=None):
    t0 = time.time()
    try:
        ev = evio.Event(d, root=root, wires_dir=wires_dir, stagea_root=stagea_root, tmpdir=tmpdir)
    except RuntimeError as x:            # the lane left no labels (glabels json "error": no live pctree, no steiner cloud): an empty event, as v1
        if str(x).startswith("glabels:"):
            return "empty", str(x)[:120]
        raise
    try:
        det = det or ev.det
        if not len(ev.pts["X"]):
            return "empty", "no points in the pctree"
        t1 = time.time()
        R = rules.label_event(ev, det)
        t2 = time.time()
        pt = points.link(ev, R)
        kps, kpt = keypoints.build(ev, R, R.table)
        t3 = time.time()
        if slim:
            z = np.load(slim)
            Xs = np.column_stack([z["pt_xraw"], z["pt_y"], z["pt_z"]]).astype(np.float64)
            if len(Xs) != len(ev.pts["X"]) or not np.allclose(Xs, ev.pts["X"], atol=1e-3):
                return "fail", f"point order differs from the slim file ({len(Xs)} vs {len(ev.pts['X'])})"
        counts = {k: int((pt["cls7"] == v).sum()) for k, v in rules.CLS7.items()}
        meta = dict(tool="truth/python/labels_v2", version=VERSION, det=det, event=ev.ev, params=R.prm, bad_source=ev.bad_source, n_bad_channels=len(ev.bad),
                    n_bad_main_entries=int(R.P.n_bad_main_entries), inputs=ev.meta_inputs, wires=os.path.basename(ev.wires), log=R.log,
                    n_points=int(len(pt["cls7"])), cls7_counts=counts, n_keypoints=int(len(kps["kp_kind"])), kp_points=int((kpt["kp_type"] > 0).sum()),
                    match_dist_p50=float(np.median(pt["match_dist"][pt["match_dist"] >= 0])) if (pt["match_dist"] >= 0).any() else -1.0,
                    match_dist_p99=float(np.percentile(pt["match_dist"][pt["match_dist"] >= 0], 99)) if (pt["match_dist"] >= 0).any() else -1.0,
                    wall_s=dict(load=round(t1 - t0, 2), rules=round(t2 - t1, 2), points_kp=round(t3 - t2, 2)), peak_rss_mb=peak_rss_mb(),
                    slim_sha256=hashlib.sha256(open(slim, "rb").read()).hexdigest() if slim else None)
        A = writer.assemble(ev, R, pt, kps, kpt, meta)
        os.makedirs(outdir, exist_ok=True)
        writer.write_npz(f"{outdir}/{ev.ev}.labels_v2.npz", A)
        return "ok", f"{meta['n_points']} points inst {R.log['n_instances']} kp {meta['n_keypoints']} ({meta['kp_points']} pts) cls7 {counts} wall {meta['wall_s']} rss {meta['peak_rss_mb']}MB"
    finally:
        ev.close()


def main():
    if "--selftest" in sys.argv:
        from labels_v2 import selftest_v2
        sys.exit(selftest_v2.main())
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    g = lambda k: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else None
    if len(args) < 2:
        print(__doc__); sys.exit(2)
    d, outdir = args[0], args[1]
    ev = os.path.basename(d.rstrip("/"))
    try:
        st, why = label(d, outdir, g("--root"), g("--wires-dir"), g("--stagea-root"), g("--slim"), g("--det"), g("--tmpdir"))
    except Exception as x:   # noqa: BLE001
        import traceback
        st, why = "fail", repr(x)[:200] + " | " + traceback.format_exc().splitlines()[-3].strip()[:200]
    print(ev, st, why, flush=True)
    sys.exit(0 if st == "ok" else 1 if st == "empty" else 2)


if __name__ == "__main__":
    main()
