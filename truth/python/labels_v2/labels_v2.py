#!/usr/bin/env python3
"""labels_v2.py -- the truth labeler v2 for one campaign event (pr_ml docs 50a / 50b secs 2.1-2.8; built in doc 54).

  labels_v2.py <event dir> <out dir> [--root R] [--wires-dir D] [--stagea-root S] [--slim F] [--det DET] [--tmpdir T] [--start-tie-order] [--kink-points] [--far-length] [--c6] [--c7]
      -> <out dir>/<event>.labels_v2.npz; prints "<event> ok|empty|fail <why>"
      --start-tie-order (default off = byte-identical; pr_ml doc 71): the instance table's start / end are taken after ordering
      each run of equal-time own deposits along the track (steiner_truth.tie_order).  The depo file stores float32 times, so a
      ms-scale (cosmic) particle's deposits tie over centimetres and the first in file order is not the first along the track.
      When set, meta params gain "start_tie_order": true.
      --kink-points (default off = byte-identical; pr_ml doc 71 sec 6): adds the members kink_x / kink_y / kink_z, kink_inst,
      kink_angle, kink_to_end -- the hard scatters of one track particle (>= 15 deg between the 30 mm chords behind and ahead on
      its true trajectory, both chords full, the deposit seen), where a segment boundary is allowed but not required.  The
      key-point table and every other member are unchanged; meta gains params "kink_points": true and n_kinks.
      --far-length (default off = byte-identical; pr_ml doc 76, Xin's Q6): a far free piece of an unresolved delta (Q16 d) gets the
      usual length test: when its extent over deposits seen in >= 2 views is >= 2 steps it is a segment -- the electron's own deposits
      class delta, its radiation class EM segment, piece type 1, its instance (still origin delta_far_piece) cls_obj delta or shower by
      charge majority; shorter far pieces stay LE blips as before.  meta params gain "far_length": true, log n_delta_far_long.
      --c6 (default off = byte-identical; pr_ml doc 82, Xin 2026-11-16: the adopted E-1 target C6): implies --far-length and changes
      the per-point class pt_cls7 only (instances, cls_obj, pieces, overlap bits, key points as without it): (a) radiation (own
      particle another e / gamma) of a Michel / delta class farther than 15 mm from the same instance's Michel / delta electron ->
      EM segment or blip by piece type; (b) an EM segment deposit of a delta electron -> delta; absorbed / swallowed EM blips -> EM
      segment.  meta params gain "c6": true and "far_length": true; log n_c6_detached, n_c6_trunk, n_c6_absorbed.
      --c7 (default off = byte-identical; pr_ml doc 86, Xin 2026-11-17: label what is visible): implies --c6.  pt_cls7 becomes the
      C7 class (hidden stubs, unresolved and buried tracks, long hidden track ends and short track pieces take the class of the
      charge that covers them, or EM blip); new members pt_c7_piece (the visible piece, -1 none), pt_c7_flag (1 folded stub,
      2 folded track, 4 hidden track now shower, 8 now EM blip, 16 zone, 32 hidden stretch of a through-going track),
      pt_c7_place_dist (mm to the nearest place, -1 none), the piece table pc_* (type in cls_obj numbering, inst, link = G4 id,
      graph, len, e0, e1, place0, place1), the place table pl_* (xyz, kinds bitmask, tag bitmask, member pairs) and the implied
      graph seg_* (see c7.py).  Every other member as with --c6; meta params gain "c7": true, the log n_c7_* counts.
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
ke_mev, vis_frac, ratio, chord_mm, run_mm, wait_ns, main_ext_mm (the delta / Michel main free piece), max_ext_mm (an EM instance's largest piece, >= 2-view extent), start_x/y/z and end_x/y/z (first / last own deposit, reco frame, mm),
g4_start_x/y/z and g4_end_x/y/z (the G4 start / end, true frame, mm; for a vertex fit).  Key points: kp_x/y/z, kp_kind, kp_inst.  meta: json
(byte-reproducible); run: json with the timings and peak RSS (the only member that differs between two runs).
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


def label(d, outdir, root=None, wires_dir=None, stagea_root=None, slim=None, det=None, tmpdir=None, params=None):
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
        R = rules.label_event(ev, det, params)
        t2 = time.time()
        pt = points.link(ev, R)
        kps, kpt = keypoints.build(ev, R, R.table)
        c7d = None
        if R.prm.get("c7"):               # --c7 (pr_ml doc 86): visible pieces and key-point places; pt["cls7"] becomes the C7 class
            from labels_v2 import c7
            c7d, c7log = c7.build(ev, R, pt, kps, det); R.log.update(c7log)
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
                    slim_sha256=hashlib.sha256(open(slim, "rb").read()).hexdigest() if slim else None)
        # timings and memory go to a separate member, so that `meta` and every array are byte-reproducible run to run (the SDCC
        # gate compares all members but `run`; two runs of one event differed only in wall_s, campaign session 2026-11-05)
        run = dict(wall_s=dict(load=round(t1 - t0, 2), rules=round(t2 - t1, 2), points_kp=round(t3 - t2, 2)), peak_rss_mb=peak_rss_mb())
        kk = keypoints.kinks(ev, R, R.table) if R.prm.get("kink_points") else None   # before meta is frozen below
        if kk is not None:
            meta["n_kinks"] = int(len(kk["kink_inst"]))
        A = writer.assemble(ev, R, pt, kps, kpt, meta)
        if kk is not None:
            A.update(kk)
        if c7d is not None:
            A.update({k: (writer.narrow(v) if np.asarray(v).dtype.kind in "iu" else np.asarray(v)) for k, v in c7d.items()})
        A["run"] = np.frombuffer(json.dumps(run, sort_keys=True).encode(), np.uint8)
        os.makedirs(outdir, exist_ok=True)
        writer.write_npz(f"{outdir}/{ev.ev}.labels_v2.npz", A)
        return "ok", f"{meta['n_points']} points inst {R.log['n_instances']} kp {meta['n_keypoints']} ({meta['kp_points']} pts) cls7 {counts} wall {run['wall_s']} rss {run['peak_rss_mb']}MB"
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
        st, why = label(d, outdir, g("--root"), g("--wires-dir"), g("--stagea-root"), g("--slim"), g("--det"), g("--tmpdir"),
                        {k: True for f, k in (("--start-tie-order", "start_tie_order"), ("--kink-points", "kink_points"), ("--far-length", "far_length"), ("--c6", "c6"), ("--c7", "c7")) if f in sys.argv} or None)
    except Exception as x:   # noqa: BLE001
        import traceback
        st, why = "fail", repr(x)[:200] + " | " + traceback.format_exc().splitlines()[-3].strip()[:200]
    print(ev, st, why, flush=True)
    sys.exit(0 if st == "ok" else 1 if st == "empty" else 2)


if __name__ == "__main__":
    main()
