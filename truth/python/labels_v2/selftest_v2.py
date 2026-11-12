#!/usr/bin/env python3
"""selftest_v2.py -- the v2 rules on synthetic events drawn on the real SBND wires (pr_ml doc 54 step 5).
  (a) a muon with an unresolved 6 mm delta (folded into the track) and a resolved 80 mm delta at 60 deg (class delta, own
      instance, attachment key point); the muon is a Track with two physical end key points
  (b) two parallel muons 4 mm apart: never resolved -> Tracks with the unresolved flag (Q17), every deposit overlapping, no key
      points; a third far muon is a Track with two end key points
  (c) a photon's 100 mm conversion electron (segment) and a 3 mm Compton blob 60 mm beyond it (blip, same shower instance,
      d_own_seg ~ 60 mm); an 8 mm proton on the segment is swallowed (shower instance, stretch 3); an 8 mm proton at the
      shower start keeps class stub (vertex guard)
  (d) a stopped mu+ with a 40 mm Michel (class Michel, key point) and a mu- decaying in flight (its electron: an EM object,
      a shower with one segment; the mu- end under its own daughter shower stays a physical end)
  (e) a neutron's photon with a 4 mm Compton blob: a G3 LE blip cluster (class blip, own instance, origin neutron, no kp)
  (f) bad channels: W channels under part of a muon masked -> not seen in W there, still seen in U / V, still free
  (g) every deposit past one plane's wires (pr_ml doc 58 sec 11.1): no entry in that plane (None), seen in the other two, labelled
Run: python3 selftest_v2.py [wires file]; exit 0 when every check passes, prints FAIL lines otherwise."""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import truth_table as TT                                   # noqa: E402
from labels_v2 import evio, rules, points, keypoints        # noqa: E402
from labels_v2.rules import OBJ, CLS7                       # noqa: E402
from labels_v2.keypoints import KP                          # noqa: E402

WIRES = "/home/xqian/toolkit-dev/wire-cell-data/sbnd-wires-geometry-v0206.json.bz2"
SPEED = 1.563; XORIG = -2020.5; XPLANE = -2014.5
INPUTS = dict(depo_toff_us=0.0, reco_toff_us=0.0, speed_mmus=SPEED, tie_order=True)
TJ = dict(time_offset_us=0.0, speed_mmus=SPEED, extra=dict(real_assoc=dict(extra_smear_tran_induction=0.6, extra_smear_tran_collection=0.6), extra_smear_long=1000.0))
C = 299.792458
FAILS = []


def check(cond, msg):
    if not cond:
        FAILS.append(msg); print("FAIL", msg)
    else:
        print("ok  ", msg)


class Builder:
    def __init__(self):
        self.rows = []; self.dep = []          # dep: (tid, P3, T0, q, e)

    def particle(self, tid, pdg, mother, process, start, end, ke_mev, t_start=0.0, q_per_step=3000.0, step=0.5, deposits=True, blob=None):
        start = np.asarray(start, float); end = np.asarray(end, float)
        L = np.linalg.norm(end - start)
        self.rows.append(dict(tid=tid, pdg=pdg, mother=mother, process=process, start=start, end=end, ke=ke_mev, t=t_start))
        if deposits:
            if blob is not None:                                   # a small blob of n deposits within radius r around start
                n, r = blob
                rng = np.random.default_rng(tid)
                P = start + rng.normal(0, r / 2, (n, 3))
                self.dep.append((tid, P, np.full(n, t_start), np.full(n, q_per_step), np.full(n, q_per_step / 2.4e4)))
            else:
                n = max(int(L / step), 2)
                s = np.linspace(0, L, n)
                P = start + np.outer(s / max(L, 1e-9), end - start)
                self.dep.append((tid, P, t_start + s / C, np.full(n, q_per_step), np.full(n, q_per_step / 2.4e4)))

    def event(self, bad=None, every=3):
        pf = np.zeros((len(self.rows), len(TT.PF_COLUMNS))); c = TT.pf_col
        for i, r in enumerate(self.rows):
            pf[i, c("nu_row")] = -1; pf[i, c("trackid")] = r["tid"]; pf[i, c("parent_trackid")] = r["mother"]; pf[i, c("mother_trackid")] = r["mother"]
            pf[i, c("pdg")] = r["pdg"]; pf[i, c("process")] = TT.process_code(r["process"]); pf[i, c("KE")] = r["ke"] / 1e3; pf[i, c("E")] = r["ke"] / 1e3
            pf[i, [c("start_x"), c("start_y"), c("start_z")]] = r["start"] / 10.0; pf[i, [c("end_x"), c("end_y"), c("end_z")]] = r["end"] / 10.0
            pf[i, c("start_t")] = r["t"]; pf[i, c("end_t")] = r["t"] + np.linalg.norm(r["end"] - r["start"]) / C
            pf[i, c("orig_trackid")] = r["tid"]; pf[i, c("pf_node")] = -1; pf[i, c("cls")] = -1; pf[i, c("merged_into")] = -1
        tids = np.concatenate([np.full(len(P), t) for t, P, _, _, _ in self.dep]); P3 = np.concatenate([P for _, P, _, _, _ in self.dep])
        T0 = np.concatenate([t for _, _, t, _, _ in self.dep]); q = np.concatenate([x for _, _, _, x, _ in self.dep]); e = np.concatenate([x for _, _, _, _, x in self.dep])
        n = len(tids)
        t_drift = (P3[:, 0] - XORIG) / SPEED * 1000.0                   # ns: reco x = xorig + t * speed recovers the true x
        D = np.column_stack([t_drift, -q, np.full(n, XPLANE), P3[:, 1], P3[:, 2], np.full(n, 1.0), np.full(n, 1.0)])
        I = np.column_stack([np.arange(n), np.zeros(n), np.zeros(n), np.full(n, -1)]).astype(np.int64)
        truth = dict(truth_pf=pf, sed_ids=np.column_stack([tids, tids]).astype(np.int32), sed_energy=e.astype(np.float32), gen_labels=["generator"], truth_nu=np.zeros((0, len(TT.NU_COLUMNS))), meta={})
        sel = np.arange(0, n, every)
        pts = dict(X=P3[sel].copy(), tid=tids[sel], orig_tid=tids[sel])
        ev = evio.Event.synthetic("sbnd", WIRES, D, I, P3, T0, truth, pts, INPUTS, TJ, bad=bad)
        return ev


def run(ev, params=None):
    R = rules.label_event(ev, "sbnd", params)
    pt = points.link(ev, R)
    kps, kpt = keypoints.build(ev, R, R.table)
    return R, pt, kps, kpt


def inst_of(R, tid):
    """compact instance id of the instance rooted at G4 id tid (-1 none)"""
    T = R.table; m = np.flatnonzero(T["inst_g4_id"] == tid)
    return int(m[0]) if len(m) else -1


def dep_of(ev, tid):
    return np.flatnonzero(ev.otid == tid)


def maj(a):
    v, c = np.unique(a, return_counts=True); return int(v[np.argmax(c)])


def test_a():
    print("--- (a) muon with an unresolved and a resolved delta")
    B = Builder()
    B.particle(1, 13, 0, "primary", (-1600, 0, 1000), (-1400, 0, 1400), 2000.0)                 # tilted in x: not isochronous
    B.particle(2, 11, 1, "muIoni", (-1550, 0, 1100), (-1547, 0, 1106), 1.0, t_start=100 / C, q_per_step=2500.0)
    d = np.array([0.3, np.sin(np.radians(60)), np.cos(np.radians(60))]); d /= np.linalg.norm(d)
    B.particle(3, 11, 1, "muIoni", (-1500, 0, 1200), np.array([-1500, 0, 1200]) + 80 * d, 20.0, t_start=200 / C, q_per_step=2500.0)
    ev = B.event(); R, pt, kps, kpt = run(ev)
    i1, i2, i3 = inst_of(R, 1), inst_of(R, 2), inst_of(R, 3)
    check(i1 >= 0 and R.table["inst_cls_obj"][i1] == OBJ["track"], "muon is a Track")
    check(i2 < 0, "the 6 mm delta has no instance (folded into the muon)")
    d2 = dep_of(ev, 2); check(maj(R.inst[d2]) == i1 and maj(R.cls7[d2]) == CLS7["track"] and (R.overlap[d2] & rules.OV_FOLDED).all(), "its deposits: muon instance, class track, folded bit")
    check(i3 >= 0 and R.table["inst_cls_obj"][i3] == OBJ["delta"] and (R.table["inst_flags"][i3] & rules.FL_RESOLVED), "the 80 mm delta is resolved, class delta, own instance")
    d3 = dep_of(ev, 3); check(maj(R.cls7[d3]) == CLS7["delta"] and maj(R.inst[d3]) == i3, "its deposits carry class delta and its instance")
    check(R.table["inst_parent_inst"][i3] == i1 and R.table["inst_edge"][i3] == rules.EDGE["delta"], "delta -> parent edge to the muon")
    kk = kps["kp_kind"].tolist()
    check(kk.count(KP["track_end"]) == 2 and kk.count(KP["delta_attach"]) == 1, f"key points: 2 track ends + 1 delta attachment (got {kk})")
    att = np.column_stack([kps["kp_x"], kps["kp_y"], kps["kp_z"]])[np.array(kk) == KP["delta_attach"]][0]
    check(np.linalg.norm(att - np.array([-1500, 0, 1200])) < 1.5, "the attachment point is the delta's start")
    # points follow their deposits
    p3 = ev.pts["orig_tid"] == 3; check((pt["cls7"][p3] == CLS7["delta"]).mean() > 0.9 and (kpt["kp_type"][p3] > 0).any(), "points of the delta: class delta, some within 10 mm of the attachment")
    check(pt["match_dist"].max() < 1.0, "every synthetic point matched to its own deposit")
    # the muon's not-free deposits are those under the deltas' rings only
    d1 = dep_of(ev, 1); nf = ~R.free[d1]
    check(0 < nf.mean() < 0.3, f"muon deposits not free: {nf.mean():.2f} (near the deltas only)")


def test_b():
    print("--- (b) two parallel muons 4 mm apart, one far muon")
    B = Builder()
    B.particle(1, 13, 0, "primary", (-1600, 0, 1000), (-1400, 0, 1400), 2000.0)
    B.particle(2, 13, 0, "primary", (-1600, 4, 1000), (-1400, 4, 1400), 2000.0)
    B.particle(3, 13, 0, "primary", (-1600, 300, 1000), (-1400, 300, 1400), 2000.0)
    ev = B.event(); R, pt, kps, kpt = run(ev)
    for t in (1, 2):
        i = inst_of(R, t); d = dep_of(ev, t)
        check(R.table["inst_cls_obj"][i] == OBJ["track"] and (R.table["inst_flags"][i] & rules.FL_UNRESOLVED), f"muon {t}: Track with the unresolved flag (Q17: long, never resolved)")
        check((~R.free[d]).mean() > 0.95, f"muon {t}: {100 * (~R.free[d]).mean():.0f} % of deposits overlapping")
        check(maj(R.cls7[d]) == CLS7["track"], f"muon {t}: deposits class track")
    i3 = inst_of(R, 3); check(R.table["inst_cls_obj"][i3] == OBJ["track"] and R.free[dep_of(ev, 3)].all(), "far muon: Track, every deposit free")
    kk = kps["kp_kind"].tolist(); check(kk.count(KP["track_end"]) == 2, f"key points: the far muon's 2 ends only among track ends (got {kk})")


def test_c():
    print("--- (c) a shower: segment + far blip, a swallowed stub and a vertex stub")
    B = Builder()
    B.particle(1, 22, 0, "primary", (-1500, 500, 1900), (-1500, 500, 2000), 300.0, deposits=False)
    B.particle(2, 11, 1, "conv", (-1500, 500, 2000), (-1450, 500, 2100), 250.0, t_start=100 / C, q_per_step=3000.0)
    B.particle(3, 11, 1, "compt", (-1420, 500, 2160), (-1420, 500, 2163), 2.0, t_start=260 / C, q_per_step=2000.0, blob=(10, 1.5))
    B.particle(4, 2212, 0, "primary", (-1475, 500, 2050), (-1471, 500, 2058), 30.0, q_per_step=6000.0)
    B.particle(5, 2212, 0, "primary", (-1501, 500, 2001), (-1497, 500, 2009), 30.0, q_per_step=6000.0)
    ev = B.event(); R, pt, kps, kpt = run(ev)
    i1 = inst_of(R, 1)
    check(i1 >= 0 and R.table["inst_cls_obj"][i1] == OBJ["shower"] and R.table["inst_nseg"][i1] == 1, "the photon's object is a shower with one segment")
    d2, d3 = dep_of(ev, 2), dep_of(ev, 3)
    check(maj(R.piece_type[d2]) == 1 and maj(R.cls7[d2]) == CLS7["segment"] and maj(R.inst[d2]) == i1, "conversion electron: segment, class EM segment, shower instance")
    check(maj(R.piece_type[d3]) == 2 and maj(R.cls7[d3]) == CLS7["blip"] and maj(R.inst[d3]) == i1, "Compton blob: blip, class EM blip, the same shower instance")
    check(50 < np.median(R.d_own_seg[d3]) < 70, f"blip distance to its own segment {np.median(R.d_own_seg[d3]):.0f} mm")
    d4 = dep_of(ev, 4); i4 = inst_of(R, 4)
    check(i4 < 0 and maj(R.inst[d4]) == i1 and maj(R.stretch[d4]) == rules.STR_SWALLOWED and maj(R.cls_obj[d4]) == OBJ["shower"], "proton on the segment: swallowed into the shower (stretch 3)")
    d5 = dep_of(ev, 5); i5 = inst_of(R, 5)
    check(i5 >= 0 and R.table["inst_cls_obj"][i5] == OBJ["stub"] and (R.table["inst_flags"][i5] & rules.FL_VERTEX_STUB), "proton at the shower start: stays a stub, vertex-stub flag")
    kk = kps["kp_kind"].tolist(); check(kk.count(KP["shower_start"]) == 1, f"one shower-start key point (got {kk})")
    ss = np.column_stack([kps["kp_x"], kps["kp_y"], kps["kp_z"]])[np.array(kk) == KP["shower_start"]][0]
    check(np.linalg.norm(ss - np.array([-1500, 500, 2000])) < 2.5, "the shower start is the conversion electron's first deposit")


def test_d():
    print("--- (d) a Michel and a decay in flight")
    B = Builder()
    B.particle(1, -13, 0, "primary", (-1600, -500, 3000), (-1500, -500, 3300), 100.0)
    B.particle(2, -11, 1, "Decay", (-1500, -500, 3300), (-1480, -460, 3300), 40.0, t_start=300 / C + 2000.0, q_per_step=2500.0)
    B.particle(3, 13, 0, "primary", (-1300, -800, 3600), (-1200, -800, 3900), 2000.0)       # apart in z and drift time from the mu+
    B.particle(4, 11, 3, "Decay", (-1200, -800, 3900), (-1180, -760, 3900), 500.0, t_start=300 / C, q_per_step=2500.0)
    ev = B.event(); R, pt, kps, kpt = run(ev)
    i2 = inst_of(R, 2); d2 = dep_of(ev, 2)
    check(i2 >= 0 and R.table["inst_cls_obj"][i2] == OBJ["michel"] and maj(R.cls7[d2]) == CLS7["michel"], "mu+ decay electron after 2 us: Michel")
    check(abs(R.table["inst_wait_ns"][i2] - 2000.0) < 1.0, f"wait {R.table['inst_wait_ns'][i2]:.0f} ns")
    i4 = inst_of(R, 4); d4 = dep_of(ev, 4)
    check(i4 >= 0 and R.table["inst_kind"][i4] == 3 and R.table["inst_cls_obj"][i4] == OBJ["shower"] and maj(R.cls7[d4]) == CLS7["segment"], "mu- decay in flight: an EM object, shower with a segment")
    kk = kps["kp_kind"].tolist(); check(kk.count(KP["michel_start"]) == 1 and kk.count(KP["shower_start"]) == 1 and kk.count(KP["track_end"]) == 4, f"key points: Michel start, shower start, 4 track ends (got {kk})")
    i3 = inst_of(R, 3); check(not (R.table["inst_flags"][i3] & rules.FL_HIDDEN_END), "the mu- end under its own decay-in-flight shower is the decay vertex, not a hidden end (doc 54 sec 4)")
    i1 = inst_of(R, 1); check(R.table["inst_cls_obj"][i1] == OBJ["track"] and R.table["inst_parent_inst"][i2] == i1 and R.table["inst_edge"][i2] == rules.EDGE["michel"], "Michel -> parent edge to its muon")


def test_e():
    print("--- (e) a neutron's photon with a Compton blob")
    B = Builder()
    B.particle(1, 2112, 0, "primary", (-1500, 800, 3900), (-1500, 800, 4000), 10.0, deposits=False)
    B.particle(2, 22, 1, "nCapture", (-1500, 800, 4000), (-1500, 800, 4000), 6.0, deposits=False)
    B.particle(3, 11, 2, "compt", (-1500, 800, 4000), (-1500, 800, 4004), 1.0, q_per_step=2500.0, blob=(10, 2.0))
    B.particle(4, 13, 0, "primary", (-1600, 1500, 1000), (-1400, 1500, 1400), 2000.0)        # a control track
    ev = B.event(); R, pt, kps, kpt = run(ev)
    i2 = inst_of(R, 2); d3 = dep_of(ev, 3)
    check(i2 >= 0 and R.table["inst_kind"][i2] == 3 and R.table["inst_cls_obj"][i2] == OBJ["le"], "the photon roots a G3 LE blip cluster")
    check(maj(R.cls7[d3]) == CLS7["blip"] and maj(R.inst[d3]) == i2 and maj(R.piece_type[d3]) == 2, "its deposits: class blip, its instance")
    check(R.table["origins"][R.table["inst_origin"][i2]] == "neutron", "origin neutron")
    check(R.table["inst_nseg"][i2] == 0 and KP["shower_start"] not in kps["kp_kind"].tolist(), "no segment, no shower-start key point")


def test_f():
    print("--- (f) bad channels")
    B = Builder()
    B.particle(1, 13, 0, "primary", (-1600, 1200, 1000), (-1400, 1200, 1400), 2000.0)
    ev0 = B.event(); g = ev0.geo[0][0][2]
    zs = np.arange(1100.0, 1200.0, 1.0); s = (zs - g["cen0"][2]) * g["pdir"][2] / g["pitch"]
    chans = np.unique(g["chan"][np.clip(np.rint(s).astype(int), 0, len(g["chan"]) - 1)])
    bad = {int(c): [(0, 100000)] for c in chans}
    ev = B.event(bad=bad); R, pt, kps, kpt = run(ev)
    d = dep_of(ev, 1); inside = (ev.P3[d, 2] >= 1102) & (ev.P3[d, 2] <= 1198); outside = (ev.P3[d, 2] < 1095) | (ev.P3[d, 2] > 1205)
    check((~R.P.seen[d[inside], 2]).all() and R.P.seen[d[outside], 2].all(), "deposits under the masked W channels are not seen in W, the others are")
    check(R.P.seen[d[inside], 0].all() and R.P.seen[d[inside], 1].all(), "they are still seen in U and V")
    check(R.free[d].all() and R.table["inst_cls_obj"][inst_of(R, 1)] == OBJ["track"], "still free, still a Track")
    check(R.P.n_bad_main_entries == int(inside.sum()) or abs(R.P.n_bad_main_entries - inside.sum()) <= 20, f"bad main entries counted: {R.P.n_bad_main_entries} vs {int(inside.sum())} deposits inside")


def test_g():
    print("--- (g) an event that one plane does not see (every deposit past that plane's wires)")
    # pr_ml doc 58 sec 11.1 (6e4, fdvd nc_b32052_e46 / nc_b32025_e1): every deposit past the last collection strip of its CRP, inside
    # both induction planes; the face's entries were appended with zero rows and Projection._totals failed on ch.min().  Here: a blob
    # just before the first SBND W wire (z < 0), inside U and V.
    B = Builder()
    B.particle(1, 11, 0, "primary", (-1500, -1000, -25), (-1500, -1000, -25), 1.0, blob=(20, 3.0))
    ev = B.event(); g = ev.geo[0][0]
    s = {pi: np.einsum("ij,j->i", ev.P3 - g[pi]["cen0"], g[pi]["pdir"]) / g[pi]["pitch"] for pi in range(3)}
    check(ev.placed.all() and (s[2] < -1).all() and all(((s[p] > 1) & (s[p] < len(g[p]["chan"]) - 2)).all() for p in (0, 1)),
          "setup: every deposit placed, past the W wires, inside U and V")
    try:
        R, pt, kps, kpt = run(ev)
    except Exception as x:   # noqa: BLE001
        check(False, f"labels an event that W does not see ({x!r:.80})"); return
    check(R.P.ENT[2] is None and R.P.ENT[0] is not None and R.P.ENT[1] is not None, "no W entry (None), U and V entries")
    check(not ev.inside.any() and not R.P.seen.any(), "outside the face's active box (W's extent): seen in no view, as the rule says")
    check(len(pt["cls7"]) == len(ev.pts["X"]), "every point labelled")


def test_h():
    print("--- (h) start / end of a particle whose deposits tie in time (a ms-scale cosmic, float32 times)")
    # pr_ml doc 71: PDHD / PDVD cosmic tracks, about 50 deposits over 16 mm share one stored time; the stored start was the first of the
    # run in file order, more than 10 mm past the true first point in a third of the tracks.  Here: a 300 mm muon whose times are
    # quantised to runs of 60 mm, the deposits stored last-first so that file order is the wrong way round inside every run.
    B = Builder()
    a, b = np.array([-1500.0, -1000.0, 100.0]), np.array([-1500.0, -1000.0, 400.0])
    B.particle(1, 13, 0, "primary", a, b, 1000.0)
    tid, P, T, q, e = B.dep[-1]; s = np.linalg.norm(P - a, axis=1); T = T[0] + np.floor(s / 60.0) * 60.0 / C
    B.dep[-1] = (tid, P[::-1].copy(), T[::-1].copy(), q, e)
    ev = B.event(); own = dep_of(ev, 1); first = own[np.argmin(np.linalg.norm(ev.P3[own] - a, axis=1))]; last = own[np.argmin(np.linalg.norm(ev.P3[own] - b, axis=1))]
    out = {}
    for nm, prm in (("off", None), ("on", dict(start_tie_order=True))):
        R = run(ev, prm)[0]; i = inst_of(R, 1); t = R.table
        out[nm] = (np.array([t[f"inst_start_{c}"][i] for c in "xyz"]), np.array([t[f"inst_end_{c}"][i] for c in "xyz"]), R.prm)
    check(np.linalg.norm(out["off"][0] - ev.Xreco[first]) > 30 and "start_tie_order" not in out["off"][2], "off: the stored start is the file-order first of the tied run (over 30 mm past the true first), no new key in params")
    check(np.linalg.norm(out["on"][0] - ev.Xreco[first]) < 1.0, "on: the stored start is the particle's first deposit")
    check(np.linalg.norm(out["on"][1] - ev.Xreco[last]) < 1.0, "on: the stored end is the particle's last deposit")


def test_i():
    print("--- (i) hard scatter of one particle (kink points, pr_ml doc 71 sec 6)")
    # one muon: 200 mm along z, then 200 mm at 30 deg (a kink in the middle, both chords full); a second muon bends 30 deg only 10 mm
    # before its end (no point with both 30 mm chords full sees 15 deg: no kink); a third is straight.  Deposits stored in time order.
    B = Builder(); th = np.radians(30.0); u = np.array([0.0, np.sin(th), np.cos(th)])

    def bent(tid, a, l1, l2):
        m = a + np.array([0.0, 0.0, l1]); b = m + l2 * u
        B.particle(tid, 13, 0, "primary", a, b, 1000.0)
        _, P, T, q, e = B.dep[-1]; n = len(P); s = np.linspace(0, l1 + l2, n)
        P = np.where((s <= l1)[:, None], a + np.outer(s, [0.0, 0.0, 1.0]), m + np.outer(s - l1, u)); B.dep[-1] = (tid, P, T[0] + s / C, q, e)
        return m
    m1 = bent(1, np.array([-1500.0, -1000.0, 100.0]), 200.0, 200.0)
    bent(2, np.array([-1300.0, -500.0, 100.0]), 200.0, 10.0)
    B.particle(3, 13, 0, "primary", np.array([-1100.0, 0.0, 100.0]), np.array([-1100.0, 0.0, 500.0]), 1000.0)
    ev = B.event(); R, pt, kps, kpt = run(ev)
    K = keypoints.kinks(ev, R, R.table); i1 = inst_of(R, 1)
    check(len(K["kink_inst"]) == 1 and int(K["kink_inst"][0]) == i1, "one kink, on the bent muon (none on the straight one, none where the bend is 10 mm from the end)")
    if len(K["kink_inst"]):
        p = np.array([K["kink_x"][0], K["kink_y"][0], K["kink_z"][0]])
        check(np.linalg.norm(p - m1) < 3.0 and abs(float(K["kink_angle"][0]) - 30.0) < 3.0, "the kink is at the bend (within 3 mm), angle 30 deg")
        check(abs(float(K["kink_to_end"][0]) - 200.0) < 3.0, "distance to the trajectory's end 200 mm")
    check(not any(k.startswith("kink") for k in kps) and "kink_points" not in R.prm, "off: the key-point table and params carry nothing new")
    # the deposit at the bend is not seen (FD-VD nc_b11009_e27, a 52 deg pion kink whose own deposit sat in an unseen pixel while 91 % of the
    # track was seen): the kink goes to the nearest seen own deposit within 10 mm along the trajectory; none within 10 mm: no kink
    own = dep_of(ev, 1); dm = np.linalg.norm(ev.P3[own] - m1, axis=1); keep = R.P.seen.copy()
    R.P.seen[own[dm < 4.0]] = False; K = keypoints.kinks(ev, R, R.table)
    check(len(K["kink_inst"]) == 1 and 3.0 < np.linalg.norm(np.array([K["kink_x"][0], K["kink_y"][0], K["kink_z"][0]]) - m1) < 10.0, "bend deposit unseen over 4 mm: the kink moves to the nearest seen own deposit")
    R.P.seen[own[dm < 15.0]] = False; K = keypoints.kinks(ev, R, R.table)
    check(len(K["kink_inst"]) == 0, "nothing seen within 10 mm of the bend: no kink")
    R.P.seen[:] = keep


def main():
    global WIRES
    if len(sys.argv) > 1 and not sys.argv[1].startswith("--"):
        WIRES = sys.argv[1]
    for t in (test_a, test_b, test_c, test_d, test_e, test_f, test_g, test_h, test_i):
        t()
    print("selftest_v2:", "PASS" if not FAILS else f"{len(FAILS)} FAIL")
    return 0 if not FAILS else 1


if __name__ == "__main__":
    sys.exit(main())
