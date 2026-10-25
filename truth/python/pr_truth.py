#!/usr/bin/env python3
"""pr_ml doc 23 sec 11 (WP0.2 "R" labels): MC-truth labels of the pattern-recognition chain's reconstructed objects
(doc 15 sec 2.4) and, sec 13, truth beside every fit point.

Inputs of one event:
  --pctree   the PR job's point-tree dump with the clouds of the clus visitor PRTruthExport (pr_truth stage):
             pr_segments / pr_seg_points = each PR segment's associated cluster points (p_index), plus the per-point
             "truth" cloud of every blob (TruthPointLabeler, doc 19 T6) and the blob charge (scalar cloud)
  --calib    the same job's PrDisplayDump file (pr_display stage, calib-pr-evt<N>.json): segments with their fit points,
             showers, vertices and the neutrino vertex, with the same ids (seg id = cluster_id*1000 + segment id)
  --depos / --truth / --wires / offsets   as steiner_truth.py (true trajectories, true neutrino vertex)

Reconstructed objects, three levels:
  segment       one PR segment (its associate_points)
  particle      a shower (the union of the segments carrying its shower_id) or a segment in no shower: the objects
                of the Wire-Cell particle-flow tree
  interaction   every segment of a neutrino candidate
Truth units: the PF node of the per-point labels (pf_node: the doc 15 sec 2.3 merging of TruthPointLabeler) for
segments and particles, the interaction row (iid) for interactions; the track key (tid) is reported beside them.
Ghost (tid -1), shared (-2) and unlabelled (-3) points count in an object's size but are never a match.

Per object (count and charge, charge = blob charge / blob points): the best unit, purity = shared / object,
completeness = shared / unit (whole event; the "bundle" variant counts the unit's points in the clusters holding the
candidate's segments only), IoU, and purity_lab = shared / labelled points of the object (ghosts left out; the "clean"
selection of the direction, track / shower, PID and fit-point summaries is purity_lab >= 0.8, since a third of an
FD-HD object can be ghost points).  Matching (doc 15 sec 2.4): SPINE purity >= 0.5 and completeness >= 0.5; Pandora
purity >= 0.5, completeness >= 0.1, >= 5 shared points; PQ with one-to-one IoU >= 0.5.  Instances with fewer than
--min-points (10, as SPICE) points are left out of efficiencies and PQ (stated in the output).

Derived labels: segment true pdg (of its best G4 key; the PF node's dominant pdg is kept as pf_pdg, since the merging
folds nuclear fragments and scattered protons into a node) and class (the matched PF node's majority class),
dirsign_true (the sign of the true direction of travel along the fit-point order, from the per-point dir of the
matched key's points; doc 23 sec 11 compares it with the reconstructed dirsign), per-vertex is_nu_vertex /
is_vertex (a true interaction vertex / particle start within 1 cm), and the distance of the reconstructed neutrino
vertex to the vertex of the interaction the candidate matched.  Interactions are the truth_nu rows: "nu" (a neutrino),
"gun" (a single-particle generator row), "corsika" (cosmic rays: SBND's CORSIKA, the ProtoDUNEs' cosmicgenerator; a
candidate matched to it is reported, never scored as a vertex) or "radiological" (doc 24: the ProtoDUNE ar39 / ar42 /
kr85 ... rows, never scored); a vertex is the median start of the row's visible primaries in the
reconstructed frame (steiner_truth's nu_vertices), else the generator vertex.

Doc 24 (the labels doc 23 sec 11 deferred; new output keys only, the doc 23 keys are unchanged):
  fit points   boundary = a true particle boundary within 1 cm (doc 15 sec 2.4): d_boundary = distance to the nearest
               start, end or kink (a direction change > 15 deg between the chords 3 cm behind and ahead, local maximum
               within +-3 cm: the C++ DepoPolyline::keypoints rule) of any trajectory unit >= 1 cm long; d_boundary_own =
               the same restricted to the matched key's own units
  PF parent    parent of every PF node: the node of the first G4 ancestor (mother chain of the node's root key) whose
               points belong to another node, -1 at the interaction (the merged truth_pf is not in the dumps); per
               segment true_parent_node
  edges        per candidate, every pair of segments that share a PR vertex or are both >= --min-points: adjacent,
               reco_same_shower (same shower_id), same_pf, same_shower_true (same PF node of class shower),
               same_interaction (the two best nodes' iid), parent_ab / parent_ba (node of a is the parent of node of b /
               the reverse)

Fit points (sec 13): every fit point of every segment with x, y, z, dQ, dx, residual range, reduced chi2 as the PR
chain wrote them, plus the distance to the matched key's trajectory (d_true, mm, every particle of the key searched)
and the true dE/dx (MeV/cm) of the nearest associated point of that key.  dQ/dx (e/cm) and dE/dx are exported side by
side; no recombination model is applied here.

Truth donor (doc 24 sec 6, PDVD -beam): the gapfill retile rebuilds the beam bundle's blobs after beam_deghost and
the rebuilt blobs carry no truth cloud.  --truth-donor <tree> (the cascade tree labelled by TruthPointLabeler with the
clustering job's settings: wct-pr-perevt.jsonnet beam_deghost_truth_dump) fills ONLY live blob nodes without a truth
cloud, from the donor blob with the same node_matches_blob keys (wpid, slice_index_min, u / v / w wire index
min / max): the donor rows as they are when the 3d points are identical, else the row of the nearest donor point of
that blob.  Nodes with no donor blob (the retile's gap-fill blobs) stay unlabelled (-3).  The retile write-back also
stores no blob charge (it lives in the Steiner graph), so a filled node with charge <= 0 takes the donor blob's charge
as its weight (the same cell's solved charge).  Counts in "truth_donor".
Without --truth-donor the outputs are those of the version before it.

  pr_truth.py --pctree <pr dump> --calib <calib-pr-evt*.json> --depos <sim-depos.tar.bz2> --truth <truth.npz>
              --wires <wires> --depo-toff-us T --reco-toff-us R --speed-mmus V --event NAME --out DIR [--truth-donor TREE]
"""
import argparse, json, os, sys
import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from steiner_truth import load_tree, per_node, build_trajectories, KeyTraj, true_vertices  # noqa: E402
import truth_table as T  # noqa: E402

CM = 10.0
TRACK, SHOWER, MICHEL, DELTA, LE = 0, 1, 2, 3, 4
# doc 24: generator labels of cosmic rays (SBND CORSIKA, ProtoDUNE cosmicgenerator) and of the ProtoDUNE radiologicals
COSMIC_GEN = {"corsika", "cosmicgenerator"}
RADIO_GEN = {"ar39", "ar42", "kr85", "k42fromar42", "rn222", "po210", "k40", "co60", "u238", "th232", "ar42fromk42"}
SCORED = ("nu", "gun")       # interaction kinds whose vertex is scored


# ---------------------------------------------------------------------------- inputs
BLOB_KEYS = ("wpid", "slice_index_min", "u_wire_index_min", "u_wire_index_max", "v_wire_index_min", "v_wire_index_max",
             "w_wire_index_min", "w_wire_index_max")   # TruthPointLabel.cxx node_matches_blob


def donor_truth(tree, donor_fn):
    """truth clouds for the live blob nodes of `tree` that have none, from the donor tree (see the module doc)."""
    def keyed(t):
        S, P = per_node(t, "scalar"), per_node(t, "3d")
        out = {}
        for nd, sc in S.items():
            if nd in P and all(k in sc for k in BLOB_KEYS):
                out.setdefault(tuple(int(sc[k][0]) for k in BLOB_KEYS), []).append(nd)
        return out, S, P
    dt = load_tree(donor_fn)
    DT = per_node(dt, "truth")
    dkey, DS, DP = keyed(dt)
    _, S, P = keyed(tree)
    TR = per_node(tree, "truth")
    st = {k: 0 for k in ("nodes", "points", "q", "exact_nodes", "exact_points", "exact_q", "nearest_nodes", "nearest_points",
                         "nearest_q", "none_nodes", "none_points", "none_q", "ambiguous_nodes", "q_from_donor_nodes",
                         "q_from_donor")}
    fill = {}
    for nd, sc in S.items():
        if nd in TR or nd not in P or not all(k in sc for k in BLOB_KEYS):
            continue
        X = np.column_stack([P[nd]["x"], P[nd]["y"], P[nd]["z"]]).astype(float)
        q = float(sc["charge"][0]) if "charge" in sc else 0.0
        st["nodes"] += 1; st["points"] += len(X); st["q"] += q
        cand = [d for d in dkey.get(tuple(int(sc[k][0]) for k in BLOB_KEYS), []) if d in DT]
        if not cand:
            st["none_nodes"] += 1; st["none_points"] += len(X); st["none_q"] += q
            continue
        if len(cand) > 1:
            st["ambiguous_nodes"] += 1
        best = None
        for d in cand:
            Y = np.column_stack([DP[d]["x"], DP[d]["y"], DP[d]["z"]]).astype(float)
            if len(Y) == len(X) and np.array_equal(Y, X):
                best = (d, None); break
            if best is None:
                best = (d, cKDTree(Y).query(X)[1])
        d, idx = best
        dq = float(DS[d]["charge"][0]) if "charge" in DS[d] else 0.0
        if not q > 0:
            st["q_from_donor_nodes"] += 1; st["q_from_donor"] += dq; q = dq
        if idx is None:
            fill[nd] = {k: np.asarray(v) for k, v in DT[d].items()}
            fill[nd]["_q_blob"] = q
            st["exact_nodes"] += 1; st["exact_points"] += len(X); st["exact_q"] += q
        else:
            fill[nd] = {k: np.asarray(v)[idx] for k, v in DT[d].items()}
            fill[nd]["_q_blob"] = q
            st["nearest_nodes"] += 1; st["nearest_points"] += len(X); st["nearest_q"] += q
    return fill, st


def event_points(tree, donor=None):
    """every live point of the dump in cluster / node order: raw xyz, truth columns, charge; cluster offsets.
    donor: extra truth clouds for nodes without one (donor_truth)."""
    P3, TR, SC = per_node(tree, "3d"), per_node(tree, "truth"), per_node(tree, "scalar")
    if donor:
        TR = {**TR, **{nd: v for nd, v in donor.items() if nd not in TR}}
    co = tree["cluster_of"]
    nodes = sorted(P3)
    cols = {k: [] for k in ("x", "y", "z", "tid", "pf", "cls", "iid", "pdg", "dx", "dy", "dz", "dedx", "q", "cl")}
    first = {}
    n = 0
    for nd in nodes:
        p = P3[nd]
        m = len(p["x"])
        cl = int(co[nd])
        first.setdefault(cl, n)
        tr = TR.get(nd)
        g = (lambda k, dflt: np.asarray(tr[k]) if tr is not None else np.full(m, dflt))
        sc = SC.get(nd)
        qb = float(sc["charge"][0]) / max(int(sc["npoints"][0]), 1) if sc is not None else 0.0
        if tr is not None and "_q_blob" in tr and not qb > 0:   # --truth-donor: the donor blob's charge
            qb = float(tr["_q_blob"]) / max(int(sc["npoints"][0]) if sc is not None else m, 1)
        for k, v in (("x", p["x"]), ("y", p["y"]), ("z", p["z"]), ("tid", g("tid", -3)), ("pf", g("pf_node", -1)),
                     ("cls", g("cls", -1)), ("iid", g("iid", -1)), ("pdg", g("pdg", 0)), ("dx", g("dir_x", 0.0)),
                     ("dy", g("dir_y", 0.0)), ("dz", g("dir_z", 0.0)), ("dedx", g("dedx_true", 0.0)),
                     ("q", np.full(m, qb)), ("cl", np.full(m, cl))):
            cols[k].append(np.asarray(v))
        n += m
    out = {k: np.concatenate(v) if v else np.zeros(0) for k, v in cols.items()}
    for k in ("tid", "pf", "cls", "iid", "pdg", "cl"):
        out[k] = out[k].astype(np.int64)
    out["X"] = np.column_stack([out["x"], out["y"], out["z"]]).astype(float)
    out["D"] = np.column_stack([out["dx"], out["dy"], out["dz"]]).astype(float)
    out["first"] = first
    return out


def segment_rows(tree, ev):
    """seg (nu_index, seg_id) -> global point rows, from pr_seg_points; join checks."""
    SP, SS = per_node(tree, "pr_seg_points"), per_node(tree, "pr_segments")
    rows, chk = {}, {"points": 0, "src_index": 0, "src_coord_exact": 0, "xyz_mismatch": 0, "segments": 0, "no_cloud": 0}
    for cl, sp in SP.items():
        g = ev["first"][cl] + sp["p_index"].astype(np.int64)
        R = np.column_stack([sp["x_raw"], sp["y_raw"], sp["z_raw"]]).astype(float)
        chk["points"] += len(g)
        chk["src_index"] += int((np.asarray(sp["src"]) == 0).sum())
        # stale stored index: the nearest cluster point, exact when it sits on the cloud row (d_join 0)
        chk["src_coord_exact"] += int(((np.asarray(sp["src"]) == 1) & (np.asarray(sp["d_join"]) == 0)).sum())
        chk["xyz_mismatch"] += int((np.abs(ev["X"][g] - R).max(axis=1) > 0).sum())
        key = np.column_stack([sp["nu_index"].astype(np.int64), sp["seg_id"].astype(np.int64)])
        for (nu, sid) in np.unique(key, axis=0):
            rows[(int(nu), int(sid))] = g[(key[:, 0] == nu) & (key[:, 1] == sid)]
    for cl, ss in SS.items():
        for nu, sid, nc in zip(ss["nu_index"], ss["seg_id"], ss["n_cloud"]):
            chk["segments"] += 1
            if int(nc) == 0:
                chk["no_cloud"] += 1
                rows.setdefault((int(nu), int(sid)), np.zeros(0, dtype=np.int64))
    return rows, chk


def calib_blocks(fn):
    c = json.load(open(fn))
    if c.get("candidates"):
        return [(int(b["nu_index"]), b) for b in c["candidates"]]
    return [(0, c)]


def interaction_rows(truth, trajs):
    """per truth_nu row: kind ("nu" for a neutrino, "gun" for a non-CORSIKA generator row without one, "corsika") and
    the vertex in the reconstructed frame: the median start of the row's visible primaries (as steiner_truth's
    nu_vertices), else the generator vertex (cm -> mm, flagged)."""
    if not truth:
        return {}
    z = np.load(truth)
    nu, pf = z["truth_nu"], z["truth_pf"]
    c, cn = T.pf_col, T.nu_col
    labels = T.u8_to_strings(z["gen_labels"]) if "gen_labels" in z.files else []
    out = {}
    for r in range(len(nu)):
        sel = pf[:, c("nu_row")] == r
        gl = sorted(set(int(v) for v in pf[sel, c("gen_label")]))
        names = [labels[g] for g in gl if 0 <= g < len(labels)]
        pdg = int(nu[r, cn("pdg")])
        kind = ("nu" if abs(pdg) in (12, 14, 16) else
                "corsika" if names and all(n in COSMIC_GEN for n in names) else
                "radiological" if names and all(n in RADIO_GEN for n in names) else "gun")
        prim = pf[sel & (pf[:, c("mother_trackid")] == 0)]
        S = [trajs[int(t)].start for t in prim[:, c("trackid")].astype(int) if int(t) in trajs]
        if S:
            v, src = np.median(np.array(S), axis=0), "primaries"
        else:
            v, src = np.array([nu[r, cn("vtx_x")], nu[r, cn("vtx_y")], nu[r, cn("vtx_z")]]) * CM, "generator"
        out[r] = dict(kind=kind, pdg=pdg, vertex=v, vertex_src=src, edep=float(nu[r, cn("edep")]))
    return out


# ---------------------------------------------------------------------------- matching
def match(rows, unit, w, tot_unit):
    """best unit of an object (rows), its shared weight, purity, completeness, IoU; unit < 0 never matches."""
    W = float(w[rows].sum())
    u = unit[rows]
    ok = u >= 0
    if not ok.any() or W <= 0:
        return -1, 0.0, 0.0, 0.0, 0.0, W
    keys, inv = np.unique(u[ok], return_inverse=True)
    sh = np.bincount(inv, weights=w[rows][ok])
    b = int(np.argmax(sh))                      # ties: the smaller unit (np.unique order)
    best, s = int(keys[b]), float(sh[b])
    t = tot_unit.get(best, s)
    return best, s, s / W, s / t if t > 0 else 0.0, s / (W + t - s) if (W + t - s) > 0 else 0.0, W


def totals(unit, w, sel=None):
    m = unit >= 0 if sel is None else (unit >= 0) & sel
    keys, inv = np.unique(unit[m], return_inverse=True)
    return dict(zip(keys.tolist(), np.bincount(inv, weights=w[m]).tolist()))


def label_objects(objs, unit, ev, tot_event, tot_bundle, min_points):
    """per object: count- and charge-weighted matches."""
    cnt = np.ones(len(unit))
    out = []
    for o in objs:
        r = o["rows"]
        row = dict(o["meta"], npts=int(len(r)), q=float(ev["q"][r].sum()),
                   ghost_frac=float((ev["tid"][r] == -1).mean()) if len(r) else 0.0)
        lab = r[unit[r] >= 0]
        for wn, w in (("n", cnt), ("q", ev["q"])):
            b, s, p, c, iou, W = match(r, unit, w, tot_event[wn])
            row.update({f"best_{wn}": b, f"shared_{wn}": s, f"purity_{wn}": p, f"completeness_{wn}": c, f"iou_{wn}": iou})
            Wl = float(w[lab].sum())
            row[f"purity_lab_{wn}"] = s / Wl if Wl > 0 else 0.0     # among the labelled points (ghosts left out)
            tb = tot_bundle[wn].get(b, 0.0)
            row[f"completeness_bundle_{wn}"] = s / tb if tb > 0 else 0.0
        row["spine"] = int(row["purity_n"] >= 0.5 and row["completeness_n"] >= 0.5)
        row["pandora"] = int(row["purity_n"] >= 0.5 and row["completeness_n"] >= 0.1 and row["shared_n"] >= 5)
        row["small"] = int(len(r) < min_points)
        out.append(row)
    return out


def pq(objs_rows, unit, truth_units, min_points):
    """one-to-one IoU >= 0.5 matching (count), pimm/SPINE PQ: reco objects and truth units with >= min_points."""
    reco = [r for r in objs_rows if len(r) >= min_points]
    tu = {k: n for k, n in truth_units.items() if n >= min_points}
    pairs = []
    for i, r in enumerate(reco):
        u = unit[r]
        for k in np.unique(u[u >= 0]):
            if int(k) not in tu:
                continue
            s = float((u == k).sum())
            iou = s / (len(r) + tu[int(k)] - s)
            if iou > 0.5:
                pairs.append((iou, i, int(k)))
    used_r, used_t, tp = set(), set(), []
    for iou, i, k in sorted(pairs, reverse=True):
        if i in used_r or k in used_t:
            continue
        used_r.add(i); used_t.add(k); tp.append(iou)
    TP, FP, FN = len(tp), len(reco) - len(tp), len(tu) - len(tp)
    sq = float(np.mean(tp)) if tp else 0.0
    rq = TP / (TP + 0.5 * FP + 0.5 * FN) if (TP + FP + FN) else 0.0
    return {"TP": TP, "FP": FP, "FN": FN, "SQ": sq, "RQ": rq, "PQ": sq * rq}


# ---------------------------------------------------------------------------- doc 24: boundaries, PF parents, edges
def keypoints(pl, window=30.0, angle_min=np.radians(15.0)):
    """start / end and the kinks of one Polyline, the C++ DepoPolyline::keypoints rule (mm)."""
    V, s = pl.V, pl.s
    n = len(V)
    ends = [V[0]] + ([V[-1]] if n > 1 else [])
    if n < 3 or pl.length <= 0:
        return ends, []
    i = np.arange(1, n - 1)
    j = np.clip(np.searchsorted(s, s[i] - window, side="right") - 1, 0, None)
    k = np.clip(np.searchsorted(s, s[i] + window, side="left"), None, n - 1)
    ok = (s[i] - s[j] >= 0.5 * window) & (s[k] - s[i] >= 0.5 * window)
    b = V[i] - V[j]
    a = V[k] - V[i]
    nb, na = np.linalg.norm(b, axis=1), np.linalg.norm(a, axis=1)
    ok &= (nb > 0) & (na > 0)
    ang = np.zeros(n)
    c = np.einsum("ij,ij->i", a, b) / np.maximum(na * nb, 1e-12)
    ang[i[ok]] = np.arccos(np.clip(c[ok], -1.0, 1.0))
    kinks = []
    for m in np.flatnonzero(ang >= angle_min):
        if m == 0 or m == n - 1:
            continue
        lo = np.searchsorted(s, s[m] - window, side="left")
        hi = np.searchsorted(s, s[m] + window, side="right")
        w = np.arange(max(lo, 1), min(hi, n - 1))
        w = w[w != m]
        if np.any(ang[w] > ang[m]) or np.any((ang[w] == ang[m]) & (w < m)):
            continue
        kinks.append(V[m])
    return ends, kinks


def boundary_index(trajs, lmin=CM):
    """k-d trees of the boundary points (start / end / kinks of every unit >= lmin) of all keys and per key."""
    allp, own = [], {}
    for key, tr in trajs.items():
        pts = []
        for u in tr.units:
            if u.length < lmin:
                continue
            e, kk = keypoints(u)
            pts += e + kk
        if pts:
            own[key] = np.array(pts)
            allp.append(own[key])
    kd = cKDTree(np.concatenate(allp)) if allp else None
    return kd, own


def pf_parents(ev, truth):
    """PF node -> parent node (-1 at the interaction).  The node's root key is the key none of whose G4 ancestors is
    another key of the same node (most points first on a tie); its parent is the node of the first ancestor that
    carries points of another node."""
    m = (ev["pf"] >= 0) & (ev["tid"] >= 0)
    if not truth or not m.any():
        return {}
    pairs, cnt = np.unique(np.column_stack([ev["tid"][m], ev["pf"][m]]), axis=0, return_counts=True)
    key2node, best = {}, {}
    for (t, k), c in zip(pairs.tolist(), cnt.tolist()):
        if c > best.get(t, 0):
            key2node[t], best[t] = k, c
    z = np.load(truth)["truth_pf"]
    c_ = T.pf_col
    mother = {int(t): int(mo) for t, mo in zip(z[:, c_("trackid")], z[:, c_("mother_trackid")])}

    def ancestors(t):
        out, seen = [], set()
        t = mother.get(t, 0)
        while t > 0 and t not in seen and len(out) < 10000:
            out.append(t); seen.add(t); t = mother.get(t, 0)
        return out
    keys_of = {}
    for t, k in key2node.items():
        keys_of.setdefault(k, []).append(t)
    parent = {}
    for k, ks in keys_of.items():
        kset = set(ks)
        roots = [t for t in ks if not any(a in kset for a in ancestors(t))] or ks
        root = max(roots, key=lambda t: (best[t], -t))
        parent[k] = -1
        for a in ancestors(root):
            n2 = key2node.get(a)
            if n2 is not None and n2 != k:
                parent[k] = n2
                break
    return parent


# ---------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pctree", required=True)
    ap.add_argument("--calib", required=True)
    ap.add_argument("--depos", required=True)
    ap.add_argument("--truth", default="")
    ap.add_argument("--wires", required=True)
    ap.add_argument("--depo-toff-us", type=float, required=True)
    ap.add_argument("--reco-toff-us", type=float, default=0.0)
    ap.add_argument("--speed-mmus", type=float, required=True)
    ap.add_argument("--min-points", type=int, default=10)
    ap.add_argument("--event", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--truth-donor", default="", help="labelled cascade tree for blobs without a truth cloud (doc 24 sec 6)")
    ap.add_argument("--tie-order", action="store_true", help="order equal-time depos along the track (pr_ml doc 25 sec 7; default off)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    truth = a.truth or None
    res = {"event": a.event, "inputs": {k: v for k, v in vars(a).items() if k != "truth_donor" or v}}

    tree = load_tree(a.pctree)
    donor = None
    if a.truth_donor:
        donor, res["truth_donor"] = donor_truth(tree, a.truth_donor)
    ev = event_points(tree, donor)
    seg_rows, chk = segment_rows(tree, ev)
    res["join"] = chk
    blocks = calib_blocks(a.calib)
    trajs, tinfo = build_trajectories(a.depos, truth, a.wires, a.depo_toff_us, a.reco_toff_us, a.speed_mmus, tie_order_on=a.tie_order)
    res.update(tinfo)
    ptrajs = {}
    for k, tr in trajs.items():
        for o, u in zip(tr.origs, tr.units):
            ptrajs[o] = KeyTraj(o, [u], [o], [0.0])
    irows = interaction_rows(truth, trajs)
    nuv = [d["vertex"] for d in irows.values() if d["kind"] in SCORED]
    vtx = true_vertices(ptrajs, 1 * CM, 1 * CM)
    res["interactions"] = {int(r): dict(kind=d["kind"], pdg=d["pdg"], vertex_mm=np.round(d["vertex"], 2).tolist(),
                                        vertex_src=d["vertex_src"], edep=d["edep"]) for r, d in irows.items()}

    # PF nodes: dominant key, its pdg / primary flag from the truth table, majority class and interaction
    pf_pdg, pf_prim = {}, {}
    if truth:
        z = np.load(truth)["truth_pf"]
        c = T.pf_col
        tid2row = {int(t): i for i, t in enumerate(z[:, c("trackid")])}
    pf_info = {}
    for k in np.unique(ev["pf"][ev["pf"] >= 0]):
        m = ev["pf"] == k
        keys, cnt = np.unique(ev["tid"][m], return_counts=True)
        key = int(keys[np.argmax(cnt)])
        cl_, cc = np.unique(ev["cls"][m], return_counts=True)
        ii, ic = np.unique(ev["iid"][m], return_counts=True)
        pdg, prim = 0, 0
        if truth and key in tid2row:
            pdg = int(z[tid2row[key], c("pdg")]); prim = int(z[tid2row[key], c("mother_trackid")] == 0)
        pf_info[int(k)] = dict(key=key, pdg=pdg, primary=prim, cls=int(cl_[np.argmax(cc)]), iid=int(ii[np.argmax(ic)]),
                               npts=int(m.sum()), q=float(ev["q"][m].sum()))
    parent_of = pf_parents(ev, truth)            # doc 24
    bkd, bown = boundary_index(trajs)
    res["n_interactions_visible"] = int(len(np.unique(ev["iid"][ev["iid"] >= 0])))
    res["interaction_points"] = {int(k): int(n) for k, n in zip(*np.unique(ev["iid"][ev["iid"] >= 0], return_counts=True))}

    cnt = np.ones(len(ev["tid"]))
    tot_event = {"pf": {"n": totals(ev["pf"], cnt), "q": totals(ev["pf"], ev["q"])},
                 "tid": {"n": totals(ev["tid"], cnt), "q": totals(ev["tid"], ev["q"])},
                 "iid": {"n": totals(ev["iid"], cnt), "q": totals(ev["iid"], ev["q"])}}

    seg_out, par_out, int_out, vtx_out, fit_out, edge_out = [], [], [], [], [], []
    res["candidates"] = []
    for nu, blk in blocks:
        segs = blk.get("segments", [])
        sid = [int(s["id"]) for s in segs]
        rows_of = {s: seg_rows.get((nu, s), np.zeros(0, dtype=np.int64)) for s in sid}
        allr = np.unique(np.concatenate([r for r in rows_of.values()] + [np.zeros(0, dtype=np.int64)]))
        bundle_cl = np.isin(ev["cl"], np.unique(ev["cl"][allr])) if len(allr) else np.zeros(len(ev["cl"]), bool)
        tot_bundle = {u: {"n": totals(ev[u], cnt, bundle_cl), "q": totals(ev[u], ev["q"], bundle_cl)} for u in ("pf", "tid", "iid")}
        cand = {"nu_index": nu, "n_segments": len(segs), "n_segments_in_export": int(sum((nu, s) in seg_rows for s in sid)),
                "n_export_not_in_calib": int(sum(1 for (n2, s2) in seg_rows if n2 == nu and s2 not in set(sid))),
                "n_points_assigned": int(len(allr))}

        # segments
        sobjs = [{"rows": rows_of[int(s["id"])], "meta": dict(nu_index=nu, seg_id=int(s["id"]), flag_shower=int(bool(s["flag_shower"])),
                  reco_pdg=int(s["particle_id"]), dirsign=int(s["dirsign"]), shower_id=int(s["shower_id"]),
                  length_cm=float(s.get("length", 0.0)), n_fit=len(s.get("points", [])))} for s in segs]
        srows = label_objects(sobjs, ev["pf"], ev, tot_event["pf"], tot_bundle["pf"], a.min_points)
        tl = label_objects(sobjs, ev["tid"], ev, tot_event["tid"], tot_bundle["tid"], a.min_points)
        for s, r, t2, o in zip(segs, srows, tl, sobjs):
            pi = pf_info.get(r["best_n"], {})
            # pdg of the best G4 key (a PF node can absorb other particles: nuclear fragments, scattered protons,
            # so its dominant pdg is not the segment's); class and primary flag of the PF node
            bt = t2["best_n"]
            tpdg = int(z[tid2row[bt], c("pdg")]) if truth and bt in tid2row else 0
            r.update(true_pdg=tpdg, pf_pdg=pi.get("pdg", 0), true_cls=pi.get("cls", -1), true_primary=pi.get("primary", 0),
                     best_tid=bt, purity_tid=t2["purity_n"], completeness_tid=t2["completeness_n"],
                     true_parent_node=parent_of.get(r["best_n"], -1) if r["best_n"] >= 0 else -2)
            # dirsign_true: true travel direction along the fit-point order, from the matched key's points
            pts = s.get("points", [])
            r["dirsign_true"] = 0
            if len(pts) >= 2 and len(o["rows"]):
                A = np.array([pts[0]["x"], pts[0]["y"], pts[0]["z"]]) * CM
                B = np.array([pts[-1]["x"], pts[-1]["y"], pts[-1]["z"]]) * CM
                rr = o["rows"][ev["tid"][o["rows"]] == r["best_tid"]]
                if len(rr) and np.linalg.norm(B - A) > 0:
                    v = float(np.einsum("ij,j->", ev["D"][rr], (B - A) / np.linalg.norm(B - A)))
                    r["dirsign_true"] = int(np.sign(v))
            # fit points (sec 13): truth beside the PR chain's own fit outputs
            key = r["best_tid"]
            if pts:
                F = np.array([[p["x"], p["y"], p["z"]] for p in pts]) * CM
                dtr = np.full(len(F), np.nan)
                if key in trajs:
                    dtr = trajs[key].nearest(F)[0]
                dedx = np.full(len(F), np.nan)
                rr = o["rows"][ev["tid"][o["rows"]] == key] if len(o["rows"]) else o["rows"]
                if len(rr):
                    _, j = cKDTree(ev["X"][rr]).query(F)
                    dedx = ev["dedx"][rr[j]]
                num = (lambda v: float("nan") if v is None else float(v))   # the JSON writes NaN as null
                Fq = np.where(np.isfinite(F), F, 0.0)                       # doc 24: distance to the true boundaries
                dba = bkd.query(Fq)[0] if bkd is not None else np.full(len(F), np.inf)
                dbo = cKDTree(bown[key]).query(Fq)[0] if key in bown else np.full(len(F), np.inf)
                bad = ~np.isfinite(F).all(axis=1)
                dba[bad] = np.nan; dbo[bad] = np.nan
                for i, p in enumerate(pts):
                    fit_out.append(dict(nu_index=nu, seg_id=r["seg_id"], i=i, x=num(p["x"]), y=num(p["y"]), z=num(p["z"]),
                                        dQ=num(p["dQ"]), dx=num(p["dx"]), rr=num(p.get("rr", -1.0)),
                                        reduced_chi2=num(p.get("reduced_chi2", 0.0)),
                                        best_tid=key, true_cls=r["true_cls"], true_pdg=r["true_pdg"], purity=r["purity_lab_n"],
                                        d_true=float(dtr[i]), dedx_true=float(dedx[i]),
                                        d_boundary=float(dba[i]), d_boundary_own=float(dbo[i]), boundary=int(dba[i] <= CM)))
        seg_out += srows
        # doc 24: segment-pair edges (adjacent pairs, and every pair of segments with >= min_points)
        vset = [{int(s.get("start_vertex_id", -1)), int(s.get("end_vertex_id", -1))} - {-1} for s in segs]
        for ia in range(len(segs)):
            ra = srows[ia]
            for ib in range(ia + 1, len(segs)):
                rb = srows[ib]
                adj = bool(vset[ia] & vset[ib])
                if not adj and (ra["small"] or rb["small"]):
                    continue
                na, nb = ra["best_n"], rb["best_n"]
                both = na >= 0 and nb >= 0
                ia_, ib_ = pf_info.get(na, {}).get("iid", -1), pf_info.get(nb, {}).get("iid", -1)
                edge_out.append(dict(nu_index=nu, seg_a=ra["seg_id"], seg_b=rb["seg_id"], adjacent=int(adj),
                                     small_a=ra["small"], small_b=rb["small"],
                                     clean=int(not ra["small"] and not rb["small"] and ra["purity_lab_n"] >= 0.8 and rb["purity_lab_n"] >= 0.8),
                                     reco_same_shower=int(ra["shower_id"] >= 0 and ra["shower_id"] == rb["shower_id"]),
                                     same_pf=int(both and na == nb),
                                     same_shower_true=int(both and na == nb and pf_info.get(na, {}).get("cls") == SHOWER),
                                     same_interaction=int(both and ia_ >= 0 and ia_ == ib_),
                                     parent_ab=int(both and na != nb and parent_of.get(nb) == na),
                                     parent_ba=int(both and na != nb and parent_of.get(na) == nb)))

        # particles: showers (union of their segments) and the segments in no shower
        showers = {int(sh["id"]): sh for sh in blk.get("showers", [])}
        members = {}
        for s in segs:
            if int(s["shower_id"]) >= 0:
                members.setdefault(int(s["shower_id"]), []).append(int(s["id"]))
        pobjs = []
        for shid, ms in sorted(members.items()):
            sh = showers.get(shid, {})
            pobjs.append({"rows": np.unique(np.concatenate([rows_of[m] for m in ms])),
                          "meta": dict(nu_index=nu, obj_id=shid, kind=1, reco_pdg=int(sh.get("particle_id", 0)), n_seg=len(ms),
                                       kine_best=float("nan") if sh.get("kine_best") is None else float(sh["kine_best"]))})   # JSON null = NaN
        for s in segs:
            if int(s["shower_id"]) < 0:
                pobjs.append({"rows": rows_of[int(s["id"])], "meta": dict(nu_index=nu, obj_id=int(s["id"]), kind=0,
                              reco_pdg=int(s["particle_id"]), n_seg=1, kine_best=0.0)})
        prows = label_objects(pobjs, ev["pf"], ev, tot_event["pf"], tot_bundle["pf"], a.min_points)
        for r in prows:
            pi = pf_info.get(r["best_n"], {})
            r.update(true_pdg=pi.get("pdg", 0), true_cls=pi.get("cls", -1), true_primary=pi.get("primary", 0))
        par_out += prows
        # truth-side efficiency and PQ, particles vs PF nodes of the matched interaction with >= min_points in the bundle
        irow = label_objects([{"rows": allr, "meta": dict(nu_index=nu)}], ev["iid"], ev, tot_event["iid"], tot_bundle["iid"], a.min_points)[0]
        int_out.append(irow)
        nu_iid = irow["best_n"]
        tu = {k: int(n) for k, n in tot_bundle["pf"]["n"].items() if pf_info.get(k, {}).get("iid") == nu_iid}
        big = {k for k, n in tu.items() if n >= a.min_points}
        sp = {r["best_n"] for r in prows if r["spine"]}
        pa = {r["best_n"] for r in prows if r["pandora"]}
        cand["truth_particles"] = len(big)
        cand["efficiency_spine"] = len(big & sp) / len(big) if big else None
        cand["efficiency_pandora"] = len(big & pa) / len(big) if big else None
        cand["pq"] = pq([o["rows"] for o in pobjs], ev["pf"], tu, a.min_points)
        cand["interaction"] = {k: irow[k] for k in ("best_n", "purity_n", "purity_lab_n", "completeness_n", "completeness_bundle_n", "purity_q", "completeness_q")}
        cand["interaction"]["kind"] = irows.get(nu_iid, {}).get("kind", "none")

        # vertices: right-vertex flags and the neutrino vertex distance
        NV = np.array(nuv) if nuv else np.zeros((0, 3))
        TVX = np.array([p for _, p in vtx]) if vtx else np.zeros((0, 3))
        for v in blk.get("vertices", []):
            f = v.get("fit", {})
            P = np.array([f.get("x", 0.0), f.get("y", 0.0), f.get("z", 0.0)]) * CM
            dn = float(np.min(np.linalg.norm(NV - P, axis=1))) if len(NV) else np.inf
            dv = float(np.min(np.linalg.norm(TVX - P, axis=1))) if len(TVX) else np.inf
            vtx_out.append(dict(nu_index=nu, vtx_id=int(v["id"]), is_main=int(bool(v.get("is_main"))), degree=int(v.get("degree", 0)),
                                d_nu=dn, d_vertex=dv, is_nu_vertex=int(dn <= CM), is_vertex=int(dv <= CM)))
        mv = blk.get("main_vertex")
        if mv and nu_iid in irows and irows[nu_iid]["kind"] in SCORED:
            # the vertex of the interaction the candidate matched (doc 23 sec 11: never a CORSIKA row)
            P = np.array([mv["x"], mv["y"], mv["z"]]) * CM
            TV = irows[nu_iid]["vertex"]
            cand["nu_vertex"] = {"reco_mm": np.round(P, 2).tolist(), "true_mm": np.round(TV, 2).tolist(),
                                 "true_src": irows[nu_iid]["vertex_src"], "dist_mm": float(np.linalg.norm(P - TV)),
                                 "dxyz_mm": np.round(P - TV, 2).tolist()}
        res["candidates"].append(cand)

    def summ(rows, extra=()):
        big = [r for r in rows if not r["small"]]
        out = {"n": len(rows), "n_ge_min": len(big)}
        if big:
            out.update(spine=float(np.mean([r["spine"] for r in big])), pandora=float(np.mean([r["pandora"] for r in big])),
                       purity_q_median=float(np.median([r["purity_q"] for r in big])),
                       completeness_q_median=float(np.median([r["completeness_q"] for r in big])),
                       ghost_frac_median=float(np.median([r["ghost_frac"] for r in big])))
        return out
    res["segments"] = summ(seg_out)
    res["particles"] = summ(par_out)
    # track / shower and direction of the segments with a clean match
    clean = [r for r in seg_out if not r["small"] and r["purity_lab_n"] >= 0.8 and r["true_cls"] >= 0]
    ts = np.zeros((2, 2), dtype=int)          # rows: true track / EM; columns: reco track / shower
    for r in clean:
        ts[int(r["true_cls"] != TRACK), r["flag_shower"]] += 1
    res["track_shower"] = {"true_track_reco_track": int(ts[0, 0]), "true_track_reco_shower": int(ts[0, 1]),
                           "true_em_reco_track": int(ts[1, 0]), "true_em_reco_shower": int(ts[1, 1])}
    dd = [r for r in clean if r["true_cls"] == TRACK and r["length_cm"] >= 5 and r["dirsign_true"] != 0]
    res["direction"] = {"n": len(dd), "reco_dirsign_0": int(sum(r["dirsign"] == 0 for r in dd)),
                        "agree": int(sum(r["dirsign"] == r["dirsign_true"] for r in dd)),
                        "opposite": int(sum(r["dirsign"] == -r["dirsign_true"] for r in dd))}
    # doc 24: fit-point boundaries and the segment-pair edges
    if fit_out:
        db = np.array([f["d_boundary"] for f in fit_out])
        res["boundary"] = {"n_fit": len(fit_out), "boundary": int(sum(f["boundary"] for f in fit_out)),
                           "d_boundary_median_mm": float(np.nanmedian(db)) if np.isfinite(db).any() else None}
    E = edge_out
    ce = [e for e in E if e["clean"]]
    res["edges"] = {"n": len(E), "adjacent": int(sum(e["adjacent"] for e in E)), "n_clean": len(ce),
                    "same_pf": int(sum(e["same_pf"] for e in ce)), "same_shower_true": int(sum(e["same_shower_true"] for e in ce)),
                    "reco_same_shower": int(sum(e["reco_same_shower"] for e in ce)),
                    "shower_both": int(sum(e["reco_same_shower"] and e["same_shower_true"] for e in ce)),
                    "same_interaction": int(sum(e["same_interaction"] for e in ce)),
                    "parent": int(sum(e["parent_ab"] or e["parent_ba"] for e in ce)),
                    "parent_adjacent": int(sum((e["parent_ab"] or e["parent_ba"]) and e["adjacent"] for e in ce))}
    pid = {}
    for r in clean:
        if r["true_cls"] == TRACK and not r["flag_shower"] and abs(r["true_pdg"]) in (13, 211, 321, 2212):
            k = f"{abs(r['true_pdg'])}->{abs(r['reco_pdg'])}"
            pid[k] = pid.get(k, 0) + 1
    res["pid_tracks"] = dict(sorted(pid.items()))
    fits = [f for f in fit_out if f["true_cls"] == TRACK and f["purity"] >= 0.8 and np.isfinite(f["d_true"])]
    if fits:
        d = np.array([f["d_true"] for f in fits])
        res["fit_points"] = {"n": len(fit_out), "n_track_clean": len(fits), "d_true_median_mm": float(np.median(d)),
                             "d_true_le_5mm": float((d <= 5).mean()), "d_true_le_15mm": float((d <= 15).mean())}
        for pdg in (13, 2212, 211):
            g = [f for f in fits if abs(f["true_pdg"]) == pdg and f["dx"] > 0 and f["dedx_true"] > 0]
            if len(g) >= 20:
                ratio = np.array([f["dQ"] / f["dx"] / f["dedx_true"] for f in g])
                res["fit_points"][f"dqdx_over_dedx_{pdg}"] = {"n": len(g), "median": float(np.median(ratio)),
                                                              "p16": float(np.percentile(ratio, 16)), "p84": float(np.percentile(ratio, 84))}

    def to_npz(rows, prefix, out):
        if not rows:
            return
        for k in rows[0]:
            out[prefix + k] = np.array([r[k] for r in rows])
    arr = {}
    for rows, pre in ((seg_out, "seg_"), (par_out, "par_"), (int_out, "int_"), (vtx_out, "vtx_"), (fit_out, "fit_"), (edge_out, "edge_")):
        to_npz(rows, pre, arr)
    np.savez_compressed(os.path.join(a.out, a.event + "_rlabels.npz"), **arr)
    json.dump(res, open(os.path.join(a.out, a.event + ".json"), "w"), indent=1, default=float)
    c0 = res["candidates"][0] if res["candidates"] else {}
    print(a.event, "join", chk, "segs", res["segments"], "particles", res["particles"],
          "eff", c0.get("efficiency_spine"), c0.get("efficiency_pandora"), "PQ", c0.get("pq", {}).get("PQ"),
          "nu_vtx_mm", c0.get("nu_vertex", {}).get("dist_mm"), "dir", res["direction"], "fit", res.get("fit_points", {}).get("d_true_median_mm"))


if __name__ == "__main__":
    main()
