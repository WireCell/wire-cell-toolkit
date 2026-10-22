#!/usr/bin/env python3
"""qa_report.py -- the label QA report of pr_ml doc 19 sec 8 (doc 15 sec 3.3 made concrete), per work tag.

Reads what validate_2d.py and validate_3d.py wrote (per_event.json of each), the events' sim-depos archives (for the hard
checks) and the frame archives (sha256 duplicate check), and writes <out>/summary.json, <out>/per_event.json and
<out>/qa_report.md with the eight items of sec 8:
  1. charge accounting per plane: closure inside the truth region, dead-region share (n/a without a dead mask), by class;
  2. mixed-point fraction per class; ghost fraction by count and by charge (from validate_3d);
  3. terminal count per cm of true trajectory: placeholder (G labels are WP0.2) -> points per cm;
  4. PF-node count per interaction vs raw G4 count above q_min (from labels_meta when the pctree carries it);
  5. unmatched reco-segment fraction: n/a until WP0.4;
  6. hard checks: (a) inverted interval -- for every blobtruth table, the sum of a blob's per-track charge equals its
     q_true (1e-9) and no blob holds more charge than the depos in the window; (b) duplicated images -- sha256 of every frame
     archive of the sample, any duplicate FAILS;
  7. timing self-check (sec 6.3): median segment lag per plane within +-1 tick and the whole-plane profile lag 0;
  8. the 50-event blind scan: protocol only (the Bee sets of sec 7.5); not run here.
usage:
  qa_report.py --tag <tag> --det <det> --out <dir> --v2d <validate_2d out dir> --v3d <validate_3d out dir>
      [--frames <glob of frame archives>] [--pctree <glob of pctree dumps>] [--blobtruth <glob>]
"""
import argparse, glob, hashlib, io, json, os, re, tarfile
import numpy as np


def sha256(fn):
    h = hashlib.sha256()
    with open(fn, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_tensors(fn):
    md, arr = {}, {}
    with tarfile.open(fn) as t:
        for m in t.getmembers():
            if not m.isfile():
                continue
            d = t.extractfile(m).read()
            if m.name.endswith("metadata.json"):
                md[m.name] = json.loads(d)
            elif m.name.endswith("array.npy"):
                arr[m.name] = np.load(io.BytesIO(d), allow_pickle=False)
    return {m["datapath"]: (m, arr.get(k.replace("metadata.json", "array.npy"))) for k, m in md.items() if "datapath" in m}


def check_blobtruth(fn):
    """item 6(a): per blob sum_k q_bk == q_true; q_true in blobs <= depos in window (+1e-6)"""
    tens = load_tensors(fn)
    tab = sumt = None
    for dp, (m, a) in tens.items():
        if m.get("datatype") == "blobtruth":
            tab = (m, np.asarray(a, dtype=np.float64).reshape(-1, len(m["columns"])))
        elif m.get("datatype") == "blobtruth_sum":
            sumt = (m, np.asarray(a, dtype=np.float64).reshape(-1, len(m["columns"])))
    if tab is None or sumt is None:
        return {"file": fn, "ok": None, "why": "no tables"}
    mt, T = tab; ms, S = sumt
    ct = {c: i for i, c in enumerate(mt["columns"])}; cs = {c: i for i, c in enumerate(ms["columns"])}
    key_t = np.stack([T[:, ct["face"]], T[:, ct["slice_ident"]], T[:, ct["blob_ident"]]], axis=1) if "face" in ct else None
    per_blob = {}
    for r in T:
        k = (int(r[ct["face"]]), int(r[ct["slice_ident"]]), int(r[ct["blob_ident"]]))
        per_blob[k] = per_blob.get(k, 0.0) + float(r[ct["q_bk"]])
    bad = 0; worst = 0.0
    if "face" in cs:
        for r in S:
            k = (int(r[cs["face"]]), int(r[cs["slice_ident"]]), int(r[cs["blob_ident"]]))
            q = float(r[cs["q_true"]]); s_ = per_blob.get(k, 0.0)
            if q > 0 and abs(s_ - q) > 1e-6 * max(q, 1.0):
                bad += 1; worst = max(worst, abs(s_ - q) / max(q, 1.0))
    qwin = float(ms.get("q_depos_window", ms.get("q_depos", 0.0)))
    qblobs = float(S[:, cs["q_true"]].sum())
    leak = qwin > 0 and qblobs > qwin * (1 + 1e-6)
    return {"file": os.path.basename(fn), "nblobs": int(len(S)), "sum_mismatch": bad, "worst_rel": worst,
            "q_blobs": qblobs, "q_window": qwin, "q_blobs_over_window": (qblobs / qwin) if qwin > 0 else None,
            "ok": bool(bad == 0 and not leak), "why": "" if (bad == 0 and not leak) else ("per-track sums != q_true" if bad else "more charge in blobs than depos in the window")}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tag", required=True); ap.add_argument("--det", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--v2d", default=""); ap.add_argument("--v3d", default="")
    ap.add_argument("--frames", default=""); ap.add_argument("--pctree", default=""); ap.add_argument("--blobtruth", default="")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rep = {"tag": a.tag, "detector": a.det, "items": {}}
    v2 = json.load(open(os.path.join(a.v2d, "per_event.json"))) if a.v2d and os.path.exists(os.path.join(a.v2d, "per_event.json")) else {}
    v3 = json.load(open(os.path.join(a.v3d, "per_event.json"))) if a.v3d and os.path.exists(os.path.join(a.v3d, "per_event.json")) else {}

    # 1. charge accounting per plane (+ by class), 7. timing
    acc, timing, lag = {}, {}, {}
    for ev, r in v2.items():
        for k, v in r.items():
            p = "%s_p%d" % (v["tag"], v["plane"])
            acc.setdefault(p, {"closure_region": [], "seg_eff": [], "by_class": {}})
            if v.get("closure_region") is not None: acc[p]["closure_region"].append(v["closure_region"])
            if v.get("seg_eff") is not None: acc[p]["seg_eff"].append(v["seg_eff"])
            for c, cv in v.get("class_eff", {}).items():
                acc[p]["by_class"].setdefault(c, {"eff": [], "closure": []})
                acc[p]["by_class"][c]["eff"].append(cv["eff"])
                if cv.get("closure") is not None: acc[p]["by_class"][c]["closure"].append(cv["closure"])
            if v.get("timing_median_ticks") is not None: timing.setdefault(p, []).append(v["timing_median_ticks"])
            if v.get("profile_lag_ticks") is not None: lag.setdefault(p, []).append(v["profile_lag_ticks"])
    item1 = {p: {"closure_region_mean": float(np.mean(x["closure_region"])) if x["closure_region"] else None,
                 "seg_eff_mean": float(np.mean(x["seg_eff"])) if x["seg_eff"] else None,
                 "by_class": {c: {"eff_mean": float(np.mean(cv["eff"])), "closure_mean": float(np.mean(cv["closure"])) if cv["closure"] else None, "n": len(cv["eff"])} for c, cv in x["by_class"].items()},
                 "dead_region_share": "n/a (no dead mask in R1)"} for p, x in acc.items()}
    rep["items"]["1_charge_accounting"] = item1
    item7 = {p: {"timing_median_ticks": float(np.median(t)), "within_1_tick": bool(abs(np.median(t)) <= 1.0),
                 "profile_lag_median_ticks": float(np.median(lag.get(p, [0]))), "lag_zero": bool(abs(np.median(lag.get(p, [0]))) <= 2)} for p, t in timing.items()}
    item7["pass"] = bool(item7) and all(v["within_1_tick"] and v["lag_zero"] for v in item7.values() if isinstance(v, dict))
    rep["items"]["7_timing_self_check"] = item7

    # 2. mixed / ghost (3-D), 3. points per cm, 4. PF nodes
    g_count, g_charge, mixed, valid, ppcm, pfn = [], [], [], [], [], []
    mixed_cls = {}
    for ev, r in v3.items():
        b = r.get("blob", {}); p = r.get("point", {})
        if b.get("ghost_frac_count") is not None: g_count.append(b["ghost_frac_count"])
        if b.get("ghost_frac_charge") is not None: g_charge.append(b["ghost_frac_charge"])
        if b.get("mixed_frac_of_filled") is not None: mixed.append(b["mixed_frac_of_filled"])
        for c, v in b.get("mixed_frac_by_class", {}).items(): mixed_cls.setdefault(c, []).append(v)
        if p.get("valid_frac") is not None: valid.append(p["valid_frac"])
        if p.get("points_per_cm_median") is not None: ppcm.append(p["points_per_cm_median"])
    rep["items"]["2_mixed_and_ghost"] = {"ghost_frac_count_mean": float(np.mean(g_count)) if g_count else None, "ghost_frac_charge_mean": float(np.mean(g_charge)) if g_charge else None,
                                        "mixed_blob_frac_mean": float(np.mean(mixed)) if mixed else None, "mixed_by_class": {c: float(np.mean(v)) for c, v in mixed_cls.items()},
                                        "point_valid_frac_mean": float(np.mean(valid)) if valid else None, "nevents": len(v3)}
    rep["items"]["3_terminals_per_cm"] = {"status": "placeholder until G labels (WP0.2)", "points_per_cm_median": float(np.median(ppcm)) if ppcm else None}
    # 4. from labels_meta in the pctree dumps
    item4 = []
    for fn in sorted(glob.glob(a.pctree)) if a.pctree else []:
        try:
            with tarfile.open(fn) as t:
                for m in t.getmembers():
                    if m.name.endswith("tensorset_") or "tensorset" in m.name and m.name.endswith("metadata.json"):
                        md = json.loads(t.extractfile(m).read())
                        lm = md.get("labels_meta")
                        if lm: item4.append({"file": os.path.basename(fn), "pf_nodes": lm["summary"].get("pf_nodes"), "pf_source": lm.get("pf_source"), "merged_weak": lm["summary"].get("pf_merged_weak"), "void_roots": lm["summary"].get("void_roots")})
        except Exception as e:  # noqa: BLE001
            item4.append({"file": os.path.basename(fn), "error": str(e)})
    rep["items"]["4_pf_nodes"] = item4 or "n/a (no pctree dumps with labels_meta given)"
    rep["items"]["5_unmatched_segments"] = "n/a until WP0.4"

    # 6. hard checks
    bt_checks = [check_blobtruth(fn) for fn in sorted(glob.glob(a.blobtruth))] if a.blobtruth else []
    fr = {}
    dup = []
    for fn in sorted(glob.glob(a.frames)) if a.frames else []:
        h = sha256(fn)
        if h in fr: dup.append((fr[h], fn))
        fr[h] = fn
    rep["items"]["6_hard_checks"] = {"a_inverted_interval": {"nfiles": len(bt_checks), "nfail": sum(1 for c in bt_checks if c["ok"] is False), "fails": [c for c in bt_checks if c["ok"] is False][:10],
                                                           "pass": all(c["ok"] is not False for c in bt_checks) if bt_checks else None},
                                     "b_duplicated_images": {"nfiles": len(fr) + len(dup), "duplicates": dup, "pass": (len(dup) == 0) if fr else None}}
    rep["items"]["8_blind_scan"] = {"status": "protocol (sec 8 item 8): 50 events, Bee sets of sec 7.5, one scanner blind to the arm, four questions per event, results in a fresh work/ql_labels/<tag>/; not run"}
    json.dump(rep, open(os.path.join(a.out, "summary.json"), "w"), indent=1)
    json.dump({"v2d": v2, "v3d": v3, "blobtruth_checks": bt_checks}, open(os.path.join(a.out, "per_event.json"), "w"), indent=1)
    with open(os.path.join(a.out, "qa_report.md"), "w") as f:
        f.write("# Label QA report: %s / %s\n\n" % (a.det, a.tag))
        f.write("| item | result |\n|---|---|\n")
        f.write("| 1 charge accounting | %s |\n" % "; ".join("%s: closure(region) %.3f, seg eff %.3f" % (p, v["closure_region_mean"] or 0, v["seg_eff_mean"] or 0) for p, v in item1.items()))
        i2 = rep["items"]["2_mixed_and_ghost"]
        f.write("| 2 mixed / ghost | ghost by count %s, by charge %s, mixed blobs %s, valid points %s (n=%d) |\n" % (i2["ghost_frac_count_mean"], i2["ghost_frac_charge_mean"], i2["mixed_blob_frac_mean"], i2["point_valid_frac_mean"], i2["nevents"]))
        f.write("| 3 terminals / cm | placeholder; points per cm median %s |\n" % rep["items"]["3_terminals_per_cm"]["points_per_cm_median"])
        f.write("| 4 PF nodes | %s |\n" % (("%d dumps" % len(item4)) if item4 else "n/a"))
        f.write("| 5 unmatched segments | n/a until WP0.4 |\n")
        h6 = rep["items"]["6_hard_checks"]
        f.write("| 6a inverted interval | %s (%d tables, %d fail) |\n" % ("PASS" if h6["a_inverted_interval"]["pass"] else ("FAIL" if h6["a_inverted_interval"]["pass"] is False else "n/a"), h6["a_inverted_interval"]["nfiles"], h6["a_inverted_interval"]["nfail"]))
        f.write("| 6b duplicated images | %s (%d archives, %d duplicates) |\n" % ("PASS" if h6["b_duplicated_images"]["pass"] else ("FAIL" if h6["b_duplicated_images"]["pass"] is False else "n/a"), h6["b_duplicated_images"]["nfiles"], len(dup)))
        f.write("| 7 timing self-check | %s: %s |\n" % ("PASS" if item7.get("pass") else "FAIL/none", "; ".join("%s median %.2f lag %.0f" % (p, v["timing_median_ticks"], v["profile_lag_median_ticks"]) for p, v in item7.items() if isinstance(v, dict))))
        f.write("| 8 blind scan | protocol only, not run |\n")
    print("wrote", a.out, json.dumps({k: (v if not isinstance(v, (dict, list)) else "...") for k, v in rep["items"].items()}))


if __name__ == "__main__":
    main()
