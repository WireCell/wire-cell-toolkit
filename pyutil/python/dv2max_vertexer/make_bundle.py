"""Build the model bundle (bundle.json + pass1.pt + pass2.pt + ranker_fold<k>.txt + isotonic_fold<k>.json) from a
DV2-MAX run directory.  A build tool, not imported at run time; it needs pyyaml and torch (e.g. the sbnd-vertexer
venv).  The ranker text and the isotonic thresholds come from the deterministic fold re-fit (see README.md).

usage: python make_bundle.py --run <sbnd-vertexer>/dv2max/runs/r7_cloud --refit <dir with ranker_fold0.txt,
       isotonic_fold0.json, feature_names.json, refit2_report.json> --out <wct>/cfg/dv2max/r7_cloud [--fold 0]
"""
import argparse
import datetime
import hashlib
import json
import shutil
from pathlib import Path

FORMAT = "dv2max-bundle-v1"
X_SHIFT_CM = 0.31        # r7/prep_r7_cache.py X_SHIFT_CM (reco inputs moved toward the cathode; labels never)
LINK_CM = 3.0            # r7/prep_r7_cache.py LINK_CM (3 cm connected components = cluster ids)


def sha256(path, chunk=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def main():
    import torch
    import yaml
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--refit", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--fold", type=int, default=0)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    run, refit, out = Path(a.run), Path(a.refit), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    k = int(a.fold)
    fam = f"f{k}"
    cfg = yaml.safe_load(open(run / "config.yaml"))
    files, source = {}, {}
    for which in ("pass1", "pass2"):
        src = run / "models" / which / fam / f"seed{a.seed}" / "final.pt"
        st = torch.load(src, map_location="cpu", weights_only=False)
        sd = st["model"]
        bad = [n for n, v in sd.items() if v.is_floating_point() and v.dtype != torch.float32]
        if bad:
            raise ValueError(f"{src}: non-fp32 tensors {bad[:3]}")
        dst = out / f"{which}.pt"
        torch.save(sd, dst)
        files[which] = {"path": dst.name, "sha256": sha256(dst), "bytes": dst.stat().st_size}
        meta = st.get("meta", {})
        source[which] = {"final_pt": str(src), "final_pt_sha256": sha256(src),
                         **{c: meta.get(c) for c in ("train_sig", "train_folds", "steps", "final", "swa_n_averaged",
                                                     "compute_dtype", "backend", "n_params")}}
    for name, fn in (("ranker", f"ranker_fold{k}.txt"), ("isotonic", f"isotonic_fold{k}.json")):
        shutil.copyfile(refit / fn, out / fn)
        files[name] = {"path": fn, "sha256": sha256(out / fn), "bytes": (out / fn).stat().st_size}
    names = json.loads((refit / "feature_names.json").read_text())
    rep = json.loads((refit / "refit2_report.json").read_text()) if (refit / "refit2_report.json").exists() else {}
    p1, p2, rc = cfg["pass1"], cfg["pass2"], cfg["ranker"]
    man = {
        "format": FORMAT,
        "model": cfg["name"],
        "family": fam,
        "family_folds": [k],
        "seed": a.seed,
        "trained_on_folds": source["pass1"]["train_folds"],
        "x_shift_cm": X_SHIFT_CM,
        "link_cm": LINK_CM,
        "frame": "models work in the label frame = raw reco frame with x -= sign(x) * x_shift_cm; DV2_Vertex returns "
                 "x += sign(x) * x_shift_cm (raw reco frame)",
        "input": "one DL-vertex call: cloud x, y, z (cm), q = dQ * q_scale + q_offset; first n_vertex_rows points = "
                 "PR-graph vertex rows (ranker candidates 'cand'); traditional vertex = tag_nu",
        "inference": {
            "device": "cpu", "backend": "torch (dv2max_vertexer.sparse)", "precision": "fp32, no autocast",
            "min_points": cfg["data"]["min_points"],
            "pass1": {"voxel_cm": p1["voxel"], "views": p1["infer"]["views"], "shifts": p1["infer"]["shifts"],
                      "decode": {c: p1["infer"][c] for c in p1["infer"] if c not in ("views", "shifts")}},
            "pass2": {"voxel_cm": p2["voxel"], "views": p2["infer"]["views"], "shifts": p2["infer"]["shifts"],
                      "decode": {c: p2["infer"][c] for c in p2["infer"] if c not in ("views", "shifts")},
                      "crop": {c: p2[c] for c in ("crop_cm", "upstream_frac", "frame_radius", "heat_input",
                                                  "mode_input", "mode_sigma")}},
            "ranker": {c: rc[c] for c in ("model", "pool_sources", "expand_fit", "dense_radius", "anchors",
                                          "feature_set", "nms_cm", "top_k_out", "calibrate")},
        },
        "ranker": {"fold": k, "kind": "LightGBM LambdaMART (lambdarank), Booster text", "n_features": len(names),
                   "feature_names": names, "trained_on": f"union rows of target events of folds != {k}, "
                                                         f"variant row_family({k}, fold)",
                   "refit_report": rep},
        "files": files,
        "source": {"run": str(run), **source},
        "config": {"name": cfg["name"], "inputs": cfg["inputs"], "pass1": p1, "pass2": p2, "ranker": rc,
                   "data": {"min_points": cfg["data"]["min_points"], "samples": cfg["data"]["samples"]}},
        "created": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "created_by": "dv2max_vertexer/make_bundle.py",
    }
    (out / "bundle.json").write_text(json.dumps(man, indent=1, default=str) + "\n")
    print(json.dumps({n: f for n, f in files.items()}, indent=1))


if __name__ == "__main__":
    main()
