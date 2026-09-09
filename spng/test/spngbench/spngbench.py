#!/usr/bin/env python3
"""
spngbench: SPNG-vs-OSP compute-performance benchmark driver.

This drives the factored benchmark config (spngbench.jsonnet) which runs exactly
one signal-processing chain (OSP or SPNG) from depos through drift + detsim to
signals.  Timing is harvested from the wire-cell log two ways:

  1. Log-timestamp *phase* intervals: startup, configuration, time-to-first-
     execute, execution, time-to-last-execute, finalize.  (See analyze_phases.)
  2. Per-node *Timer* accounting: wall-sec and core-sec for every DFP node,
     parsed from "Timer:" lines and rolled up into per-node mean+/-stdev plus
     forward-vs-rest SP category sums.  (See analyze_nodes / rollup_*.)

The first-stage report is a per-config JSON (schema "spngbench-config/1").  Two
such reports (OSP and SPNG at the same core/device point) can be combined into a
comparison JSON (schema "spngbench-compare/1").

Core-count controls:
  - Wire-Cell node parallelism: TbbFlow + TbbDataFlowGraph.max_threads=wc_cores.
  - Torch intra-op threads:      OMP_NUM_THREADS=torch_cores (env).
  - Device:                      -A device=cpu|gpu|gpu0|gpu1 (jsonnet TLA).

This module is import-friendly (functions have no CLI side effects) and also has
a small argparse CLI for single-config runs.  The grid scan lives in the driver
that imports this (see the epic spng-fra).
"""

import os
import re
import sys
import json
import time
import socket
import argparse
import subprocess
from pathlib import Path
from statistics import mean, stdev as _sstdev

HERE = Path(__file__).resolve().parent
JSONNET = HERE / "spngbench.jsonnet"
DEFAULT_MODEL = "/nfs/data/1/calcuttj/toolkit_testing/legacy_roiuniter_090326.ts"

# ---------------------------------------------------------------------------
# Node classification.  Rule-based and intentionally extensible: to add a new
# forward-inference node type or a new SP implementation, extend these sets /
# prefixes.  Everything not matched is "other" (simulation, drift, I/O, fanning).
# ---------------------------------------------------------------------------
FORWARD_CLASSES = {
    "WireCell::Pytorch::DNNROIFinding",   # OSP forward (DDNROI)
    "WireCell::SPNG::TensorForward",      # SPNG forward (ROIUNITER)
}
SP_OTHER_CLASSES = {
    "WireCell::SigProc::OmnibusSigProc",  # OSP decon+ROI (non-forward SP)
}
SP_OTHER_PREFIXES = (
    "WireCell::SPNG::",                   # all SPNG nodes except the forward one
)


def classify(cls):
    """Return the timing category of a node class: forward | sp_other | other."""
    if cls in FORWARD_CLASSES:
        return "forward"
    if cls in SP_OTHER_CLASSES:
        return "sp_other"
    for pre in SP_OTHER_PREFIXES:
        if cls.startswith(pre):
            return "sp_other"
    return "other"


# ---------------------------------------------------------------------------
# Log parsing: timestamps and phase markers.
# ---------------------------------------------------------------------------
# Timestamps look like "[18:26:27.076] ...".  spdlog occasionally drops the
# leading '[' under thread contention, so it is optional here.  Continuation
# lines of multi-line log entries carry no timestamp and are skipped.
_TS_RE = re.compile(r"^\[?(\d{2}):(\d{2}):(\d{2})\.(\d{3})\]")

# Phase markers.  Each is (key, [needles], keep_last); any needle matches.  The
# needles are chosen to be engine-agnostic: "executing app:" and the finalize/
# config lines come from the main logger under both Pgrapher and TbbFlow, while
# exec_end is Pgrapher's "graph execution complete" or TbbFlow's "totals:".
_PHASE_MARKERS = [
    ("config_start",    ["configuring component:"],   False),
    ("config_end",      ["configured component:"],     True),
    ("exec_start",      ["executing app:"],            False),
    ("exec_end",        ["graph execution complete", "totals: wall="], False),
    ("finalize_start",  ["finalizing component:"],     False),
]

# Pgrapher per-node "Timer:" payload (mirrors wirecell.util.logtimes; inlined so
# this driver is self-contained and needs no wire-cell-python import).
_TIMER_RE = re.compile(
    r'Timer:\s*'
    r'(?P<wall>[-+0-9.eE]+)\s*wall-sec,\s*'
    r'(?P<core>[-+0-9.eE]+)\s*core-sec:\s*'
    r'\((?P<cls>[^)]*)\)\s*'
    r'"(?P<inst>[^"]*)"'
)
_TIMER_TOTAL_RE = re.compile(
    r'Timer:\s*Total\s*'
    r'(?P<wall>[-+0-9.eE]+)\s*wall-sec,\s*'
    r'(?P<core>[-+0-9.eE]+)\s*core-sec'
)
# TbbFlow (TbbDataFlowGraph, summary>=1) per-node timing.  wall is total across
# calls in MILLISECONDS; core is in seconds; class is in [..] not (..).
_TBB_RE = re.compile(
    r'calls=(?P<calls>\d+)\s+time=(?P<wall_ms>[-+0-9.eE]+)\s+mean=\S+\s+max=\S+\s+'
    r'\[wall-ms\]\s+core=(?P<core>[-+0-9.eE]+)\s+\[s\]\s+'
    r'\[(?P<cls>[^\]]*)\]\s+"(?P<inst>[^"]*)"'
)
_TBB_TOTAL_RE = re.compile(
    r'totals:\s*wall=(?P<wall>[-+0-9.eE]+)\s*s,\s*core=(?P<core>[-+0-9.eE]+)\s*s'
)


def _parse_ts(line):
    m = _TS_RE.match(line)
    if not m:
        return None
    h, mi, s, ms = (int(x) for x in m.groups())
    return h * 3600 + mi * 60 + s + ms / 1000.0


def analyze_phases(lines):
    """
    Extract phase-boundary timestamps and derived interval durations.

    Robust to spdlog reordering: first/last are min/max over all timestamped
    lines, and each marker uses the timestamp on its own line (event time), not
    file position.  Returns None if no timestamps were found.
    """
    all_ts = []
    marker = {}
    for line in lines:
        ts = _parse_ts(line)
        if ts is None:
            continue
        all_ts.append(ts)
        for key, needles, keep_last in _PHASE_MARKERS:
            if (keep_last or key not in marker) and any(n in line for n in needles):
                marker[key] = ts
    if not all_ts:
        return None

    first = min(all_ts)
    last = max(all_ts)
    rel = {k: round(v - first, 3) for k, v in marker.items()}
    rel["first"] = 0.0
    rel["last"] = round(last - first, 3)

    def diff(a, b):
        if a in rel and b in rel:
            d = rel[b] - rel[a]
            return round(d, 3) if d >= 0 else round(d, 3)
        return None

    durations = {
        "startup":       diff("first", "config_start"),
        "config":        diff("config_start", "config_end"),
        "first_execute": diff("config_end", "exec_start"),
        "execution":     diff("exec_start", "exec_end"),
        "last_execute":  diff("exec_end", "finalize_start"),
        "finalize":      diff("finalize_start", "last"),
        "total":         diff("first", "last"),
    }
    return {"markers": rel, "durations": durations}


def analyze_nodes(lines):
    """Return (nodes, timer_total) parsed from Timer log lines.

    nodes: list of {instance, class, wall-sec, core-sec}.
    timer_total: {wall, core} from the 'Timer: Total' line, or None.
    """
    nodes = []
    total = None
    for line in lines:
        # Totals (either engine).
        mt = _TIMER_TOTAL_RE.search(line) or _TBB_TOTAL_RE.search(line)
        if mt:
            total = {"wall": float(mt.group("wall")), "core": float(mt.group("core"))}
            continue
        # Pgrapher per-node.
        m = _TIMER_RE.search(line)
        if m:
            nodes.append({
                "instance": m.group("inst"),
                "class": m.group("cls"),
                "wall-sec": float(m.group("wall")),
                "core-sec": float(m.group("core")),
            })
            continue
        # TbbFlow per-node (wall in ms -> s).
        m = _TBB_RE.search(line)
        if m:
            nodes.append({
                "instance": m.group("inst"),
                "class": m.group("cls"),
                "wall-sec": float(m.group("wall_ms")) / 1000.0,
                "core-sec": float(m.group("core")),
                "calls": int(m.group("calls")),
            })
    return nodes, total


def analyze_log(logfile):
    """Parse one wire-cell log into {phases, nodes, timer_total}."""
    with open(logfile, errors="replace") as fp:
        lines = fp.readlines()
    phases = analyze_phases(lines)
    nodes, total = analyze_nodes(lines)
    return {"phases": phases, "nodes": nodes, "timer_total": total}


# ---------------------------------------------------------------------------
# Running one wire-cell job.
# ---------------------------------------------------------------------------
_OOM_RE = re.compile(r"out of memory|CUDA error: out of memory|CUDA_ERROR_OUT_OF_MEMORY",
                     re.IGNORECASE)


def build_cmd(which, device, wc_cores, input, model_file, output, logfile,
              detname="pdhd", engine="TbbFlow", verbosity=0):
    return [
        "wire-cell",
        "-c", str(JSONNET),
        "-l", str(logfile),
        "-L", "debug",
        "-A", f"input={input}",
        "-A", f"model_file={model_file}",
        "-A", f"output={output}",
        "-A", f"which={which}",
        "-A", f"detname={detname}",
        "-A", f"device={device}",
        "-A", f"engine={engine}",
        "-A", f"wc_cores={wc_cores}",
        "-A", f"verbosity={verbosity}",
    ]


def run_wirecell(which, device, wc_cores, torch_cores, input, outdir,
                 model_file=DEFAULT_MODEL, detname="pdhd", engine="TbbFlow",
                 tag="", verbosity=0, extra_env=None, dry_run=False):
    """
    Run one wire-cell job and classify the outcome.

    Returns a dict describing the run; the log is left on disk for analyze_log.
    outcome is one of: ok | oom | error | skipped(dry).
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    stem = tag or f"{which}-{device}-wc{wc_cores}-omp{torch_cores}"
    logfile = outdir / f"{stem}.log"
    output = outdir / f"{stem}.npz"

    cmd = build_cmd(which, device, wc_cores, input, model_file, output, logfile,
                    detname=detname, engine=engine, verbosity=verbosity)

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(torch_cores)
    if extra_env:
        env.update(extra_env)

    result = {
        "which": which, "device": device, "wc_cores": wc_cores,
        "torch_cores": torch_cores, "engine": engine, "input": str(input),
        "detname": detname, "logfile": str(logfile), "output": str(output),
        "cmd": cmd, "stem": stem,
    }

    if dry_run:
        result["outcome"] = "skipped"
        result["returncode"] = None
        result["wall_clock"] = 0.0
        return result

    t0 = time.time()
    proc = subprocess.run(cmd, env=env, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT)
    result["wall_clock"] = round(time.time() - t0, 3)
    result["returncode"] = proc.returncode
    captured = proc.stdout.decode(errors="replace") if proc.stdout else ""

    # Scan both captured stdio and the log file for OOM markers.
    oom = bool(_OOM_RE.search(captured))
    if not oom and logfile.exists():
        try:
            with open(logfile, errors="replace") as fp:
                oom = bool(_OOM_RE.search(fp.read()))
        except OSError:
            pass

    if proc.returncode == 0:
        result["outcome"] = "ok"
    elif oom:
        result["outcome"] = "oom"
    else:
        result["outcome"] = "error"
    if oom:
        result["oom"] = True
    # Keep a short tail of captured output for post-mortem on failures.
    if result["outcome"] != "ok":
        result["stdio_tail"] = "\n".join(captured.splitlines()[-30:])
    return result


# ---------------------------------------------------------------------------
# Rollups.
# ---------------------------------------------------------------------------
def _stat(values):
    n = len(values)
    if n == 0:
        return {"mean": None, "stdev": None, "n": 0, "min": None, "max": None}
    m = mean(values)
    sd = _sstdev(values) if n > 1 else 0.0
    return {
        "mean": round(m, 4), "stdev": round(sd, 4), "n": n,
        "min": round(min(values), 4), "max": round(max(values), 4),
    }


def rollup_nodes(analyses):
    """Per-node mean+/-stdev over a list of analyze_log() results.

    Keyed by "class\\tinstance" because instance names (e.g. "tpc0") are reused
    across classes.
    """
    acc = {}
    for a in analyses:
        for nd in a["nodes"]:
            key = f"{nd['class']}\t{nd['instance']}"
            slot = acc.setdefault(key, {"class": nd["class"], "instance": nd["instance"],
                                        "wall": [], "core": []})
            slot["wall"].append(nd["wall-sec"])
            slot["core"].append(nd["core-sec"])
    out = {}
    for key, slot in acc.items():
        out[key] = {
            "class": slot["class"],
            "instance": slot["instance"],
            "category": classify(slot["class"]),
            "wall": _stat(slot["wall"]),
            "core": _stat(slot["core"]),
        }
    return out


def _cat_sums(nodes):
    sums = {c: {"wall": 0.0, "core": 0.0} for c in ("forward", "sp_other", "other")}
    for nd in nodes:
        c = classify(nd["class"])
        sums[c]["wall"] += nd["wall-sec"]
        sums[c]["core"] += nd["core-sec"]
    return sums


def rollup_categories(analyses):
    """Sum forward / sp_other / other per run, then mean+/-stdev across runs.

    Also derive sp_total (forward + sp_other) and all_total (from Timer: Total).
    """
    per = [_cat_sums(a["nodes"]) for a in analyses]
    out = {}
    for cat in ("forward", "sp_other", "other"):
        out[cat] = {
            "wall": _stat([p[cat]["wall"] for p in per]),
            "core": _stat([p[cat]["core"] for p in per]),
        }
    out["sp_total"] = {
        "wall": _stat([p["forward"]["wall"] + p["sp_other"]["wall"] for p in per]),
        "core": _stat([p["forward"]["core"] + p["sp_other"]["core"] for p in per]),
    }
    tot_wall = [a["timer_total"]["wall"] for a in analyses if a.get("timer_total")]
    tot_core = [a["timer_total"]["core"] for a in analyses if a.get("timer_total")]
    out["all_total"] = {"wall": _stat(tot_wall), "core": _stat(tot_core)}
    return out


def rollup_phases(analyses):
    keys = ["startup", "config", "first_execute", "execution",
            "last_execute", "finalize", "total"]
    out = {}
    for k in keys:
        vals = [a["phases"]["durations"][k]
                for a in analyses
                if a.get("phases") and a["phases"]["durations"].get(k) is not None]
        out[k] = _stat(vals) if vals else None
    return out


def _ngpu():
    try:
        out = subprocess.run(["nvidia-smi", "-L"], stdout=subprocess.PIPE,
                             stderr=subprocess.DEVNULL, timeout=10)
        return len([l for l in out.stdout.decode().splitlines() if l.strip()])
    except Exception:
        return 0


def make_report(which, device, wc_cores, torch_cores, runs, engine="TbbFlow",
                detname="pdhd"):
    """Build a per-config report (schema spngbench-config/1) from run dicts.

    Each run dict is a run_wirecell() result; ok runs are analyzed and rolled up.
    """
    ok = [r for r in runs if r.get("outcome") == "ok"]
    analyses = [analyze_log(r["logfile"]) for r in ok]

    if any(r.get("outcome") == "oom" for r in runs):
        overall = "oom" if not ok else "partial-oom"
    elif ok and len(ok) == len(runs):
        overall = "ok"
    elif ok:
        overall = "partial"
    else:
        overall = "error"

    report = {
        "schema": "spngbench-config/1",
        "meta": {
            "which": which, "device": device,
            "wc_cores": wc_cores, "torch_cores": torch_cores,
            "engine": engine, "detname": detname,
            "host": socket.gethostname(), "ngpu": _ngpu(),
            "nrun": len(runs), "nok": len(ok),
            "created": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "cmd": runs[0]["cmd"] if runs else None,
        },
        "outcome": overall,
        "runs": [{
            "input": r.get("input"), "outcome": r.get("outcome"),
            "returncode": r.get("returncode"), "wall_clock": r.get("wall_clock"),
            "logfile": r.get("logfile"), "stem": r.get("stem"),
            **({"stdio_tail": r["stdio_tail"]} if r.get("stdio_tail") else {}),
        } for r in runs],
        "phases": rollup_phases(analyses) if analyses else None,
        "summary": {
            "per_node": rollup_nodes(analyses),
            "categories": rollup_categories(analyses),
        } if analyses else None,
    }
    return report


def bench_config(which, device, wc_cores, torch_cores, inputs, outdir,
                 repeats=1, model_file=DEFAULT_MODEL, detname="pdhd",
                 engine="TbbFlow", verbosity=0, dry_run=False):
    """Run one (which, device, wc_cores, torch_cores) config over inputs x repeats.

    Writes report-<stem>.json into outdir and returns the report dict.
    """
    if isinstance(inputs, (str, Path)):
        inputs = [inputs]
    runs = []
    for ii, inp in enumerate(inputs):
        for rr in range(repeats):
            suffix = ""
            if len(inputs) > 1:
                suffix += f"-in{ii}"
            if repeats > 1:
                suffix += f"-rep{rr}"
            tag = f"{which}-{device}-wc{wc_cores}-omp{torch_cores}{suffix}"
            runs.append(run_wirecell(
                which, device, wc_cores, torch_cores, inp, outdir,
                model_file=model_file, detname=detname, engine=engine,
                tag=tag, verbosity=verbosity, dry_run=dry_run))

    report = make_report(which, device, wc_cores, torch_cores, runs,
                         engine=engine, detname=detname)
    stem = f"{which}-{device}-wc{wc_cores}-omp{torch_cores}"
    rpath = Path(outdir) / f"report-{stem}.json"
    rpath.write_text(json.dumps(report, indent=2))
    report["_path"] = str(rpath)
    return report


def combine(osp_report, spng_report, outdir=None, stem=None):
    """Combine an OSP and a SPNG per-config report into a comparison JSON.

    schema spngbench-compare/1.
    """
    def cat(rep, name, field="wall"):
        try:
            return rep["summary"]["categories"][name][field]["mean"]
        except (TypeError, KeyError):
            return None

    def ratio(a, b):
        if a and b:
            return round(a / b, 3)
        return None

    comparison = {}
    for name in ("forward", "sp_other", "sp_total", "all_total"):
        o = cat(osp_report, name)
        s = cat(spng_report, name)
        comparison[name] = {
            "osp_wall_mean": o, "spng_wall_mean": s,
            "spng_over_osp": ratio(s, o),
        }

    meta_src = osp_report["meta"]
    combined = {
        "schema": "spngbench-compare/1",
        "meta": {
            "device": meta_src["device"],
            "wc_cores": meta_src["wc_cores"],
            "torch_cores": meta_src["torch_cores"],
            "detname": meta_src["detname"],
            "host": meta_src["host"], "ngpu": meta_src["ngpu"],
            "created": time.strftime("%Y-%m-%dT%H:%M:%S"),
        },
        "outcome": {"osp": osp_report["outcome"], "spng": spng_report["outcome"]},
        "comparison": comparison,
        "osp": osp_report,
        "spng": spng_report,
    }
    if outdir:
        stem = stem or f"{meta_src['device']}-wc{meta_src['wc_cores']}-omp{meta_src['torch_cores']}"
        cpath = Path(outdir) / f"compare-{stem}.json"
        cpath.write_text(json.dumps(combined, indent=2))
        combined["_path"] = str(cpath)
    return combined


# ---------------------------------------------------------------------------
# CLI (single-config; grid scan lives in the higher-level driver).
# ---------------------------------------------------------------------------
def _add_common(p):
    p.add_argument("--device", default="cpu", help="cpu|gpu|gpu0|gpu1 (default cpu)")
    p.add_argument("--wc-cores", type=int, default=1, help="Wire-Cell TBB max_threads")
    p.add_argument("--torch-cores", type=int, default=1, help="OMP_NUM_THREADS for torch")
    p.add_argument("--input", default=None, help="depo .npz (repeatable)", action="append")
    p.add_argument("--outdir", default="spngbench-pdhd", help="work/output directory")
    p.add_argument("--repeats", type=int, default=1)
    p.add_argument("--model-file", default=DEFAULT_MODEL)
    p.add_argument("--detname", default="pdhd")
    p.add_argument("--engine", default="TbbFlow", help="TbbFlow (parallel) or Pgrapher (serial)")
    p.add_argument("--dry-run", action="store_true")


def _default_inputs(args):
    if args.input:
        return args.input
    # Fall back to the standard muon depos relative to the toolkit checkout.
    # HERE = <toolkit>/spng/test/spngbench, so parents[2] = <toolkit>.
    cand = HERE.parents[2] / "test" / "data" / "muon-depos.npz"
    return [str(cand)]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    pr = sub.add_parser("run", help="run one SP chain (osp or spng) at one config")
    pr.add_argument("which", choices=["osp", "spng"])
    _add_common(pr)

    pp = sub.add_parser("pair", help="run both osp and spng at one config and combine")
    _add_common(pp)

    pa = sub.add_parser("analyze", help="(re)analyze an existing log into a report fragment")
    pa.add_argument("logfile")

    args = ap.parse_args(argv)

    if args.cmd == "analyze":
        print(json.dumps(analyze_log(args.logfile), indent=2))
        return 0

    inputs = _default_inputs(args)

    if args.cmd == "run":
        rep = bench_config(args.which, args.device, args.wc_cores, args.torch_cores,
                           inputs, args.outdir, repeats=args.repeats,
                           model_file=args.model_file, detname=args.detname,
                           engine=args.engine, dry_run=args.dry_run)
        print(f"outcome={rep['outcome']} wrote {rep['_path']}")
        return 0

    if args.cmd == "pair":
        osp = bench_config("osp", args.device, args.wc_cores, args.torch_cores,
                           inputs, args.outdir, repeats=args.repeats,
                           model_file=args.model_file, detname=args.detname,
                           engine=args.engine, dry_run=args.dry_run)
        spng = bench_config("spng", args.device, args.wc_cores, args.torch_cores,
                            inputs, args.outdir, repeats=args.repeats,
                            model_file=args.model_file, detname=args.detname,
                            engine=args.engine, dry_run=args.dry_run)
        comb = combine(osp, spng, outdir=args.outdir)
        print(f"osp={osp['outcome']} spng={spng['outcome']} wrote {comb['_path']}")
        return 0

    return 1


if __name__ == "__main__":
    sys.exit(main())
