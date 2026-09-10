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
import shutil
import argparse
import subprocess
from pathlib import Path
from statistics import mean, stdev as _sstdev

HERE = Path(__file__).resolve().parent
JSONNET = HERE / "spngbench.jsonnet"
DEFAULT_MODEL = "/nfs/data/1/calcuttj/toolkit_testing/legacy_roiuniter_090326.ts"


def log(msg):
    """Progress message to stderr (stdout stays reserved for JSON output)."""
    print(msg, file=sys.stderr, flush=True)

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


def build_cmd(stage, device, wc_cores, input, model_file, output, logfile,
              detname="pdhd", engine="TbbFlow", gpu_scheme="none", ngpu=1, napa=1,
              verbosity=0):
    return [
        "wire-cell",
        "-c", str(JSONNET),
        "-l", str(logfile),
        "-L", "debug",
        "-A", f"input={input}",
        "-A", f"model_file={model_file}",
        "-A", f"output={output}",
        "-A", f"stage={stage}",
        "-A", f"detname={detname}",
        "-A", f"device={device}",
        "-A", f"gpu_scheme={gpu_scheme}",
        "-A", f"ngpu={ngpu}",
        "-A", f"napa={napa}",
        "-A", f"engine={engine}",
        "-A", f"wc_cores={wc_cores}",
        "-A", f"verbosity={verbosity}",
    ]


def run_wirecell(stage, device, wc_cores, torch_cores, input, outdir,
                 model_file=DEFAULT_MODEL, detname="pdhd", engine="TbbFlow",
                 gpu_scheme="none", ngpu=1, napa=1, tag="", output=None,
                 verbosity=0, extra_env=None, dry_run=False):
    """
    Run one wire-cell job (stage=sim|osp|spng) and classify the outcome.

    Returns a dict describing the run; the log is left on disk for analyze_log.
    outcome is one of: ok | oom | error | skipped(dry).
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    stem = tag or f"{stage}-{device}-wc{wc_cores}-omp{torch_cores}"
    logfile = outdir / f"{stem}.log"
    output = Path(output) if output else outdir / f"{stem}.npz"

    cmd = build_cmd(stage, device, wc_cores, input, model_file, output, logfile,
                    detname=detname, engine=engine, gpu_scheme=gpu_scheme,
                    ngpu=ngpu, napa=napa, verbosity=verbosity)

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(torch_cores)
    if extra_env:
        env.update(extra_env)

    result = {
        "stage": stage, "device": device, "wc_cores": wc_cores,
        "torch_cores": torch_cores, "engine": engine, "input": str(input),
        "gpu_scheme": gpu_scheme, "ngpu": ngpu, "napa": napa,
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


def config_stem(stage, device, wc_cores, torch_cores, gpu_scheme="none", ngpu=1):
    """Filename stem for one config; encodes the GPU scheme when sharded."""
    s = f"{stage}-{device}-wc{wc_cores}-omp{torch_cores}"
    if gpu_scheme and gpu_scheme != "none":
        s += f"-{gpu_scheme}{ngpu}"
    return s


def make_report(stage, device, wc_cores, torch_cores, runs, engine="TbbFlow",
                detname="pdhd", gpu_scheme="none", ngpu=1, napa=1):
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
            "stage": stage, "device": device,
            "wc_cores": wc_cores, "torch_cores": torch_cores,
            "gpu_scheme": gpu_scheme, "ngpu": ngpu, "napa": napa,
            "engine": engine, "detname": detname,
            "host": socket.gethostname(), "host_ngpu": _ngpu(),
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


def ensure_adc(depo_input, outdir, model_file=DEFAULT_MODEL, detname="pdhd",
               engine="Pgrapher", napa=1, dry_run=False):
    """Run the sim stage (depos -> ADC frame file(s)) once for a depo input, cached.

    Returns the ADC file base path.  For napa>1 the sim writes one ADC file per
    APA (base-tpc<N>.npz); the base path is what OSP/SPNG are given and the
    jsonnet re-derives the per-APA names.  The sim is CPU-only and shared by OSP
    and SPNG, so it runs once per depo input and is reused across the grid.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    stem = "adc-" + Path(depo_input).stem
    adc = outdir / f"{stem}.npz"
    # Existence probe: base name for napa==1, first per-APA file for napa>1.
    probe = adc if napa == 1 else outdir / f"{stem}-tpc0.npz"
    if probe.exists() and probe.stat().st_size > 0:
        return str(adc)
    if dry_run:
        return str(adc)
    r = run_wirecell("sim", "cpu", 1, 1, depo_input, outdir,
                     model_file=model_file, detname=detname, engine=engine,
                     napa=napa, tag=stem, output=str(adc))
    if r["outcome"] != "ok":
        log(f"WARNING: sim failed for {depo_input}: {r['outcome']}")
    return str(adc)


def bench_config(stage, device, wc_cores, torch_cores, inputs, outdir,
                 repeats=1, model_file=DEFAULT_MODEL, detname="pdhd",
                 engine="TbbFlow", gpu_scheme="none", ngpu=1, napa=1,
                 verbosity=0, dry_run=False):
    """Run one config (stage osp|spng) over inputs x repeats.

    `inputs` are ADC frame files (already simmed).  Writes report-<stem>.json
    and returns the report dict.
    """
    if isinstance(inputs, (str, Path)):
        inputs = [inputs]
    base = config_stem(stage, device, wc_cores, torch_cores, gpu_scheme, ngpu)
    runs = []
    for ii, inp in enumerate(inputs):
        for rr in range(repeats):
            suffix = ""
            if len(inputs) > 1:
                suffix += f"-in{ii}"
            if repeats > 1:
                suffix += f"-rep{rr}"
            runs.append(run_wirecell(
                stage, device, wc_cores, torch_cores, inp, outdir,
                model_file=model_file, detname=detname, engine=engine,
                gpu_scheme=gpu_scheme, ngpu=ngpu, napa=napa,
                tag=base + suffix, verbosity=verbosity, dry_run=dry_run))

    report = make_report(stage, device, wc_cores, torch_cores, runs,
                         engine=engine, detname=detname,
                         gpu_scheme=gpu_scheme, ngpu=ngpu, napa=napa)
    rpath = Path(outdir) / f"report-{base}.json"
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

    m = spng_report["meta"]     # SPNG carries the GPU scheme; OSP is single-device
    combined = {
        "schema": "spngbench-compare/1",
        "meta": {
            "device": osp_report["meta"]["device"],
            "wc_cores": m["wc_cores"],
            "torch_cores": m["torch_cores"],
            "gpu_scheme": m.get("gpu_scheme", "none"), "ngpu": m.get("ngpu", 1),
            "detname": m["detname"],
            "host": m["host"], "host_ngpu": m.get("host_ngpu"),
            "created": time.strftime("%Y-%m-%dT%H:%M:%S"),
        },
        "outcome": {"osp": osp_report["outcome"], "spng": spng_report["outcome"]},
        "comparison": comparison,
        "osp": osp_report,
        "spng": spng_report,
    }
    if outdir:
        scheme_tag = ("-" + m["gpu_scheme"] + str(m["ngpu"])) if m.get("gpu_scheme", "none") != "none" else ""
        stem = stem or f"{m['device']}-wc{m['wc_cores']}-omp{m['torch_cores']}{scheme_tag}"
        cpath = Path(outdir) / f"compare-{stem}.json"
        cpath.write_text(json.dumps(combined, indent=2))
        combined["_path"] = str(cpath)
    return combined


# ---------------------------------------------------------------------------
# Configuration and grid enumeration.
# ---------------------------------------------------------------------------
def default_config():
    """The default benchmark configuration.

    Everything the grid needs lives here so a committed YAML/JSON file can
    reproduce or extend a run.  Designed to grow: add inputs, bump core/GPU
    ranges, or add detectors without touching the driver code.
    """
    return {
        "detname": "pdhd",
        "model_file": DEFAULT_MODEL,
        "inputs": None,                 # None -> the standard muon-depos.npz
        "outdir": "spngbench-pdhd",
        "repeats": 1,
        # Number of APAs (per-APA pipelines) each job covers.  Fixed for a grid
        # run (not a scanned axis); clamped to the detector's physical APA count.
        # More APAs give more independent pipelines that multiple wire-cell cores
        # can run concurrently, so the wc_cores range is capped at napa.
        "napa": 1,
        "which": ["osp", "spng"],
        "engine": "TbbFlow",
        "hw": {"max_hyperthreads": 64},  # cap: skip cpu cells with wc*torch > this
        "cpu_grid": {
            "wc_cores": [1, 2, 3, 4],
            "torch_cores": [1, 2, 4, 8, 16, 32],
        },
        "gpu_grid": {
            "device": "gpu",
            "wc_cores": [1, 2, 3, 4],
            "torch_cores": 1,           # torch always 1 core on GPU
        },
        "shard_grid": {                 # multi-GPU sharding of the SPNG graph
            "device": "gpu0",           # base device (OSP + unmatched nodes)
            "schemes": ["transverse", "longitudinal"],
            "ngpu": [2],                # GPU counts to try
            "wc_cores": [4],
            "torch_cores": 1,
        },
    }


# Physical APA counts for known detectors (used to clamp napa and the wc_cores
# range; the jsonnet clamps the actual pipeline count regardless).
DETECTOR_MAX_APA = {"pdhd": 4, "pdsp": 6, "pdvd": 8, "fdhd": 150, "fdvd": 160}


def clamp_napa(cfg):
    """Clamp cfg['napa'] to the detector's physical APA count (best effort)."""
    napa = max(1, int(cfg.get("napa", 1)))
    phys = DETECTOR_MAX_APA.get(cfg.get("detname", "pdhd"))
    return min(napa, phys) if phys else napa


def _merge(base, over):
    out = dict(base)
    for k, v in (over or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _merge(out[k], v)
        else:
            out[k] = v
    return out


def load_config(path):
    """Load a JSON or YAML config file, merged over default_config()."""
    text = Path(path).read_text()
    if str(path).endswith((".yaml", ".yml")):
        import yaml                     # optional dependency
        over = yaml.safe_load(text)
    else:
        over = json.loads(text)
    return _merge(default_config(), over)


def _cell(device, wc, torch, gpu_scheme="none", ngpu=1):
    return {"device": device, "wc_cores": wc, "torch_cores": torch,
            "gpu_scheme": gpu_scheme, "ngpu": ngpu}


def _wc_range(cfg, wc_sel, default):
    """Wire-cell core counts to scan, capped at napa (independent APA pipelines)."""
    napa = clamp_napa(cfg)
    return [w for w in (wc_sel or default) if w <= napa]


def cpu_cells(cfg, wc_sel=None, torch_sel=None):
    """Enumerate CPU grid cells (device=cpu), skipping wc*torch > max_hyperthreads."""
    g = cfg["cpu_grid"]
    cap = cfg["hw"]["max_hyperthreads"]
    wcs = _wc_range(cfg, wc_sel, g["wc_cores"])
    tcs = torch_sel or g["torch_cores"]
    return [_cell("cpu", wc, tc) for wc in wcs for tc in tcs if wc * tc <= cap]


def gpu_cells(cfg, wc_sel=None):
    """Enumerate single-GPU grid cells (all GPU-capable nodes on device, torch=1)."""
    g = cfg["gpu_grid"]
    wcs = _wc_range(cfg, wc_sel, g["wc_cores"])
    return [_cell(g.get("device", "gpu"), wc, g.get("torch_cores", 1)) for wc in wcs]


def shard_cells(cfg, wc_sel=None):
    """Enumerate multi-GPU sharded cells: (scheme x ngpu x wc_cores)."""
    g = cfg["shard_grid"]
    wcs = _wc_range(cfg, wc_sel, g["wc_cores"])
    return [_cell(g.get("device", "gpu0"), wc, g.get("torch_cores", 1),
                  gpu_scheme=scheme, ngpu=ng)
            for scheme in g["schemes"] for ng in g["ngpu"] for wc in wcs]


def select_cells(cfg, mode, wc_sel=None, torch_sel=None):
    """Return the list of cells for the selected mode: cpu|gpu|shard|all."""
    cells = []
    if mode in ("cpu", "all"):
        cells += cpu_cells(cfg, wc_sel, torch_sel)
    if mode in ("gpu", "all"):
        cells += gpu_cells(cfg, wc_sel)
    if mode in ("shard", "two-gpu", "all"):
        cells += shard_cells(cfg, wc_sel)
    return cells


# ---------------------------------------------------------------------------
# Grid execution.
# ---------------------------------------------------------------------------
def _cat_mean(rep, name, field="wall"):
    try:
        return rep["summary"]["categories"][name][field]["mean"]
    except (TypeError, KeyError):
        return None


def _cell_summary(cell, reports, compare):
    """Compact per-cell summary for the grid index."""
    s = {k: cell[k] for k in ("device", "wc_cores", "torch_cores", "gpu_scheme", "ngpu")}
    for w, rep in reports.items():
        s[w] = {
            "outcome": rep["outcome"],
            "forward_wall": _cat_mean(rep, "forward"),
            "sp_rest_wall": _cat_mean(rep, "sp_other"),
            "sp_total_wall": _cat_mean(rep, "sp_total"),
            "all_total_wall": _cat_mean(rep, "all_total"),
            "report_path": rep.get("_path"),
        }
    if compare:
        s["compare_path"] = compare.get("_path")
        s["spng_over_osp_sp_total"] = compare["comparison"]["sp_total"]["spng_over_osp"]
    return s


def run_grid(cfg, cells, which_list, dry_run=False):
    """Run each grid cell (osp/spng as selected), writing per-cell compare JSONs
    and an incrementally-updated grid index.  GPU OOM and missing GPUs are
    recorded, never fatal.
    """
    outdir = Path(cfg["outdir"])
    outdir.mkdir(parents=True, exist_ok=True)
    depos = cfg["inputs"] or [str(HERE.parents[2] / "test" / "data" / "muon-depos.npz")]
    host_ngpu = _ngpu()
    napa = clamp_napa(cfg)

    # Run the sim once per depo input; OSP and SPNG share the ADC frame file(s).
    log(f"sim: producing ADC frames for {len(depos)} depo input(s), napa={napa}")
    adc_inputs = [ensure_adc(d, str(outdir), model_file=cfg["model_file"],
                             detname=cfg["detname"], napa=napa, dry_run=dry_run)
                  for d in depos]

    index = {
        "schema": "spngbench-grid/1",
        "meta": {
            "detname": cfg["detname"], "host": socket.gethostname(), "host_ngpu": host_ngpu,
            "engine": cfg["engine"], "repeats": cfg["repeats"], "napa": napa,
            "depos": depos, "adc_inputs": adc_inputs, "which": which_list,
            "created": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "ncells": len(cells),
        },
        "cells": [],
    }
    index_path = outdir / "grid-index.json"

    def flush():
        index_path.write_text(json.dumps(index, indent=2))

    for i, cell in enumerate(cells):
        device = cell["device"]
        scheme, ngpu = cell.get("gpu_scheme", "none"), cell.get("ngpu", 1)
        stag = f"-{scheme}{ngpu}" if scheme != "none" else ""
        log(f"[cell {i+1}/{len(cells)}] {device} wc{cell['wc_cores']} "
            f"omp{cell['torch_cores']}{stag} which={which_list}")

        # A sharded cell needing more GPUs than present is recorded, not run.
        needs = ngpu if scheme != "none" else (1 if device.startswith("gpu") else 0)
        if needs > host_ngpu:
            index["cells"].append({**{k: cell[k] for k in
                                      ("device", "wc_cores", "torch_cores", "gpu_scheme", "ngpu")},
                                   "skipped": True,
                                   "reason": f"needs {needs} GPUs, have {host_ngpu}"})
            flush()
            continue

        reports = {}
        for w in which_list:
            # OSP is a single-device chain; only SPNG carries the GPU scheme.
            w_scheme = scheme if w == "spng" else "none"
            reports[w] = bench_config(
                w, device, cell["wc_cores"], cell["torch_cores"], adc_inputs, str(outdir),
                repeats=cfg["repeats"], model_file=cfg["model_file"],
                detname=cfg["detname"], engine=cfg["engine"],
                gpu_scheme=w_scheme, ngpu=ngpu, napa=napa, dry_run=dry_run)
        compare = None
        if "osp" in reports and "spng" in reports:
            compare = combine(reports["osp"], reports["spng"], outdir=str(outdir))
        index["cells"].append(_cell_summary(cell, reports, compare))
        flush()

    log(f"wrote grid index {index_path} ({len(index['cells'])} cells)")
    index["_path"] = str(index_path)
    return index


# ---------------------------------------------------------------------------
# Human-readable reporting (Markdown / HTML / LaTeX from a benchmark JSON).
# ---------------------------------------------------------------------------
# Report categories, finer than the rollup's forward/sp_other/other: separate
# I/O from compute, and within compute separate DNN (the neural-net forward)
# from the rest of signal processing.
_IO_HINTS = ("::Sio::", "FrameFile", "DepoFile", "TensorFile")
CAT_ORDER = ["dnn", "sp_nondnn", "io", "other"]
CAT_LABEL = {"dnn": "DNN forward", "sp_nondnn": "SP (non-DNN)",
             "io": "I/O", "other": "other"}
CAT_COLOR = {"dnn": "#d95f5f", "sp_nondnn": "#5f8fbf", "io": "#7fb37f",
             "other": "#c8c8c8"}


def report_category(cls):
    if cls in FORWARD_CLASSES:
        return "dnn"
    if any(h in cls for h in _IO_HINTS):
        return "io"
    if classify(cls) == "sp_other":
        return "sp_nondnn"
    return "other"


def _short(cls):
    return cls.split("::")[-1]


def per_node_list(report):
    """Flat list of timed nodes from a report's per-node rollup."""
    pn = (report.get("summary") or {}).get("per_node") or {}
    out = []
    for v in pn.values():
        out.append({
            "class": v["class"], "short": _short(v["class"]),
            "instance": v["instance"], "cat": report_category(v["class"]),
            "wall": v["wall"]["mean"] or 0.0, "wall_sd": v["wall"]["stdev"] or 0.0,
            "core": v["core"]["mean"] or 0.0,
        })
    return sorted(out, key=lambda n: -n["wall"])


def cat_sums(nodes):
    d = {c: 0.0 for c in CAT_ORDER}
    for n in nodes:
        d[n["cat"]] += n["wall"]
    return d


def load_bench(path):
    """Load a benchmark JSON and return {'meta':..., 'reports': {stage: report}}.

    Accepts compare/1 (osp+spng), config/1 (one stage), or grid/1 (uses the
    first cell that has a compare file).
    """
    j = json.loads(Path(path).read_text())
    schema = j.get("schema", "")
    if schema.startswith("spngbench-compare"):
        return {"meta": j["meta"], "reports": {"osp": j["osp"], "spng": j["spng"]}}
    if schema.startswith("spngbench-config"):
        return {"meta": j["meta"], "reports": {j["meta"]["stage"]: j}}
    if schema.startswith("spngbench-grid"):
        base = Path(path).parent
        for c in j.get("cells", []):
            cp = c.get("compare_path")
            if cp and Path(cp).exists():
                return load_bench(cp)
            if cp and (base / Path(cp).name).exists():
                return load_bench(base / Path(cp).name)
        raise SystemExit("grid index has no usable compare file; point report at a compare-*.json")
    raise SystemExit(f"unrecognized benchmark JSON schema: {schema!r}")


def _fmt(x, nd=2):
    return "-" if x is None else f"{x:.{nd}f}"


# ---- figures (PNG, shared by all three output formats) ----
def make_figures(reports, figdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figdir = Path(figdir)
    figdir.mkdir(parents=True, exist_ok=True)
    figs = {}
    stages = [s for s in ("osp", "spng") if s in reports]

    # (1) Stacked category composition, OSP vs SPNG.
    sums = {s: cat_sums(per_node_list(reports[s])) for s in stages}
    fig, ax = plt.subplots(figsize=(5, 4))
    bottoms = {s: 0.0 for s in stages}
    x = range(len(stages))
    for cat in CAT_ORDER:
        vals = [sums[s][cat] for s in stages]
        ax.bar(x, vals, bottom=[bottoms[s] for s in stages],
               color=CAT_COLOR[cat], label=CAT_LABEL[cat], width=0.6)
        for s in stages:
            bottoms[s] += sums[s][cat]
    ax.set_xticks(list(x))
    ax.set_xticklabels([s.upper() for s in stages])
    ax.set_ylabel("wall time [s]")
    ax.set_title("Time composition by category")
    ax.legend(fontsize=8)
    fig.tight_layout()
    p = figdir / "categories.png"
    fig.savefig(p, dpi=120)
    plt.close(fig)
    figs["categories"] = p.name

    # (2) Per-stage top-node horizontal bars, colored by category.
    for s in stages:
        nodes = per_node_list(reports[s])[:14]
        if not nodes:
            continue
        fig, ax = plt.subplots(figsize=(6.5, max(2.5, 0.36 * len(nodes))))
        y = range(len(nodes))
        ax.barh(list(y), [n["wall"] for n in nodes],
                color=[CAT_COLOR[n["cat"]] for n in nodes],
                xerr=[n["wall_sd"] for n in nodes], error_kw=dict(lw=0.6))
        ax.set_yticks(list(y))
        ax.set_yticklabels([f'{n["short"]}:{n["instance"]}'[:34] for n in nodes], fontsize=7)
        ax.invert_yaxis()
        ax.set_xlabel("wall time [s]")
        ax.set_title(f"{s.upper()} — top nodes by wall time")
        fig.tight_layout()
        p = figdir / f"nodes_{s}.png"
        fig.savefig(p, dpi=120)
        plt.close(fig)
        figs[f"nodes_{s}"] = p.name

    return figs


# ---- annotated flow graphs (GraphViz dot -> PNG) ----
def _tlas_from_cmd(cmd):
    jsonnet, tlas = None, {}
    i = 0
    while i < len(cmd):
        if cmd[i] == "-c" and i + 1 < len(cmd):
            jsonnet = cmd[i + 1]; i += 2; continue
        if cmd[i] == "-A" and i + 1 < len(cmd):
            k, _, v = cmd[i + 1].partition("="); tlas[k] = v; i += 2; continue
        i += 1
    return jsonnet, tlas


def _dot_node_ids(dot_text):
    """Collect node identifiers referenced by edges in a dot file."""
    ids = set()
    for line in dot_text.splitlines():
        if "->" not in line:
            continue
        for side in line.split("->"):
            side = side.strip()
            m = re.match(r'\s*("(?:[^"]*)"|[A-Za-z0-9_.]+)', side)
            if m:
                ids.add(m.group(1).strip('"'))
    return ids


def render_flow_graph(report, out_png, device_hint=None):
    """Render this report's flow graph annotated with per-node wall/core time.

    Runs `wcpy pgraph dotify` for the topology, appends node statements with
    timing labels/colors, and renders with `dot`.  Returns the PNG basename or
    None if any tool is missing.
    """
    if not shutil.which("dot") or not shutil.which("wcpy"):
        return None
    cmd = report.get("meta", {}).get("cmd")
    if not cmd:
        return None
    jsonnet, tlas = _tlas_from_cmd(cmd)
    if not jsonnet:
        return None

    out_png = Path(out_png)
    base_dot = out_png.with_suffix(".base.dot")
    ann_dot = out_png.with_suffix(".dot")

    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = ""
    env.setdefault("MPLCONFIGDIR", "/tmp/mpl-spngbench")
    dargs = ["wcpy", "pgraph", "dotify", "--no-services", "--no-params"]
    for k, v in tlas.items():
        dargs += ["-A", f"{k}={v}"]
    dargs += [jsonnet, str(base_dot)]
    try:
        r = subprocess.run(dargs, env=env, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, timeout=120)
        if r.returncode != 0 or not base_dot.exists():
            return None
    except Exception:
        return None

    dot_text = base_dot.read_text()
    ids = _dot_node_ids(dot_text)
    # Map node timing by dot id "<Type>_<instance>".  dotify uses the factory
    # type name, which for SPNG nodes carries an "SPNG" prefix the demangled
    # Timer class short lacks (e.g. dot "SPNGKernelConvolve" vs "KernelConvolve"),
    # so index by both spellings.
    lut = {}
    for n in per_node_list(report):
        lut[f'{n["short"]}_{n["instance"]}'] = n
        lut[f'SPNG{n["short"]}_{n["instance"]}'] = n

    stmts = []
    for nid in ids:
        n = lut.get(nid)
        if not n:
            continue
        cat = n["cat"]
        label = f'{n["short"]}\\n{n["instance"]}\\nwall={n["wall"]:.2f}s core={n["core"]:.2f}s'
        stmts.append(f'  "{nid}" [style=filled, fillcolor="{CAT_COLOR[cat]}", '
                     f'label="{label}"];')
    if stmts:
        idx = dot_text.rfind("}")
        dot_text = dot_text[:idx] + "\n" + "\n".join(stmts) + "\n}\n"
    ann_dot.write_text(dot_text)

    try:
        r = subprocess.run(["dot", "-Tpng", "-o", str(out_png), str(ann_dot)],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=120)
        if r.returncode != 0 or not out_png.exists():
            return None
    except Exception:
        return None
    return out_png.name


# ---- format emitters ----
def _node_rows(report, n=20):
    return [(x["short"], x["instance"], CAT_LABEL[x["cat"]], x["wall"], x["wall_sd"], x["core"])
            for x in per_node_list(report)[:n]]


def _summary_context(bench):
    meta, reports = bench["meta"], bench["reports"]
    ctx = {"meta": meta, "stages": [s for s in ("osp", "spng") if s in reports]}
    ctx["cat"] = {s: cat_sums(per_node_list(reports[s])) for s in ctx["stages"]}
    ctx["total"] = {s: sum(ctx["cat"][s].values()) for s in ctx["stages"]}
    return ctx


def emit_markdown(bench, ctx, figs, graphs):
    m, R = bench["meta"], bench["reports"]
    L = []
    L.append(f"# spngbench summary — {m.get('detname','?')}\n")
    L.append(f"host **{m.get('host','?')}**, device **{m.get('device','?')}**, "
             f"wc_cores **{m.get('wc_cores','?')}**, torch_cores **{m.get('torch_cores','?')}**, "
             f"gpu_scheme **{m.get('gpu_scheme','none')}** (ngpu {m.get('ngpu',1)})\n")

    L.append("\n## OSP vs SPNG\n")
    L.append("| category | " + " | ".join(s.upper() + " [s]" for s in ctx["stages"]) + " |")
    L.append("|---|" + "---|" * len(ctx["stages"]))
    for cat in CAT_ORDER:
        L.append(f"| {CAT_LABEL[cat]} | " +
                 " | ".join(_fmt(ctx['cat'][s][cat]) for s in ctx["stages"]) + " |")
    L.append(f"| **total** | " + " | ".join(f"**{_fmt(ctx['total'][s])}**" for s in ctx["stages"]) + " |")
    if "osp" in ctx["stages"] and "spng" in ctx["stages"]:
        o, s = ctx["total"]["osp"], ctx["total"]["spng"]
        L.append(f"\nSPNG/OSP total wall ratio: **{_fmt(s / o if o else None, 2)}×**\n")
    L.append(f"\n![category composition]({figs['categories']})\n")

    L.append("\n## Per-node timing\n")
    for s in ctx["stages"]:
        L.append(f"\n### {s.upper()}\n")
        if f"nodes_{s}" in figs:
            L.append(f"![{s} nodes]({figs[f'nodes_{s}']})\n")
        L.append("\n| node | instance | category | wall [s] | ±sd | core [s] |")
        L.append("|---|---|---|--:|--:|--:|")
        for sh, inst, cl, w, sd, co in _node_rows(R[s]):
            L.append(f"| {sh} | {inst} | {cl} | {_fmt(w)} | {_fmt(sd)} | {_fmt(co)} |")

    if graphs:
        L.append("\n## Flow graphs (per-node CPU/GPU timing)\n")
        for s in ctx["stages"]:
            if graphs.get(s):
                L.append(f"\n### {s.upper()}\n")
                L.append(f"![{s} flow graph]({graphs[s]})\n")
    return "\n".join(L) + "\n"


def emit_html(bench, ctx, figs, graphs):
    m, R = bench["meta"], bench["reports"]
    def esc(x):
        return str(x).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    H = ["<!doctype html><meta charset='utf-8'><title>spngbench summary</title>",
         "<style>body{font-family:sans-serif;margin:2em;max-width:60em}"
         "table{border-collapse:collapse}td,th{border:1px solid #ccc;padding:3px 8px}"
         "th{background:#eee}img{max-width:100%}</style>"]
    H.append(f"<h1>spngbench summary — {esc(m.get('detname','?'))}</h1>")
    H.append(f"<p>host <b>{esc(m.get('host','?'))}</b>, device <b>{esc(m.get('device','?'))}</b>, "
             f"wc_cores <b>{m.get('wc_cores','?')}</b>, torch_cores <b>{m.get('torch_cores','?')}</b>, "
             f"gpu_scheme <b>{esc(m.get('gpu_scheme','none'))}</b> (ngpu {m.get('ngpu',1)})</p>")
    H.append("<h2>OSP vs SPNG</h2><table><tr><th>category</th>" +
             "".join(f"<th>{s.upper()} [s]</th>" for s in ctx["stages"]) + "</tr>")
    for cat in CAT_ORDER:
        H.append("<tr><td>" + CAT_LABEL[cat] + "</td>" +
                 "".join(f"<td>{_fmt(ctx['cat'][s][cat])}</td>" for s in ctx["stages"]) + "</tr>")
    H.append("<tr><td><b>total</b></td>" +
             "".join(f"<td><b>{_fmt(ctx['total'][s])}</b></td>" for s in ctx["stages"]) + "</tr></table>")
    H.append(f"<p><img src='{figs['categories']}'></p>")
    H.append("<h2>Per-node timing</h2>")
    for s in ctx["stages"]:
        H.append(f"<h3>{s.upper()}</h3>")
        if f"nodes_{s}" in figs:
            H.append(f"<p><img src='{figs[f'nodes_{s}']}'></p>")
        H.append("<table><tr><th>node</th><th>instance</th><th>category</th>"
                 "<th>wall [s]</th><th>±sd</th><th>core [s]</th></tr>")
        for sh, inst, cl, w, sd, co in _node_rows(R[s]):
            H.append(f"<tr><td>{esc(sh)}</td><td>{esc(inst)}</td><td>{cl}</td>"
                     f"<td>{_fmt(w)}</td><td>{_fmt(sd)}</td><td>{_fmt(co)}</td></tr>")
        H.append("</table>")
    if graphs:
        H.append("<h2>Flow graphs (per-node CPU/GPU timing)</h2>")
        for s in ctx["stages"]:
            if graphs.get(s):
                H.append(f"<h3>{s.upper()}</h3><p><img src='{graphs[s]}'></p>")
    return "\n".join(H) + "\n"


def emit_latex(bench, ctx, figs, graphs, fragment=False, figpre=""):
    m, R = bench["meta"], bench["reports"]
    def esc(x):
        return str(x).replace("_", r"\_").replace("&", r"\&").replace("%", r"\%")
    def fig(name):
        return figpre + name
    ncol = len(ctx["stages"])
    T = []
    if not fragment:
        T += [r"\documentclass{article}",
              r"\usepackage{lmodern}\usepackage{graphicx}\usepackage{booktabs}\usepackage[margin=1in]{geometry}",
              r"\begin{document}"]
    sec = r"\subsection*{" if fragment else r"\section*{"
    T.append(sec + "spngbench summary --- " + esc(m.get("detname", "?")) + "}")
    T.append(f"host \\texttt{{{esc(m.get('host','?'))}}}, device \\texttt{{{esc(m.get('device','?'))}}}, "
             f"wc\\_cores {m.get('wc_cores','?')}, torch\\_cores {m.get('torch_cores','?')}, "
             f"gpu\\_scheme \\texttt{{{esc(m.get('gpu_scheme','none'))}}} (ngpu {m.get('ngpu',1)}).")
    sub = r"\subsubsection*{" if fragment else r"\subsection*{"
    T.append(sub + "OSP vs SPNG}")
    T.append(r"\begin{tabular}{l" + "r" * ncol + "}\\toprule")
    T.append("category & " + " & ".join(s.upper() for s in ctx["stages"]) + r" \\\midrule")
    for cat in CAT_ORDER:
        T.append(f"{CAT_LABEL[cat]} & " +
                 " & ".join(_fmt(ctx['cat'][s][cat]) for s in ctx["stages"]) + r" \\")
    T.append(r"\midrule total & " + " & ".join(_fmt(ctx['total'][s]) for s in ctx["stages"]) + r" \\\bottomrule")
    T.append(r"\end{tabular}")
    T.append(r"\begin{center}\includegraphics[width=0.6\textwidth]{" + fig(figs["categories"]) + r"}\end{center}")
    T.append(sub + "Per-node timing}")
    for s in ctx["stages"]:
        T.append(r"\paragraph{" + s.upper() + "}")
        if f"nodes_{s}" in figs:
            T.append(r"\begin{center}\includegraphics[width=0.8\textwidth]{" + fig(figs[f"nodes_{s}"]) + r"}\end{center}")
        T.append(r"\begin{tabular}{lllrrr}\toprule")
        T.append(r"node & instance & category & wall [s] & sd & core [s] \\\midrule")
        for sh, inst, cl, w, sd, co in _node_rows(R[s], n=14):
            T.append(f"{esc(sh)} & {esc(inst)} & {cl.replace('(non-DNN)','')} & "
                     f"{_fmt(w)} & {_fmt(sd)} & {_fmt(co)} " + r"\\")
        T.append(r"\bottomrule\end{tabular}")
    if graphs:
        T.append(sub + "Flow graphs}")
        for s in ctx["stages"]:
            if graphs.get(s):
                T.append(r"\paragraph{" + s.upper() + "}")
                T.append(r"\begin{center}\includegraphics[width=\textwidth]{" + fig(graphs[s]) + r"}\end{center}")
    if not fragment:
        T.append(r"\end{document}")
    return "\n".join(T) + "\n"


def report_command(bench_path, outdir, formats=("md", "html", "tex"), with_graphs=True,
                   tex_fragment=False, figpre=""):
    """Write a Markdown/HTML/LaTeX summary directory from a benchmark JSON.

    tex_fragment: emit summary.tex as a \\input-able fragment (no preamble);
    figpre: prefix for figure paths in the tex fragment (for inclusion from a
    parent document in a different directory).
    """
    bench = load_bench(bench_path)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    figs = make_figures(bench["reports"], outdir)
    graphs = {}
    if with_graphs:
        for s in bench["reports"]:
            g = render_flow_graph(bench["reports"][s], outdir / f"graph_{s}.png")
            if g:
                graphs[s] = g
    ctx = _summary_context(bench)

    written = []
    if "md" in formats:
        p = outdir / "summary.md"; p.write_text(emit_markdown(bench, ctx, figs, graphs)); written.append(p)
    if "html" in formats:
        p = outdir / "summary.html"; p.write_text(emit_html(bench, ctx, figs, graphs)); written.append(p)
    if "tex" in formats:
        p = outdir / "summary.tex"
        p.write_text(emit_latex(bench, ctx, figs, graphs, fragment=tex_fragment, figpre=figpre))
        written.append(p)
    return {"outdir": str(outdir), "figures": figs, "graphs": graphs,
            "written": [str(p) for p in written]}


# ---------------------------------------------------------------------------
# Grid-level report (across all grid points from a grid-index.json).
# ---------------------------------------------------------------------------
def _report_cat4(report):
    """Return {DNN, sp_other, other, total} wall seconds from a per-node report."""
    d = {"DNN": 0.0, "sp_other": 0.0, "other": 0.0, "total": 0.0}
    for n in per_node_list(report):
        d["total"] += n["wall"]
        if n["cat"] == "dnn":
            d["DNN"] += n["wall"]
        elif n["cat"] == "sp_nondnn":
            d["sp_other"] += n["wall"]
        else:
            d["other"] += n["wall"]
    return d


SUBPIX = ["DNN", "sp_other", "other", "total"]          # the four categories
GPU_MODE_LABEL = {0: "CPU only", 1: "CPU + 1 GPU", 2: "CPU + 2 GPU", 3: "CPU + 3 GPU",
                  4: "CPU + 4 GPU"}


def _gpu_count(meta):
    """Number of GPUs a config uses: shard ngpu, or 1 for a single-GPU device, else 0."""
    if meta.get("gpu_scheme", "none") != "none":
        return int(meta.get("ngpu", 1))
    if str(meta.get("device", "")).startswith("gpu"):
        return 1
    return 0


def grid_cells(base, compare_path=None):
    """Reconstruct the full grid from every compare-*.json in a directory.

    The grid-index.json only records the last grid run, but each run leaves a
    compare-<config>.json per cell, so we glob those to recover all device modes.
    Returns a list of records with config, gpu count, per-stage outcome, and the
    per-category times, plus the compare file path for per-point summaries.
    """
    base = Path(base)
    cells = []
    for cf in sorted(base.glob("compare-*.json")):
        try:
            j = json.loads(cf.read_text())
        except Exception:
            continue
        if not j.get("schema", "").startswith("spngbench-compare"):
            continue
        m = j["meta"]
        rec = {
            "path": str(cf), "meta": m, "gpu": _gpu_count(m),
            "wc": m["wc_cores"], "torch": m["torch_cores"],
            "scheme": m.get("gpu_scheme", "none"), "ngpu": int(m.get("ngpu", 1)),
            "outcome": {"osp": j["osp"]["outcome"], "spng": j["spng"]["outcome"]},
            "cat": {},
        }
        for s in ("osp", "spng"):
            rec["cat"][s] = _report_cat4(j[s]) if j[s].get("outcome") == "ok" else None
        cells.append(rec)
    return cells


def _mode_axes(cells, gpu):
    """Axes for one device mode: X = wc cores, Y = torch cores (or scheme for shard)."""
    ms = [c for c in cells if c["gpu"] == gpu]
    xs = sorted({c["wc"] for c in ms})
    if gpu >= 2:
        ys = sorted({c["scheme"] for c in ms})
        ylabels = ys
        ykey = lambda c: c["scheme"]
    else:
        ys = sorted({c["torch"] for c in ms})
        ylabels = [f"omp{t}" for t in ys]
        ykey = lambda c: c["torch"]
    return ms, xs, ys, ylabels, ykey


def make_grid_matrix(cells, figdir):
    """Outer product of matrix views: (device modes) x (DNN, sp_other, other, total).

    One composite figure; rows = device modes present, columns = categories.  In
    each view every cell is a grid point split into two sub-pixels: left = OSP,
    right = SPNG, coloured by that category's wall time.  Colour scale is per
    category (shared across device modes and OSP/SPNG), log-scaled to span the
    CPU-to-GPU dynamic range.  Not-tested cells are black; crashed sub-pixels
    white.  Returns {"matrix": name} or {}.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import LogNorm
    try:
        cmap = matplotlib.colormaps["viridis"]
    except AttributeError:
        import matplotlib.cm as cm
        cmap = cm.get_cmap("viridis")

    gpus = sorted({c["gpu"] for c in cells})
    if not gpus:
        return {}
    nrows, ncols = len(gpus), len(SUBPIX)

    # Per-category log norm shared across device modes and both stages.
    norms = {}
    for cat in SUBPIX:
        vals = [c["cat"][s][cat] for c in cells for s in ("osp", "spng")
                if c["cat"][s] and c["cat"][s][cat] > 0]
        if vals:
            norms[cat] = LogNorm(vmin=max(min(vals), max(vals) * 1e-3), vmax=max(vals))

    fig, axs = plt.subplots(nrows, ncols, squeeze=False,
                            figsize=(2.6 * ncols + 1, 2.4 * nrows + 1.2))
    for gi, gpu in enumerate(gpus):
        ms, xs, ys, ylabels, ykey = _mode_axes(cells, gpu)
        xidx = {x: i for i, x in enumerate(xs)}
        yidx = {y: i for i, y in enumerate(ys)}
        loc = {(yidx[ykey(c)], xidx[c["wc"]]): c for c in ms}
        for cj, cat in enumerate(SUBPIX):
            ax = axs[gi][cj]
            norm = norms.get(cat)
            for yi in range(len(ys)):
                for xi in range(len(xs)):
                    x0, y0 = xi, yi
                    c = loc.get((yi, xi))
                    if c is None:
                        ax.add_patch(Rectangle((x0, y0), 1, 1, facecolor="black", edgecolor="none"))
                    else:
                        for half, s in ((0, "osp"), (1, "spng")):
                            d = c["cat"][s]
                            if d is None:
                                ax.add_patch(Rectangle((x0 + half * 0.5, y0), 0.5, 1,
                                                       facecolor="white", edgecolor="red", lw=0.3))
                            else:
                                v = max(d[cat], norm.vmin) if norm else 0
                                ax.add_patch(Rectangle((x0 + half * 0.5, y0), 0.5, 1,
                                                       facecolor=cmap(norm(v)) if norm else "black",
                                                       edgecolor="none"))
                    ax.add_patch(Rectangle((x0, y0), 1, 1, fill=False, edgecolor="#999", lw=0.5))
            ax.set_xlim(0, max(1, len(xs))); ax.set_ylim(0, max(1, len(ys))); ax.invert_yaxis()
            ax.set_xticks([i + 0.5 for i in range(len(xs))]); ax.set_xticklabels(xs, fontsize=7)
            ax.set_yticks([i + 0.5 for i in range(len(ys))]); ax.set_yticklabels(ylabels, fontsize=7)
            ax.tick_params(length=0)
            if gi == 0:
                ax.set_title(cat, fontsize=10)
            if cj == 0:
                ax.set_ylabel(GPU_MODE_LABEL.get(gpu, f"{gpu} GPU"), fontsize=9)
            if gi == nrows - 1:
                ax.set_xlabel("wc cores", fontsize=8)

    # One horizontal colour bar per category column (its own log scale).
    fig.subplots_adjust(left=0.10, right=0.98, top=0.90, bottom=0.16, hspace=0.35, wspace=0.25)
    for cj, cat in enumerate(SUBPIX):
        if cat not in norms:
            continue
        x0 = 0.10 + cj * (0.88 / ncols) + 0.02
        cax = fig.add_axes([x0, 0.05, 0.88 / ncols - 0.05, 0.015])
        cb = fig.colorbar(ScalarMappable(norm=norms[cat], cmap=cmap), cax=cax,
                          orientation="horizontal")
        cb.ax.tick_params(labelsize=6)
        cax.set_title(f"{cat} [s]", fontsize=7)
    fig.suptitle("wall time — each cell: left=OSP, right=SPNG   (black=not tested, white=crashed)",
                 fontsize=10)

    p = Path(figdir) / "matrix.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    return {"matrix": p.name}


def make_grid_bars(cells, figdir):
    """Bar charts of higher-level trends across all grid points."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figs = {}
    figdir = Path(figdir)
    ok = [c for c in cells if c["cat"]["osp"] or c["cat"]["spng"]]
    ok = sorted(ok, key=lambda c: (c["gpu"], c["torch"], c["scheme"], c["wc"]))
    if not ok:
        return figs

    def lab(c):
        base = GPU_MODE_LABEL.get(c["gpu"], f"{c['gpu']}gpu").replace("CPU + ", "").replace("CPU only", "cpu")
        conf = c["scheme"][:4] if c["gpu"] >= 2 else f"omp{c['torch']}"
        return f"{base}/{conf}/wc{c['wc']}"

    labels = [lab(c) for c in ok]
    osp_t = [c["cat"]["osp"]["total"] if c["cat"]["osp"] else 0.0 for c in ok]
    spng_t = [c["cat"]["spng"]["total"] if c["cat"]["spng"] else 0.0 for c in ok]
    x = range(len(labels))
    fig, ax = plt.subplots(figsize=(max(6, 0.3 * len(labels) + 2), 4))
    ax.bar([i - 0.2 for i in x], osp_t, width=0.4, label="OSP", color="#5f8fbf")
    ax.bar([i + 0.2 for i in x], spng_t, width=0.4, label="SPNG", color="#d95f5f")
    ax.set_xticks(list(x)); ax.set_xticklabels(labels, rotation=90, fontsize=6)
    ax.set_ylabel("total wall time [s]"); ax.set_yscale("log"); ax.legend()
    ax.set_title("Total wall time per grid point (log)")
    fig.tight_layout()
    p = figdir / "bars_total.png"; fig.savefig(p, dpi=120); plt.close(fig)
    figs["total"] = p.name

    rl, ratios = [], []
    for c in ok:
        if c["cat"]["osp"] and c["cat"]["spng"] and c["cat"]["osp"]["total"]:
            rl.append(lab(c)); ratios.append(c["cat"]["spng"]["total"] / c["cat"]["osp"]["total"])
    if rl:
        fig, ax = plt.subplots(figsize=(max(6, 0.3 * len(rl) + 2), 3.5))
        ax.bar(range(len(rl)), ratios, color="#7a5fbf")
        ax.axhline(1.0, color="k", lw=0.6)
        ax.set_xticks(range(len(rl))); ax.set_xticklabels(rl, rotation=90, fontsize=6)
        ax.set_ylabel("SPNG / OSP total wall"); ax.set_title("SPNG-over-OSP ratio per grid point")
        fig.tight_layout()
        p = figdir / "bars_ratio.png"; fig.savefig(p, dpi=120); plt.close(fig)
        figs["ratio"] = p.name
    return figs


def _grid_overview_rows(cells):
    out = []
    for c in sorted(cells, key=lambda c: (c["gpu"], c["torch"], c["scheme"], c["wc"])):
        conf = c["scheme"] if c["gpu"] >= 2 else f"omp{c['torch']}"
        def cell4(d):
            return tuple(_fmt(d[k]) for k in SUBPIX) if d else ("crash",) * 4
        out.append((GPU_MODE_LABEL.get(c["gpu"], f"{c['gpu']}gpu"), conf, c["wc"],
                    cell4(c["cat"]["osp"]), cell4(c["cat"]["spng"])))
    return out


def render_stage_graph(tlas, out_png):
    """Render a plain data-flow graph for one stage via `wcpy pgraph dotify`.

    Uses --no-services --no-params so service components and config parameters
    are omitted.  Returns the PNG basename or None if a tool is missing.
    """
    if not shutil.which("dot") or not shutil.which("wcpy"):
        return None
    out_png = Path(out_png)
    base_dot = out_png.with_suffix(".dot")
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = ""
    env.setdefault("MPLCONFIGDIR", "/tmp/mpl-spngbench")
    dargs = ["wcpy", "pgraph", "dotify", "--no-services", "--no-params"]
    for k, v in tlas.items():
        dargs += ["-A", f"{k}={v}"]
    dargs += [str(JSONNET), str(base_dot)]
    try:
        r = subprocess.run(dargs, env=env, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, timeout=120)
        if r.returncode != 0 or not base_dot.exists():
            return None
        r = subprocess.run(["dot", "-Tpng", "-o", str(out_png), str(base_dot)],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=120)
        if r.returncode != 0 or not out_png.exists():
            return None
    except Exception:
        return None
    return out_png.name


def make_job_graphs(grid, cells, outdir):
    """Render sim/osp/spng data-flow graphs and gather the driving input files.

    Returns (job_graphs, inputs) where job_graphs maps stage -> png basename and
    inputs is the list of depo files used to drive the benchmark.
    """
    depos = grid["meta"].get("depos") or []
    adc = grid["meta"].get("adc_inputs") or []
    # Recover model_file / detname from a representative per-point command.
    tlas0 = {}
    if cells:
        try:
            b = load_bench(cells[0]["path"])
            _, tlas0 = _tlas_from_cmd(b["reports"]["osp"]["meta"].get("cmd", []))
        except Exception:
            pass
    common = dict(model_file=tlas0.get("model_file", DEFAULT_MODEL), output="signals.npz",
                  detname=tlas0.get("detname", "pdhd"), device="cpu", gpu_scheme="none",
                  ngpu=1, napa=grid["meta"].get("napa", tlas0.get("napa", 1)),
                  engine="Pgrapher", wc_cores=1, verbosity=0)
    stage_input = {"sim": (depos[0] if depos else "depos.npz"),
                   "osp": (adc[0] if adc else "adc.npz"),
                   "spng": (adc[0] if adc else "adc.npz")}
    job_graphs = {}
    for st in ("sim", "osp", "spng"):
        g = render_stage_graph(dict(common, stage=st, input=stage_input[st]),
                               Path(outdir) / f"job_{st}.png")
        if g:
            job_graphs[st] = g
    return job_graphs, depos


JOB_DESC = {
    "sim": "**sim** — depos → drift → detsim → ADC frame file.  Run once per input; "
           "its ADC output feeds both OSP and SPNG.",
    "osp": "**osp** — ADC frame file → OmnibusSigProc → DNN-ROI → signal frame file.",
    "spng": "**spng** — ADC frame file → SPNG (decon / filters / ROI-uniter DNN) → "
            "signal frame file.",
}


def grid_report_command(grid_path, outdir, formats=("md", "html", "tex"),
                        point_graphs=False):
    """Produce a cross-grid report from a grid-index.json.

    Renders OSP/SPNG sub-pixel matrices and trend bars, a grid overview table,
    per-grid-point summaries (linked in md/html, included in latex), and a TOC.
    """
    grid = json.loads(Path(grid_path).read_text())
    if not grid.get("schema", "").startswith("spngbench-grid"):
        raise SystemExit("grid-report expects a grid-index.json (schema spngbench-grid/*)")
    base = Path(grid_path).parent
    outdir = Path(outdir)
    (outdir / "points").mkdir(parents=True, exist_ok=True)

    # The grid-index records only the last run; reconstruct the full grid (all
    # device modes) from every compare-*.json in the same directory.
    cells = grid_cells(base)
    mats = make_grid_matrix(cells, outdir)
    bars = make_grid_bars(cells, outdir)
    job_graphs, inputs = make_job_graphs(grid, cells, outdir)

    # Per-point summaries (one sub-directory each), for every compare file.
    points = []
    for c in sorted(cells, key=lambda c: (c["gpu"], c["torch"], c["scheme"], c["wc"])):
        conf = c["scheme"] if c["gpu"] >= 2 else f"omp{c['torch']}"
        rlab = f"{GPU_MODE_LABEL.get(c['gpu'], str(c['gpu']) + 'gpu')} {conf}"
        stem = f"{c['gpu']}gpu-{conf}-wc{c['wc']}".replace(" ", "_")
        pdir = outdir / "points" / stem
        report_command(c["path"], str(pdir), formats=formats,
                       with_graphs=point_graphs, tex_fragment=True,
                       figpre=f"points/{stem}/")
        ratio = None
        if c["cat"]["osp"] and c["cat"]["spng"] and c["cat"]["osp"]["total"]:
            ratio = c["cat"]["spng"]["total"] / c["cat"]["osp"]["total"]
        points.append((stem, rlab, c["wc"], ratio))

    ctx = {"grid": grid, "mats": mats, "bars": bars,
           "overview": _grid_overview_rows(cells), "points": points,
           "ncells": len(cells), "jobs": job_graphs, "inputs": inputs}
    written = []
    if "md" in formats:
        p = outdir / "grid-summary.md"; p.write_text(_emit_grid_md(ctx)); written.append(p)
    if "html" in formats:
        p = outdir / "grid-summary.html"; p.write_text(_emit_grid_html(ctx)); written.append(p)
    if "tex" in formats:
        p = outdir / "grid-summary.tex"; p.write_text(_emit_grid_tex(ctx)); written.append(p)
    return {"outdir": str(outdir), "written": [str(p) for p in written],
            "npoints": len(points), "matrices": mats, "bars": bars}


_MATRIX_BLURB = ("The outer product of device modes (rows) and node categories "
                 "(columns).  Each cell is a grid point split into two sub-pixels: "
                 "**left = OSP, right = SPNG**, coloured by that category's wall time "
                 "on a per-column log scale.  Not-tested cells are black; crashed "
                 "sub-pixels white.")


def _emit_grid_md(ctx):
    m = ctx["grid"]["meta"]
    L = [f"# spngbench grid summary — {m.get('detname','?')}\n",
         f"host **{m.get('host','?')}** ({m.get('host_ngpu','?')} GPU), engine "
         f"**{m.get('engine','?')}**, {ctx['ncells']} grid points.\n",
         "\n## Contents\n",
         "- [Jobs](#jobs)",
         "- [Grid matrix](#grid-matrix)",
         "- [Trends](#trends)",
         "- [Overview table](#overview-table)",
         "- Per-grid-point summaries:"]
    for stem, rlab, w, ratio in ctx["points"]:
        L.append(f"  - [{rlab} / wc{w}](points/{stem}/summary.md) (SPNG/OSP {_fmt(ratio,2)}×)")

    L.append("\n## Jobs\n")
    L.append("The benchmark factors the work into three wire-cell jobs.  Data-flow "
             "graphs below are from `wcpy pgraph dotify --no-services --no-params` "
             "(service components and config parameters omitted).\n")
    L.append("\n**Input depo file(s):**\n")
    for f in (ctx["inputs"] or ["(none recorded)"]):
        L.append(f"- `{f}`")
    for st in ("sim", "osp", "spng"):
        L.append(f"\n### {st}\n")
        L.append(JOB_DESC[st] + "\n")
        if ctx["jobs"].get(st):
            L.append(f"\n![{st} graph]({ctx['jobs'][st]})\n")

    L.append("\n## Grid matrix\n")
    L.append(_MATRIX_BLURB.replace("**", "**") + "\n")
    if ctx["mats"].get("matrix"):
        L.append(f"\n![grid matrix]({ctx['mats']['matrix']})\n")

    L.append("\n## Trends\n")
    for k in ("total", "ratio"):
        if ctx["bars"].get(k):
            L.append(f"\n![{k}]({ctx['bars'][k]})\n")

    L.append("\n## Overview table\n")
    L.append("| mode | config | wc | OSP DNN/spo/oth/tot | SPNG DNN/spo/oth/tot |")
    L.append("|---|---|--:|---|---|")
    for mode, conf, w, o, s in ctx["overview"]:
        L.append(f"| {mode} | {conf} | {w} | {'/'.join(map(str,o))} | {'/'.join(map(str,s))} |")
    return "\n".join(L) + "\n"


def _emit_grid_html(ctx):
    m = ctx["grid"]["meta"]
    def esc(x): return str(x).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    H = ["<!doctype html><meta charset='utf-8'><title>spngbench grid summary</title>",
         "<style>body{font-family:sans-serif;margin:2em;max-width:72em}"
         "table{border-collapse:collapse}td,th{border:1px solid #ccc;padding:2px 6px}"
         "th{background:#eee}img{max-width:100%}</style>",
         f"<h1>spngbench grid summary — {esc(m.get('detname','?'))}</h1>",
         f"<p>host <b>{esc(m.get('host','?'))}</b> ({m.get('host_ngpu','?')} GPU), engine "
         f"<b>{esc(m.get('engine','?'))}</b>, {ctx['ncells']} grid points.</p>",
         "<h2 id='toc'>Contents</h2><ul>",
         "<li><a href='#jobs'>Jobs</a></li>",
         "<li><a href='#matrix'>Grid matrix</a></li>",
         "<li><a href='#trends'>Trends</a></li>",
         "<li><a href='#overview'>Overview table</a></li>",
         "<li>Per-grid-point summaries:<ul>"]
    for stem, rlab, w, ratio in ctx["points"]:
        H.append(f"<li><a href='points/{stem}/summary.html'>{esc(rlab)} / wc{w}</a> "
                 f"(SPNG/OSP {_fmt(ratio,2)}×)</li>")
    H.append("</ul></li></ul>")
    H.append("<h2 id='jobs'>Jobs</h2>")
    H.append("<p>The benchmark factors the work into three wire-cell jobs.  Graphs are "
             "from <code>wcpy pgraph dotify --no-services --no-params</code> (service "
             "components and config parameters omitted).</p>")
    H.append("<p><b>Input depo file(s):</b></p><ul>")
    for f in (ctx["inputs"] or ["(none recorded)"]):
        H.append(f"<li><code>{esc(f)}</code></li>")
    H.append("</ul>")
    for st in ("sim", "osp", "spng"):
        H.append(f"<h3>{st}</h3><p>{esc(JOB_DESC[st]).replace('**','')}</p>")
        if ctx["jobs"].get(st):
            H.append(f"<p><img src='{ctx['jobs'][st]}'></p>")
    H.append("<h2 id='matrix'>Grid matrix</h2>")
    H.append("<p>" + _MATRIX_BLURB.replace("**", "") + "</p>")
    if ctx["mats"].get("matrix"):
        H.append(f"<p><img src='{ctx['mats']['matrix']}'></p>")
    H.append("<h2 id='trends'>Trends</h2>")
    for k in ("total", "ratio"):
        if ctx["bars"].get(k):
            H.append(f"<p><img src='{ctx['bars'][k]}'></p>")
    H.append("<h2 id='overview'>Overview table</h2>")
    H.append("<table><tr><th>mode</th><th>config</th><th>wc</th>"
             "<th>OSP DNN/spo/oth/tot</th><th>SPNG DNN/spo/oth/tot</th></tr>")
    for mode, conf, w, o, s in ctx["overview"]:
        H.append(f"<tr><td>{esc(mode)}</td><td>{esc(str(conf))}</td><td>{w}</td>"
                 f"<td>{'/'.join(map(str,o))}</td><td>{'/'.join(map(str,s))}</td></tr>")
    H.append("</table>")
    return "\n".join(H) + "\n"


def _emit_grid_tex(ctx):
    m = ctx["grid"]["meta"]
    def esc(x): return str(x).replace("_", r"\_").replace("&", r"\&").replace("%", r"\%")
    T = [r"\documentclass{article}",
         r"\usepackage{lmodern}\usepackage{graphicx}\usepackage{booktabs}\usepackage[margin=1in]{geometry}",
         r"\usepackage{hyperref}",
         r"\begin{document}",
         r"\title{spngbench grid summary --- " + esc(m.get("detname", "?")) + "}",
         r"\author{}\date{}\maketitle",
         r"\tableofcontents\newpage",
         r"\section{Jobs}",
         r"The benchmark factors the work into three wire-cell jobs.  Data-flow graphs "
         r"are from \texttt{wcpy pgraph dotify -{}-no-services -{}-no-params} (service "
         r"components and config parameters omitted).",
         r"\paragraph{Input depo file(s):}~\\"]
    for f in (ctx["inputs"] or ["(none recorded)"]):
        T.append(r"\texttt{" + esc(f) + r"}\\")
    for st in ("sim", "osp", "spng"):
        T.append(r"\subsection*{" + st + "}")
        T.append(esc(JOB_DESC[st].replace("**", "")))
        if ctx["jobs"].get(st):
            T.append(r"\begin{center}\includegraphics[width=\textwidth,height=0.55\textheight,"
                     r"keepaspectratio]{" + ctx["jobs"][st] + r"}\end{center}")
    T += [r"\clearpage",
         r"\section{Grid matrix}",
         "The outer product of device modes (rows) and node categories (columns). "
         "Each cell is a grid point split into two sub-pixels: left = OSP, right = SPNG, "
         "coloured by that category's wall time on a per-column log scale.  "
         "Not-tested cells are black; crashed sub-pixels white."]
    if ctx["mats"].get("matrix"):
        T.append(r"\begin{center}\includegraphics[width=\textwidth]{" + ctx["mats"]["matrix"] + r"}\end{center}")
    T.append(r"\section{Trends}")
    for k in ("total", "ratio"):
        if ctx["bars"].get(k):
            T.append(r"\begin{center}\includegraphics[width=\textwidth]{" + ctx["bars"][k] + r"}\end{center}")
    T.append(r"\section{Overview table}")
    T.append(r"\small\begin{tabular}{lllll}\toprule")
    T.append(r"mode & config & wc & OSP D/s/o/t & SPNG D/s/o/t \\\midrule")
    for mode, conf, w, o, s in ctx["overview"]:
        T.append(f"{esc(mode)} & {esc(str(conf))} & {w} & {esc('/'.join(map(str,o)))} & "
                 f"{esc('/'.join(map(str,s)))} " + r"\\")
    T.append(r"\bottomrule\end{tabular}\normalsize")
    T.append(r"\section{Per-grid-point summaries}")
    for stem, rlab, w, ratio in ctx["points"]:
        T.append(r"\clearpage")
        T.append(r"\input{points/" + stem + "/summary.tex}")
    T.append(r"\end{document}")
    return "\n".join(T) + "\n"


# ---------------------------------------------------------------------------
# CLI.
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
    p.add_argument("--gpu-scheme", default="none", choices=["none", "transverse", "longitudinal"],
                   help="SPNG GPU sharding scheme (default none)")
    p.add_argument("--ngpu", type=int, default=1, help="number of GPUs the scheme may use")
    p.add_argument("--napa", type=int, default=1,
                   help="number of APAs / per-APA pipelines (clamped to detector max)")
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

    pg = sub.add_parser("grid", help="scan a grid of (device, wc_cores, torch_cores) configs")
    pg.add_argument("--config", default=None, help="JSON/YAML config file (merged over defaults)")
    pg.add_argument("--mode", choices=["cpu", "gpu", "shard", "two-gpu", "all"], default="cpu",
                    help="which cell family to run (default cpu); 'shard' = multi-GPU")
    pg.add_argument("--which", choices=["osp", "spng", "both"], default="both")
    pg.add_argument("--wc-cores", type=int, nargs="+", default=None,
                    help="subset of wire-cell core counts")
    pg.add_argument("--torch-cores", type=int, nargs="+", default=None,
                    help="subset of torch core counts (cpu mode)")
    pg.add_argument("--napa", type=int, default=None,
                    help="number of APAs / per-APA pipelines (fixed, not scanned)")
    pg.add_argument("--outdir", default=None)
    pg.add_argument("--repeats", type=int, default=None)
    pg.add_argument("--input", action="append", default=None)
    pg.add_argument("--model-file", default=None)
    pg.add_argument("--detname", default=None)
    pg.add_argument("--dry-run", action="store_true")

    psc = sub.add_parser("show-config", help="print the effective config (defaults + --config)")
    psc.add_argument("--config", default=None)

    prep = sub.add_parser("report", help="write a Markdown/HTML/LaTeX summary from a benchmark JSON")
    prep.add_argument("bench_json", help="a compare-*.json, report-*.json, or grid-index.json")
    prep.add_argument("-o", "--outdir", default="spngbench-report", help="output directory")
    prep.add_argument("--formats", default="md,html,tex",
                      help="comma list of md,html,tex (default all)")
    prep.add_argument("--no-graphs", action="store_true", help="skip the flow-graph figures")

    pgr = sub.add_parser("grid-report", help="cross-grid summary from a grid-index.json")
    pgr.add_argument("grid_json", help="a grid-index.json")
    pgr.add_argument("-o", "--outdir", default="spngbench-grid-report", help="output directory")
    pgr.add_argument("--formats", default="md,html,tex", help="comma list of md,html,tex")
    pgr.add_argument("--point-graphs", action="store_true",
                     help="also render flow graphs in each per-point summary (slow)")

    args = ap.parse_args(argv)

    if args.cmd == "analyze":
        print(json.dumps(analyze_log(args.logfile), indent=2))
        return 0

    if args.cmd == "show-config":
        cfg = load_config(args.config) if args.config else default_config()
        print(json.dumps(cfg, indent=2))
        return 0

    if args.cmd == "report":
        formats = tuple(f.strip() for f in args.formats.split(",") if f.strip())
        res = report_command(args.bench_json, args.outdir, formats=formats,
                             with_graphs=not args.no_graphs)
        log(f"wrote {len(res['written'])} summary file(s), "
            f"{len(res['figures'])} figure(s), {len(res['graphs'])} flow graph(s)")
        for w in res["written"]:
            print(w)
        return 0

    if args.cmd == "grid-report":
        formats = tuple(f.strip() for f in args.formats.split(",") if f.strip())
        res = grid_report_command(args.grid_json, args.outdir, formats=formats,
                                  point_graphs=args.point_graphs)
        log(f"wrote grid summary + {res['npoints']} per-point summaries; "
            f"matrices={list(res['matrices'])} bars={list(res['bars'])}")
        for w in res["written"]:
            print(w)
        return 0

    if args.cmd == "grid":
        cfg = load_config(args.config) if args.config else default_config()
        # CLI overrides win over the config file.
        for key, val in (("outdir", args.outdir), ("repeats", args.repeats),
                         ("model_file", args.model_file), ("detname", args.detname),
                         ("napa", args.napa)):
            if val is not None:
                cfg[key] = val
        if args.input:
            cfg["inputs"] = args.input
        which_list = ["osp", "spng"] if args.which == "both" else [args.which]
        cells = select_cells(cfg, args.mode, wc_sel=args.wc_cores, torch_sel=args.torch_cores)
        if not cells:
            log("no cells selected")
            return 1
        log(f"grid: mode={args.mode} which={which_list} cells={len(cells)} "
            f"outdir={cfg['outdir']}")
        idx = run_grid(cfg, cells, which_list, dry_run=args.dry_run)
        # Compact human-facing summary to stdout.
        for c in idx["cells"]:
            stag = f" {c['gpu_scheme']}{c['ngpu']}" if c.get("gpu_scheme", "none") != "none" else ""
            label = f"  {c['device']} wc{c['wc_cores']} omp{c['torch_cores']}{stag}"
            if c.get("skipped"):
                print(f"{label}: SKIPPED ({c['reason']})")
            else:
                osp = c.get("osp", {}); spng = c.get("spng", {})
                print(f"{label}: osp={osp.get('outcome')}({osp.get('sp_total_wall')}s) "
                      f"spng={spng.get('outcome')}({spng.get('sp_total_wall')}s) "
                      f"spng/osp={c.get('spng_over_osp_sp_total')}")
        print(f"grid index: {idx['_path']}")
        return 0

    depos = _default_inputs(args)

    # run/pair: sim each depo input to an ADC file once, then run the SP chain(s).
    napa = clamp_napa({"detname": args.detname, "napa": args.napa})
    adc_inputs = [ensure_adc(d, args.outdir, model_file=args.model_file,
                             detname=args.detname, napa=napa, dry_run=args.dry_run)
                  for d in depos]

    if args.cmd == "run":
        scheme = args.gpu_scheme if args.which == "spng" else "none"
        rep = bench_config(args.which, args.device, args.wc_cores, args.torch_cores,
                           adc_inputs, args.outdir, repeats=args.repeats,
                           model_file=args.model_file, detname=args.detname,
                           engine=args.engine, gpu_scheme=scheme, ngpu=args.ngpu,
                           napa=napa, dry_run=args.dry_run)
        print(f"outcome={rep['outcome']} wrote {rep['_path']}")
        return 0

    if args.cmd == "pair":
        osp = bench_config("osp", args.device, args.wc_cores, args.torch_cores,
                           adc_inputs, args.outdir, repeats=args.repeats,
                           model_file=args.model_file, detname=args.detname,
                           engine=args.engine, napa=napa, dry_run=args.dry_run)
        spng = bench_config("spng", args.device, args.wc_cores, args.torch_cores,
                            adc_inputs, args.outdir, repeats=args.repeats,
                            model_file=args.model_file, detname=args.detname,
                            engine=args.engine, gpu_scheme=args.gpu_scheme, ngpu=args.ngpu,
                            napa=napa, dry_run=args.dry_run)
        comb = combine(osp, spng, outdir=args.outdir)
        print(f"osp={osp['outcome']} spng={spng['outcome']} wrote {comb['_path']}")
        return 0

    return 1


if __name__ == "__main__":
    sys.exit(main())
