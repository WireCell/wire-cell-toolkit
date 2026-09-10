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
              detname="pdhd", engine="TbbFlow", gpu_scheme="none", ngpu=1, verbosity=0):
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
        "-A", f"engine={engine}",
        "-A", f"wc_cores={wc_cores}",
        "-A", f"verbosity={verbosity}",
    ]


def run_wirecell(stage, device, wc_cores, torch_cores, input, outdir,
                 model_file=DEFAULT_MODEL, detname="pdhd", engine="TbbFlow",
                 gpu_scheme="none", ngpu=1, tag="", output=None,
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
                    ngpu=ngpu, verbosity=verbosity)

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(torch_cores)
    if extra_env:
        env.update(extra_env)

    result = {
        "stage": stage, "device": device, "wc_cores": wc_cores,
        "torch_cores": torch_cores, "engine": engine, "input": str(input),
        "gpu_scheme": gpu_scheme, "ngpu": ngpu,
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
                detname="pdhd", gpu_scheme="none", ngpu=1):
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
            "gpu_scheme": gpu_scheme, "ngpu": ngpu,
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
               engine="Pgrapher", dry_run=False):
    """Run the sim stage (depos -> ADC frame file) once for a depo input, cached.

    Returns the ADC file path.  The sim is CPU-only and shared by OSP and SPNG,
    so it is run once per depo input and reused across the whole grid.  Later the
    returned ADC file can be replaced by one made from real detector data.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    stem = "adc-" + Path(depo_input).stem
    adc = outdir / f"{stem}.npz"
    if adc.exists() and adc.stat().st_size > 0:
        return str(adc)
    if dry_run:
        return str(adc)
    r = run_wirecell("sim", "cpu", 1, 1, depo_input, outdir,
                     model_file=model_file, detname=detname, engine=engine,
                     tag=stem, output=str(adc))
    if r["outcome"] != "ok":
        log(f"WARNING: sim failed for {depo_input}: {r['outcome']}")
    return str(adc)


def bench_config(stage, device, wc_cores, torch_cores, inputs, outdir,
                 repeats=1, model_file=DEFAULT_MODEL, detname="pdhd",
                 engine="TbbFlow", gpu_scheme="none", ngpu=1,
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
                gpu_scheme=gpu_scheme, ngpu=ngpu,
                tag=base + suffix, verbosity=verbosity, dry_run=dry_run))

    report = make_report(stage, device, wc_cores, torch_cores, runs,
                         engine=engine, detname=detname,
                         gpu_scheme=gpu_scheme, ngpu=ngpu)
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


def cpu_cells(cfg, wc_sel=None, torch_sel=None):
    """Enumerate CPU grid cells (device=cpu), skipping wc*torch > max_hyperthreads."""
    g = cfg["cpu_grid"]
    cap = cfg["hw"]["max_hyperthreads"]
    wcs = wc_sel or g["wc_cores"]
    tcs = torch_sel or g["torch_cores"]
    return [_cell("cpu", wc, tc) for wc in wcs for tc in tcs if wc * tc <= cap]


def gpu_cells(cfg, wc_sel=None):
    """Enumerate single-GPU grid cells (all GPU-capable nodes on device, torch=1)."""
    g = cfg["gpu_grid"]
    wcs = wc_sel or g["wc_cores"]
    return [_cell(g.get("device", "gpu"), wc, g.get("torch_cores", 1)) for wc in wcs]


def shard_cells(cfg, wc_sel=None):
    """Enumerate multi-GPU sharded cells: (scheme x ngpu x wc_cores)."""
    g = cfg["shard_grid"]
    wcs = wc_sel or g["wc_cores"]
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

    # Run the sim once per depo input; OSP and SPNG share the ADC frame file(s).
    log(f"sim: producing ADC frames for {len(depos)} depo input(s)")
    adc_inputs = [ensure_adc(d, str(outdir), model_file=cfg["model_file"],
                             detname=cfg["detname"], dry_run=dry_run) for d in depos]

    index = {
        "schema": "spngbench-grid/1",
        "meta": {
            "detname": cfg["detname"], "host": socket.gethostname(), "host_ngpu": host_ngpu,
            "engine": cfg["engine"], "repeats": cfg["repeats"],
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
                gpu_scheme=w_scheme, ngpu=ngpu, dry_run=dry_run)
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
              r"\usepackage{graphicx}\usepackage{booktabs}\usepackage[margin=1in]{geometry}",
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
def _cell_row(cell):
    """(sort_key, label) for a cell's compute-config row of the grid matrix."""
    if cell.get("gpu_scheme", "none") != "none":
        return ((2, cell["ngpu"], cell["gpu_scheme"]), f"GPU×{cell['ngpu']} {cell['gpu_scheme'][:5]}")
    if str(cell.get("device", "")).startswith("gpu"):
        return ((1, 1, ""), "GPU×1")
    return ((0, cell["torch_cores"], ""), f"CPU omp{cell['torch_cores']}")


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


SUBPIX = ["DNN", "sp_other", "other", "total"]   # 2x2: [DNN sp_other; other total]


def _grid_matrix_data(grid, base):
    """Assemble the grid matrix: axes and per-cell category times / status.

    Returns (rows, cols, data) where data[stage][(ri,ci)] is a dict with keys
    'status' (ok|crash|absent) and the four category times when ok.  Cells not
    present in the grid are 'absent'; cells present but crashed are 'crash'.
    """
    cells = grid.get("cells", [])
    cols = sorted({c["wc_cores"] for c in cells})
    rowset = {}
    for c in cells:
        k, lab = _cell_row(c)
        rowset[k] = lab
    rows = sorted(rowset)
    row_labels = [rowset[k] for k in rows]
    col_index = {w: i for i, w in enumerate(cols)}
    row_index = {k: i for i, k in enumerate(rows)}

    stages = grid["meta"].get("which", ["osp", "spng"])
    data = {s: {} for s in stages}
    for c in cells:
        rk, _ = _cell_row(c)
        ri, ci = row_index[rk], col_index[c["wc_cores"]]
        skipped = c.get("skipped")
        cpath = c.get("compare_path")
        loaded = None
        if not skipped and cpath:
            for cand in (Path(cpath), base / Path(cpath).name):
                if cand.exists():
                    loaded = load_bench(cand)
                    break
        for s in stages:
            key = (ri, ci)
            if skipped:
                data[s][key] = {"status": "absent"}
            elif loaded and s in loaded["reports"]:
                oc = c.get(s, {}).get("outcome")
                if oc == "ok":
                    d = _report_cat4(loaded["reports"][s])
                    d["status"] = "ok"
                    data[s][key] = d
                else:
                    data[s][key] = {"status": "crash"}
            else:
                # No compare loaded: fall back to the cell summary numbers.
                cs = c.get(s, {})
                if cs.get("outcome") == "ok" and cs.get("all_total_wall") is not None:
                    dnn = cs.get("forward_wall") or 0.0
                    spo = cs.get("sp_rest_wall") or 0.0
                    tot = cs.get("all_total_wall") or 0.0
                    data[s][key] = {"status": "ok", "DNN": dnn, "sp_other": spo,
                                    "other": max(0.0, tot - dnn - spo), "total": tot}
                else:
                    data[s][key] = {"status": "crash" if cs else "absent"}
    return rows, row_labels, cols, data


def make_grid_matrices(grid, base, figdir):
    """Render the OSP and SPNG sub-pixel heatmap matrices on a shared scale."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    import matplotlib.cm as cm
    from matplotlib.colors import Normalize
    try:
        cmap = matplotlib.colormaps["viridis"]
    except AttributeError:
        cmap = cm.get_cmap("viridis")

    rows, row_labels, cols, data = _grid_matrix_data(grid, base)
    stages = [s for s in ("osp", "spng") if s in data]
    R, C = len(rows), len(cols)
    if R == 0 or C == 0:
        return {}, {"rows": row_labels, "cols": cols}

    vmax = 0.0
    for s in stages:
        for cell in data[s].values():
            if cell.get("status") == "ok":
                vmax = max(vmax, cell["total"])
    vmax = vmax or 1.0
    norm = Normalize(vmin=0.0, vmax=vmax)
    # 2x2 sub-pixel offsets: DNN(0,0) sp_other(0,1) other(1,0) total(1,1)
    off = {"DNN": (0, 0), "sp_other": (1, 0), "other": (0, 1), "total": (1, 1)}

    figs = {}
    figdir = Path(figdir)
    for s in stages:
        fig, ax = plt.subplots(figsize=(max(4, 1.1 * C + 2), max(3, 0.6 * R + 1.5)))
        for (ri, ci), cell in data[s].items():
            x0, y0 = ci, ri
            st = cell.get("status")
            if st == "absent":
                ax.add_patch(Rectangle((x0, y0), 1, 1, facecolor="black", edgecolor="none"))
            elif st == "crash":
                ax.add_patch(Rectangle((x0, y0), 1, 1, facecolor="white", edgecolor="red", lw=0.5))
            elif st == "ok":
                for cat, (dx, dy) in off.items():
                    ax.add_patch(Rectangle((x0 + dx * 0.5, y0 + dy * 0.5), 0.5, 0.5,
                                           facecolor=cmap(norm(cell[cat])), edgecolor="none"))
            ax.add_patch(Rectangle((x0, y0), 1, 1, fill=False, edgecolor="#888", lw=0.8))
        ax.set_xlim(0, C); ax.set_ylim(0, R); ax.invert_yaxis()
        ax.set_xticks([c + 0.5 for c in range(C)]); ax.set_xticklabels(cols)
        ax.set_yticks([r + 0.5 for r in range(R)]); ax.set_yticklabels(row_labels, fontsize=8)
        ax.set_xlabel("wire-cell cores (wc_cores)")
        ax.set_title(f"{s.upper()} wall time  [sub-pixels: DNN | sp_other / other | total]")
        ax.tick_params(length=0)
        fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, label="wall time [s]", shrink=0.8)
        fig.tight_layout()
        p = figdir / f"matrix_{s}.png"
        fig.savefig(p, dpi=130)
        plt.close(fig)
        figs[s] = p.name
    return figs, {"rows": row_labels, "cols": cols}


def make_grid_bars(grid, base, figdir):
    """Bar charts of higher-level trends across the grid."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows, row_labels, cols, data = _grid_matrix_data(grid, base)
    stages = [s for s in ("osp", "spng") if s in data]
    figs = {}
    figdir = Path(figdir)

    # Total wall per cell, OSP vs SPNG, as grouped bars over all grid points.
    labels, osp_t, spng_t = [], [], []
    for ri, rlab in enumerate(row_labels):
        for ci, w in enumerate(cols):
            key = (ri, ci)
            od = data.get("osp", {}).get(key, {})
            sd = data.get("spng", {}).get(key, {})
            if od.get("status") == "ok" or sd.get("status") == "ok":
                labels.append(f"{rlab}/wc{w}")
                osp_t.append(od.get("total", 0.0) if od.get("status") == "ok" else 0.0)
                spng_t.append(sd.get("total", 0.0) if sd.get("status") == "ok" else 0.0)
    if labels:
        x = range(len(labels))
        fig, ax = plt.subplots(figsize=(max(6, 0.32 * len(labels) + 2), 4))
        ax.bar([i - 0.2 for i in x], osp_t, width=0.4, label="OSP", color="#5f8fbf")
        ax.bar([i + 0.2 for i in x], spng_t, width=0.4, label="SPNG", color="#d95f5f")
        ax.set_xticks(list(x)); ax.set_xticklabels(labels, rotation=90, fontsize=6)
        ax.set_ylabel("total wall time [s]"); ax.legend()
        ax.set_title("Total wall time per grid point")
        fig.tight_layout()
        p = figdir / "bars_total.png"; fig.savefig(p, dpi=120); plt.close(fig)
        figs["total"] = p.name

    # SPNG/OSP total ratio per grid point.
    rlabels, ratios = [], []
    for ri, rlab in enumerate(row_labels):
        for ci, w in enumerate(cols):
            od = data.get("osp", {}).get((ri, ci), {})
            sd = data.get("spng", {}).get((ri, ci), {})
            if od.get("status") == "ok" and sd.get("status") == "ok" and od.get("total"):
                rlabels.append(f"{rlab}/wc{w}"); ratios.append(sd["total"] / od["total"])
    if rlabels:
        fig, ax = plt.subplots(figsize=(max(6, 0.32 * len(rlabels) + 2), 3.5))
        ax.bar(range(len(rlabels)), ratios, color="#7a5fbf")
        ax.axhline(1.0, color="k", lw=0.6)
        ax.set_xticks(range(len(rlabels))); ax.set_xticklabels(rlabels, rotation=90, fontsize=6)
        ax.set_ylabel("SPNG / OSP total wall"); ax.set_title("SPNG-over-OSP ratio per grid point")
        fig.tight_layout()
        p = figdir / "bars_ratio.png"; fig.savefig(p, dpi=120); plt.close(fig)
        figs["ratio"] = p.name
    return figs


def _grid_overview_rows(grid, base):
    rows, row_labels, cols, data = _grid_matrix_data(grid, base)
    out = []
    for ri, rlab in enumerate(row_labels):
        for ci, w in enumerate(cols):
            od = data.get("osp", {}).get((ri, ci), {})
            sd = data.get("spng", {}).get((ri, ci), {})
            if od.get("status") == "absent" and sd.get("status") == "absent":
                continue
            def cell4(d):
                if d.get("status") != "ok":
                    return (d.get("status", "-"),) * 4
                return (_fmt(d["DNN"]), _fmt(d["sp_other"]), _fmt(d["other"]), _fmt(d["total"]))
            out.append((rlab, w, cell4(od), cell4(sd)))
    return out


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

    mats, axes = make_grid_matrices(grid, base, outdir)
    bars = make_grid_bars(grid, base, outdir)

    # Per-point summaries (one sub-directory each).
    points = []
    for c in grid.get("cells", []):
        if c.get("skipped"):
            continue
        cp = c.get("compare_path")
        cand = None
        for x in ((Path(cp),) if cp else ()) + ((base / Path(cp).name,) if cp else ()):
            if x.exists():
                cand = x; break
        if not cand:
            continue
        _, rlab = _cell_row(c)
        stem = f"{rlab}-wc{c['wc_cores']}".replace(" ", "_").replace("×", "x")
        pdir = outdir / "points" / stem
        report_command(str(cand), str(pdir), formats=formats,
                       with_graphs=point_graphs, tex_fragment=True,
                       figpre=f"points/{stem}/")
        points.append((stem, rlab, c["wc_cores"], c.get("spng_over_osp_sp_total")))

    ctx = {"grid": grid, "mats": mats, "bars": bars, "axes": axes,
           "overview": _grid_overview_rows(grid, base), "points": points}
    written = []
    if "md" in formats:
        p = outdir / "grid-summary.md"; p.write_text(_emit_grid_md(ctx)); written.append(p)
    if "html" in formats:
        p = outdir / "grid-summary.html"; p.write_text(_emit_grid_html(ctx)); written.append(p)
    if "tex" in formats:
        p = outdir / "grid-summary.tex"; p.write_text(_emit_grid_tex(ctx)); written.append(p)
    return {"outdir": str(outdir), "written": [str(p) for p in written],
            "npoints": len(points), "matrices": mats, "bars": bars}


def _emit_grid_md(ctx):
    g = ctx["grid"]; m = g["meta"]
    L = [f"# spngbench grid summary — {m.get('detname','?')}\n",
         f"host **{m.get('host','?')}** ({m.get('host_ngpu','?')} GPU), engine **{m.get('engine','?')}**, "
         f"{m.get('ncells','?')} cells, {len(ctx['points'])} run.\n"]
    L.append("\n## Contents\n")
    L.append("- [Grid matrices](#grid-matrices)")
    L.append("- [Trends](#trends)")
    L.append("- [Overview table](#overview-table)")
    L.append("- Per-grid-point summaries:")
    for stem, rlab, w, ratio in ctx["points"]:
        L.append(f"  - [{rlab} / wc{w}](points/{stem}/summary.md) (SPNG/OSP {_fmt(ratio,2)}×)")

    L.append("\n## Grid matrices\n")
    L.append("Each cell is a grid point; 2×2 sub-pixels are **DNN** (top-left), "
             "**sp_other** (top-right), **other** (bottom-left), **total** (bottom-right), "
             "on a shared time colour scale.  Black = not tested, white = crashed.\n")
    for s in ("osp", "spng"):
        if ctx["mats"].get(s):
            L.append(f"\n![{s} matrix]({ctx['mats'][s]})\n")

    L.append("\n## Trends\n")
    for k in ("total", "ratio"):
        if ctx["bars"].get(k):
            L.append(f"\n![{k}]({ctx['bars'][k]})\n")

    L.append("\n## Overview table\n")
    L.append("| config | wc | OSP DNN/spo/oth/tot | SPNG DNN/spo/oth/tot |")
    L.append("|---|--:|---|---|")
    for rlab, w, o, s in ctx["overview"]:
        L.append(f"| {rlab} | {w} | {'/'.join(map(str,o))} | {'/'.join(map(str,s))} |")
    return "\n".join(L) + "\n"


def _emit_grid_html(ctx):
    g = ctx["grid"]; m = g["meta"]
    def esc(x): return str(x).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    H = ["<!doctype html><meta charset='utf-8'><title>spngbench grid summary</title>",
         "<style>body{font-family:sans-serif;margin:2em;max-width:70em}"
         "table{border-collapse:collapse}td,th{border:1px solid #ccc;padding:2px 6px}"
         "th{background:#eee}img{max-width:100%}</style>",
         f"<h1>spngbench grid summary — {esc(m.get('detname','?'))}</h1>",
         f"<p>host <b>{esc(m.get('host','?'))}</b> ({m.get('host_ngpu','?')} GPU), engine "
         f"<b>{esc(m.get('engine','?'))}</b>, {m.get('ncells','?')} cells, {len(ctx['points'])} run.</p>",
         "<h2 id='toc'>Contents</h2><ul>",
         "<li><a href='#matrices'>Grid matrices</a></li>",
         "<li><a href='#trends'>Trends</a></li>",
         "<li><a href='#overview'>Overview table</a></li>",
         "<li>Per-grid-point summaries:<ul>"]
    for stem, rlab, w, ratio in ctx["points"]:
        H.append(f"<li><a href='points/{stem}/summary.html'>{esc(rlab)} / wc{w}</a> "
                 f"(SPNG/OSP {_fmt(ratio,2)}×)</li>")
    H.append("</ul></li></ul>")
    H.append("<h2 id='matrices'>Grid matrices</h2>")
    H.append("<p>Each cell is a grid point; 2×2 sub-pixels are DNN (top-left), sp_other "
             "(top-right), other (bottom-left), total (bottom-right), shared colour scale. "
             "Black = not tested, white = crashed.</p>")
    for s in ("osp", "spng"):
        if ctx["mats"].get(s):
            H.append(f"<p><img src='{ctx['mats'][s]}'></p>")
    H.append("<h2 id='trends'>Trends</h2>")
    for k in ("total", "ratio"):
        if ctx["bars"].get(k):
            H.append(f"<p><img src='{ctx['bars'][k]}'></p>")
    H.append("<h2 id='overview'>Overview table</h2>")
    H.append("<table><tr><th>config</th><th>wc</th><th>OSP DNN/spo/oth/tot</th>"
             "<th>SPNG DNN/spo/oth/tot</th></tr>")
    for rlab, w, o, s in ctx["overview"]:
        H.append(f"<tr><td>{esc(rlab)}</td><td>{w}</td><td>{'/'.join(map(str,o))}</td>"
                 f"<td>{'/'.join(map(str,s))}</td></tr>")
    H.append("</table>")
    return "\n".join(H) + "\n"


def _emit_grid_tex(ctx):
    g = ctx["grid"]; m = g["meta"]
    def esc(x): return str(x).replace("_", r"\_").replace("&", r"\&").replace("%", r"\%")
    T = [r"\documentclass{article}",
         r"\usepackage{graphicx}\usepackage{booktabs}\usepackage[margin=1in]{geometry}",
         r"\usepackage{hyperref}",
         r"\begin{document}",
         r"\title{spngbench grid summary --- " + esc(m.get("detname", "?")) + "}",
         r"\author{}\date{}\maketitle",
         r"\tableofcontents\newpage",
         r"\section{Grid matrices}",
         r"Each cell is a grid point; 2$\times$2 sub-pixels are DNN, sp\_other, other, total "
         r"on a shared time colour scale.  Black = not tested, white = crashed."]
    for s in ("osp", "spng"):
        if ctx["mats"].get(s):
            T.append(r"\begin{center}\includegraphics[width=0.8\textwidth]{" + ctx["mats"][s] + r"}\end{center}")
    T.append(r"\section{Trends}")
    for k in ("total", "ratio"):
        if ctx["bars"].get(k):
            T.append(r"\begin{center}\includegraphics[width=\textwidth]{" + ctx["bars"][k] + r"}\end{center}")
    T.append(r"\section{Overview table}")
    T.append(r"\small\begin{tabular}{llll}\toprule")
    T.append(r"config & wc & OSP D/s/o/t & SPNG D/s/o/t \\\midrule")
    for rlab, w, o, s in ctx["overview"]:
        T.append(f"{esc(rlab)} & {w} & {esc('/'.join(map(str,o)))} & {esc('/'.join(map(str,s)))} " + r"\\")
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
                         ("model_file", args.model_file), ("detname", args.detname)):
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
    adc_inputs = [ensure_adc(d, args.outdir, model_file=args.model_file,
                             detname=args.detname, dry_run=args.dry_run) for d in depos]

    if args.cmd == "run":
        scheme = args.gpu_scheme if args.which == "spng" else "none"
        rep = bench_config(args.which, args.device, args.wc_cores, args.torch_cores,
                           adc_inputs, args.outdir, repeats=args.repeats,
                           model_file=args.model_file, detname=args.detname,
                           engine=args.engine, gpu_scheme=scheme, ngpu=args.ngpu,
                           dry_run=args.dry_run)
        print(f"outcome={rep['outcome']} wrote {rep['_path']}")
        return 0

    if args.cmd == "pair":
        osp = bench_config("osp", args.device, args.wc_cores, args.torch_cores,
                           adc_inputs, args.outdir, repeats=args.repeats,
                           model_file=args.model_file, detname=args.detname,
                           engine=args.engine, dry_run=args.dry_run)
        spng = bench_config("spng", args.device, args.wc_cores, args.torch_cores,
                            adc_inputs, args.outdir, repeats=args.repeats,
                            model_file=args.model_file, detname=args.detname,
                            engine=args.engine, gpu_scheme=args.gpu_scheme, ngpu=args.ngpu,
                            dry_run=args.dry_run)
        comb = combine(osp, spng, outdir=args.outdir)
        print(f"osp={osp['outcome']} spng={spng['outcome']} wrote {comb['_path']}")
        return 0

    return 1


if __name__ == "__main__":
    sys.exit(main())
