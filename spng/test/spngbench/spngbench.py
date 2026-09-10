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
import threading
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


def analyze_timeline(path):
    """Parse a TbbFlow node timeline JSON into (nodes, total).

    Each node carries per-execution [start,end] intervals; wall-sec is the sum of
    interval lengths (== the old Timer wall total), core-sec from core_sum.
    """
    j = json.loads(Path(path).read_text())
    nodes = []
    tot_wall = tot_core = 0.0
    for n in j.get("nodes", []):
        ivals = n.get("intervals", [])
        wall = sum(b - a for a, b in ivals) if ivals else n.get("wall_sum", 0.0)
        core = n.get("core_sum", 0.0)
        tot_wall += wall
        tot_core += core
        nodes.append({
            "instance": n.get("instance", ""), "class": n.get("class", ""),
            "wall-sec": wall, "core-sec": core, "calls": n.get("calls", 0),
            "intervals": ivals,
        })
    return nodes, {"wall": tot_wall, "core": tot_core}


def analyze_log(logfile):
    """Parse one wire-cell run into {phases, nodes, timer_total, timeline}.

    Per-node timing comes from the TbbFlow node timeline JSON when present (it
    carries per-execution intervals); otherwise it falls back to parsing the
    log's Timer/summary lines.  Phase intervals always come from log timestamps.
    """
    with open(logfile, errors="replace") as fp:
        lines = fp.readlines()
    phases = analyze_phases(lines)
    tl = Path(str(logfile)[:-4] + ".timeline.json") if str(logfile).endswith(".log") \
        else Path(str(logfile) + ".timeline.json")
    if tl.exists():
        nodes, total = analyze_timeline(tl)
        return {"phases": phases, "nodes": nodes, "timer_total": total,
                "timeline": str(tl)}
    nodes, total = analyze_nodes(lines)
    return {"phases": phases, "nodes": nodes, "timer_total": total, "timeline": None}


# ---------------------------------------------------------------------------
# Running one wire-cell job.
# ---------------------------------------------------------------------------
_OOM_RE = re.compile(r"out of memory|CUDA error: out of memory|CUDA_ERROR_OUT_OF_MEMORY",
                     re.IGNORECASE)


def build_cmd(stage, device, wc_cores, input, model_file, output, logfile,
              detname="pdhd", engine="TbbFlow", gpu_scheme="none", ngpu=1, napa=1,
              apa=-1, timeline="", verbosity=0):
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
        "-A", f"apa={apa}",
        "-A", f"engine={engine}",
        "-A", f"wc_cores={wc_cores}",
        "-A", f"timeline={timeline}",
        "-A", f"verbosity={verbosity}",
    ]


# ---------------------------------------------------------------------------
# Memory sampler: samples the child process's CPU RSS and per-process GPU VRAM
# at a fixed rate while wire-cell runs, timestamped with CLOCK_REALTIME so the
# samples align with the TbbFlow node timeline.
# ---------------------------------------------------------------------------
def _read_rss_bytes(pid):
    try:
        with open(f"/proc/{pid}/status") as fp:
            for line in fp:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1]) * 1024   # kB -> bytes
    except OSError:
        return None
    return None


def _nvml_init():
    """Return a pynvml module with handles, or None (fall back to nvidia-smi)."""
    try:
        import pynvml
        pynvml.nvmlInit()
        return pynvml
    except Exception:
        return None


def _vram_bytes_nvml(nvml, pid):
    total = 0
    try:
        for i in range(nvml.nvmlDeviceGetCount()):
            h = nvml.nvmlDeviceGetHandleByIndex(i)
            for p in nvml.nvmlDeviceGetComputeRunningProcesses(h):
                if p.pid == pid and p.usedGpuMemory:
                    total += int(p.usedGpuMemory)
    except Exception:
        pass
    return total


def _vram_bytes_smi(pid):
    try:
        out = subprocess.run(
            ["nvidia-smi", "--query-compute-apps=pid,used_memory",
             "--format=csv,noheader,nounits"],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, timeout=3).stdout.decode()
        total = 0
        for line in out.splitlines():
            parts = [x.strip() for x in line.split(",")]
            if len(parts) >= 2 and parts[0].isdigit() and int(parts[0]) == pid:
                total += int(parts[1]) * 1024 * 1024   # MiB -> bytes
        return total
    except Exception:
        return 0


class MemorySampler(threading.Thread):
    """Sample a PID's RSS (+ per-process VRAM when gpu=True) at `rate` Hz."""
    def __init__(self, pid, rate=20.0, gpu=False):
        super().__init__(daemon=True)
        self.pid = pid
        self.dt = 1.0 / rate
        self.rate = rate
        self.gpu = gpu
        self.samples = []                 # list of [t, rss_bytes, vram_bytes]
        self._halt = threading.Event()
        self._nvml = _nvml_init() if gpu else None

    def _vram(self):
        if not self.gpu:
            return 0
        if self._nvml:
            return _vram_bytes_nvml(self._nvml, self.pid)
        return _vram_bytes_smi(self.pid)

    def run(self):
        while not self._halt.is_set():
            t = time.time()
            rss = _read_rss_bytes(self.pid)
            if rss is None:
                break                     # process gone
            self.samples.append([t, rss, self._vram()])
            rest = self.dt - (time.time() - t)
            if rest > 0:
                self._halt.wait(rest)

    def stop(self):
        self._halt.set()
        self.join(timeout=3)

    def result(self):
        rss = [s[1] for s in self.samples]
        vram = [s[2] for s in self.samples]
        return {
            "clock": "CLOCK_REALTIME", "rate_hz": self.rate,
            "gpu": self.gpu, "vram_tool": ("pynvml" if self._nvml else "nvidia-smi"),
            "nsamples": len(self.samples),
            "peak_rss": max(rss, default=0), "peak_vram": max(vram, default=0),
            "samples": self.samples,
        }


def _vram_map_smi(pids):
    """One nvidia-smi call -> {pid: used_vram_bytes} for pids in the given set."""
    try:
        out = subprocess.run(
            ["nvidia-smi", "--query-compute-apps=pid,used_memory",
             "--format=csv,noheader,nounits"],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, timeout=3).stdout.decode()
        m = {}
        for line in out.splitlines():
            parts = [x.strip() for x in line.split(",")]
            if len(parts) >= 2 and parts[0].isdigit():
                pid = int(parts[0])
                if pid in pids:
                    m[pid] = m.get(pid, 0) + int(parts[1]) * 1024 * 1024
        return m
    except Exception:
        return {}


class GroupMemorySampler(threading.Thread):
    """Sample the SUMMED RSS (+ summed per-process VRAM) across a set of PIDs at
    `rate` Hz.  For the process-parallel axis: the peak of the summed footprint is
    the group's simultaneous RAM/VRAM demand (the number that decides whether N
    jobs co-fit on one host / one shared GPU)."""
    def __init__(self, pids, rate=20.0, gpu=False):
        super().__init__(daemon=True)
        self.pids = list(pids)
        self.dt = 1.0 / rate
        self.rate = rate
        self.gpu = gpu
        self.samples = []                 # [t, sum_rss_bytes, sum_vram_bytes]
        self.peak_rss = {p: 0 for p in self.pids}
        self.peak_vram = {p: 0 for p in self.pids}
        self._halt = threading.Event()
        self._nvml = _nvml_init() if gpu else None

    def _vram_map(self):
        if not self.gpu:
            return {}
        if self._nvml:
            return {p: _vram_bytes_nvml(self._nvml, p) for p in self.pids}
        return _vram_map_smi(set(self.pids))

    def run(self):
        while not self._halt.is_set():
            t = time.time()
            srss = 0
            alive = False
            for p in self.pids:
                r = _read_rss_bytes(p)
                if r is not None:
                    alive = True
                    srss += r
                    if r > self.peak_rss[p]:
                        self.peak_rss[p] = r
            vm = self._vram_map()
            svram = sum(vm.values())
            for p, v in vm.items():
                if v > self.peak_vram.get(p, 0):
                    self.peak_vram[p] = v
            if not alive:
                break                     # all processes gone
            self.samples.append([t, srss, svram])
            rest = self.dt - (time.time() - t)
            if rest > 0:
                self._halt.wait(rest)

    def stop(self):
        self._halt.set()
        self.join(timeout=3)

    def result(self):
        rss = [s[1] for s in self.samples]
        vram = [s[2] for s in self.samples]
        return {
            "clock": "CLOCK_REALTIME", "rate_hz": self.rate, "gpu": self.gpu,
            "vram_tool": ("pynvml" if self._nvml else "nvidia-smi"),
            "npids": len(self.pids), "nsamples": len(self.samples),
            "peak_rss_sum": max(rss, default=0),
            "peak_vram_sum": max(vram, default=0),
            "per_pid_peak_rss": self.peak_rss,
            "per_pid_peak_vram": self.peak_vram,
            "samples": self.samples,
        }


def run_wirecell(stage, device, wc_cores, torch_cores, input, outdir,
                 model_file=DEFAULT_MODEL, detname="pdhd", engine="TbbFlow",
                 gpu_scheme="none", ngpu=1, napa=1, tag="", output=None,
                 mem_rate=20.0, verbosity=0, extra_env=None, dry_run=False):
    """
    Run one wire-cell job (stage=sim|osp|spng) and classify the outcome.

    Collects a TbbFlow node timeline JSON (per-execution intervals) and a memory
    profile JSON (RSS + per-process VRAM, sampled at mem_rate Hz) alongside the
    log.  outcome is one of: ok | oom | error | skipped(dry).
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    stem = tag or f"{stage}-{device}-wc{wc_cores}-omp{torch_cores}"
    logfile = outdir / f"{stem}.log"
    output = Path(output) if output else outdir / f"{stem}.npz"
    # The node timeline is a TbbFlow feature; skip it for Pgrapher (e.g. sim).
    timeline = outdir / f"{stem}.timeline.json" if engine == "TbbFlow" else None
    memprofile = outdir / f"{stem}.memprofile.json"

    cmd = build_cmd(stage, device, wc_cores, input, model_file, output, logfile,
                    detname=detname, engine=engine, gpu_scheme=gpu_scheme,
                    ngpu=ngpu, napa=napa, timeline=(str(timeline) if timeline else ""),
                    verbosity=verbosity)

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(torch_cores)
    if extra_env:
        env.update(extra_env)

    result = {
        "stage": stage, "device": device, "wc_cores": wc_cores,
        "torch_cores": torch_cores, "engine": engine, "input": str(input),
        "gpu_scheme": gpu_scheme, "ngpu": ngpu, "napa": napa,
        "detname": detname, "logfile": str(logfile), "output": str(output),
        "timeline": str(timeline) if timeline else None,
        "memprofile": str(memprofile), "cmd": cmd, "stem": stem,
    }

    if dry_run:
        result["outcome"] = "skipped"
        result["returncode"] = None
        result["wall_clock"] = 0.0
        return result

    gpu_run = str(device).startswith("gpu") or gpu_scheme != "none"
    t0 = time.time()
    proc = subprocess.Popen(cmd, env=env, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)
    sampler = MemorySampler(proc.pid, rate=mem_rate, gpu=gpu_run)
    sampler.start()
    captured = proc.stdout.read().decode(errors="replace") if proc.stdout else ""
    proc.wait()
    sampler.stop()
    result["wall_clock"] = round(time.time() - t0, 3)
    result["returncode"] = proc.returncode
    try:
        memprofile.write_text(json.dumps(sampler.result()))
    except OSError:
        pass

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
# Process-parallel axis: launch nproc concurrent single-APA jobs (one per APA),
# each wc1/omp1 on a shared CPU or one shared GPU, and measure the aggregate
# wall time and the summed RAM/VRAM footprint.  This models running one job per
# APA in parallel (as in production) and, crucially, the VRAM pressure of N jobs
# sharing one GPU -- the number that decides how many co-fit (e.g. RTX-4090 24GB
# vs L40S 48GB).
# ---------------------------------------------------------------------------
def run_proc_group(which, device, nproc, adc_base, outdir, model_file=DEFAULT_MODEL,
                   detname="pdhd", engine="TbbFlow", napa_phys=4, mem_rate=20.0,
                   verbosity=0, dry_run=False):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    stem = f"proc-{which}-{device}-n{nproc}"
    gpu_run = str(device).startswith("gpu")

    def adc_for(k):
        cand = str(Path(str(adc_base)).with_suffix("")) + f"-tpc{k}.npz"
        return cand if Path(cand).exists() else str(adc_base)

    base_meta = {
        "which": which, "device": device, "nproc": nproc, "wc_cores": 1,
        "torch_cores": 1, "gpu_scheme": "none", "ngpu": 1, "napa_phys": napa_phys,
        "detname": detname, "engine": engine, "host": socket.gethostname(),
        "host_ngpu": _ngpu(), "created": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }

    if dry_run:
        rec = {"schema": "spngbench-proc/1", "meta": base_meta, "outcome": "skipped",
               "wall_clock": 0.0, "procs": [], "memory": None}
        path = outdir / f"{stem}.json"
        path.write_text(json.dumps(rec, indent=2))
        rec["_path"] = str(path)
        return rec

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"

    launched = []                         # (apa, popen, outfh, logfile, timeline, output, cmd)
    t0 = time.time()
    for k in range(nproc):
        logfile = outdir / f"{stem}-apa{k}.log"
        output = outdir / f"{stem}-apa{k}.npz"
        timeline = outdir / f"{stem}-apa{k}.timeline.json" if engine == "TbbFlow" else None
        cmd = build_cmd(which, device, 1, adc_for(k), model_file, output, logfile,
                        detname=detname, engine=engine, gpu_scheme="none", ngpu=1,
                        napa=napa_phys, apa=k,
                        timeline=(str(timeline) if timeline else ""), verbosity=verbosity)
        outfh = open(outdir / f"{stem}-apa{k}.out", "wb")
        p = subprocess.Popen(cmd, env=env, stdout=outfh, stderr=subprocess.STDOUT)
        launched.append((k, p, outfh, logfile, timeline, output, cmd))

    sampler = GroupMemorySampler([p.pid for _, p, *_ in launched],
                                 rate=mem_rate, gpu=gpu_run)
    sampler.start()
    for _, p, outfh, *_ in launched:
        p.wait()
        outfh.close()
    sampler.stop()
    wall = round(time.time() - t0, 3)

    procs = []
    n_ok = 0
    any_oom = False
    for k, p, _outfh, logfile, timeline, output, cmd in launched:
        oom = False
        if logfile.exists():
            try:
                with open(logfile, errors="replace") as fp:
                    oom = bool(_OOM_RE.search(fp.read()))
            except OSError:
                pass
        oc = "ok" if p.returncode == 0 else ("oom" if oom else "error")
        n_ok += (oc == "ok")
        any_oom = any_oom or oom
        procs.append({
            "apa": k, "returncode": p.returncode, "outcome": oc,
            "input": adc_for(k), "logfile": str(logfile),
            "timeline": str(timeline) if timeline else None, "output": str(output),
        })
    outcome = ("ok" if n_ok == nproc else
               "oom" if any_oom else
               "partial" if n_ok else "error")

    memres = sampler.result()
    memprofile = outdir / f"{stem}.memprofile.json"
    try:
        memprofile.write_text(json.dumps(memres))
    except OSError:
        pass

    rec = {
        "schema": "spngbench-proc/1",
        "meta": {**base_meta, "cmd": launched[0][6] if launched else None},
        "outcome": outcome,
        "wall_clock": wall,
        "procs": procs,
        "memory": {
            "peak_rss_sum": memres["peak_rss_sum"],
            "peak_vram_sum": memres["peak_vram_sum"],
            "per_pid_peak_rss": memres["per_pid_peak_rss"],
            "per_pid_peak_vram": memres["per_pid_peak_vram"],
            "vram_tool": memres["vram_tool"], "nsamples": memres["nsamples"],
        },
        "memprofile": str(memprofile),
    }
    path = outdir / f"{stem}.json"
    path.write_text(json.dumps(rec, indent=2))
    rec["_path"] = str(path)
    return rec


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


# Map the fine-grained node category to the 3-way memory-blame category.
_BLAME_CAT = {"dnn": "dnn", "sp_nondnn": "sp_other", "io": "other", "other": "other"}


def load_memprofile(run):
    p = run.get("memprofile") if isinstance(run, dict) else run
    if p and Path(p).exists():
        try:
            return json.loads(Path(p).read_text())
        except Exception:
            return None
    return None


def memory_blame(analysis, mp):
    """Correlate node intervals with memory samples -> peak RSS/VRAM per category.

    A sample is attributed to a category if its time lies within any execution
    interval of a node in that category.  Categories overlap under parallelism,
    so per-category peaks are not additive.  Requires the timeline (intervals);
    returns only the overall peaks if intervals are absent.
    """
    samples = (mp or {}).get("samples", [])
    out = {"peak_rss": (mp or {}).get("peak_rss", 0),
           "peak_vram": (mp or {}).get("peak_vram", 0), "by_category": {}}
    if not samples:
        return out
    cats = {"dnn": [], "sp_other": [], "other": []}
    for nd in analysis.get("nodes", []):
        key = _BLAME_CAT.get(report_category(nd["class"]), "other")
        cats[key].extend(nd.get("intervals", []))
    for key, ivals in cats.items():
        ivals = sorted(ivals)
        rss = vram = 0
        for t, r, v in samples:
            if any(a <= t <= b for a, b in ivals):
                rss = max(rss, r); vram = max(vram, v)
        out["by_category"][key] = {"peak_rss": rss, "peak_vram": vram}
    return out


def memory_summary(runs, analyses):
    """Aggregate memory blame over ok runs (peak = max across runs)."""
    blames = []
    for r, a in zip(runs, analyses):
        mp = load_memprofile(r)
        if mp:
            blames.append(memory_blame(a, mp))
    if not blames:
        return None
    agg = {"peak_rss": max(b["peak_rss"] for b in blames),
           "peak_vram": max(b["peak_vram"] for b in blames),
           "by_category": {}}
    for key in ("dnn", "sp_other", "other"):
        rss = [b["by_category"].get(key, {}).get("peak_rss", 0) for b in blames]
        vram = [b["by_category"].get(key, {}).get("peak_vram", 0) for b in blames]
        agg["by_category"][key] = {"peak_rss": max(rss, default=0),
                                   "peak_vram": max(vram, default=0)}
    return agg


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
            "timeline": r.get("timeline"), "memprofile": r.get("memprofile"),
            **({"stdio_tail": r["stdio_tail"]} if r.get("stdio_tail") else {}),
        } for r in runs],
        "phases": rollup_phases(analyses) if analyses else None,
        "summary": {
            "per_node": rollup_nodes(analyses),
            "categories": rollup_categories(analyses),
        } if analyses else None,
        "memory": memory_summary(ok, analyses) if analyses else None,
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
        "proc_grid": {                  # process-level parallelism: nproc concurrent
            "devices": ["cpu", "gpu0"], # single-APA jobs (wc1/omp1), CPU-only or 1 GPU
            "nproc": [1, 2, 3, 4],      # (never 2-GPU); nproc clamped to physical APAs
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


def proc_cells(cfg, nproc_sel=None):
    """Enumerate process-parallel cells: nproc concurrent single-APA jobs on each
    configured device (cpu / one GPU), wc1/omp1.  nproc is clamped to the
    detector's physical APA count."""
    g = cfg.get("proc_grid", {})
    phys = DETECTOR_MAX_APA.get(cfg.get("detname", "pdhd")) or max(
        (nproc_sel or g.get("nproc", [1])), default=1)
    nprocs = [n for n in (nproc_sel or g.get("nproc", [1, 2, 3, 4])) if 1 <= n <= phys]
    devs = g.get("devices", ["cpu", "gpu0"])
    return [{"device": d, "wc_cores": 1, "torch_cores": 1, "gpu_scheme": "none",
             "ngpu": 1, "nproc": n} for d in devs for n in nprocs]


def select_cells(cfg, mode, wc_sel=None, torch_sel=None, nproc_sel=None):
    """Return the list of cells for the selected mode: cpu|gpu|shard|proc|all."""
    cells = []
    if mode in ("cpu", "all"):
        cells += cpu_cells(cfg, wc_sel, torch_sel)
    if mode in ("gpu", "all"):
        cells += gpu_cells(cfg, wc_sel)
    if mode in ("shard", "two-gpu", "all"):
        cells += shard_cells(cfg, wc_sel)
    if mode in ("proc", "all"):
        cells += proc_cells(cfg, nproc_sel)
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

    # The process-parallel axis addresses individual APAs by per-APA ADC file, so
    # the sim must produce at least max(nproc) per-APA files.
    max_nproc = max((c["nproc"] for c in cells if c.get("nproc")), default=0)
    sim_napa = max(napa, max_nproc)

    # Run the sim once per depo input; OSP and SPNG share the ADC frame file(s).
    log(f"sim: producing ADC frames for {len(depos)} depo input(s), napa={sim_napa}")
    adc_inputs = [ensure_adc(d, str(outdir), model_file=cfg["model_file"],
                             detname=cfg["detname"], napa=sim_napa, dry_run=dry_run)
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
        "proc": [],                     # process-parallel records (schema proc/1)
    }
    index_path = outdir / "grid-index.json"

    def flush():
        index_path.write_text(json.dumps(index, indent=2))

    phys_napa = DETECTOR_MAX_APA.get(cfg["detname"]) or sim_napa

    for i, cell in enumerate(cells):
        device = cell["device"]

        # Process-parallel cell: nproc concurrent single-APA jobs (per which).
        if cell.get("nproc"):
            nproc = cell["nproc"]
            if device.startswith("gpu") and host_ngpu < 1:
                index["proc"].append({"device": device, "nproc": nproc,
                                      "skipped": True, "reason": "needs 1 GPU, have 0"})
                flush()
                continue
            log(f"[cell {i+1}/{len(cells)}] proc {device} nproc{nproc} which={which_list}")
            for w in which_list:
                rec = run_proc_group(w, device, nproc, adc_inputs[0], str(outdir),
                                     model_file=cfg["model_file"], detname=cfg["detname"],
                                     engine=cfg["engine"], napa_phys=phys_napa,
                                     dry_run=dry_run)
                mem = rec.get("memory") or {}
                index["proc"].append({
                    "which": w, "device": device, "nproc": nproc,
                    "outcome": rec["outcome"], "wall_clock": rec["wall_clock"],
                    "peak_rss_sum": mem.get("peak_rss_sum"),
                    "peak_vram_sum": mem.get("peak_vram_sum"),
                    "path": rec.get("_path"),
                })
            flush()
            continue

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

    if ctx.get("memfigs"):
        L.append("\n## Memory\n")
        L.append("Peak RSS / VRAM, and per-category peak (a category's peak is the "
                 "max sample while any node of that category is active; overlapping, "
                 "not additive):\n")
        L.append("\n| stage | peak RSS [MB] | peak VRAM [MB] | DNN | sp_other | other (RSS MB) |")
        L.append("|---|--:|--:|--:|--:|--:|")
        for s in ctx["stages"]:
            mem = (R[s].get("memory") or {})
            bc = mem.get("by_category", {})
            def mb(x): return _fmt((x or 0) / 1e6, 0)
            L.append(f"| {s.upper()} | {mb(mem.get('peak_rss'))} | {mb(mem.get('peak_vram'))} | "
                     f"{mb(bc.get('dnn',{}).get('peak_rss'))} | {mb(bc.get('sp_other',{}).get('peak_rss'))} | "
                     f"{mb(bc.get('other',{}).get('peak_rss'))} |")
        for s in ctx["stages"]:
            if ctx["memfigs"].get(s):
                L.append(f"\n### {s.upper()} memory + node activity\n")
                L.append(f"![{s} memory]({ctx['memfigs'][s]})\n")

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
    if ctx.get("memfigs"):
        H.append("<h2>Memory</h2>")
        H.append("<table><tr><th>stage</th><th>peak RSS [MB]</th><th>peak VRAM [MB]</th>"
                 "<th>DNN</th><th>sp_other</th><th>other (RSS MB)</th></tr>")
        for s in ctx["stages"]:
            mem = (R[s].get("memory") or {}); bc = mem.get("by_category", {})
            def mb(x): return _fmt((x or 0) / 1e6, 0)
            H.append(f"<tr><td>{s.upper()}</td><td>{mb(mem.get('peak_rss'))}</td>"
                     f"<td>{mb(mem.get('peak_vram'))}</td><td>{mb(bc.get('dnn',{}).get('peak_rss'))}</td>"
                     f"<td>{mb(bc.get('sp_other',{}).get('peak_rss'))}</td>"
                     f"<td>{mb(bc.get('other',{}).get('peak_rss'))}</td></tr>")
        H.append("</table>")
        for s in ctx["stages"]:
            if ctx["memfigs"].get(s):
                H.append(f"<h3>{s.upper()} memory + node activity</h3><p><img src='{ctx['memfigs'][s]}'></p>")

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
    if ctx.get("memfigs"):
        T.append(sub + "Memory}")
        for s in ctx["stages"]:
            if ctx["memfigs"].get(s):
                T.append(r"\paragraph{" + s.upper() + "}")
                T.append(r"\begin{center}\includegraphics[width=\textwidth,height=0.5\textheight,"
                         r"keepaspectratio]{" + fig(ctx["memfigs"][s]) + r"}\end{center}")

    if graphs:
        T.append(sub + "Flow graphs}")
        for s in ctx["stages"]:
            if graphs.get(s):
                T.append(r"\paragraph{" + s.upper() + "}")
                T.append(r"\begin{center}\includegraphics[width=\textwidth]{" + fig(graphs[s]) + r"}\end{center}")
    if not fragment:
        T.append(r"\end{document}")
    return "\n".join(T) + "\n"


def _resolve_json(path, base):
    if not path:
        return None
    for cand in (Path(path), Path(base) / Path(path).name):
        if cand.exists():
            try:
                return json.loads(cand.read_text())
            except Exception:
                return None
    return None


def make_memory_figures(bench, base, figdir):
    """Per-stage memory-vs-time profile + node-activity Gantt (needs memprofile
    and timeline JSON, i.e. runs made with memory profiling on).  Returns
    {stage: png} for whichever stages have data."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figdir = Path(figdir)
    figs = {}
    for stage in [s for s in ("osp", "spng") if s in bench["reports"]]:
        rep = bench["reports"][stage]
        runs = rep.get("runs") or []
        if not runs:
            continue
        mp = _resolve_json(runs[0].get("memprofile"), base)
        if not mp or not mp.get("samples"):
            continue
        tl = _resolve_json(runs[0].get("timeline"), base)
        samples = mp["samples"]
        t0 = samples[0][0]
        ts = [s[0] - t0 for s in samples]
        rss = [s[1] / 1e6 for s in samples]
        vram = [s[2] / 1e6 for s in samples]
        has_v = bool(mp.get("gpu")) and max(vram) > 0

        nodes = [n for n in (tl or {}).get("nodes", []) if n.get("intervals")]
        nodes.sort(key=lambda n: min(a for a, b in n["intervals"]))
        nrows = len(nodes)
        gh = max(1.2, min(9.0, 0.14 * nrows))
        fig, (axp, axg) = plt.subplots(
            2, 1, figsize=(9, 2.6 + gh), sharex=True,
            gridspec_kw={"height_ratios": [2.4, gh]})
        axp.plot(ts, rss, color="#5f8fbf", lw=1.2, label="RSS")
        axp.set_ylabel("RSS [MB]", color="#5f8fbf")
        axp.set_title(f"{stage.upper()} — memory vs time and node activity")
        # X+Y grid on the profile; X grid on the timeline, both sharing X so their
        # vertical grid lines align.
        axp.grid(True, which="major", axis="both", alpha=0.3, linestyle=":")
        if has_v:
            axv = axp.twinx()
            axv.plot(ts, vram, color="#d95f5f", lw=1.2, label="VRAM")
            axv.set_ylabel("VRAM [MB]", color="#d95f5f")
        # Draw the intervals and remember, per node type (class), the first one so
        # we can label the trace with the node-type name.
        first_of_class = {}
        for i, n in enumerate(nodes):
            cls = n["class"]
            col = CAT_COLOR[report_category(cls)]
            for a, b in n["intervals"]:
                axg.hlines(i, a - t0, b - t0, color=col, lw=2.2)
            start = min(a for a, b in n["intervals"]) - t0
            if cls not in first_of_class or start < first_of_class[cls][1]:
                first_of_class[cls] = (i, start)
        for cls, (i, start) in first_of_class.items():
            axg.text(start, i - 0.4, _short(cls), fontsize=5, va="bottom", ha="left",
                     color=CAT_COLOR[report_category(cls)], clip_on=True)
        axg.set_ylim(-1, max(1, nrows)); axg.invert_yaxis()
        axg.set_yticks([]); axg.set_ylabel(f"{nrows} nodes")
        axg.set_xlabel("time since run start [s]")
        axg.grid(True, which="major", axis="x", alpha=0.3, linestyle=":")
        # category legend
        from matplotlib.patches import Patch
        axg.legend(handles=[Patch(color=CAT_COLOR[c], label=CAT_LABEL[c])
                            for c in ("dnn", "sp_nondnn", "io", "other")],
                   fontsize=6, ncol=4, loc="upper right")
        fig.tight_layout()
        p = figdir / f"mem_{stage}.png"
        fig.savefig(p, dpi=120)
        plt.close(fig)
        figs[stage] = p.name
    return figs


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
    base = Path(bench_path).parent

    figs = make_figures(bench["reports"], outdir)
    memfigs = make_memory_figures(bench, base, outdir)
    graphs = {}
    if with_graphs:
        for s in bench["reports"]:
            g = render_flow_graph(bench["reports"][s], outdir / f"graph_{s}.png")
            if g:
                graphs[s] = g
    ctx = _summary_context(bench)
    ctx["memfigs"] = memfigs

    written = []
    if "md" in formats:
        p = outdir / "summary.md"; p.write_text(emit_markdown(bench, ctx, figs, graphs)); written.append(p)
    if "html" in formats:
        p = outdir / "summary.html"; p.write_text(emit_html(bench, ctx, figs, graphs)); written.append(p)
    if "tex" in formats:
        p = outdir / "summary.tex"
        p.write_text(emit_latex(bench, ctx, figs, graphs, fragment=tex_fragment, figpre=figpre))
        written.append(p)
    return {"outdir": str(outdir), "figures": figs, "graphs": graphs, "memfigs": memfigs,
            "written": [str(p) for p in written]}


# ---------------------------------------------------------------------------
# Grid-level report (across all grid points from a grid-index.json).
# ---------------------------------------------------------------------------
def _report_cat4(report):
    """Return {DNN, sp_other, other, total} wall seconds from a per-node report.

    NOTE: these are SUMS of per-node wall spans -- a work/utilization measure, not
    job latency.  Under TBB concurrency each node's wall span inflates (GPU
    serialization + CPU contention), so this sum RISES with wc-cores even as the
    real elapsed latency falls.  For latency use _report_latency() instead.
    """
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


def _report_latency(report):
    """Elapsed job latency (seconds) for a stage report -- the meaningful 'how
    long did it take' number, distinct from the sum-of-node-walls work measure.

    Prefers the timeline execution span (max end - min start over all node
    intervals), then the log 'execution' phase, then the whole-process
    wall_clock.  Returns None if none are available.
    """
    # 1) Timeline span: the true SP execution wall (excludes startup/model load).
    best = None
    for r in (report.get("runs") or []):
        tl = r.get("timeline")
        if tl and Path(tl).exists():
            try:
                j = json.loads(Path(tl).read_text())
                ivs = [(a, b) for n in j.get("nodes", [])
                       for a, b in n.get("intervals", []) if b > a]
                if ivs:
                    span = max(b for _, b in ivs) - min(a for a, _ in ivs)
                    best = span if best is None else min(best, span)
            except Exception:
                pass
    if best is not None:
        return best
    # 2) 'execution' phase from the log timestamps (engine-agnostic).
    try:
        ex = report["phases"]["execution"]["mean"]
        if ex is not None:
            return ex
    except (TypeError, KeyError):
        pass
    # 3) Whole-process wall clock (includes startup + model load).
    wc = [r.get("wall_clock") for r in (report.get("runs") or []) if r.get("wall_clock")]
    return min(wc) if wc else None


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
            "cat": {}, "lat": {}, "mem": {},
        }
        for s in ("osp", "spng"):
            ok = j[s].get("outcome") == "ok"
            rec["cat"][s] = _report_cat4(j[s]) if ok else None
            rec["lat"][s] = _report_latency(j[s]) if ok else None
            rec["mem"][s] = j[s].get("memory")   # None for pre-mem-profiling runs
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
    x = range(len(labels))

    # (1) Elapsed job latency -- the meaningful "how long did it take" metric.
    #     This is what shows real speedup: it falls as wc-cores rise.
    osp_l = [c["lat"]["osp"] or 0.0 for c in ok]
    spng_l = [c["lat"]["spng"] or 0.0 for c in ok]
    fig, ax = plt.subplots(figsize=(max(6, 0.3 * len(labels) + 2), 4))
    ax.bar([i - 0.2 for i in x], osp_l, width=0.4, label="OSP", color="#5f8fbf")
    ax.bar([i + 0.2 for i in x], spng_l, width=0.4, label="SPNG", color="#d95f5f")
    ax.set_xticks(list(x)); ax.set_xticklabels(labels, rotation=90, fontsize=6)
    ax.set_ylabel("elapsed wall time [s]"); ax.set_yscale("log"); ax.legend()
    ax.set_title("Elapsed execution wall time per grid point (log) -- lower is faster")
    fig.tight_layout()
    p = figdir / "bars_latency.png"; fig.savefig(p, dpi=120); plt.close(fig)
    figs["latency"] = p.name

    # (2) Sum of per-node wall spans -- a WORK/utilization measure (NOT latency);
    #     it rises with wc-cores as concurrent nodes' spans inflate.  Kept for
    #     comparison but clearly distinguished from latency above.
    osp_t = [c["cat"]["osp"]["total"] if c["cat"]["osp"] else 0.0 for c in ok]
    spng_t = [c["cat"]["spng"]["total"] if c["cat"]["spng"] else 0.0 for c in ok]
    fig, ax = plt.subplots(figsize=(max(6, 0.3 * len(labels) + 2), 4))
    ax.bar([i - 0.2 for i in x], osp_t, width=0.4, label="OSP", color="#5f8fbf")
    ax.bar([i + 0.2 for i in x], spng_t, width=0.4, label="SPNG", color="#d95f5f")
    ax.set_xticks(list(x)); ax.set_xticklabels(labels, rotation=90, fontsize=6)
    ax.set_ylabel("sum of per-node wall [s]"); ax.set_yscale("log"); ax.legend()
    ax.set_title("Total node-work per grid point (sum of per-node wall; utilization, not latency)")
    fig.tight_layout()
    p = figdir / "bars_total.png"; fig.savefig(p, dpi=120); plt.close(fig)
    figs["total"] = p.name

    rl, ratios = [], []
    for c in ok:
        if c["lat"]["osp"] and c["lat"]["spng"]:
            rl.append(lab(c)); ratios.append(c["lat"]["spng"] / c["lat"]["osp"])
    if rl:
        fig, ax = plt.subplots(figsize=(max(6, 0.3 * len(rl) + 2), 3.5))
        ax.bar(range(len(rl)), ratios, color="#7a5fbf")
        ax.axhline(1.0, color="k", lw=0.6)
        ax.set_xticks(range(len(rl))); ax.set_xticklabels(rl, rotation=90, fontsize=6)
        ax.set_ylabel("SPNG / OSP elapsed wall"); ax.set_title("SPNG-over-OSP elapsed-latency ratio per grid point")
        fig.tight_layout()
        p = figdir / "bars_ratio.png"; fig.savefig(p, dpi=120); plt.close(fig)
        figs["ratio"] = p.name
    return figs


def make_grid_membars(cells, figdir):
    """Peak-RSS (and peak-VRAM) per grid point, OSP vs SPNG, split into the
    DNN / sp_other / other categories (per-category peak; overlapping)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    figs = {}
    figdir = Path(figdir)
    ok = [c for c in cells if (c.get("mem", {}).get("osp") or c.get("mem", {}).get("spng"))]
    ok = sorted(ok, key=lambda c: (c["gpu"], c["torch"], c["scheme"], c["wc"]))
    if not ok:
        return figs

    def lab(c):
        b = GPU_MODE_LABEL.get(c["gpu"], f"{c['gpu']}gpu").replace("CPU + ", "").replace("CPU only", "cpu")
        return f"{b}/{c['scheme'][:4] if c['gpu']>=2 else 'omp'+str(c['torch'])}/wc{c['wc']}"

    labels = [lab(c) for c in ok]
    cats = [("dnn", "dnn"), ("sp_other", "sp_nondnn"), ("other", "other")]

    def val(c, s, catkey, which):
        m = (c.get("mem", {}).get(s) or {}).get("by_category", {}).get(catkey, {})
        return (m.get(which) or 0) / 1e6   # MB

    for which, title, fname in (("peak_rss", "Peak RSS", "peakram"),
                                ("peak_vram", "Peak VRAM", "peakvram")):
        allvals = [val(c, s, ck, which) for c in ok for s in ("osp", "spng") for ck, _ in cats]
        if not any(v > 0 for v in allvals):
            continue
        n = len(ok)
        w = 0.13
        # Extra whitespace between grid points: space the case centers by `pitch`
        # (> the ~0.78 cluster width) so each point's two triplets stay grouped and
        # clearly separated from its neighbours.
        pitch = 1.5
        x = [i * pitch for i in range(n)]
        offs = {"osp": [-3 * w, -2 * w, -1 * w], "spng": [1 * w, 2 * w, 3 * w]}
        fig, ax = plt.subplots(figsize=(max(6, 0.42 * pitch * n + 2), 4.2))
        for s in ("osp", "spng"):
            for k, (ck, colkey) in enumerate(cats):
                ax.bar([xi + offs[s][k] for xi in x], [val(c, s, ck, which) for c in ok],
                       width=w, color=CAT_COLOR[colkey],
                       hatch=("" if s == "osp" else "//"), edgecolor="white", linewidth=0.2)
        ax.set_xticks(x); ax.set_xticklabels(labels, rotation=90, fontsize=6)
        ax.set_xlim(-pitch * 0.5, (n - 0.5) * pitch)
        ax.set_ylabel(f"{title} [MB]")
        ax.set_title(f"{title} per grid point (left triplet = OSP, right hatched = SPNG)")
        ax.legend(handles=[Patch(color=CAT_COLOR[c], label=CAT_LABEL[c]) for _, c in cats]
                  + [Patch(facecolor="#ccc", hatch="//", label="SPNG (hatched)")],
                  fontsize=7, ncol=4)
        fig.tight_layout()
        p = figdir / f"bars_{fname}.png"; fig.savefig(p, dpi=120); plt.close(fig)
        figs[fname] = p.name
    return figs


# ---------------------------------------------------------------------------
# Process-parallel axis reporting: read the proc-*.json records and chart the
# aggregate wall time and the summed RAM/VRAM footprint vs the number of
# simultaneous single-APA jobs, per stage and device.
# ---------------------------------------------------------------------------
def proc_records(base):
    """Load every proc-*.json (schema spngbench-proc/1) in a directory."""
    base = Path(base)
    recs = []
    for pf in sorted(base.glob("proc-*.json")):
        try:
            j = json.loads(pf.read_text())
        except Exception:
            continue
        if not j.get("schema", "").startswith("spngbench-proc"):
            continue
        recs.append(j)
    return recs


def _proc_series(recs):
    """Group proc records into {(which, device): [(nproc, wall, rss, vram, outcome)...]}."""
    ser = {}
    for j in recs:
        m = j["meta"]; mem = j.get("memory") or {}
        key = (m["which"], m["device"])
        ser.setdefault(key, []).append((
            m["nproc"], j.get("wall_clock"),
            mem.get("peak_rss_sum") or 0, mem.get("peak_vram_sum") or 0,
            j.get("outcome")))
    for key in ser:
        ser[key] = sorted(ser[key])
    return ser


def make_proc_figures(recs, figdir):
    """Wall / peak-RAM / peak-VRAM vs number of simultaneous single-APA jobs."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figs = {}
    figdir = Path(figdir)
    ser = _proc_series(recs)
    if not ser:
        return figs

    STY = {"osp": "-o", "spng": "-s"}
    COL = {"cpu": "#5f8fbf", "gpu0": "#d95f5f", "gpu": "#d95f5f"}

    def lab(which, dev):
        return f"{which.upper()} / {dev}"

    def plot(idx, ylabel, title, fname, scale=1.0, gpu_only=False):
        fig, ax = plt.subplots(figsize=(6, 4))
        drew = False
        for (which, dev), pts in sorted(ser.items()):
            if gpu_only and not dev.startswith("gpu"):
                continue
            xs = [p[0] for p in pts]
            ys = [p[idx] / scale if p[idx] else None for p in pts]
            if not any(v for v in ys):
                continue
            ax.plot(xs, ys, STY.get(which, "-o"), color=COL.get(dev, "#7a5fbf"),
                    label=lab(which, dev))
            drew = True
        if not drew:
            plt.close(fig)
            return
        ax.set_xlabel("number of simultaneous single-APA jobs (nproc)")
        ax.set_ylabel(ylabel); ax.set_title(title)
        xs_all = sorted({p[0] for pts in ser.values() for p in pts})
        ax.set_xticks(xs_all)
        ax.grid(True, alpha=0.3); ax.legend(fontsize=8)
        fig.tight_layout()
        p = figdir / fname; fig.savefig(p, dpi=120); plt.close(fig)
        figs[fname.split(".")[0].replace("proc_", "")] = p.name

    plot(1, "elapsed wall time [s]",
         "Process-parallel: wall time vs simultaneous jobs", "proc_wall.png")
    plot(2, "summed peak RSS [GB]",
         "Process-parallel: total RAM vs simultaneous jobs", "proc_ram.png", scale=1e9)
    plot(3, "summed peak VRAM [GB]",
         "Process-parallel: total VRAM (shared GPU) vs simultaneous jobs",
         "proc_vram.png", scale=1e9, gpu_only=True)
    return figs


def _proc_table_rows(recs):
    """Rows (which, device, nproc, outcome, wall, RAM_GB, VRAM_GB) for the table."""
    rows = []
    for (which, dev), pts in sorted(_proc_series(recs).items()):
        for nproc, wall, rss, vram, outcome in pts:
            rows.append((which, dev, nproc, outcome,
                         wall, (rss or 0) / 1e9, (vram or 0) / 1e9))
    return rows


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
    bars.update(make_grid_membars(cells, outdir))
    job_graphs, inputs = make_job_graphs(grid, cells, outdir)

    # Process-parallel axis (proc-*.json), if any were run into this directory.
    precs = proc_records(base)
    proc_figs = make_proc_figures(precs, outdir) if precs else {}
    proc_rows = _proc_table_rows(precs) if precs else []

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
           "ncells": len(cells), "jobs": job_graphs, "inputs": inputs,
           "proc_figs": proc_figs, "proc_rows": proc_rows}
    written = []
    if "md" in formats:
        p = outdir / "grid-summary.md"; p.write_text(_emit_grid_md(ctx)); written.append(p)
    if "html" in formats:
        p = outdir / "grid-summary.html"; p.write_text(_emit_grid_html(ctx)); written.append(p)
    if "tex" in formats:
        p = outdir / "grid-summary.tex"; p.write_text(_emit_grid_tex(ctx)); written.append(p)
    return {"outdir": str(outdir), "written": [str(p) for p in written],
            "npoints": len(points), "matrices": mats, "bars": bars,
            "proc_figs": proc_figs}


_MATRIX_BLURB = (
    "The outer product of device modes (rows) and node categories (columns).  "
    "Each cell is a grid point split into two sub-pixels: **left = OSP, "
    "right = SPNG**, coloured by that category's wall time on a per-column log "
    "scale.  Not-tested cells are black; crashed sub-pixels white.\n\n"
    "The category columns sum per-node wall time by kind: **DNN** = the neural-net "
    "forward (`Pytorch::DNNROIFinding` for OSP, `SPNG::TensorForward` for SPNG); "
    "**sp_other** = all other signal-processing nodes (decon, filters, ROI); "
    "**other** = everything else, chiefly frame I/O; **total** = the sum of all "
    "three.  (These are node-work sums, not latency — see Trends.)\n\n"
    "Within each cell the axes are the intra-job parallelism: **x = wire-cell "
    "cores** (`TbbDataFlowGraph.max_threads`, the number of graph nodes that may "
    "run at once) and **y = omp** (`OMP_NUM_THREADS`, torch intra-op threads) for "
    "CPU / single-GPU rows, or the shard scheme for multi-GPU rows.")


def _blurb_html(text):
    """Render the light markdown used in blurbs (**bold**, `code`, blank-line
    paragraphs) as HTML paragraphs."""
    import re as _re
    paras = text.split("\n\n")
    out = []
    for para in paras:
        p = (para.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))
        p = _re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", p)
        p = _re.sub(r"`(.+?)`", r"<code>\1</code>", p)
        out.append(f"<p>{p}</p>")
    return "".join(out)


def _blurb_tex(text):
    """Render the light markdown used in blurbs as LaTeX (\\textbf, \\texttt,
    paragraph breaks)."""
    import re as _re
    paras = text.split("\n\n")
    out = []
    for para in paras:
        p = para.replace("_", r"\_").replace("&", r"\&").replace("%", r"\%")
        p = _re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", p)
        p = _re.sub(r"`(.+?)`", r"\\texttt{\1}", p)
        out.append(p)
    return "\n\n".join(out)


def _emit_grid_md(ctx):
    m = ctx["grid"]["meta"]
    L = [f"# spngbench grid summary — {m.get('detname','?')}\n",
         f"host **{m.get('host','?')}** ({m.get('host_ngpu','?')} GPU), engine "
         f"**{m.get('engine','?')}**, {ctx['ncells']} grid points.\n",
         "\n## Contents\n",
         "- [Jobs](#jobs)",
         "- [Grid matrix](#grid-matrix)",
         "- [Trends](#trends)",
         *(["- [Process parallelism](#process-parallelism)"] if ctx.get("proc_rows") else []),
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
    L.append("*`latency` is elapsed execution wall time (job duration) — the metric "
             "that shows real speedup: it falls as wire-cell cores rise. `total` is the "
             "sum of per-node wall spans, a work/utilization measure that RISES with "
             "cores as concurrent nodes' spans inflate; it is not latency.*\n")
    for k in ("latency", "total", "ratio", "peakram", "peakvram"):
        if ctx["bars"].get(k):
            L.append(f"\n![{k}]({ctx['bars'][k]})\n")

    if ctx.get("proc_rows"):
        L.append("\n## Process parallelism\n")
        L.append("A separate axis: instead of one job over N APAs, run *N independent "
                 "single-APA jobs at once* (wc1/omp1), on the CPU or **one shared GPU**. "
                 "Wall time measures throughput at fixed per-job cost; the summed "
                 "RAM/VRAM is the simultaneous footprint that decides how many jobs "
                 "co-fit on a host or a shared GPU (e.g. 24 GB RTX-4090 vs 48 GB L40S).\n")
        for k in ("wall", "ram", "vram"):
            if ctx["proc_figs"].get(k):
                L.append(f"\n![proc {k}]({ctx['proc_figs'][k]})\n")
        L.append("\n| stage | device | nproc | outcome | wall [s] | RAM [GB] | VRAM [GB] |")
        L.append("|---|---|--:|---|--:|--:|--:|")
        for which, dev, nproc, outcome, wall, ram, vram in ctx["proc_rows"]:
            L.append(f"| {which} | {dev} | {nproc} | {outcome} | {_fmt(wall)} | "
                     f"{_fmt(ram)} | {_fmt(vram)} |")

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
         *(["<li><a href='#proc'>Process parallelism</a></li>"] if ctx.get("proc_rows") else []),
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
    H.append(_blurb_html(_MATRIX_BLURB))
    if ctx["mats"].get("matrix"):
        H.append(f"<p><img src='{ctx['mats']['matrix']}'></p>")
    H.append("<h2 id='trends'>Trends</h2>")
    H.append("<p><em><code>latency</code> is elapsed execution wall time (job duration) "
             "&mdash; the metric that shows real speedup: it falls as wire-cell cores "
             "rise. <code>total</code> is the sum of per-node wall spans, a "
             "work/utilization measure that RISES with cores as concurrent nodes' spans "
             "inflate; it is not latency.</em></p>")
    for k in ("latency", "total", "ratio", "peakram", "peakvram"):
        if ctx["bars"].get(k):
            H.append(f"<p><img src='{ctx['bars'][k]}'></p>")
    if ctx.get("proc_rows"):
        H.append("<h2 id='proc'>Process parallelism</h2>")
        H.append("<p>A separate axis: run <em>N independent single-APA jobs at once</em> "
                 "(wc1/omp1) on the CPU or <b>one shared GPU</b>.  Wall time measures "
                 "throughput at fixed per-job cost; the summed RAM/VRAM is the "
                 "simultaneous footprint that decides how many jobs co-fit on a host or a "
                 "shared GPU (e.g. 24 GB RTX-4090 vs 48 GB L40S).</p>")
        for k in ("wall", "ram", "vram"):
            if ctx["proc_figs"].get(k):
                H.append(f"<p><img src='{ctx['proc_figs'][k]}'></p>")
        H.append("<table><tr><th>stage</th><th>device</th><th>nproc</th><th>outcome</th>"
                 "<th>wall [s]</th><th>RAM [GB]</th><th>VRAM [GB]</th></tr>")
        for which, dev, nproc, outcome, wall, ram, vram in ctx["proc_rows"]:
            H.append(f"<tr><td>{esc(which)}</td><td>{esc(dev)}</td><td>{nproc}</td>"
                     f"<td>{esc(outcome)}</td><td>{_fmt(wall)}</td><td>{_fmt(ram)}</td>"
                     f"<td>{_fmt(vram)}</td></tr>")
        H.append("</table>")
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
         _blurb_tex(_MATRIX_BLURB)]
    if ctx["mats"].get("matrix"):
        T.append(r"\begin{center}\includegraphics[width=\textwidth]{" + ctx["mats"]["matrix"] + r"}\end{center}")
    T.append(r"\section{Trends}")
    T.append(r"\emph{\texttt{latency} is elapsed execution wall time (job duration) --- "
             r"the metric that shows real speedup: it falls as wire-cell cores rise. "
             r"\texttt{total} is the sum of per-node wall spans, a work/utilization "
             r"measure that RISES with cores as concurrent nodes' spans inflate; it is "
             r"not latency.}" + "\n")
    for k in ("latency", "total", "ratio", "peakram", "peakvram"):
        if ctx["bars"].get(k):
            T.append(r"\begin{center}\includegraphics[width=\textwidth]{" + ctx["bars"][k] + r"}\end{center}")
    if ctx.get("proc_rows"):
        T.append(r"\clearpage\section{Process parallelism}")
        T.append(r"A separate axis: run \emph{N independent single-APA jobs at once} "
                 r"(wc1/omp1) on the CPU or \textbf{one shared GPU}.  Wall time measures "
                 r"throughput at fixed per-job cost; the summed RAM/VRAM is the "
                 r"simultaneous footprint that decides how many jobs co-fit on a host or "
                 r"a shared GPU (e.g.\ 24\,GB RTX-4090 vs 48\,GB L40S).")
        for k in ("wall", "ram", "vram"):
            if ctx["proc_figs"].get(k):
                T.append(r"\begin{center}\includegraphics[width=0.8\textwidth]{"
                         + ctx["proc_figs"][k] + r"}\end{center}")
        T.append(r"\small\begin{tabular}{llrlrrr}\toprule")
        T.append(r"stage & device & nproc & outcome & wall [s] & RAM [GB] & VRAM [GB] \\\midrule")
        for which, dev, nproc, outcome, wall, ram, vram in ctx["proc_rows"]:
            T.append(f"{esc(which)} & {esc(dev)} & {nproc} & {esc(outcome)} & "
                     f"{_fmt(wall)} & {_fmt(ram)} & {_fmt(vram)} " + r"\\")
        T.append(r"\bottomrule\end{tabular}\normalsize")
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
    pg.add_argument("--mode", choices=["cpu", "gpu", "shard", "two-gpu", "proc", "all"],
                    default="cpu",
                    help="which cell family to run (default cpu); 'shard' = multi-GPU, "
                         "'proc' = process-parallel (N single-APA jobs at once)")
    pg.add_argument("--which", choices=["osp", "spng", "both"], default="both")
    pg.add_argument("--wc-cores", type=int, nargs="+", default=None,
                    help="subset of wire-cell core counts")
    pg.add_argument("--torch-cores", type=int, nargs="+", default=None,
                    help="subset of torch core counts (cpu mode)")
    pg.add_argument("--nproc", type=int, nargs="+", default=None,
                    help="subset of process counts (proc mode)")
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
        cells = select_cells(cfg, args.mode, wc_sel=args.wc_cores,
                             torch_sel=args.torch_cores, nproc_sel=args.nproc)
        if not cells:
            log("no cells selected")
            return 1
        log(f"grid: mode={args.mode} which={which_list} cells={len(cells)} "
            f"outdir={cfg['outdir']}")
        idx = run_grid(cfg, cells, which_list, dry_run=args.dry_run)
        # Compact human-facing summary to stdout.
        for c in idx.get("proc", []):
            if c.get("skipped"):
                print(f"  proc {c['device']} n{c['nproc']}: SKIPPED ({c['reason']})")
                continue
            rss = c.get("peak_rss_sum"); vram = c.get("peak_vram_sum")
            rss_s = f"{rss/1e9:.2f}GB" if rss else "-"
            vram_s = f"{vram/1e9:.2f}GB" if vram else "-"
            print(f"  proc {c['which']} {c['device']} n{c['nproc']}: "
                  f"{c['outcome']} wall={c['wall_clock']}s RAM={rss_s} VRAM={vram_s}")
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
