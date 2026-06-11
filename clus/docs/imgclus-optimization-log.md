# PDHD / PDVD imaging + clustering optimization log

Running log of the memory-footprint and wall-time optimization work that
follows the baseline measured in [imgclus-resource-profile.md](imgclus-resource-profile.md).

**Ground rule:** every change applied here must produce *functionally
identical* outputs, verified by content-hash A/B comparison on a fixed event
set. Ideas that would change results (even benignly) are **not** applied;
they are collected in the "Result-changing ideas" section at the end for
review.

## A/B methodology

- Harness: `wcp-porting-img/abtest/` (`run_events.sh`, `ab_compare.sh`,
  `hash_archive.py`, `timecmd.py`). Snapshots of all per-event output
  archives + wall/peak-RSS metadata land in `abtest/snap/<label>/`.
- Output archives embed wall-clock timestamps (tar member `mtime=time(0)` in
  `util/inc/WireCellUtil/custard/custard.hpp`; zip likewise), so raw archive
  bytes are never reproducible. Comparison hashes **member payloads**
  (sha256 over name+payload, members sorted by name).
- Peak RSS via `getrusage(RUSAGE_CHILDREN)` (exact high-water mark over the
  full child tree); wall via wrapper timer. Events run 6-way concurrent on a
  64-core host, `OMP_NUM_THREADS=MKL_NUM_THREADS=1`.

### Fixed event set

| Tag | Detector | Run / Evt | Why |
|---|---|---|---|
| hd-typ  | PDHD | 027409 / 0  | typical |
| hd-typ2 | PDHD | 027980 / 3  | second typical |
| hd-busy | PDHD | 028084 / 18 | busy tail (981 s / 7.2 GB in profile) |
| hd-max  | PDHD | 027305 / 0  | extreme (1501 s / 10.1 GB in profile) |
| vd-typ  | PDVD | 039349 / 0  | typical |
| vd-busy | PDVD | 039252 / 5  | PDVD tail |

A 22-event spot-check list spanning every profiled run is in
`abtest/spot_events.txt`; content hashes of the pre-optimization outputs are
stashed in `abtest/snap/preexist/hashes.txt`.

## Baseline (this work, current master add080f5, snapshot `base1`)

6-way concurrent, wall = stage wall-clock (s), RSS = exact peak over the
process tree (MB).

| Event | img wall | img RSS | clus wall | clus RSS |
|---|---|---|---|---|
| hd-typ  027409/0  | 66   | 703   | 28  | 704 |
| hd-typ2 027980/3  | 70   | 691   | 31  | 769 |
| hd-busy 028084/18 | 656  | 7180  | 338 | 2833 |
| hd-max  027305/0  | 1028 | 10113 | 530 | 3284 |
| vd-typ  039349/0  | 140  | 441   | 26  | 745 |
| vd-busy 039252/5  | 327  | 1918  | 145 | 1879 |

## Determinism findings

- `base1` vs `base2` (two back-to-back full reruns, same binary, same
  inputs): **all 138 output archives byte-identical at member level** for
  all 6 events, both stages, including the all-APA clustering outputs. The
  imaging+clustering chain is fully deterministic on PDHD/PDVD, so the
  byte-identity gate applies to every output file with no exclusions.
- Content hashes of outputs that pre-dated this work (stash
  `snap/preexist`) partially disagree with `base1` (notably masked-fork
  imaging archives). Since run-to-run determinism is proven, those
  disagreements reflect input/config state during the original profile
  batch (e.g. PDVD `sp-frames` symlinks re-pointed mid-batch), not
  nondeterminism. The preexist stash is therefore NOT used as a comparison
  baseline; the spot-check baseline is regenerated with current code.

## Optimization entries

### 1. Per-anode sequential imaging (`run_img_evt.sh -P`) — script-level, no C++ change

One `wire-cell` process per anode instead of all anodes in one process.
Per-anode configs are pre-compiled in parallel with `wcsonnet` (also removes
the single large gojsonnet compile from the critical path); the per-anode
processes then run sequentially. Implemented in
`wcp-porting-img/{pdhd,pdvd}/run_img_evt.sh` (wcp-porting-validation repo).

**A/B verdict: PASS** — all imaging archives and all downstream clustering
outputs byte-identical to `base1` on all 6 events (snapshot `peranode`).

| Event | img wall (s) | img RSS (MB) |
|---|---|---|
| hd-typ  | 66 → **50**  | 702 → 657 |
| hd-typ2 | 70 → **53**  | 691 → 663 |
| hd-busy | 656 → **595** | 7180 → 7152 |
| hd-max  | 1028 → **996** | 10113 → 10090 |
| vd-typ  | 140 → **50** | 440 → 396 |
| vd-busy | 327 → **183** | 1918 → 1876 |

Take-aways:
- **Wall improves everywhere** (PDVD typical 2.8x: the 8-anode jsonnet
  compile + geometry init dominated). Clustering unchanged.
- **The RSS tail is NOT a sum-of-anodes effect.** On hd-max the peak barely
  moved: one anode alone (anode0, 809 s of the 996 s) holds the ~10 GB.
  Per-anode wall split on hd-max: 809/135/18/31 s; on hd-busy: 106/12/465/8 s.
  The tail memory therefore lives *inside* a single anode's imaging pipeline
  (cluster graph + per-pass rebuilds), and is attacked by the in-process
  work below, not by process splitting.
- `-P` is opt-in; default behavior of the scripts is unchanged.

## Result-changing ideas (NOT applied — for review)

Collected as encountered; see end of file.
