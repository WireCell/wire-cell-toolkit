# TrackFitting `update_association`: an exact competitor index for the exclusion scan

Status: **byte-identical** (perf change, always on; gates below). It changes only how the competing segments are found when
`fit_exclusion` is on; every keep / drop decision is the one the full scan makes. No parameter, no config change.

## Repro

```bash
# unit tests (3 new cases in the existing file; new symbols: link against this build's own libraries, see "Verification")
./build/clus/wcdoctest-clus -tc="*update_association*"
# the slow event (pr_ml doc 87b; wcp-porting-validation pr_ml/docs/87_results/scripts):
#   run_d87_ex_event.sh <tag> <stored lane event dir> <pin dir>        one stored lane event on a given library pin
#   d87_ex_cmp.py <run A> <run B>                                      outputs compared (archives by member hash, ROOT by value)
```

## Symptom

On the pr_ml doc 87b tuning-event fit lane, Wire-Cell's own pattern recognition (arm `tradF`) took 4.4-4.5 h on one FD-HD event
(`nue_b11094_e0`), against 17 s for the median event. With the `spd_solver` knob on, the solve was no longer the cause.

## Root cause

A CPU profile of the event (16,344 s under the profiler; 4,084,163 samples at 250 Hz) put **79.3 %** of all samples in the
three competitor loops of `TrackFitting::update_association`, one per plane:
- In each `form_map_graph` pass, for each segment, every 2-D pixel it claims is tested. If the segment's own distance
  (`exclusion_closest_2d_dis`) is 0.3 cm or more, every other segment of the cluster is asked in turn, one k-d tree query each,
  until one is as close or closer (`<=`, early break).
- The cost is pixels × segments per pass. The passes repeat because the vertex-examination steps call `do_multi_tracking` again
  and again. The event's main cluster has 414 segments.
- Per-line profile view and the cumulative table: pr_ml doc 87b, 2026-11-23 12:03 EST note (`87_results/read/prof/`).

## Why it hid

The scan is exact and cheap on ordinary clusters (tens of segments). Only a cluster with hundreds of segments, revisited by many
passes, shows it, and then only as time.

## Fix

One index per `form_map_graph` call (RAII scope; only when `fit_exclusion` is on and the debug dump is off):
- **A grid per (plane, face, apa):** 1 cm cells; each occupied cell lists the segments with a stored 2-D point in it (the points
  of each segment's chosen cloud, `fit` else `main`, as `exclusion_closest_2d_dis` chooses).
- **Query:** the cells covering the segment's own distance, enlarged by × (1 + 1e-9) + 1e-6 mm, give the candidates. Each
  candidate is asked with the unchanged `exclusion_closest_2d_dis`, with `<=` and the early break, in `all_segments` order.
- **Kept current:** a segment's fit cloud is replaced after its own loop (`set_fit_associate_vec`). So each `update_association`
  call refreshes every segment's chosen cloud (object and size), and a grid re-indexes exactly the segments that changed before
  it is read. The grid holds the clouds it indexed alive, so a freed cloud's address cannot come back as the same cloud.
- **The full scan still runs when:** the own distance is 200 mm or more (the 1e9 sentinel among them); the chosen clouds'
  `wpid_params` differ; `all_segments` is not the vector the index was built on; a 2-D tree is not in owned storage; the debug
  dump is on.
- **Two read-only helpers on `DynamicPointCloud`:** `project_2d` (the projection `get_closest_2d_dis` makes) and `kd2d_find`
  (the existing tree, or nullptr; `kd2d()` would create an empty one).

**Why the decisions are identical:** a competitor the grid leaves out has every stored point outside the enlarged radius, so its
distance from the pixel exceeds the own distance and it could not have dropped the pixel. Every candidate is decided by the
unchanged function, in the same order, with the same tie rule.

No new dependency; no plugin-to-plugin link.

## Verification

**Unit tests** (`clus/test/doctest_update_association.cxx`, 3 new cases):
- 40 random segments × 726 cells: decisions with and without the index are equal; the index was used; some cells are kept and some
  dropped.
- An exact tie still drops the cell. Mixed `wpid_params` fall back to the full scan.
- With the 4 existing cases, 7 cases and 74 assertions pass. The whole clus suite passes, **553 of 553**.
- **Negative control:** halving the search radius makes the random case fail.
- The test calls new symbols, so `wcdoctest-clus` must link against this build's libraries ahead of an older installed
  `local/lib` (waf's test link finds the installed library first); the run above linked by hand against the build tree.

**Byte-identical gates** (2026-11-23; side A = the libraries before the change, side B = after; same binaries, only
`LD_LIBRARY_PATH` differs; scripts and logs in wcp-porting-validation `pr_ml/docs/87_results/read/ex_gates/`):

| gate | labels | result |
|---|---|---|
| uBooNE qlport two-gate | `d87exA_ub` / `d87exB_ub` | zips **35 / 35** content-identical; tagger identical on **35 / 35** |
| PDHD / PDVD `abtest` clus-only | `d87exA_pd` / `d87exB_pd` | **OVERALL: PASS** |
| SBND PR chain, 60 events (`fit_exclusion` on in production) | `work-d87ex{A,B}-mcp60` | **120 / 120** archives byte-identical |
| pr_ml dev PR lane null, 8 events (`fit_exclusion` on) | `n2_pin4ex` / `n2_pin5ex` | **8 / 8** equal by value |
| the slow event, `nue_b11094_e0` | before (unprofiled and profiled) / after | `mabc-pr.zip`, `pr-pctree.tar.gz` equal by member hash; `tracking-pr.root` 1,360 branches, 0 different |

SBND, the lane and the slow event run the new path; uBooNE and PDHD / PDVD check that the exclusion-off path is untouched.

## Performance

`nue_b11094_e0`, one thread, the same `cfg.json` and TrackFitting file (`spd_solver` 1):

| run | wall | peak RSS |
|---|---|---|
| before, unprofiled | 16,317 s | 6.43 GB |
| before, under the CPU profiler | 16,344 s | 6.54 GB |
| **after** | **3,225 s** | 6.59 GB |

- **5.06× faster:** 13,092 s saved. The profile put 12,912 s in the scan.
- **Peak RSS +161 MB (+2.5 %)** against the unprofiled run, for the grid's cells.
- The two "before" runs differ in log level (debug with the solver census, against info under the profiler) and agree within 27 s.
- Details: pr_ml doc 87b, 2026-11-23 notes 12:03 (profile), 12:09 (registration), 12:35 (build), 12:53 (gates) and 13:35 (this
  event).

## Limits

- The full scan remains for the fallback cases above.
- Only one slow event was timed. On ordinary events the scan was already cheap, and the gates show no change in their outputs.
