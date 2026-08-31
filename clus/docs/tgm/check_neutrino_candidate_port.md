# Porting `check_neutrino_candidate()` into `TaggerCheckTGM` (two-TPC generalization)

Port of the prototype's Dijkstra path-topology neutrino veto,
`WCPPID::ToyFiducial::check_neutrino_candidate`
(`prototype_base/pid/src/Cosmic_tagger.h` lines 1677–1985), into
`clus/src/TaggerCheckTGM.cxx`, generalized from uBooNE's single TPC to the
multi-volume (SBND two-TPC) toolkit geometry.  This closes the v1 gap where
every in-beam-window bundle was conservatively exempt from the TGM tag: the
prototype gates those tags on `!check_neutrino_candidate(...)`, which was
unported, so v1 never tagged them (see the file-header note in
`TaggerCheckTGM.cxx` and `tgm_implementation_comparison.md`, flash-type row).

Behavior ships behind a new config knob **`check_neutrino_candidate`
(default `false`)**.  OFF = v1 behavior, proven byte-identical below.

## Repro

```bash
cd /nfs/data/1/xqian/toolkit-dev/toolkit
wcbuild     # ./wcb build --notests -p && ./wcb install --notests -p
./build/clus/wcdoctest-clus

# Motivating event: SBND MCP2025C reco1, run 18255 subrun 1 event 284349,
# bundle at flash 1.555 us (in the 0.2-2.2 us beam window).  Anode->downstream
# TPC0 crosser, CASE-A geometry (dbg: pe1 (-201.2,127.5,350.9) x-wall exit,
# pe2 (-87.1,54.4,500.9) z-wall exit, mid_inside true, len 202.1/426.3 cm,
# 202.1 > 0.45*426.3 = 191.8), yet TGM=false in v1 because of the beam guard.
SB=/nfs/data/1/xqian/toolkit-dev/wcp-porting-img/sbnd/sbnd_xin
cd $SB && SBND_INPUT_DIR=$PWD/input_files_reco1/extracted-mcp2025c-10evt \
  SBND_WORK_ROOT=$PWD/work-mcp10 ./run_nusel_evt.sh data -nucand 1
# (idx 1 = evt 284349; drop -nucand for the legacy conservative behavior)
```

## What the prototype function does

Given the two extreme points `wcp1`, `wcp2` of a TGM-candidate pair, walk the
in-cluster shortest path between them and declare "neutrino candidate" (=
suppress the TGM tag) when the path shows non-muon topology:

1. **Path** (proto 1678–1722): `Create_graph(ct_point_cloud)`;
   `dijkstra_shortest_paths(wcp1)`; `cal_shortest_path(wcp2)`.  Resample
   twice: `path_wcps_vec` keeps points > 0.5 cm apart (kink search);
   `path_wcps_vec1` caps spacing at 1 cm by interpolation (gap search).  NB
   the prototype re-reads `back()` after each interpolation push — the
   resulting shrinking-step pattern is reproduced verbatim.
2. **Gap veto** (proto 1725–1830): first decide `flag_2view_check` — start
   true, drop to false when the chord is quasi-parallel to a wire direction
   (U/V < 10°, W < 5° in the prolonged-signal drift-angle construction) or
   quasi-isochronous (|90°−∠(drift, chord)| < 5°).  Then per path sample:
   the counters reset when ≥2 planes (1 when `!flag_2view_check`) have CT
   points within 1 cm, or the point is in a dead region; otherwise accumulate
   `num_nth` (consecutive in-FV samples), `min_dis` (closest approach to an
   endpoint), `num_bad` (`!is_good_point`).  Fire when `num_nth > 7 &&
   min_dis < 25 cm && num_bad > 7`.
3. **Kink veto** (proto 1834–1983): slide along `path_wcps_vec` with ±5-point
   local directions (`dir1/dir2`), ±14-point PCA directions (`dir3/dir4` via
   `calc_PCA_dir`), and ±(3..12]-point mean directions (`dir5/dir6`);
   `cut1` = how many of the three direction pairs break by > 25°, `cut2` =
   how many difference vectors are > 5° off isochronous.  A run of
   `cut1>=3 && cut2>=2` samples fires when the legs to the endpoints open
   ≥ 25–60° (three-tier condition, precedence preserved; 1-sample tier needs
   > 5 cm legs) **and** the kink apex is inside the FV and not in a dead
   region (dead apex tolerated only when the opening ≥ 45°).

## Call sites (prototype ↔ toolkit)

The prototype guards its beam-flash (`main_flash->get_type()==2`) branches
with `!check_neutrino_candidate && temp_length > 0.45*length_limit`.  The
toolkit uses the beam-window test on `cluster_t0` in place of the flash type
(v1 design).  With the knob ON the three sites become prototype-faithful:

| Site | Prototype | Toolkit knob ON | Toolkit knob OFF (v1) |
|---|---|---|---|
| CASE A, `flag_check` true | 1421–1437: beam → veto + length arbitrate; else tag | same, on `in_beam_window` | in-beam → never tag |
| CASE A, `!flag_check_again` | 1466–1477: veto + length, **no beam gate** | same, all clusters | in-beam → never tag; out-of-beam → length only |
| CASE B, both ends explained | 1564–1581: beam → veto + length on the `at(0)` extremes; else tag | same | in-beam → never tag |

Note the middle row: with the knob ON the veto also arbitrates for
out-of-beam clusters, exactly as in the prototype.  On the 10-event
validation sample this changed no out-of-beam verdict.

The `flag_2view_check` parameter is `true` at every prototype call site (the
explicit `true` for `i==0&&k==1` equals the header default) — kept as the
default of the ported method.

## Reused toolkit ports (nothing re-invented)

| Prototype facility | Toolkit equivalent used |
|---|---|
| `Create_graph(ct_point_cloud)` + `dijkstra_shortest_paths` + `cal_shortest_path` | `cluster.graph_algorithms("ctpc", m_dv, m_pcts).shortest_path(i1, i2)` with `get_closest_point_index` / `indices_to_points` |
| `ToyCTPointCloud::get_closest_points(p, r, plane)` | `Grouping::get_closest_points(p_raw, r, apa, face, pind)` |
| `ToyCTPointCloud::is_good_point` (0.6 cm, ch_range 1, allowed_bad 1) | `Grouping::is_good_point(p_raw, apa, face)` — identical defaults |
| `ToyFiducial::inside_dead_region` | `FiducialUtils::inside_dead_region(p_raw, apa, face)` (minimal_views 2, the established STM translation) |
| `inside_fiducial_volume(p, offset_x)` | the component's own `inside_fv()` (IFiducial + `fv_tolerance` margins, offset_x = 0 post-switch_scope) — same test check_tgm itself uses |
| `PR3DCluster::calc_PCA_dir` | `Facade::calc_pca_dir` |
| hard-coded U/V/W directions (proto 1736–1739) | `Facade::compute_wireplane_params` per-(apa,face) wire directions |
| raw coordinates (single volume) | `IPCTransformSet::backward(p, cluster_t0, face, apa)` per point — the `FiducialUtils::check_signal_processing` pattern |

## Two-TPC generalizations (single-TPC prototype → multi-volume toolkit)

1. **Coordinate spaces.** The path, FV tests and kink geometry run in the
   T0-corrected default scope (what `check_tgm` already uses; `offset_x = 0`).
   Every CT-point-cloud / dead-region query resolves the sample's
   `(apa, face)` via `DetectorVolumes::contained_by` and backward-transforms
   into that volume's raw coordinates.
2. **Cathode gap.** A two-TPC crosser path passes through x ≈ 0 where no
   volume contains the samples.  Such points are treated like dead regions
   (counter reset): no hits can exist there, so counting them as "bad" would
   spuriously fire the gap veto on every genuine cathode crosser.  A kink
   apex outside every volume likewise counts as "in dead region" (suppresses
   the veto — conservative toward the prototype's cosmic verdict).
3. **Wire directions.** The prolonged-signal test evaluates the per-(apa,face)
   wire frames of *both* endpoints when they differ; either frame being
   quasi-parallel drops the 2-view requirement.  (For SBND the two faces'
   {U,V} direction sets coincide up to the swap, and the test is symmetric
   under it, so this matters only for future detectors.)
4. **Drift direction.** `drift_dir = (1,0,0)` is kept: every use is of the
   form |90° − ∠(drift, v)| or |dir.x|, both invariant under the drift-sign
   flip between the two volumes (same argument as `drift_dir_abs` in
   `check_tgm`).

## Deliberate deviations from the prototype text

- `dir3/dir5` (and `dir4/dir6`) are computed **once after** the ±14-point
  gather loop instead of being recomputed on every `j` iteration; only the
  final all-points value survives in the prototype, so the results are
  identical (efficiency hoist, noted per the review criteria).
- A degenerate/empty shortest path returns `false` (no veto) instead of
  crashing; unreachable-destination paths from `GraphAlgorithms` are short
  and handled naturally by the resampling.
- No other logic, constant, or precedence was changed; the three-tier kink
  condition is transcribed with the prototype's `&&`/`||` precedence
  (`G1 || (G2 && legs>10cm) || (G3 && legs>15cm)`, each
  `Gn = (open>θn && sym>5.5°) || open>60°`).

## Config

- C++: `TaggerCheckTGM` gains mixin `NeedPCTS` (`pc_transforms`, already
  emitted by the jsonnet builder) and knob `check_neutrino_candidate`
  (default `false`).
- `cfg/pgrapher/common/clus.jsonnet` `tagger_check_tgm(...,
  check_neutrino_candidate=false)` — key-suppression idiom, key omitted when
  off ⇒ byte-identical pre-port compiled config.
- `cfg/pgrapher/experiment/sbnd/clus.jsonnet` threads
  `tgm_neutrino_candidate=false` through `clus_pr` / `pr(...)`.
- Runner (`wcp-porting-img`): `sbnd_xin/wct-pr-perevt.jsonnet` TLA
  `tgm_neutrino_candidate=false`; `run_nusel_evt.sh` flag `-nucand` (or env
  `SBND_TGM_NUCAND=1`).  Default off.

## Verification

**Byte-identical, knob OFF** (gate labels; scratch =
`~/tmp/claude-25225/.../scratchpad`):

- Compiled-config proof: `wct-pr-perevt.jsonnet` compiled with knob off
  against the pre-change cfg tree — **byte-identical** (`cfg_old_off.json` ==
  `cfg_new_off.json`); with knob on the `check_neutrino_candidate` key
  appears (`cfg_new_on.json`).
- Output gate: all 10 MCP2025C events re-run knob-off with the new library,
  `mabc-pr.zip` member-content hashes (`abtest/hash_archive.py`) vs the
  pre-change production outputs `work-mcp10/nusel_evt*/mabc-pr.zip`:
  **10/10 PASS** (label `scratchpad/tgm_gate`, e.g. evt284349
  `5a043a01822b…`).
- `./build/clus/wcdoctest-clus`: 41/41 pass.
- Freshness proof done (`local/lib/libWireCellClus.so` newer than source).

**Knob ON effect** (label `scratchpad/tgm_on10`): across the 10-event sample
exactly one verdict changes — run 18255 evt 284349 cluster 11 (the 1.555 µs
in-beam bundle, the motivating anode→downstream TPC0 crosser):
`TGM=false → TGM=true` (veto returns "not a neutrino candidate";
202.1 cm > 0.45 × 426.3 cm).  Its table label flips STM → TGM.  No
out-of-beam verdict changes anywhere in the sample.

**Status: knob-off path byte-identical (gate PASS above).  Knob-on is a
behavior change — enable via `-nucand` and revalidate against hand scans
before adopting as a production default.**  The veto's protective direction
(a true neutrino with through-going-looking extremes being spared) is not
exercised by this 10-event sample; a nueCC-enriched scan is the natural
follow-up before default-ON.
