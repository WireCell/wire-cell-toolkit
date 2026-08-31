# STM Tagger Review: Prototype vs Toolkit

**Scope:** Logic fidelity, bug detection, efficiency, determinism, multi-APA/face correctness.
**Prototype source:** `prototype_base/pid/src/ToyFiducial.cxx` (functions at the line ranges cited below).
**Toolkit source:** `clus/src/TaggerCheckSTM.cxx` + `clus/src/FiducialUtils.cxx` + `clus/src/Clustering_Util.cxx` (FCCheckResult).
**Date reviewed:** 2026-04-10.

---

## Functional-fidelity table

| Prototype function | Toolkit location | Logic | Bugs | Efficiency | Determinism | Multi-APA |
|---|---|---|---|---|---|---|
| `check_stm` (entry) | `check_stm_conditions`; `visit()` | ✅ clean visit() | ✅ STM flag set | ✅ unified run_pass | ✅ | ⚠️ drift from first wpid (TGM anode: ✅ fixed) |
| `inside_fiducial_volume` | `FiducialUtils::inside_fiducial_volume` | ⚠️ tolerance unused | ✅ | ✅ | ✅ | ✅ (world coords) |
| `check_dead_volume` | `FiducialUtils::check_dead_volume` | ✅ | ✅ | ✅ | ✅ | ✅ (per-point contained_by) |
| `check_signal_processing` | `FiducialUtils::check_signal_processing` | ✅ | ✅ via cluster_fc_check | ✅ | ✅ | ✅ (per-point contained_by) |
| `search_other_tracks` | `TaggerCheckSTM::search_other_tracks` | ✅ | ✅ | see §7 | ✅ | ✅ (per-point wpid_array) |
| `check_other_tracks` | `TaggerCheckSTM::check_other_tracks` | ✅ | ✅ | ✅ | ✅ | ✅ sign-agnostic |
| `find_first_kink` | `TaggerCheckSTM::find_first_kink` | ✅ shorted-Y restored | ✅ | ✅ | ✅ | ✅ sign-agnostic |
| `detect_proton` | `TaggerCheckSTM::detect_proton` | ✅ | ✅ (fixed -1 return) | ✅ | ✅ | ✅ (per-point paf) |
| `eval_stm` | `TaggerCheckSTM::eval_stm_core` | ✅ | ✅ (fixed -1 return) | ✅ L[] hoisted (A11) | ✅ | N/A (uses world-L) |
| `inside_dead_region` | `FiducialUtils::inside_dead_region` | ⚠️ different algorithm (see §10) | ✅ | ✅ | ✅ | ✅ (apa/face explicit) |
| `AddDeadRegion` | *(not ported)* | ✅ (handled by PC-tree) | — | — | ✅ improvement | ✅ |
| `check_other_clusters` | `TaggerCheckSTM::check_other_clusters` | ✅ | ✅ | ✅ | ✅ | ✅ |
| `check_full_detector_dead` | *(not ported)* | ⚠️ missing event_type bit 6 | — | — | — | — |

Legend: ✅ good, ⚠️ concern/partial, ❌ bug.

---

## 1. `check_stm` / `check_stm_conditions` (entry point)

### Prototype
`prototype_base/pid/src/ToyFiducial.cxx:405–1022`

Entry point. Two rounds of endpoint detection (steiner with `flag_cosmic=true`, then `=false`), each using `inside_fiducial_volume` + `VHoughTrans` + `check_signal_processing` + `check_dead_volume` per extreme point. Then `first_wcp/last_wcp` selection, forward pass (rough_path → tracking × 2 → `find_first_kink` → TGM check → `eval_stm` ladder → `search_other_tracks` → `check_other_tracks` → `detect_proton`), and if `flag_double_end`, backward pass with same structure.

### Toolkit
- Full pipeline: `TaggerCheckSTM::check_stm_conditions` at `clus/src/TaggerCheckSTM.cxx:2135–2697`.
- Entry in the visitor pattern: `TaggerCheckSTM::visit()` at `clus/src/TaggerCheckSTM.cxx:79–492`.

### Issues

**BUG-1 (CRITICAL): `visit()` is a development stub — pipeline not active.**
`visit()` (lines 79–492) overrides the raw `wcpts()` with 48 hard-coded UBoone test-event points (lines 168–218) and calls `break_segment` at a hard-coded position (lines 361–390). The real `check_stm_conditions()` call is commented out at line 440. The component currently does not perform any STM/TGM classification; it only exercises the tracking infrastructure.

**Fix:** Remove the 48-point `wcpts()` override and the hard-coded `break_segment`. Replace the stub `visit()` body with a call to `check_stm_conditions(cluster, associated_clusters)` and handle its return value (see BUG-2).

**BUG-2 (CRITICAL): `Flags::STM` never set.**
`check_stm_conditions` returns `true` when the cluster passes the STM test, but never calls `cluster.set_flag(Flags::STM)`. The TGM branch sets `Flags::TGM` correctly (lines 2334, 2348, 2539). STM is silently dropped.

**Fix:** When `visit()` (or `check_stm_conditions` itself) detects STM, add `cluster.set_flag(Flags::STM)`.

**LOGIC-1: Endpoint selection delegated to `cluster_fc_check`.**
The prototype's endpoint logic (lines 427–644) is replicated in `Facade::cluster_fc_check` (`Clustering_Util.cxx:74–287`). That function correctly calls `check_signal_processing` and `check_dead_volume` per endpoint (lines 162–179) and is multi-APA aware via `dv->contained_by(p1)`. ✅

**LOGIC-2: Backward pass is missing the one-outside + dead-volume TGM branch.**
Prototype forward pass (lines 777–789): `else if (!inside_fiducial_volume(pts.front()) && left_L < 3cm)` → check `check_dead_volume` → possibly tag TGM. The backward pass (lines 914–931) has no equivalent `else if` branch. Toolkit forward (2339–2352) has this branch; toolkit backward (2512–2541) also omits it. This asymmetry **matches the prototype exactly** — preserved intentionally.

**LOGIC-3: Forward `left_L < 5cm && dQdx < 1.8` condition.**
Prototype forward (line 803) has `left_L < 5*units::cm && (dQ/(...))/50e3 < 1.8` in the short-track reset. Prototype backward (line 941-943) is missing this condition. Toolkit forward (line 2369) has it; toolkit backward (2553-2555) also omits it. Matches prototype — preserved. ✅

**MULTI-APA-1: `drift_dir` in `check_stm_conditions` derived from `*wpids.begin()`.**
Lines 2167–2174 pick the first `WirePlaneId` from `cluster.grouping()->wpids()` to compute `drift_dir`. For single-face clusters this is correct. For clusters spanning multiple faces (e.g. DUNE opposite-drift faces) the drift direction may be wrong for half the track's points.

**Fix:** Thread per-point `drift_dir` through the TGM anode check (lines 2325–2330 and 2530–2535), using `m_dv->face_dirx(wpid)` for the wpid closest to each endpoint.

**DETERMINISM-1: `std::map<Cluster*, std::vector<Cluster*>> main_to_associated` in `visit()` (line 105).**
Pointer-keyed map. With a single key in practice, the ordering is benign, but replace with `std::map<int, std::vector<Cluster*>>` keyed by `Cluster::ident()` when the stub is cleaned up.

**EFFICIENCY-1: `check_stm_conditions` forward/backward duplication.**
Lines 2222–2441 (forward) and 2444–2625 (backward) are structurally identical — same L/Q accumulation, same TGM test, same `eval_stm` ladder, same `search_other_tracks` / `check_other_tracks` / `detect_proton` sequence. Factor into a private lambda or helper taking `(first, last)` and `flag_is_backward`.

**EFFICIENCY-2: `eval_stm` called up to 4× per pass, each rebuilding `L[]`.**
The ladder at lines 2385–2427 (forward) calls `eval_stm` up to 4 times for the same `adjusted_segment`. Each call reiterates all fit points to build `L[i]` and `dQ_dx[i]` (lines 1514–1526). Hoist the `L[i]` / `dQ_dx[i]` construction outside the ladder and pass them in (or cache them on the segment).

---

## 2. `inside_fiducial_volume`

### Prototype
`prototype_base/pid/src/ToyFiducial.cxx:1756–1818`

Two modes: flat polygon check (`boundary_xy/xz`) when `tolerance_vec==NULL`, or space-charge-boundary (SCB) arrays when a tolerance vector is given. The SCB mode mutates and reverts per-slice arrays.

### Toolkit
`clus/src/FiducialUtils.cxx:72–77`

```cpp
bool FiducialUtils::inside_fiducial_volume(const Point& p,
                                           const std::vector<double>& tolerance_vec) const
{
    // currently tolerance vector is not used ...
    return m_sd.fiducial->contained(p);
}
```

### Issues

**LOGIC-1: Tolerance vector not used.** The prototype's SCB correction mode (for space-charge boundary awareness) is not implemented — the `tolerance_vec` parameter is accepted but silently ignored. The prototype uses this in several callers (e.g. the NuE tagger). For STM, the prototype always passes `NULL` (no tolerance), so this has no impact on STM-specific behavior. The comment acknowledges the gap.

No further issues. `m_sd.fiducial->contained(p)` is the correct primary fiducial check for world-frame points. ✅

---

## 3. `check_dead_volume`

### Prototype
`prototype_base/pid/src/ToyFiducial.cxx:1716–1752`

Walks from `p` along `dir` in `step` increments while inside FV. Per-step: count dead via `inside_dead_region(temp_p)`. Early-exit when 4 live points found. Final: dead_frac > 0.81 → false, else true.

### Toolkit
`clus/src/FiducialUtils.cxx:80–124`

Same structure. Per-step: uses `m_sd.dv->contained_by(temp_p)` to get `(apa, face)`, then `transform->backward(temp_p,...)` to get raw coordinates, then calls `inside_dead_region(temp_p_raw, apa, face)`. Default `cut_value=4`, `cut_ratio=0.81`. ✅

Multi-APA correct: per-step `contained_by` call gives the right `(apa, face)` as the path crosses TPC boundaries. ✅

Minor note: `m_sd.dv->contained_by` returns `wpid.apa() < 0 || wpid.face() < 0` for points outside any TPC volume; those are counted as dead (`num_points_dead++` at line 100). In the prototype, the walk exits when `inside_fiducial_volume` fails (line 1720–1725) — the two paths handle out-of-volume points differently. In practice both resolve consistently because the while condition (`inside_fiducial_volume`) should prevent out-of-TPC points from being reached. ✅

---

## 4. `check_signal_processing`

### Prototype
`prototype_base/pid/src/ToyFiducial.cxx:1668–1714`

Walks from `p` along `dir`. Per step: queries `ct_point_cloud.get_closest_points(temp_p, 1.2cm, plane)` for U/V/W. Point is "dead" (from signal-processing perspective) if ANY plane has a nearby hit OR `inside_dead_region(temp_p)`. Early-exit when 5 live points. Final: dead_frac > 0.8 → false, else true.

### Toolkit
`clus/src/FiducialUtils.cxx:128–175`

Same structure, with `get_closest_points(temp_p_raw, 1.2cm, apa, face, plane)` and per-point APA/face from `contained_by`. Default `cut_value=5`, `cut_ratio=0.8`. Multi-APA correct. ✅

**BUG-1: `check_signal_processing` is never called from `TaggerCheckSTM.cxx`.**
The prototype invokes it inside `check_stm` (lines 486, 592) during endpoint selection for each extreme point whose direction is near-collinear with a wire plane. In the toolkit, this call is made correctly inside `cluster_fc_check` (`Clustering_Util.cxx:163`), which is called by `check_stm_conditions`. ✅ (The concern flagged during planning was based on searching TaggerCheckSTM.cxx only; the function IS wired via `cluster_fc_check`.)

No action needed here.

---

## 5. `search_other_tracks` (PR3DCluster method in prototype)

### Prototype
`prototype_base/pid/src/PR3DCluster_pattern_recognition.h:2–301`

Builds secondary track segments from remaining steiner-graph points not well-described by the primary track. Uses Voronoi + Kruskal MST + connected components, then re-fits sub-segments.

### Toolkit
`TaggerCheckSTM::search_other_tracks` at `clus/src/TaggerCheckSTM.cxx:1715–2011`

**Multi-APA:** Per-point `wpid_array[i].apa()/.face()` (lines 1746, 1948). `get_closest_2d_point_info(p, plane, face, apa)` (1767, 1954). `transform->backward(p, cluster_t0, face, apa)` (1783, 1964). `get_closest_dead_chs(p_raw, 1, apa, face, plane)` (1786–1788, 1967–1973). All per-point. ✅

**Containers (determinism):** Internal adjacency lists keyed by `size_t` indices (lines 1804, 1825, 1876). No pointer-keyed maps. ✅

---

## 6. `check_other_tracks`

### Prototype
`prototype_base/pid/src/ToyFiducial.cxx:319–371`

Loops fit tracks from index 1 (secondary tracks). Uses `track_length`, `track_medium_dQ_dx`, `get_track_length(2)` for straightness. Several (length, dQdx) thresholds and a `ntracks ≥ 3` count reject.

### Toolkit
`TaggerCheckSTM::check_other_tracks` at `clus/src/TaggerCheckSTM.cxx:2059–2126`

Logic matches prototype. ✅

**MULTI-APA-1 (BUG): `drift_dir_abs(1,0,0)` hardcoded at line 2064.**
The angle `dir1.angle(drift_dir_abs)` (line 2099) is used to skip tracks nearly perpendicular to drift (threshold ±7.5°). For faces with `face_dirx = -1`, the angle is unaffected (it uses `fabs(angle-90) < 7.5`), so the magnitude check is direction-sign-agnostic. ✅ Actually this one is fine — since the check is `fabs(angle - 90) < 7.5` and `angle` is computed via `acos` which is always [0,π], using `(1,0,0)` vs `(-1,0,0)` gives the same result because the angle between a vector and its negation is supplementary. **No bug here** — the absolute-value check makes the sign irrelevant. ✅

---

## 7. `find_first_kink`

### Prototype
`prototype_base/pid/src/ToyFiducial.cxx:1024–1266`

Two sweeps over the fit path:
- **Sweep 1** (lines 1104–1175): require `refl>20 && ave>10`. Inside FV. Several rejection cuts including `inside_dead_region`. "Shorted Y" guard (W-wire 7130–7269): checks V `ch_mcell_set_map` dead channels; if dead, skip or apply stricter cuts. Charge balance cut (sum_fQ/sum_bQ > 0.6, or extended condition). Returns `max_numbers[i]`.
- **Sweep 2** (lines 1177–1261): looser (`ave>15`). "Shorted Y" V check (1190–1199): V `ch_mcell_set_map` → skip. Then U/V/W `ch_mcell_set_map` dead-channel check (1202–1231) → suppress if balanced. Same charge cut. Returns `max_numbers[i]`.
- Returns `fine_tracking_path.size()` when no kink found (sentinel).

### Toolkit
`TaggerCheckSTM::find_first_kink` at `clus/src/TaggerCheckSTM.cxx:887–1214`

**Sweep 1 (1020–1112):** Same angle conditions (refl>20, ave>10). FV check. Rejection cuts including `inside_dead_region(current_point_raw, paf.at(i).first, paf.at(i).second, 2)`. "Shorted Y" block commented out (1069–1078). Same charge balance cut. Returns `max_numbers[i]`.

**Sweep 2 (1114–1210):** Same angle conditions (refl>20, ave>15). FV check. Same rejection cuts. "Shorted Y" V-plane guard commented out (1144–1154). U/V/W `is_wire_dead(paf.at(i).first, paf.at(i).second, plane, round(p*+k), round(pt.at(i)))` (1157–1181). Charge balance cut with dead-channel suppression (1200). Returns `max_numbers[i]`.

Returns `static_cast<int>(fine_tracking_path.size())` when no kink (sentinel). ✅

### Issues

**MULTI-APA-1 (BUG): `drift_dir_abs(1,0,0)` hardcoded at line 926.**
Used to compute `para_angles[i]` = max(|angle(v10, drift)-90°|, |angle(v20, drift)-90°|). `para_angles` gates the 5-point smoothing window (must be > 12° to contribute). If drift is actually `(-1,0,0)` for a face, all angles are supplementary, and `|angle-90°|` would be `|180°-90°|=90°` vs the intended value. Since `acos(v·drift)` for `drift=(-1,0,0)` gives `180°-acos(v·(1,0,0))`, we have `|angle-90|` = same as before (symmetric around 90). The magnitude of `|angle-90|` is identical regardless of sign. **No bug** — `para_angles` computation is sign-agnostic. ✅

The `drift_dir_abs` is also used in the angle2 computation for j=0 (line 964–971). Same argument: since `|angle-90|` is computed, sign does not matter. ✅

**LOGIC-1: "Shorted Y" UBoone region guard missing in both sweeps (commented out at 1069–1078 and 1144–1154).**
The prototype sweep-1 block (1132–1144) applies to W wire range [7130, 7269] and additionally suppresses kinks by checking V dead-channel `ch_mcell_set_map`. The toolkit has this commented out and marked as a TODO. For UBoone this is a physics gap (may allow false kinks in the shorted-Y region); for other detectors it is irrelevant.

**Action item A7:** For UBoone, implement the shorted-Y guard using `cluster.grouping()->is_wire_dead(paf.at(i).first, paf.at(i).second, 1, round(pv.at(i)+k), round(pt.at(i)))` gated by the appropriate W-wire range. The W-wire range [7135, 7264] is UBoone-specific and should be made a configuration parameter or looked up from the geometry.

**LOGIC-2: `inside_dead_region` called with `minimal_views=2` (toolkit line 1063, 1140).**
Prototype calls `inside_dead_region(fine_tracking_path.at(i))` which uses mcell UV/UW/VW intersection (effectively any 2-plane overlap). Toolkit calls `inside_dead_region(..., 2)` meaning 2+ planes must be dead. Equivalent semantics. ✅

**LOGIC-3: Return value on error is `-1` instead of `pts.size()`.**
When `fits` is empty or FiducialUtils unavailable (lines 920–922, 895–897), the function returns `-1`. In `eval_stm` and `detect_proton`, the check `kink_num < 0 || kink_num >= num_pts` treats -1 as "no kink" — correct. In `check_stm_conditions`, `kink_num = -1` falls into neither `if (kink_num == 0)` nor `else if (kink_num > 0)`, so `first_loop_end = dx_size` and `second_loop_start = dx_size`, giving `exit_L = total, left_L = 0`. This then triggers the `left_L < 8cm && dQdx < 1.5` short-track reset, which sets `kink_num = dQ.size()`. The behaviour is accidentally correct but should be explicit.

**Fix:** Change the error-return paths to return `static_cast<int>(fits.size())` (the "no kink" sentinel) instead of `-1`, and remove the `kink_num < 0` branches in callers.

---

## 8. `detect_proton`

### Prototype
`prototype_base/pid/src/ToyFiducial.cxx:1268–1472`

Takes `(main_cluster, kink_num)`. Builds `pts/dQ/dx` from fit path. `end_p = pts[kink_num]` or `pts.back()` when `kink_num == pts.size()`. Loops secondary fit tracks: delta-ray / Michel-electron checks by proximity to `end_p` and direction. Cumulative L[] / dQ_dx[]. Peak search. KS test vs `g_muon` (ROOT TGraphs) and constant 50e3. Decision mixing ks1/ks2, ratio1/ratio2, peak dQdx.

### Toolkit
`TaggerCheckSTM::detect_proton` at `clus/src/TaggerCheckSTM.cxx:1216–1476`

Same algorithm. Uses `particle_data()->get_dEdx_function("muon")->scalar_function(...)` instead of ROOT TGraph (cleaner). KS test via `WireCell::kslike_compare`. ✅

**BUG-1: Error return is `-1` cast to `bool` → `true`.**
When FiducialUtils is unavailable (line 1219–1222):
```cpp
if (!fiducial_utils) {
    SPDLOG_LOGGER_DEBUG(s_log, "detect_proton: TaggerCheckSTM: ...");
    return -1;
}
```
`-1` cast to `bool` is `true`, meaning "proton detected" — the STM check is then rejected (since `if (!detect_proton(...)) return true` in `check_stm_conditions`). The error path silently prevents STM tagging.

**Fix:** Return `false` (no proton detected, allow STM) or add a third-state output. For robustness, return `false` (let the track proceed as STM candidate when proton determination is impossible).

**MULTI-APA:** Per-point `paf` used at line 1360 (`transform->backward(...)` with `paf.at(i).second, paf.at(i).first`). Secondary tracks' DPC closest-point queries use segment-level face/apa from the segment's fit data. ✅

---

## 9. `eval_stm`

### Prototype
`prototype_base/pid/src/ToyFiducial.cxx:1474–1653`

Builds `L[i]` (Euclidean cumulative from start), `dQ_dx[i]`. Finds `end_L` at `L[kink_num] - 0.5cm` (or `L.back()` for no-kink). Peak search in `[end_L - peak_range, end_L + 0.5cm]`. Then `com_range` window for comparison data, residual analysis. ROOT TH1F KolmogorovTest ("M"). Decision on `ks1-ks2`, `ratio1/ratio2`.

### Toolkit
`TaggerCheckSTM::eval_stm` at `clus/src/TaggerCheckSTM.cxx:1478–1685`

**L[] computation:** Euclidean distance between consecutive points (lines 1520–1526). Same as prototype (1484–1487). ✅

**kink sentinel:** Prototype `kink_num == pts.size()` → toolkit `kink_num < 0 || kink_num >= num_pts` (line 1530). Same semantics when `kink_num = pts.size()`. ✅

**KS test:** `WireCell::kslike_compare` replaces ROOT `KolmogorovTest("M")`. Both compute the KS statistic without correction factors. Should be equivalent; confirm with unit test if needed.

**Decision logic:** All thresholds match prototype exactly (lines 1648–1681 vs prototype 1618–1651). ✅

**BUG-1: Error return `-1` cast to `bool` → `true`.**
Line 1484: `return -1` when FiducialUtils unavailable. Interpreted as `true` → caller believes STM passes. Same fix as detect_proton: return `false`.

**EFFICIENCY-1: `L[]` and `dQ_dx[]` rebuilt on every call.**
`check_stm_conditions` calls `eval_stm` up to 4× on the same `adjusted_segment` with the same `kink_num` but varying `(peak_range, offset_length, com_range)`. The `L[]` and `dQ_dx[]` arrays (lines 1514–1526) are constant across calls. Extract them into a small struct computed once and passed through.

---

## 10. `inside_dead_region`

### Prototype
`prototype_base/pid/src/ToyFiducial.cxx:1820–1879`

Converts `p` → `time_slice` and U/V/W channel numbers using hard-coded MicroBooNE wire geometry (channel offsets 2400/4800, angles ±60°/0°). Looks up per-channel mcell sets in `ch_mcell_set_map`. Intersects UV, UW, VW pairs. Returns true if any intersection's mcell time range covers `time_slice`.

Effectively: returns true when at least 2 planes overlap in a dead mcell that spans the given drift time — a "2-of-3 planes dead in the same geometric cell" test.

### Toolkit
`clus/src/FiducialUtils.cxx:35–69`

```cpp
bool FiducialUtils::inside_dead_region(const Point& p_raw, const int apa, const int face, const int minimal_views) const
{
    const auto [tind_u, wind_u] = m_internal.live->convert_3Dpoint_time_ch(p_raw, apa, face, 0);
    // ... same for v, w
    if (m_internal.live->is_wire_dead(apa, face, 0, wind_u, tind_u)) dead_view_count++;
    // ... same for v, w
    return dead_view_count >= minimal_views;
}
```

**LOGIC-1: Different dead-region algorithm.**
The prototype uses mcell-level UV/UW/VW geometric intersection. The toolkit uses `is_wire_dead(apa, face, plane, wind, tind)` which checks the `WireDataCache::dead_wires[plane][wind_index]` time interval. If the toolkit's dead-wire map encodes the same extended ranges as the prototype's `AddDeadRegion` (with `dead_region_ch_ext` extension applied to both channel and time windows), the results should be equivalent.

The key question is whether `build_wire_cache()` (which populates `WireDataCache::dead_wires` from the PC-tree `dead_winds_*` datasets) applies a time-range extension equivalent to the prototype's `dead_region_ch_ext`. If the PC-tree already incorporates the extended ranges (i.e. was built with the extension applied), then the toolkit is equivalent. If not, toolkit `inside_dead_region` will under-report near cell boundaries.

**Action item A10:** Verify in `Facade_Grouping.cxx:build_wire_cache` whether the dead-wire time ranges stored in `dead_wires[plane][wind_index]` already include an extension. If not, clarify whether the `dead_winds_*` PC-tree datasets already encode the extension that the prototype's `AddDeadRegion` applies.

**LOGIC-2: `minimal_views` parameter.**
The toolkit is called with `minimal_views=2` in `find_first_kink` (lines 1063, 1140). The `check_dead_volume` call (inside `cluster_fc_check`) calls `inside_dead_region(temp_p_raw, apa, face)` using the default `minimal_views=1`. The prototype `check_dead_volume` calls `inside_dead_region(temp_p)` which uses the 2-plane intersection logic. Using `minimal_views=1` (any plane dead) is more liberal than the prototype's 2-plane requirement — it will flag more points as dead, making `check_dead_volume` more likely to return `true` (track can continue through dead region). This is a conservatism difference that slightly biases toward fewer exit candidates.

**Fix:** Change the default `minimal_views` from 1 to 2 in `FiducialUtils::check_dead_volume` calls to `inside_dead_region`, to match the prototype semantics.

---

## 11. `AddDeadRegion` / Dead-region storage

### Prototype
`prototype_base/pid/src/ToyFiducial.cxx:1882–1944`

`AddDeadRegion(mcell, time_slices, bad_u, bad_v, bad_w)`: builds `mcell_time_map[mcell] = {front-ext, back+ext}` (time range with `dead_region_ch_ext` extension) and populates `ch_mcell_set_map[ch]` for each bad plane's wires, with channel range extended by `±dead_region_ch_ext` and clamped to per-plane windows (U: 0–2399, V: 2400–4799, W: 4800–8255).

### Toolkit
No `AddDeadRegion` function exists. Dead regions are populated automatically when `Facade::Grouping::build_wire_cache(apa, face, plane)` reads the `dead_winds_<apa>_<face>_<plane>` PC-tree datasets. This is architecturally cleaner and inherently multi-APA (each `(apa, face, plane)` is stored separately). ✅

The prototype's pointer-keyed `mcell_time_map` (`std::map<SlimMergeGeomCell*, std::pair<int,int>>`) and `ch_mcell_set_map` (`std::map<int, std::set<SlimMergeGeomCell*>>`) are replaced by integer-keyed `dead_wires[plane][wire_index] → (t_min, t_max)` — a determinism improvement. ✅

See LOGIC-1 in §10 for the open question about whether the PC-tree data includes the `dead_region_ch_ext` extension.

---

## 12. Cross-cutting: Drift direction handling

The prototype hard-codes `TVector3 drift_dir(1,0,0)` in `check_stm` (line 414) and `find_first_kink` (line 1035). This was correct for single-TPC MicroBooNE.

Toolkit usage summary:

| Location | Line | Handling | Multi-face correct? |
|---|---|---|---|
| `check_stm_conditions` TGM angle check | 2325, 2530 | `drift_dir` from `m_dv->face_dirx(*wpids.begin())` | ⚠️ single-face |
| `find_first_kink` para_angles | 926 | `WireCell::Vector drift_dir_abs(1.0, 0.0, 0.0)` | ✅ sign-agnostic |
| `check_other_tracks` perpendicular skip | 2064 | `geo_point_t drift_dir_abs(1, 0, 0)` | ✅ sign-agnostic |
| `cluster_fc_check` wire-angle checks | Clustering_Util.cxx:134–158 | per-point `dv->contained_by(p1)` + `wpid_params` | ✅ |

The only true multi-APA concern is `check_stm_conditions`'s use of `*wpids.begin()` for the TGM anode check (lines 2317 and 2522): the check `pts.at(kink_idx).x() < 6*cm` and `pts.back().x() < 2*cm` implicitly assume anode-at-positive-x. For faces with `face_dirx = -1`, the anode is at the detector's boundary in the negative-x direction. This needs a per-face threshold.

**Action item A12:** For the TGM anode check, derive `anode_x_threshold` from `m_dv` geometry for the relevant `(apa, face)` rather than hardcoding the positive-x thresholds (2 cm and 6 cm). This makes the TGM check multi-APA correct.

---

## 13. Cross-cutting: `check_full_detector_dead`

### Prototype
`prototype_base/pid/src/ToyFiducial.cxx:276–317`

Scans `mcell_time_map` for dead channels with drift position in `[-10cm, 266cm]`. Returns `true` if > 8000 dead channels. Called at the driver level (`wire-cell-prod-stm-port.cxx:1368`) to set `event_type |= 1UL << 6`.

### Toolkit
No equivalent function exists. This is an event-level dead-channel fraction flag used to set a bit in the event type, but the toolkit's event classification uses `Flags` on the cluster, not an event-level int. The information is available from `cluster.grouping()->get_all_dead_chs(apa, face, plane)`.

**Action item A13:** If the event-level dead-detector flag (`event_type bit 6`) is needed downstream, implement it as a separate visitor or as a check inside `TaggerCheckSTM::visit()` that reads `cluster.grouping()->get_all_dead_chs(...)` for all `(apa, face, plane)` combinations.

---

## 14. Cross-cutting: STM pipeline wiring in `visit()`

The correct `visit()` body (once the dev stub is removed) should be:

```cpp
virtual void visit(Ensemble& ensemble) const {
    m_track_fitter.set_detector_volume(m_dv);
    m_track_fitter.set_pc_transforms(m_pcts);

    auto groupings = ensemble.with_name(m_grouping_name);
    if (groupings.empty()) return;
    auto& grouping = *groupings.at(0);

    Cluster* main_cluster = nullptr;
    std::vector<Cluster*> associated_clusters;
    for (auto* cluster : grouping.children()) {
        if (cluster->get_flag(Flags::main_cluster)) main_cluster = cluster;
        else if (cluster->get_flag(Flags::associated_cluster)) associated_clusters.push_back(cluster);
    }
    if (!main_cluster) return;

    bool is_stm = check_stm_conditions(*main_cluster, associated_clusters);
    if (is_stm && !main_cluster->get_flag(Flags::TGM)) {
        main_cluster->set_flag(Flags::STM);
    }
}
```

Note: `check_stm_conditions` already sets `Flags::TGM` internally, so only `Flags::STM` needs to be set here.

---

## Action items (ordered by severity)

| ID | Severity | File:line | Description | Status |
|---|---|---|---|---|
| A1 | CRITICAL | `TaggerCheckSTM.cxx:79–492` | Remove 48-point wcpts override and break_segment stub; wire `check_stm_conditions` into `visit()`. | **DONE** (dev stub guarded by `if (false &&...)`, clean `visit()` written) |
| A2 | CRITICAL | `TaggerCheckSTM.cxx:440` (commented) | Set `Flags::STM` when `check_stm_conditions` returns true and `Flags::TGM` not set. | **DONE** (in new `visit()` body) |
| A3 | HIGH | `TaggerCheckSTM.cxx:find_first_kink` | `find_first_kink` returns `-1` on error; change to return `static_cast<int>(fits.size())` ("no kink" sentinel). | **DONE** |
| A4 | HIGH | `TaggerCheckSTM.cxx:detect_proton` | `detect_proton` returns `-1` (→ true) on error; change to return `false`. | **DONE** |
| A5 | HIGH | `TaggerCheckSTM.cxx:eval_stm` | `eval_stm` returns `-1` (→ true) on error; change to return `false`. | **DONE** |
| A6 | MEDIUM | `FiducialUtils.cxx:106` | `check_dead_volume` calls `inside_dead_region` with default `minimal_views=1`; change to `minimal_views=2` to match prototype 2-plane semantics. | **DONE** |
| A7 | MEDIUM | `TaggerCheckSTM.cxx:find_first_kink` | Implement the shorted-wire-region guard in `find_first_kink` in a detector-agnostic way. | **DONE** — guard is controlled by new config parameter `shorted_y_w_range` (JSON int array `[w_min, w_max]`). Default is `[]` (disabled, detector-agnostic). Set to `[7135, 7264]` in UBoone jsonnet config to restore original behaviour. Member vars `m_shorted_y_w_min`/`m_shorted_y_w_max` default to -1. |
| A8 | MEDIUM | `TaggerCheckSTM.cxx:2167–2174` | For TGM anode check, derive `anode_x_threshold` per `(apa, face)` from `m_dv` geometry instead of hardcoded 2 cm / 6 cm. | **DONE** — added `dist_to_anode` lambda using `m_dv->contained_by` + `face_dirx` + `inner_bounds`. Falls back to `|x|` for out-of-volume points (preserves UBoone behaviour). Applied in both forward and backward pass. |
| A9 | MEDIUM | `TaggerCheckSTM.cxx:2222–2441, 2444–2625` | Factor the duplicated forward/backward pass body into a private helper. | **DONE** — unified into `run_pass` lambda taking `(start_wcp, end_wcp, is_forward)`. `is_forward` controls: early fits.size()≤3 exit, TGM dead-volume else-if branch, left_L>40cm guard, extra 5cm short-track condition. Forward pass calls `run_pass(first,last,true)`; backward `run_pass(last,first,false)`. |
| A10 | MEDIUM | `Facade_Grouping.cxx:build_wire_cache` | Verify that dead-wire time ranges in `WireDataCache::dead_wires` include the equivalent of prototype's `dead_region_ch_ext` extension. Document finding. | **FINDING (no code change needed):** `build_wire_cache` reads `dead_winds_a{apa}f{face}p{plane}` PC-tree datasets and stores raw `{start_x, end_x}` per wire index without applying any extension (`dead_region_ch_ext`). Whether the extension is already embedded in the PC-tree datasets depends on the data-generation pipeline. If the `dead_winds_*` arrays are written by the same upstream step that once called `AddDeadRegion`, the extension is already present in the data. If they come from raw detector bad-channel lists, the extension is absent. **Recommendation:** verify at the PC-tree production step. No code change in `build_wire_cache` itself. |
| A11 | LOW | `TaggerCheckSTM.cxx:eval_stm` | Hoist `L[]`/`dQ_dx[]` construction out of the 4× `eval_stm` call ladder. | **DONE** — split into `build_eval_arrays(segment)` → `STMEvalArrays{pts,L,dQ_dx,valid}` and `eval_stm_core(arrs, kink_num, …)`. `run_pass` calls `build_eval_arrays` once before the eval ladder, then calls `eval_stm_core` for each combination. Public `eval_stm(segment,…)` wrapper preserved for external callers. |
| A12 | LOW | `TaggerCheckSTM.cxx:visit()` (old line 105) | `std::map<Cluster*,…>` pointer-keyed map removed with the stub — no longer present. | **DONE** (map removed with stub) |
| A13 | LOW | *(new code)* | If needed downstream, implement `check_full_detector_dead` equivalent using `get_all_dead_chs(apa, face, plane)`. | open — not yet needed by any downstream consumer. |

---

## Verification

1. **Build:** `cd build && ninja wirecell-clus`.
2. **Pipeline test:** `wire-cell -c clus/test/uboone-mabc.jsonnet` on a known-good UBoone input. Only pipeline instantiating `TaggerCheckSTM`.
3. **Flag check:** Confirm `Flags::STM` and `Flags::TGM` appear on expected clusters in the bee output / SPDLOG debug output.
4. **Prototype parity:** Compare STM decisions against a prototype run on the same input. Document intentional divergences (multi-APA improvements, dead-region algorithm equivalence).
5. **Grep audits after fixes:**
   - No `return -1` inside `find_first_kink`, `detect_proton`, or `eval_stm`.
   - No hard-coded 48-point `wcpts()` array in `visit()`.
   - `Flags::STM` appears in `TaggerCheckSTM.cxx`.
