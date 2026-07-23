# FC / STM / Neutrino boundary check: two containment bugs fixed

Two independent defects in `Facade::cluster_fc_check`
(`clus/src/Clustering_Util.cxx`) — the two-round steiner-boundary containment
check shared by `TaggerCheckFC`, `TaggerCheckSTM` and `TaggerCheckNeutrino` —
made a genuinely-contained cathode-crossing cosmic read as an fiducial exiter.
Diagnosed on SBND MCP2025C reco1 run 18255 event 284657 main cluster 28 (an
out-of-time cosmic, matched flash t0 = −1051 µs, displays as a compact ~70 cm
track crossing the cathode, fully inside the FV box).

Both defects were independently fatal to the verdict (`is_fc =
exit_wcps.empty()`, so any one spurious exit forces `FC=false`), so both had
to be fixed for the cluster to become FC.

## Repro

```bash
cd /nfs/data/1/xqian/toolkit-dev/toolkit
wcbuild && ./build/clus/wcdoctest-clus
SB=/nfs/data/1/xqian/toolkit-dev/wcp-porting-img/sbnd/sbnd_xin
cd $SB && SBND_INPUT_DIR=$PWD/input_files_reco1/extracted-mcp2025c-10evt \
  SBND_WORK_ROOT=$PWD/work-mcp10 ./run_nusel_evt.sh data all
# evt 284657 main 28 (idx 2): FC 0 -> 1 after the fix.
# The env var WCT_FC_DEBUG (temporary instrumentation, not in the tree) was
# used during diagnosis to dump per-exit reasons and per-plane SP counts.
```

---

## Bug 1 — steiner boundary `x_t0cor` was uncorrected (raw drift-x)

### Symptom
`cluster_fc_check` read a steiner boundary point at x = −229 cm — 28 cm past
the TPC0 anode (−201 cm) — and flagged it as a direct fiducial-volume exit,
even though every displayed point of the cluster (and the cluster's own
`get_extreme_wcps` extremes) sits ≥ 74 cm inside the FV.

### Root cause
`switch_scope`'s `add_corrected_points` writes the T0-corrected `x_t0cor` only
into the **original** cluster's blob `3d` point clouds. The steiner cloud is
built afterward from a **retiled** copy of the cluster
(`CreateSteinerGraph` → `improve_cluster_2`, i.e. `improvecluster_2.cxx`), and
in that retiler the correction was applied in the wrong order:

```cpp
new_cluster.add_corrected_points(m_pcts, correction_name);  // reads get_cluster_t0() == 0
new_cluster.from(*orig_cluster);                            // copies the real T0 — too late
```

`add_corrected_points` builds `x_t0cor` from `get_cluster_t0()`, but the
freshly-retiled `new_cluster`'s T0 still defaults to **0** at that point —
`from()` (which copies `cluster_t0`) had not run yet. So the correction applied
a **zero drift shift** and `x_t0cor` came out equal to the raw drift-x.
`cluster_fc_check` then reads the steiner boundary as
`steiner_pc.get("x_t0cor")` (the default-scope x name) and gets that raw value.

### Why it hid
The error is exactly `v_drift × T0`, so it is **sub-centimeter for in-time
clusters** (the whole validated population) and only becomes large for
out-of-time cosmics. On cluster 28 it is `0.1563 cm/µs × 1051 µs = 164 cm`:
the steiner boundary points match real display tips in y,z to 0.0–0.2 cm, and
their x is off by 164.0 cm (measured), pushing them past the anode.

### Fix
`clus/src/improvecluster_2.cxx`: set the real T0 **before** the correction.

```cpp
new_cluster.set_cluster_t0(orig_cluster->get_cluster_t0());  // <-- added
new_cluster.add_corrected_points(m_pcts, correction_name);
new_cluster.from(*orig_cluster);
```

After the fix the steiner boundary of cluster 28 reads (−68.2, 103.1, 75.1) —
the real x_min tip, inside the box — and the spurious direct-FV exit is gone.

### Not fixed here (same latent pattern)
`clus/src/improvecluster_1.cxx` and `clus/src/retile_cluster.cxx`
(`RetileCluster`) contain the identical "correct-before-`from()`" ordering.
They are not exercised by the SBND steiner path (which uses
`improve_cluster_2`), so they are left untouched pending their own validation;
flagged here so the owner can decide.

---

## Bug 2 — SP gap check counted induction-only 2D hits

### Symptom
Even with bug 1 present, a second exit fired: `check_signal_processing`
(`clus/src/FiducialUtils.cxx`) declared the endpoint (−17, 54, 123) a
prolonged-signal exit, walking 127 steps outward and finding all "occupied".

### Root cause
A step was marked occupied by
`result_u>0 || result_v>0 || result_w>0 || dead` — **any single plane**.
Per-plane over the 127 steps: **U=6, V=127, W=0**. The collection plane W had
**zero** hits (no real charge beyond the track — the track genuinely ends), yet
the V induction plane registered 2D activity at every step (projection
coincidences in the dense raw region near the anode where the backward
transform lands). A lone induction-plane hit, which the collection plane
contradicts, is not evidence the track continues.

### Why it hid
Faithful port of the prototype
(`prototype_base/2dtoy/src/ToyFiducial.cxx:1622`, same OR). In-time,
well-separated tracks rarely have the outward direction graze foreign
induction activity; a cathode-crossing cosmic whose extreme point is a
mid-track side point (here z_max, since the track runs along x) does.

### Fix
`clus/src/FiducialUtils.cxx` `check_signal_processing`: require the **W
(collection)** plane to confirm charge (dead regions still count).

```cpp
auto result_w = m_internal.live->get_closest_points(temp_p_raw,1.2*units::cm,apa,face,2);
if (result_w.size() > 0 || inside_dead_region(temp_p_raw, apa, face)) {
    num_points_dead++;
}
```

W is the reliable indicator of real charge (collection essentially always
carries signal where charge exists), so U/V-only shadows no longer mark the
path occupied. Edge note: a real track segment sitting on a single dead W wire
(with U,V live) is now only caught by the ≥2-view `inside_dead_region` term,
not by W charge — acceptable, as that is a dead-region case by construction.

---

## Scope / reproducibility

Both changes live entirely in the **PR tail** (steiner → fiducialutils →
TGM/STM/FC/Neutrino), which runs only after `switch_scope`:

- `improvecluster_2.cxx`: the touched block runs only when
  `default_scope != raw_scope`, i.e. after a T0 correction — never in the main
  (pre-switch_scope) clustering pipeline, whose byte-identical img/clus/QL/
  matching gates are therefore unaffected.
- `check_signal_processing` is called only from `cluster_fc_check` (the FC/STM/
  Neutrino taggers).

These are **deliberate behavior changes** (no default-OFF knob, per owner
request — both are correctness fixes). Not byte-identical; the PR-tail tagger
verdicts change.

## Verification

- `./build/clus/wcdoctest-clus`: 41/41 pass. Freshness proof done.
- 10-event MCP2025C rerun, FC=1 count 3 → 6; **5 verdicts changed**, all on
  out-of-time cosmics or STM boundary cases:
  - FC 0→1: evt 284657 main 28 (t −1051 µs, the motivating cluster),
    evt 285999 main 20 (−1176 µs), evt 286241 main 12 (−742 µs).
  - label not-tagged→STM: evt 286021 main 15 (−577 µs),
    evt 286065 main 12 (+1132 µs).
- No in-time / small-|t0| verdict changed (bug 1's error is sub-cm there;
  bug 2 only bit the one grazing-direction endpoint).
