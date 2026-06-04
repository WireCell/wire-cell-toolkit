# Applying a per-TPC (Y,Z) position offset in SBND Q/L matching

**Status: design recommendation only — this document changes no code.** It answers the
question "if we want to apply the measured TPC0/TPC1 transverse offset so that QLMatching
uses the corrected position, where is the best place to inject it?"

## Background

The cathode-crossing diagnostic (`dump_cathode_diag`, see
`sbnd_xin/docs/cathode-crossing-diagnostic.md`) measured a small, data-specific transverse
gap between the TPC0 and TPC1 halves of a track where it crosses the central cathode:

```
combined transverse fit  T_yz = (Δy, Δz) = (−0.22, +1.34) cm
```

The drift-x part of the gap (~−1.5 cm in data) is **deliberately left to the existing t0 /
`flash_x_offset` correction** and is not re-calibrated here — closing it with a rigid x-shift
would also close the *physical* cathode dead-gap. So the new correction is purely transverse:
`(0, Δy, Δz)`.

The user's request has two specific constraints:

1. The corrected position must be **used inside QLMatching** (the charge→light prediction).
   In per-APA clustering the offset does not matter and must not be applied there.
2. The offset can be interpreted as a **mis-placement of one anode plane relative to the
   other**, so it is expected to "go with" the anode-referenced fiducial volume.

## What the code does today

The only position correction QLMatching applies is the drift/x term

```
flash_x_offset = sign_offset · flash_time · drift_speed      (QLMatching.cxx:662)
```

There is **no single "corrected-position" function** ("pos_correctedT" does not exist as a
chokepoint). `flash_x_offset` is recomputed per-flash and the corrected point is built inline
at the handful of sites that need it — every one shifts **x only**, passing y,z straight
through:

| site | what it builds |
|---|---|
| `build_bundles` (`QLMatching.cxx:724`) | corrected point fed to the photon model + PE-inclusion gate |
| `compute_endpoint_flags` (`:1551`, `:1652`) | drift-end flags; `:1652` is the cathode-fiducial test |
| `compute_two_boundary_flag` | the two-boundary diagnostic flag |
| `dump_cathode_diag` | the offset diagnostic |

Two facts make the injection point easy to reason about:

- **QLMatching is pure-3D + photon model.** It never reads the per-APA 2D wire
  associations. Those "wind" arrays are computed during per-APA clustering
  (`Facade_Grouping.cxx:448 convert_3Dpoint_time_ch`, stored at
  `PointTreeBuilding.cxx:347` and summarised per blob at `Facade_Blob.cxx:98`) and are
  **frozen** before matching. A Y,Z offset applied *locally, as QLMatching reads the 3D
  points* therefore **cannot desync the 2D association** — which is exactly why offsets are
  irrelevant in per-APA clustering and why this is the correct stage to apply them.
- The per-TPC **active-volume box** `run.y_lo/y_hi/z_lo/z_hi` comes from
  `m_dv->inner_bounds(wpid)` in the cryostat frame (`compute_geometry`, `:617-620`) and is
  used as the PE-inclusion gate (`:731`) and the ±y/±z boundary-flag distances (`:1689`).

## Recommendation

### 1. Create the chokepoint the offset should live in

Introduce the single corrected-position helper that the code is currently missing:

```cpp
geo_point_t QLMatching::corrected_point(const geo_point_t& p,
                                        double flash_x_offset,
                                        const ApaRun& run) const
{
    return { p.x() + flash_x_offset, p.y() + run.dy, p.z() + run.dz };
}
```

Carry the per-TPC `(dy, dz)` on `ApaRun`, set once in `compute_geometry` from a new jsonnet
knob, and replace the ~4 inline corrected-position constructions with this helper. The x term
is unchanged: the new offset is `(0, Δy, Δz)`.

**Default zero ⇒ correction OFF ⇒ existing production configs stay bit-identical** (the
project's toggleable-behavior convention).

### 2. Frame model: the offset "rides the anode"

Shift the reconstructed points **and** the anode-referenced active-volume box by the same Δ;
keep the genuinely shared / physical objects fixed.

| object | treatment | effect |
|---|---|---|
| reconstructed 3D points (PE prediction, drift-end flags) | `+ (0, Δy, Δz)` | corrected charge enters the photon model — the intended effect |
| active-volume box `y_lo/y_hi/z_lo/z_hi` | shift by `+ (Δy, Δz)` in `compute_geometry` | PE-inclusion gate and ±y/z boundary flags stay offset-invariant (the box moves with the anode) |
| `cathode_fiducial` CPA structure + PMT / `SemiAnalyticalModel` | **fixed** (cryostat frame) | the corrected point is tested / used against the fixed structure — this is where the offset has its real, intended effect |
| per-APA 2D "wind" arrays, `select_scope_fv` clustering FV | **untouched** | no 2D desync; honors "offsets don't matter in per-APA clustering" |

The active-volume box shift is what realises the user's "offset goes with the anode plane":
because the box rides with the points, the in-TPC containment gate is invariant, while the
cross-TPC relationship to the fixed PMTs (the photon prediction) and to the shared cathode
CPA does change — which is the physics the correction is meant to capture. At the measured
magnitude (~0.67 cm per side vs a multi-metre box) the box-shift-vs-fixed-box difference is
negligible in practice; shifting it is the conceptually clean choice.

### 3. Do **not** write the offset back into the point cloud

Keep the correction local to QLMatching's 3D reads. Persisting the shifted positions into the
stored point cloud is the one move that would desync the frozen "wind" arrays and corrupt any
downstream retile / Bee re-derivation that recomputes 2D associations from 3D positions
(`improvecluster_1.cxx:472`, `retile_cluster.cxx:251`). The offset is a read-time correction,
not a stored geometry change.

## Caveats to keep in mind

1. **The ±Δ/2 symmetric split was derived for the display.** The diagnostic constrains only
   the *relative* TPC0–TPC1 offset; splitting it symmetrically assumes the true position is the
   midpoint (the absolute per-TPC placement relative to the PMTs is not separately measured).
   This is negligible for the photon prediction at this magnitude, but it is an assumption,
   not a measurement.
2. **Interpretation is modest.** An anode mis-placement and a PMT mis-placement are partly
   degenerate. Treat this as an empirical position correction applied to the charge before the
   photon model, not as a claim about which piece of hardware actually moved.
3. **Matching-only ≠ display.** Because the offset never touches the stored point cloud, the
   Bee output is unchanged by this code path. The display correction remains the separate
   Bee-zip shift already documented in `sbnd_xin/docs/cathode-crossing-diagnostic.md`.

## Config shape (when implemented)

A per-TPC knob keyed by anode ident, default all-zero:

```jsonnet
pos_offset_tpc0: [0, -0.11, +0.67],   // cm; East,  x < 0
pos_offset_tpc1: [0, +0.11, -0.67],   // cm; West,  x >= 0
```

(values from the diagnostic's symmetric split of `T_yz = (−0.22, +1.34)` cm.)

## Files that would change when this is implemented

- `match/inc/WireCellMatch/QLMatching.h` — `ApaRun::dy/dz`, the `corrected_point` declaration,
  the new config members.
- `match/src/QLMatching.cxx` — parse the per-TPC knob in `configure()`, set `run.dy/dz` and
  shift the active-volume box in `compute_geometry`, route the ~4 corrected-position sites
  through `corrected_point`.
- the SBND Q/L jsonnet (the standalone `sbnd_xin` config) — expose `pos_offset_tpc{0,1}`,
  default zero.

No source files are changed by *this* document.
