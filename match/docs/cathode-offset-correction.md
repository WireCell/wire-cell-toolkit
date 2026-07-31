# Applying a per-TPC (Y,Z) position offset in SBND Q/L matching

**Status: IMPLEMENTED (scope route), DATA-ONLY.** The recommendation below was adopted as written:
the per-TPC `(Δy,Δz)` lives in the DetectorVolumes metadata (single source of truth) and is
materialized **after** QLMatching by the existing `T0Correction`/`switch_scope` rail as the
`{x_t0cor, y_cor, z_cor}` corrected scope.

**Data-only gating (since the `reality` gate).** This offset is an empirical *data* calibration:
it was measured from data cathode-crossers (data transverse ~1.4 cm vs MC ~0). MC carries no such
misalignment, so applying it to MC would inject a spurious shift. `pos_offset_on` is therefore no
longer a hard-coded boolean — it is derived from the `reality` parameter of
`cfg/pgrapher/experiment/sbnd/clus.jsonnet`: `reality='data'` → ON, `reality='sim'` (MC) → OFF.
The clus maker's callers thread `reality` through: `wcls-img-clus.jsonnet` (LArSoft production,
`std.extVar('reality')`), `sbnd_xin/wct-clus-matching-perevt.jsonnet` and
`sbnd_xin/wct-clus-matching-standalone.jsonnet` (the standalone dev chains, from the `mc|data` run
arg). The default is `reality='data'`, so the **data path stays byte-identical** to the previous
always-on state (verified by full-graph diff); only the MC path changes (offset now absent).

What was built (see *Files that would change* at the bottom for the exact sites):

1. **Common place** — `cfg/pgrapher/experiment/sbnd/clus.jsonnet` adds a `pos_offset` array to the
   `a0f0pA`/`a1f0pA` DetectorVolumes metadata blocks (gated by `pos_offset_on`, which is now
   `reality == 'data'` — see the data-only note above). The SAME `detector_volumes()` helper feeds
   both the per-APA DV given to QLMatching and the all-anode DV given to `switch_scope`, so one
   entry reaches both consumers via `m_dv->metadata(wpid)`.
2. **QLMatching reads it in (read-only for now)** — `ApaRun::dy/dz` are populated in
   `compute_geometry()` from `pos_offset` and logged, but **not yet consumed** (no inline
   `corrected_point`, no photon-model/active-volume change). They are parked for the forthcoming
   cross-TPC matching judgement.
3. **The scope change** — `clus/src/PCTransforms.cxx` `T0Correction` now shifts `y,z` by the
   per-TPC `(Δy,Δz)` (constant, independent of `cluster_t0`) in `forward()` and exactly inverts it
   in `backward()`, at both Point and Dataset level. Presence of `pos_offset` flips
   `output_scope()`/`stored_array_names()` to carry `y_cor/z_cor`; absence keeps the original
   `{x_t0cor,y,z}` scope (production / non-SBND **bit-identical**, verified by config diff). A hard
   `backward(forward(p)) == p` roundtrip assertion guards the inverse (`clus/src/clustering_test.cxx`).
4. **3D→2D safety** — confirmed: every 3D→wire site (`convert_3Dpoint_time_ch` callers, and the
   `connect_graph_relaxed`/TrackFitting/FiducialUtils paths that run post-`switch_scope`) reads
   `point3d_raw()` or back-transforms via `backward()` first, so the wire mapping is taken in the
   raw frame. Correctness there reduces exactly to the `backward()` inverse of #3.

The original design discussion is retained below unchanged for provenance.

---

**Original status: design recommendation only.** It answers the
question "if we want to apply the measured TPC0/TPC1 transverse offset so that QLMatching
uses the corrected position, where is the best place to inject it?"

## Background

The cathode-crossing diagnostic (`dump_cathode_diag`, see
`sbnd_xin/docs/14_cathode-crossing-diagnostic.md`) measured a small, data-specific transverse
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

### Premise correction: there is no "compare two APAs together" step

It is tempting to think of the chain as "judge points per-APA, then bring TPC0 and TPC1
into a common frame and compare them — and *that* is where the transverse offset enters."
**That cross-TPC comparison step does not exist in QLMatching today.** `run_one_apa`
(`QLMatching.cxx:465-478`) runs the *entire* pipeline — `build_bundles`, the LASSO
`fit_round1`/`fit_round2`, `apply_matched_t0s` — independently inside each per-APA
`ApaRun`. The multi-APA branch (`:426-446`) is **pure tree concatenation** (`merge_pct`
joins only the root opflash PCs and `take_children`); no fit or containment test ever
spans both TPCs.

So inside QLMatching the offset is **not** a TPC0-vs-TPC1 comparison. It is a per-TPC
**rigid shift of one TPC's charge relative to the *fixed* reference** — the PMTs (photon
model) and the shared cathode CPA. That is exactly the role caveat #1 (the symmetric
±Δ/2 split) already describes: an absolute per-TPC placement against the fixed frame, not
a relative measurement re-applied between the two halves.

### The three stages where the offset is (or is not) felt

| stage | what runs | does the offset enter here? | is it "free"? |
|---|---|---|---|
| **A. inside QLMatching** | per-APA matching (`run_one_apa`) | yes — at the photon model + cathode CPA, against the fixed frame | **no** — QLMatching reads points inline, not through a scope (see below) |
| **B. the scope change** | post-QLMatching `switch_scope` materializes the T0 (and now also the transverse) correction as a *view* | yes — both corrections baked into one corrected scope | one transform edit |
| **C. downstream DV/FV checks** | tagging / boundary / containment that read the cluster | yes — inherited automatically from B's default scope | **yes — one scope, zero call-site edits** (FV stays baked; small Y-Z boundary effect accepted, see below) |

Stage C is the user's main concern ("make the detector-volume / FV shift after the offset
clean and unconfusing"). The chosen answer is the simplest one: **shift the points, keep the FV
exactly as it is** — see *The fiducial volume* below.

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

### 2. Frame model: apply the shift only against the fixed external reference

The offset enters the points **only where charge is compared to a genuinely fixed external
reference** — the PMTs (photon model) and the shared cathode CPA. Everything that is
**anode-referenced** (the active-volume box, the per-TPC ±y/z boundary distances) stays in the
unshifted frame, which keeps its containment result **invariant** — exactly the user's "the
offset goes with the anode plane."

| object | treatment | effect |
|---|---|---|
| reconstructed 3D points at the PE-prediction site | `+ (0, Δy, Δz)` (via `corrected_point`) | corrected charge enters the photon model — the intended effect |
| active-volume box `y_lo/y_hi/z_lo/z_hi` + the PE-inclusion gate / ±y/z boundary flags | **fixed**; the gate reads **raw** `y,z` | offset-**invariant** (anode-referenced ⇒ box and points share the misplacement) |
| `cathode_fiducial` CPA structure + PMT / `SemiAnalyticalModel` | **fixed** (cryostat frame) | the corrected point is tested / used against the fixed structure — where the offset has its real effect |
| per-APA 2D "wind" arrays, `select_scope_fv` clustering FV | **untouched** | no 2D desync; honors "offsets don't matter in per-APA clustering" |

Invariance needs **no box-bounds shift and no `compute_geometry` change**: feeding the *fixed*
anode-frame box its *raw* `y,z` is exactly equivalent to shifting both the box and the points by
`+Δ`. So the active-volume gate simply keeps reading raw `y,z` while only the PE-prediction read
goes through `corrected_point`. (The *downstream* Stage-C FV, by contrast, reads the shifted
default scope against a baked boundary — the small, accepted Y-Z boundary effect documented
below. Stage A is exact because the raw/shifted choice is free at an inline read; Stage C accepts
the effect because exactness there would cost a wrapper.)

### 3. Do **not** write the offset back into the point cloud

Keep the correction local to QLMatching's 3D reads. Persisting the shifted positions into the
stored point cloud is the one move that would desync the frozen "wind" arrays and corrupt any
downstream retile / Bee re-derivation that recomputes 2D associations from 3D positions
(`improvecluster_1.cxx:472`, `retile_cluster.cxx:251`). The offset is a read-time correction,
not a stored geometry change.

## The scope / PCTransform machinery (the "other way", and how the FV fits)

There is already a non-destructive way to "change the scope of points" with an
old↔new correspondence: the `IPCTransform` + scope mechanism the T0 correction uses.
This is the cleanest home for the **transverse** part, and it also gives a crisp answer
to the fiducial-volume question.

### What the machinery does

`T0Correction` (`clus/src/PCTransforms.cxx:34`) is an `IPCTransform`. Its
`forward(Dataset, …)` (`:83`) reads the blob's raw `{x,y,z}` and emits corrected arrays,
but it is **additive**: `Cluster::add_corrected_points` (`Facade_Cluster.cxx:194`) writes
only the *new* arrays named by `stored_array_names()` (`:151`, just `{"x_t0cor"}`) into the
same `"3d"` PC and **leaves `x,y,z` in place**. The corrected coordinates are exposed as a
*view* — the registered `output_scope()` `{"3d",{"x_t0cor","y","z"}}` (`:145`) — not as an
overwrite. `ClusteringSwitchScope` (`clus/src/clustering_switch_scope.cxx:50`) runs this per
cluster, registers the scope, and makes it the cluster's **default scope** (`:77`) so
everything downstream (the all-APA merge steps and the Bee dump) reads corrected coordinates.

Two consequences matter here:

- **The old↔new "mapping" is index identity plus coexistence.** Because the transform is
  additive, point `i`'s raw `(x,y,z)` and corrected `(x_t0cor, …)` live side by side at the
  same index. Nothing is renumbered; both frames are queryable at once.
- **This is the sanctioned way to materialize a correction** — and it does *not* contradict
  §3 above. §3 forbids *overwriting* the stored positions (which would desync the frozen 2D
  "wind" arrays). The scope mechanism never overwrites; it adds named arrays and a view. So
  materializing the transverse offset *through a transform* is safe where writing it back into
  `x/y/z` would not be.

### Why the transverse offset fits the transform but the drift offset does not

`switch_scope` runs **after** QLMatching (`cfg/.../clus.jsonnet:276`; it consumes the
`cluster_t0` QLMatching wrote). That ordering is forced for T0: the x correction is
`sign·t0·v_drift`, so it is unknown until matching has *picked* a flash, and QLMatching
itself must try **every** candidate flash (`QLMatching.cxx:660`) with a different
`flash_x_offset` per flash. A single materialized x-scope therefore cannot feed the matching
decision.

The transverse offset is the opposite: `(Δy, Δz)` is a **per-TPC constant**, independent of
t0 and of which flash is being tested. So it *can* be materialized once, as a view, **before**
QLMatching — exactly the place T0 cannot go. The clean shape is therefore:

- **Stage A (inside QLMatching) — not free.** QLMatching reads points *inline*
  (`points.at(i).x()` at `QLMatching.cxx:724`) and applies `flash_x_offset` by hand; it does
  **not** read through a cluster scope. So the offset *inside matching* cannot ride the
  post-QLMatching default-scope flip — it needs either the `corrected_point` helper of §1 added
  inline at `:725-726`, or a `(Δy,Δz)` transform materialized **pre-QLMatching** that QLMatching
  explicitly reads (a `{"3d",{"x","y_cor","z_cor"}}` view in place of raw `{x,y,z}` at `:720`,
  with `flash_x_offset` still added to `x` on top). Pick by taste; the inline helper is the
  smaller change.
- **Stage B (the scope change) — recommended: extend `T0Correction` *in place*.** The existing
  post-QLMatching `switch_scope` already runs `T0Correction`. Rather than add a *sibling*
  transform (which would force a second scope + a second switch for no gain), have
  `T0Correction` also carry `(Δy,Δz)`: its `forward(Dataset)` shifts `y,z` and adds `y_cor,z_cor`
  to `stored_array_names()`/`output_scope()`, with the transverse part **ignoring `cluster_t0`**
  (it is a per-TPC constant, not a t0-dependent drift term). Source the per-`(apa,face)` `(Δy,Δz)`
  from `m_dv->metadata(wpid)` exactly as `time_offset`/`drift_speed` already are
  (`PCTransforms.cxx:43-47`). The corrected scope becomes `{"3d", {"x_t0cor", "y_cor", "z_cor"}}`,
  one view holding both corrections — which is what makes Stage C free.

A single source of truth for the per-`(apa,face)` `(Δy,Δz)` should feed both stages (the
`m_dv->metadata` knob), defaulting to zero ⇒ OFF ⇒ production bit-identical.

### The fiducial volume (Stage C): "turn it on after, off before" *is* the scope boundary

This is the user's main question — *"is there an easy way to turn the offset ON after QLMatching
and OFF before for the FV, avoiding a lot of code change?"* The answer is **yes, and it is the
existing T0 scope rail — not a new mechanism.**

Every Facade fiducial / detector-volume read funnels through **one chokepoint**: `sv()` →
`m_default_scope` (`Facade_Cluster.h:157`). `get_extreme_wcps`, `cluster_fc_check`
(`Clustering_Util.cxx`), `TaggerCheckNeutrino`, `TaggerCheckSTM` — none take an explicit scope;
they all read `point3d(i)` → `sv3d()` → `scoped_view(m_default_scope)`. And `switch_scope`
already calls `set_default_scope(correction_scope)` at exactly the post-QLMatching boundary
(`clustering_switch_scope.cxx:77`). So the **default scope is the on/off flag**:

- **before** `switch_scope`: default scope is raw `{x,y,z}` → offset **OFF**;
- **after** `switch_scope`: default scope is the corrected scope → offset **ON**.

Once Stage B has put `y_cor,z_cor` into that corrected scope (`{x_t0cor,y_cor,z_cor}`), **every
downstream DV/FV check inherits the transverse shift for free** — no per-call-site edits, no two
sets of FV objects, no `TranslatedFiducial` wrapper, no delta argument on `IFiducial`. You keep
one set of baked, unchanged FV objects; flipping which coordinate arrays the *default scope*
names does all the work. This is the **chosen decision: shift the points, leave the FV as it
is.**

### What "shift points, keep FV baked" costs — and why it is the right call here

The decision is not free of physics, and the doc is honest about it. In Y-Z the FV/detector
walls are **anode-referenced** (defined by the wire planes), and the offset is read as an anode
mis-placement — so a fully consistent treatment would move the walls by the same per-TPC Δ that
moves the points (walls and points are both anode-referenced ⇒ the containment result should be
**invariant** under the shift). Shifting only the points and leaving the baked FV fixed
therefore leaves the FV boundary effectively offset by **~Δ (≤ ~0.67 cm per side)** relative to
the points — and the same uniform default-scope flip also hands the shifted `y,z` to
`m_dv->contained_by` (the transform `filter()` at `PCTransforms.cxx:131`), the real anode wire
geometry.

This is a **deliberately accepted simplification, not an oversight**: the offset is sub-cm while
the FV margins are multi-cm, cathode-crossing points sit far from the transverse TPC boundaries
(so TPC identification via `contained_by` is unaffected), and the alternative buys exactness no
analysis here needs. **Exact Y-Z invariance** — if it is ever wanted — is a single, well-scoped
addition: wrap the baked FV in a per-TPC `TranslatedFiducial` that subtracts the per-TPC
`(Δy,Δz)` from each query point (keyed by TPC via `contained_by`/x-sign, leaving x untouched so
the T0 comparison stays real), sourced from the same `(Δy,Δz)` knob. That handles even merged
cathode-crossing clusters point-by-point. **It is intentionally not adopted here** — the small
boundary effect is the accepted cost of the clean "shift points only" route.

#### Per-call-site selection still applies *inside* QLMatching (Stage A)

The Stage-C default-scope flip covers only the *downstream* checks. **Inside QLMatching (Stage
A)** the frame distinction of §2 is still expressed as an explicit per-site coordinate choice,
because QLMatching reads inline rather than through a cluster scope:

| call site | frame | reads | offset effect |
|---|---|---|---|
| `SemiAnalyticalModel` PE prediction (`QLMatching.cxx:734`) | cathode/cryostat (PMTs fixed) | corrected `y+Δy,z+Δz` | **real** — intended |
| cathode CPA `m_cathode_fv->contained()` (`:1652`) | cathode (CPA fixed) | corrected `y+Δy,z+Δz` | **real** — intended |
| active-volume gate `run.y_lo/y_hi/z_lo/z_hi` (`:731`) | anode (rides the points) | original `y,z` | **none** — invariant |

Feeding **original** `y,z` to the unshifted anode-frame box is *exactly equivalent* to feeding
corrected `y,z` to a box shifted by `+Δ` — so even Stage A needs **no box-bounds shift and no
`TranslatedFiducial` wrapper**; it is just which `y,z` each of the ~3 inline sites reads. A
`TranslatedFiducial(child, Δy, Δz)` wrapper would only be worth adding if you chose *not* to
materialize the offset at all and wanted a single FV object that self-applies the delta — an
unnecessary single-use abstraction under the scope approach, where the delta already lives in
`y_cor,z_cor` and every FV stays baked and unchanged.

## Caveats to keep in mind

1. **The ±Δ/2 symmetric split was derived for the display.** The diagnostic constrains only
   the *relative* TPC0–TPC1 offset; splitting it symmetrically assumes the true position is the
   midpoint (the absolute per-TPC placement relative to the PMTs is not separately measured).
   This is negligible for the photon prediction at this magnitude, but it is an assumption,
   not a measurement.
2. **Interpretation is modest.** An anode mis-placement and a PMT mis-placement are partly
   degenerate. Treat this as an empirical position correction applied to the charge before the
   photon model, not as a claim about which piece of hardware actually moved.
3. **RESOLVED in the implementation — the transform (scope) route was chosen, so the display
   reflects the transverse shift directly.** The SBND Bee `clustering` point set now reads the same
   corrected coords as the scope (`common_corr_coords` → `{x_t0cor,y_cor,z_cor}` when on), so the
   transverse shift is materialized once, in the display too. **Do not also apply the separate
   Bee-zip transverse shift** documented in `sbnd_xin/docs/14_cathode-crossing-diagnostic.md` — that
   would double-count. The original both-paths discussion is kept below for provenance.

   **Matching-only ≠ display — unless you take the transform route.** The §1 inline
   `corrected_point` is local to QLMatching's reads and never touches the stored point cloud, so
   the Bee output is unchanged by it; the display correction then remains the separate Bee-zip
   shift documented in `sbnd_xin/docs/14_cathode-crossing-diagnostic.md`. If instead the transverse
   offset is carried by the `switch_scope` transform (the scope-machinery section above), it
   becomes the cluster's default scope (`clustering_switch_scope.cxx:77`) and the Bee dump reads
   it — so the display would then reflect the transverse shift directly and the separate Bee-zip
   step would be redundant. Choose one path for the display; don't apply both.

## Config shape (when implemented)

A per-TPC knob keyed by anode ident, default all-zero:

```jsonnet
pos_offset_tpc0: [0, -0.11, +0.67],   // cm; East,  x < 0
pos_offset_tpc1: [0, +0.11, -0.67],   // cm; West,  x >= 0
```

(values from the diagnostic's symmetric split of `T_yz = (−0.22, +1.34)` cm.)

## Files that would change when this is implemented

**Inline route (§1):**

- `match/inc/WireCellMatch/QLMatching.h` — `ApaRun::dy/dz`, the `corrected_point` declaration,
  the new config members.
- `match/src/QLMatching.cxx` — parse the per-TPC knob in `configure()`, set `run.dy/dz` in
  `compute_geometry`, route the corrected-position sites through `corrected_point` (cathode CPA +
  photon model read `y+dy, z+dz`; the active-volume gate keeps reading raw `y,z`, per the FV
  table above).
- the SBND Q/L jsonnet (the standalone `sbnd_xin` config) — expose `pos_offset_tpc{0,1}`,
  default zero.

**Scope/transform route (Stage B + C — recommended for the downstream DV/FV shift):**

- `clus/src/PCTransforms.cxx` — **extend `T0Correction` in place** (preferred over a sibling
  transform): `forward(Dataset)` also shifts `y,z` and adds `y_cor,z_cor`; the transverse part
  ignores `cluster_t0`; `stored_array_names()`/`output_scope()` grow to `{x_t0cor,y_cor,z_cor}`.
  Read the per-`(apa,face)` `(Δy,Δz)` from `m_dv->metadata(wpid)` alongside the existing
  `time_offset`/`drift_speed` (`PCTransforms.cxx:43-47`).
- the per-TPC `(Δy,Δz)` added to the **DetectorVolumes metadata** in the jsonnet (the natural home,
  since `T0Correction` already pulls its per-TPC constants from there), default zero.
- **No downstream FV/DV edits.** Because all Facade containment reads go through `sv()` →
  `m_default_scope` (`Facade_Cluster.h:157`) and the existing `switch_scope` makes the corrected
  scope the default (`clustering_switch_scope.cxx:77`), Stage C picks up the shift automatically.
- (Stage A only) if the offset is also wanted *inside* matching via a transform rather than the
  inline `corrected_point`, add a pre-QLMatching read so QLMatching reads the corrected scope at
  `QLMatching.cxx:720` instead of raw `{x,y,z}`. Otherwise the §1 inline helper covers Stage A.

All routes share one per-`(apa,face)` `(Δy,Δz)` config and default to zero ⇒ OFF ⇒ production
bit-identical. No source files are changed by *this* document.
