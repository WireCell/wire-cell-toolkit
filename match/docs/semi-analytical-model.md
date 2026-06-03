# The Semi-Analytical Photon Model

This note explains the optical light-prediction model used by `QLMatching`
(charge–light matching). It answers the questions that come up most often:

- Is the input a grid or an arbitrary position?
- What position range is accepted?
- How does it know which TPC (drift volume) an OpDet belongs to?
- What is the underlying mechanism — is it machine learning?
- Can it predict light on X-ARAPUCAs or only PMTs?
- What do the extra "efficiency" parameters mean physically?

Source: `match/src/SemiAnalyticalModel.cxx`,
`match/inc/WireCellMatch/SemiAnalyticalModel.h`. It is a dependency-free port of
`larsim/PhotonPropagation/SemiAnalyticalModel.{h,cxx}` — the same model LArSoft
uses for SBND fast optical simulation.

---

## 1. What it computes

For a single **scintillation point** (a 3-D position where charge/light is
produced), the model returns, **per optical detector**, a *visibility*: the
fraction of emitted scintillation photons that reach and would be seen at that
detector. Two visibilities are produced:

- **Direct / VUV** — 128 nm argon scintillation photons travelling straight from
  the point to the detector (`detectedDirectVisibilities`).
- **Reflected / VIS** — photons that first hit the cathode plane (coated with a
  wavelength-shifting reflector foil) and are re-emitted as visible light toward
  the detector (`detectedReflectedVisibilities`).

A visibility is a dimensionless geometric acceptance in roughly `[0, 1]`. It is
**not** a photoelectron count — converting visibility to predicted PE is done by
`QLMatching` using the efficiency parameters (§6).

---

## 2. Input: arbitrary position, not a grid

The entry points take a single continuous point:

```cpp
void detectedDirectVisibilities (std::vector<double>& vis, const Point& scintPoint) const;
void detectedReflectedVisibilities(std::vector<double>& vis, const Point& scintPoint) const;
```

`scintPoint` is **any** continuous `(x, y, z)` in **centimetres** — there is no
voxel grid and no lookup-table-by-cell. This is the defining difference from a
traditional *photon library*, which precomputes visibilities on a fixed 3-D
voxel grid and interpolates between voxels. Here every position is evaluated
analytically on the fly. (`QLMatching` calls it once per charge point of every
blob — see `QLMatching.cxx:387–391`.)

### Accepted range

There is no hard domain check inside the model; it will return a number for any
finite point. In practice the meaningful range is the active LAr volume, defined
by the injected `Geometry` block (SBND values, in cm):

| field                  | SBND value | meaning                                  |
|------------------------|-----------:|------------------------------------------|
| `cathode_x`            | `0.0`      | cathode plane sits at x = 0              |
| `active_center_y`      | `0.0`      | y centre of the active volume            |
| `active_center_z`      | `250.5`    | z centre                                 |
| `active_size_y`        | `407.465`  | full height in y (≈ ±203.7 cm)           |
| `active_size_z`        | `501.0`    | full length in z (0 … 501 cm)            |
| `vuv_absorption_length`| `2000.0`   | VUV attenuation length in LAr            |

(The injected `vuv_absorption_length = 2000 cm` is the operative value; note it
is far larger than the ~85 cm physical figure cited as the header default in
`SemiAnalyticalModel.h:59–61` — i.e. the JSON nearly disables VUV self-absorption
for SBND.)

Behaviour at the edges, rather than a thrown error:

- **VUV cutoff:** if the point-to-detector distance exceeds `MaxPDDistance`
  (config, default 1000 cm) the direct visibility is set to 0
  (`SemiAnalyticalModel.cxx:237`).
- **Table lookups** off the end of their range behave in three distinct ways,
  not one: the Gaisser–Hillas **4-parameter** angular bin is **clamped** to the
  last valid bin (`angle_bin`); the GH **border corrections** are **linearly
  extrapolated** in θ (`interpolate(..., extrapolate=true)`,
  `SemiAnalyticalModel.cxx:280–282, 290–292`); and the VIS 2-D depth×radius
  correction **returns the edge value** (`interpolate(..., extrapolate=false)`,
  `:108–110`). So none of them throw, but only the VIS correction and the GH
  parameter bin are truly clamped — the GH border term keeps extrapolating.
- **Gaisser–Hillas guard:** any non-finite, negative, or `>10` GH correction is
  forced to 0 (`SemiAnalyticalModel.cxx:303, 348`).

So feeding a point far outside the active volume yields small or zero
visibilities rather than a crash, but the parameterisation is only *valid* inside
the volume it was fit for.

---

## 3. How it distinguishes the two TPCs (APA 0 vs APA 1)

SBND has two drift volumes sharing a central cathode at `x = 0`. The model
enforces **same-TPC visibility** purely by the sign of x:

```cpp
// SemiAnalyticalModel.cxx:233  (direct) and :355 (reflected)
if ((scintPoint.x() < 0.) != (od.center.x() < 0.)) continue;   // vis stays 0
```

A scintillation point only illuminates detectors on the **same side of the
cathode** (same x-sign). This is the port of SBND's `SBNDOpticalPath_tool`: light
does not cross the opaque cathode. There is no APA index passed in — TPC identity
is derived entirely from geometry (the OpDet's `center.x` and the point's `x`).

For the reflected component, the cathode "hotspot" used as the secondary source
is placed at `plane_depth = ±|cathode_x|` with the sign taken from the point's x
(`SemiAnalyticalModel.cxx:317, 352`), so the reflected path is also kept on the
correct side.

(Separately, `QLMatching` applies its own per-TPC OpDet mask using the same
`center.x < cathode_x` test, so each per-APA flash is matched only against the
detectors on its side — see `match/docs/QL_algorithm.md`.)

---

## 4. Mechanism: analytic functions + small correction tables — no ML

This is **not** machine learning and **not** a precomputed voxel library. It is a
closed-form *semi-analytical* parameterisation. For each (point, detector) pair:

1. **Geometric acceptance** from the solid angle the detector subtends at the
   point, times an exponential VUV absorption factor:
   `visibility_geo = exp(-distance / λ_abs) · Ω / 4π`.
   - **PMTs (dome, type 1):** `Omega_Dome_Model` — an analytic dome solid-angle
     formula with a 9-bin angular parameter table (`SemiAnalyticalModel.cxx:534`).
   - **ARAPUCAs (flat rectangle, type 0):** `Rectangle_SolidAngle` — exact
     rectangular solid angle via `acos`, decomposed into sub-rectangles for
     off-axis points (`:465`, `:473`).
2. **Gaisser–Hillas correction** `GH(distance)` (`:413`) whose four parameters
   come from small fitted tables indexed by the emission angle θ
   (`GH_PARS_flat/dome`), plus a linear radial "border" correction
   `GH_border_*` that depends on the in-plane radius r of the point.
3. **Reflected light** repeats the idea in two hops: point → cathode hotspot
   (treated as a rectangular source), hotspot → detector, with its own 2-D
   (depth × radius) correction tables `VIS_correction_flat/dome`.

The only "data" are these compact fitted coefficient tables (Gaisser–Hillas
parameters, angular/radial border corrections, VIS correction grids) carried in
the JSON. Everything else is evaluated from analytic geometry. The helper math
(`fast_acos`, elliptic integrals for the disk model, 1-D/2-D interpolation) is
ported verbatim from LArSoft's `PhotonPropagationUtils`.

---

## 5. PMTs vs X-ARAPUCAs

The model **can predict light for both** and selects the formula by the OpDet
`type` field (`SemiAnalyticalModel.cxx:252–263`, `:373–384`):

- `type == 1` → **dome PMT** → `Omega_Dome_Model`, dome GH/VIS tables.
- `type == 0` → **flat (X)ARAPUCA** → `Rectangle_SolidAngle`, flat GH/VIS tables.

Both branches are fully implemented and exercised: the SBND `OpDets` table
carries 120 PMTs (`type=1`) and 192 X-ARAPUCAs (`type=0`), and the model
evaluates a visibility for whichever detector you pass.

**Important distinction:** producing a *visibility* for ARAPUCAs is supported by
the model, but the current SBND `QLMatching` only *uses* the PMT channels in the
fit (its OpDet mask defaults to `active_opdet_types = [1]`). So ARAPUCA
visibilities are computed-capable but not currently fed into the charge–light
likelihood. Enabling them is a config change, not a model change.

The header (`SemiAnalyticalModel.h:14–20`) documents what is **not** ported and
would raise/be wrong if used: lateral PDs, anode reflections, Xe / vertical-border
/ field-cage corrections, and disk PMTs (`type=2`). Only `orientation=0`
(anode/cathode-facing) detectors are supported.

---

## 6. The efficiency parameters — physical meaning

The model returns a **geometric/optical visibility** only. `QLMatching` turns it
into a predicted photoelectron count per detector with three extra factors
(`QLMatching.cxx:399`):

```cpp
pred_PE[idet] += q · QtoL · dir_vis · VUVEfficiency[idet]
               + q · QtoL · ref_vis · VISEfficiency[idet];
```

summed over every charge point `q` in the cluster. The factors mean:

- **`QtoL` (charge→light):** converts reconstructed charge into the number of
  scintillation photons emitted at the point. It folds in the LAr light yield and
  the anticorrelation between ionization charge collected and scintillation light
  produced. (SBND config sets `QtoL: 1.0`; the C++ default is `0.5`. It is an
  overall scale, degenerate with the absolute efficiency normalisation.)

- **`VUVEfficiency[idet]` (direct/VUV detection efficiency):** the per-detector
  probability of registering a photoelectron *given* a 128 nm VUV photon arrives.
  It bundles the wavelength-shifter (TPB) conversion efficiency, the photocathode
  / SiPM quantum efficiency, the light-collector (ARAPUCA / light-guide)
  collection efficiency, and any per-channel calibration scale.

- **`VISEfficiency[idet]` (reflected/VIS detection efficiency):** the same kind of
  detection efficiency but for the **reflected, wavelength-shifted (visible)**
  component coming off the cathode foils. It additionally folds in the cathode
  reflector + shifter behaviour.

Two efficiencies are needed because the **two light components are physically
different** — different wavelength (128 nm VUV vs reflected visible) and
different detection paths — so a single detector has different sensitivity to
each.

The efficiency arrays are per-OpDet (312 entries) and **not uniform** — they take
a few discrete values (e.g. VUV ∈ {0, 0.01752, 0.0392}, VIS ∈ {0.00271,
0.01264, 0.026, 0.0357}). The variation encodes detector type and PMT surface
treatment: entries with **`VUVEfficiency = 0` but nonzero `VISEfficiency`** are
detectors effectively blind to direct 128 nm VUV that register only the
reflected visible component (e.g. ARAPUCAs and uncoated PMTs). These tables are the
SBND defaults baked into `QLMatching.h`; a config may override them via the
`VUVEfficiency` / `VISEfficiency` JSON arrays.

---

## 7. Quick reference

| question                         | answer                                                        |
|----------------------------------|---------------------------------------------------------------|
| Grid or arbitrary point?         | Arbitrary continuous `(x,y,z)` in cm; evaluated analytically.  |
| Accepted range?                  | Active volume from `Geometry`; clamps/zeros outside, no crash. |
| VUV cutoff                       | distance > `MaxPDDistance` (1000 cm) → 0.                      |
| Which TPC?                       | By x-sign vs cathode (`x<0` ≠ `od.x<0` → not visible).         |
| Mechanism                        | Closed-form solid angle + Gaisser–Hillas + small fit tables.   |
| Machine learning?                | No.                                                            |
| PMTs only?                       | No — PMTs (type 1) **and** ARAPUCAs (type 0) both supported.   |
| ARAPUCAs in the QL fit?          | Capable, but default mask uses PMTs only (`active_opdet_types=[1]`). |
| Output unit                      | Dimensionless visibility (not PE).                             |
| Visibility → PE                  | `q · QtoL · vis · efficiency`, summed over charge points.      |

See also: `match/docs/QL_algorithm.md` (matching algorithm),
`match/docs/qlmatching-code.md` (config reference).
