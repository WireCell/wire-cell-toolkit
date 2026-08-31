# Deep Dive: Track Fitting and Calorimetry

**Source files:**
- `clus/src/TrackFitting.cxx` — ~8,350 lines, core implementation
- `clus/src/TrackFitting_Util.cxx` — utility functions
- `clus/inc/WireCellClus/TrackFitting.h` — `TrackFitting` class
- `clus/inc/WireCellClus/TrackFittingPresets.h` — parameter presets

---

## Role of TrackFitting

`TrackFitting` is not a WireCell component (it does not implement any
`IComponent` interface). It is a **utility class** instantiated by
`TaggerCheckNeutrino` and passed as a dependency to every function in
`PatternAlgorithms`. It encapsulates:

1. **Wire → space coordinate projection**: converting between 3D points and
   their observed wire/time signatures on each readout plane.
2. **Charge measurement with error**: reading dQ/dx from the raw wire data
   for each point on a segment.
3. **Track fitting** (linear regression in wire-time space): finding the best
   straight-line approximation to a segment.
4. **Calorimetric integration**: summing dE/dx along a path to get kinetic energy.
5. **2D charge map caching**: pre-collecting all wire charge data once per
   event to accelerate per-shower energy calculations.

---

## TrackFitting::Parameters

The `Parameters` struct holds all numerical constants that describe the
detector response:

```cpp
// clus/inc/WireCellClus/TrackFitting.h
struct Parameters {
    // Diffusion
    double DL = 6.4 * pow(units::cm,2)/units::second;  // longitudinal
    double DT = 9.8 * pow(units::cm,2)/units::second;  // transverse

    // Software filter broadening (wire pitch units)
    double col_sigma_w_T  = 0.188060 * 3*units::mm * 0.2; // collection plane
    double ind_sigma_u_T  = 0.402993 * 3*units::mm * 0.3; // U induction
    double ind_sigma_v_T  = 0.402993 * 3*units::mm * 0.5; // V induction

    // Relative and absolute charge uncertainties
    double rel_uncer_ind  = 0.075;    // 7.5% relative for induction
    double rel_uncer_col  = 0.05;     // 5.0% relative for collection
    double add_uncer_col  = 300.0;    // +300 electrons absolute for collection

    // Longitudinal filter (time dimension)
    double add_sigma_L    = 1.428249 * 0.5505*units::mm / 0.5;

    // Charge error parameters
    double rel_charge_uncer = 0.1;   // 10% relative
    double add_charge_uncer = 600;   // 600 electrons absolute

    // Fit quality thresholds
    double scaling_quality_th = 0.5;
    double scaling_ratio      = 0.05;

    // Default dQ/dx for unfit points
    double default_dQ_dx = 5000;  // electrons/cm
};
```

These parameters are loaded at startup from the file specified by
`trackfitting_config_file` in the jsonnet configuration. Different
parameter values are used for different detectors (uBooNE vs. SBND vs. DUNE).

---

## Fitting Modes

```cpp
enum class FittingType {
    Single,    // fit one segment at a time
    Multiple   // fit multiple segments jointly (for vertices)
};
```

**Single fitting** is used for most track segments. The segment's waypoints
are projected onto each wire plane, and a linear regression fits the best
line in wire-time space. The fit residuals (distance of each point from the
fitted line) are used to detect kinks (split points) and to compute an
overall quality metric.

**Multiple fitting** is used at vertices where two or more segments share
an endpoint. The fit must simultaneously constrain all segments to share
the same vertex point while individually minimizing each segment's residual.
This is solved as a constrained least-squares problem using Eigen's
iterative solvers.

---

## The Fit Structure

Each fitted point on a segment is described by a `PR::Fit`:

```cpp
// PRCommon.h (included via PRVertex.h / PRSegment.h)
struct Fit {
    WireCell::Point point;   // 3D position from the fit
    double range{0};         // cumulative arc length from segment start
    int    index{-1};        // index into the segment's waypoint array
    bool   flag_fix{false};  // if true, this point is held fixed in the fit
    // ...
    void reset();
    double distance(const WireCell::Point& p) const;
};
```

A `PR::Vertex` holds one `Fit` (its fitted position). A `PR::Segment` holds
a `std::vector<Fit>` (one per waypoint along the trajectory).

---

## Wire Coordinate System

LArTPCs have three wire planes:
- **U** (induction, ±60°)
- **V** (induction, ∓60°)
- **W** or **Y** (collection, 0°)

plus a drift time axis **X**. A 3D point `(x, y, z)` maps to:

```
U_wire = f_U(y, z)          (linear function of wire angle)
V_wire = f_V(y, z)
W_wire = f_W(y, z)          (= y/wire_pitch for vertical wires)
t      = x / v_drift
```

`TrackFitting` performs these projections via the `IAnodePlane` geometry
interface. The fitting is done in `(wire, time)` space separately per plane,
then the results are combined to give a 3D fit.

---

## Charge Measurement

For each point on a segment, `TrackFitting` reads the charge in the
surrounding wire-time region:

1. Project the 3D point onto each wire plane.
2. Look up the charge measurement `Q` at the closest wire-time bin.
3. Apply calibration corrections (electronics response, lifetime, SCE).
4. Estimate the charge error `σ_Q` using the parameters above:
   ```
   σ_Q = sqrt( (rel_uncer * Q)² + add_uncer² )
   ```
5. Compute the projected path length `Δl` in cm (accounting for track angle).
6. `dQ/dx = Q / Δl`

---

## 2D Charge Map Caching

The most expensive operation in shower energy reconstruction is reading charge
data from the wire planes for every shower segment. With hundreds of shower
points and multiple shower candidates per event, this becomes a performance
bottleneck.

The solution: pre-collect all charge data into maps **once** at the start of
shower reconstruction, then reuse them:

```cpp
// Called once per event at start of shower_clustering_with_nv()
void collect_charge_maps(TrackFitting& track_fitter);

// Resulting maps (owned by PatternAlgorithms)
ChargeMap m_charge_2d_u;  // map<CoordReadout, ChargeMeasurement>
ChargeMap m_charge_2d_v;
ChargeMap m_charge_2d_w;
WireMap   m_map_apa_ch_plane_wires;
```

Then `cal_kine_charge()` is called with the pre-collected maps:

```cpp
double cal_kine_charge(ShowerPtr shower,
    const ChargeMap& charge_2d_u,
    const ChargeMap& charge_2d_v,
    const ChargeMap& charge_2d_w,
    const WireMap& map_apa_ch_plane_wires,
    TrackFitting& track_fitter,
    IDetectorVolumes::pointer dv);
```

---

## Calorimetric Energy Reconstruction

### For tracks: range-based and dE/dx integration

Two complementary methods are available:

1. **Range method**: for stopping tracks, look up the range→energy relation
   in the `ParticleDataSet` tables. More precise but only applicable to
   contained tracks.

2. **dE/dx integration**: for all tracks,
   ```
   T_kin = ∫ (dE/dx)(x) dx
         ≈ Σ_i (dE/dx)_i · Δl_i
   ```
   where `(dE/dx)_i = (dQ/dx)_i / (recombination_factor × W_ion)`.
   The recombination factor accounts for charge that recombines before
   being collected (see `IRecombinationModel`).

### For showers: total charge integration

```
E_shower = (total collected charge) / (W_ion × avg_recombination_factor)
         = cal_kine_charge(shower, ...) / W_ion
```

The `cal_kine_charge()` function sums `dQ/dx · Δl` over all hits associated
with the shower. The correction factor for recombination is averaged over
the shower volume.

---

## TrackFittingPresets

```cpp
// clus/inc/WireCellClus/TrackFittingPresets.h
namespace TrackFittingPresets {
    TrackFitting::Parameters create_with_current_values();
    TrackFitting::Parameters create_for_uboone();
    TrackFitting::Parameters create_for_sbnd();
}
```

These factory functions return `Parameters` objects initialized to the
standard values for specific detectors. `TaggerCheckNeutrino` uses
`create_with_current_values()` as the default, which is then overridden by
the values in `trackfitting_config_file`.

---

## Segment Fit Quality (`Fit` residuals)

After fitting, each segment has a `fits` vector. The quality of the fit is
assessed by computing the mean squared residual:

```
χ² = (1/N) Σ_i |p_i - fit_line(p_i)|² / σ²_i
```

Segments with `χ²` above a threshold are candidates for splitting at the
point of maximum residual (the kink). This is used in `break_segments()`.
