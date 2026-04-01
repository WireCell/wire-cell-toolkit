# Deep Dive: Track/Shower Separation and Particle Identification

**Source files:**
- `clus/src/NeutrinoTrackShowerSep.cxx` — `separate_track_shower()`
- `clus/src/NeutrinoShowerClustering.cxx` — shower grouping algorithms
- `clus/src/NeutrinoKinematics.cxx` — energy/momentum reconstruction
- `clus/src/NeutrinoTaggerNuMu.cxx` — νμ CC tagger
- `clus/src/NeutrinoTaggerNuE.cxx` — νe CC tagger
- `clus/src/NeutrinoTaggerPi0.cxx` — π0 tagger
- `clus/src/NeutrinoTaggerCosmic.cxx` — cosmic rejection tagger
- `clus/inc/WireCellClus/NeutrinoPatternBase.h` — `PatternAlgorithms` class

---

## Overview

After `find_proto_vertex()` produces a stable PR graph, the reconstruction
must determine **what each segment is**. The two primary particle types are:

- **Track**: a charged particle (muon, proton, pion) traveling a nearly
  straight path. Identified by a relatively constant, moderate dE/dx.
- **Shower**: an electromagnetic shower (electron, photon → e+e- pair, π0 → γγ).
  Identified by rapidly increasing dE/dx (shower multiplicity) and a diffuse,
  cone-like spatial structure.

---

## Step 1 — Segment-Level Classification (`separate_track_shower`)

```cpp
void separate_track_shower(Graph& graph, Facade::Cluster& cluster);
```

This function iterates over all `PR::Segment` objects in the graph and
classifies each one by examining:

### Criterion 1: dE/dx profile shape

Using the charge measurements stored in the segment's `fits` vector
(populated by `TrackFitting`):

- A **track** has a roughly flat dE/dx with a possible Bragg peak at one end.
- A **shower** has a rapidly increasing dE/dx from start to peak (the shower
  maximum), then falling off.

The function computes a KS-test (Kolmogorov–Smirnov) statistic between the
observed dQ/dx distribution and a template shower profile. High KS statistic
→ shower; low → track.

### Criterion 2: Angular spread (PCA)

The points associated with a segment are collected and a principal component
analysis (PCA) is run:

```cpp
std::pair<geo_point_t, geo_vector_t> calc_PCA_main_axis(
    std::vector<geo_point_t>& points);
```

The ratio of the second to first PCA eigenvalue measures the transverse
spread. Tracks have a low ratio (pencil-like); showers have a high ratio
(cone-like).

### Criterion 3: Branching structure

A segment that has many short sub-branches hanging off it (found via the
local Steiner graph topology) is classified as a shower spine, since EM
showers produce many secondary tracks.

### Result

Each segment gets its flags updated:

```cpp
// PRSegment.h
enum class SegmentFlags {
    kShowerTrajectory = 1<<1,  // segment is part of a shower (trajectory-based)
    kShowerTopology   = 1<<2,  // segment is part of a shower (topology-based)
    kAvoidMuonCheck   = 1<<3,  // exclude from muon hypothesis test
};
```

---

## Step 2 — Direction Assignment (`determine_direction`)

```cpp
void determine_direction(Graph& graph, Facade::Cluster& cluster,
                          const ParticleDataSet::pointer& particle_data,
                          const IRecombinationModel::pointer& recomb_model);
```

For each segment, a `dirsign` (+1/-1/0) is assigned:

- **+1**: the segment points *away* from the primary vertex (particle
  originates at the vertex and travels outward).
- **−1**: the segment points *toward* the primary vertex (particle enters
  from outside — likely a cosmic).
- **0**: undetermined.

The direction is inferred by comparing the dE/dx at the two endpoints: a
Bragg peak (higher dE/dx) at the far end is consistent with a stopping
track, confirming +1 direction from the vertex. This uses the
`IRecombinationModel` to convert between dQ/dx (measured) and dE/dx
(physics quantity).

---

## Step 3 — Shower Clustering (`shower_clustering_with_nv`)

Individual shower-flagged segments must be grouped into coherent EM shower
objects (`PR::Shower`). This is done by:

```cpp
void shower_clustering_with_nv(
    int acc_segment_id,
    IndexedShowerSet& pi0_showers,
    ...,
    Graph& graph, VertexPtr main_vertex, ...);
```

The function:

1. Starts from the primary vertex and collects all segments flagged as shower.
2. Groups connected shower segments into `PR::Shower` objects using a greedy
   graph walk.
3. Handles π0 candidates separately:
   - `id_pi0_with_vertex()`: finds pairs of photon showers that reconstruct
     the π0 mass (134.98 MeV) when their invariant mass is computed.
   - `id_pi0_without_vertex()`: finds displaced shower pairs for cases where
     the π0 decays away from the main vertex.

The π0 identification uses the `Pi0KineFeatures` struct:

```cpp
struct Pi0KineFeatures {
    int    flag;       // 0=none, 1=with_vertex, 2=without_vertex
    double mass;       // reconstructed invariant mass
    double vtx_dis;    // distance from π0 vertex to main vertex
    double energy_1;   // shower 1 energy (from cal_kine_charge)
    double theta_1;    // shower 1 polar angle
    double phi_1;      // shower 1 azimuthal angle
    double energy_2;   // shower 2 energy
    double angle;      // opening angle between showers
    // ...
};
```

---

## Step 4 — Kinematics

### Track kinematics

For track segments, the kinetic energy is computed by integrating the
measured dQ/dx along the track path and applying the Bethe–Bloch formula
through the `IRecombinationModel`:

```
T_kinetic = ∫ (dE/dx) · dx   along segment
```

The `TrackFitting` class handles the integration, including corrections for:
- Track angle relative to wire planes (projected path length vs. true path)
- Dead wire gaps
- Space-charge distortions

### Shower kinematics (`calculate_shower_kinematics`)

```cpp
void calculate_shower_kinematics(
    IndexedShowerSet& showers,
    IndexedVertexSet& vertices_in_long_muon,
    IndexedSegmentSet& segments_in_long_muon,
    Graph& graph, TrackFitting& track_fitter,
    IDetectorVolumes::pointer dv,
    const ParticleDataSet::pointer& particle_data,
    const IRecombinationModel::pointer& recomb_model);
```

For EM showers, energy is reconstructed calorimetrically:

```
E_shower = cal_kine_charge(shower, ...) / W_ion
```

where `cal_kine_charge()` sums the ionization charge over all shower
segments (including sub-showers), and `W_ion` is the ionization energy per
electron-ion pair (23.6 eV in LAr, via the recombination model).

To avoid re-scanning all hits multiple times in a loop over showers,
the charge maps are pre-collected once:

```cpp
// NeutrinoPatternBase.h
void collect_charge_maps(TrackFitting& track_fitter);
// Populates: m_charge_2d_u, m_charge_2d_v, m_charge_2d_w
//            m_map_apa_ch_plane_wires
```

These cached maps are then reused by all `cal_kine_charge()` calls within
a single `shower_clustering_with_nv()` invocation, reducing the per-shower
cost from O(N_hits) to O(N_shower_points).

---

## Step 5 — Particle Identification

After kinematics are computed, higher-level particle ID taggers run:

### NeutrinoTaggerNuMu

```
clus/src/NeutrinoTaggerNuMu.cxx
clus/inc/WireCellClus/NeutrinoTaggerNuMu.h
```

Identifies νμ CC events: one long muon track + optional hadronic activity.
Uses the `ParticleDataSet` to look up the muon hypothesis (flat dE/dx,
minimal shower activity, possible Michel electron at the end).

### NeutrinoTaggerNuE

```
clus/src/NeutrinoTaggerNuE.cxx
clus/inc/WireCellClus/NeutrinoTaggerNuE.h
```

Identifies νe CC events: one EM shower (the electron) + optional hadronic
activity. The shower must start within a few cm of the vertex (no
conversion gap, unlike a photon from π0 decay).

### NeutrinoTaggerPi0

```
clus/src/NeutrinoTaggerPi0.cxx
clus/inc/WireCellClus/NeutrinoTaggerPi0.h
```

Identifies π0 events (two photon showers with invariant mass ≈ 134.98 MeV).
Operates on the `Pi0KineFeatures` struct filled during shower clustering.

### NeutrinoTaggerCosmic

```
clus/src/NeutrinoTaggerCosmic.cxx
clus/inc/WireCellClus/NeutrinoTaggerCosmic.h
```

Rejects cosmic muon contamination by checking:
- Does the cluster enter/exit through the top or bottom of the detector?
- Is the track's `dirsign` −1 everywhere (particle entering from outside)?
- Is the cluster outside the fiducial volume?

### NeutrinoTaggerSSM / SinglePhoton

```
clus/src/NeutrinoTaggerSSM.cxx
clus/src/NeutrinoTaggerSinglePhoton.cxx
```

Specialized taggers for sterile neutrino searches and single-photon
topologies (ν → ν + γ).

---

## The `ParticleDataSet`

```
clus/src/ParticleDataSet.cxx
clus/inc/WireCellClus/ParticleDataSet.h
```

A data table (loaded from a file at startup) that provides:
- Particle mass and mean free path for various species.
- Expected dE/dx vs. kinetic energy for muons, protons, pions, kaons.
- Range–energy relations for stopping particles.

All particle hypothesis tests use this dataset, making the PID logic
detector-agnostic (the dataset file changes per detector, not the code).

---

## `IRecombinationModel`

The `recombination_model` component (configured as
`wc.tn(ub.uBooNE_box_recomb_model)` in the jsonnet) converts between the
measured ionization charge density (dQ/dx) and the energy loss (dE/dx):

```
dE/dx = dQ/dx / (recombination_factor × W_ion)
```

The Box model (used for MicroBooNE) accounts for the fact that some
ionization electrons recombine with argon ions before being collected,
with the recombination fraction depending on the local charge density and
the electric field.

The `dQdx_scale` and `dQdx_offset` parameters in the jsonnet configuration
apply a linear calibration to the raw detector charge before passing it
through the recombination model.
