# Bee Visualization Output

**Source files:**
- `clus/inc/WireCellClus/MultiAlgBlobClustering.h` — `BeePointsConfig`, `BeePFConfig`, `ApaBeePoints`
- `clus/src/MultiAlgBlobClustering.cxx` — all `fill_bee_*` functions
- `WireCellUtil/Bee.h` — `Bee::Points`, `Bee::ParticleTree`, `Bee::Sink`

---

## What is the Bee Viewer?

Bee is a web-based 3D event display for LArTPC data. It reads ZIP files
containing per-event JSON files. `MultiAlgBlobClustering` produces these ZIP
files at configurable points in the pipeline. There are two distinct output
modes:

1. **Point-cloud dumps** (`bee_points`): raw 3D point sets, one per named
   grouping or algorithm, used for visual debugging at any pipeline stage.
2. **Particle-flow tree** (`bee_pf`): a hierarchical JSON tree representing
   the full reconstructed event (tracks + showers + kinematics), consumed
   by Bee's particle-flow panel.

---

## Point-Cloud Dumps

### Configuration (`BeePointsConfig`)

```cpp
struct BeePointsConfig {
    std::string name;        // Bee dataset name (e.g. "live", "steiner")
    std::string detector;    // geometry name (e.g. "uboone")
    std::string algorithm;   // algorithm label
    std::string pcname;      // which named PC to dump from each cluster
    std::string grouping;    // which Facade::Grouping to read ("live", "dead")
    std::string visitor;     // dump after this visitor runs (empty = at end)
    std::vector<std::string> coords; // coordinate names in the PC (e.g. ["x","y","z"])
    bool individual;         // true = one file per APA/face; false = global
    int  filter{1};          // 1=normal, 0=off, -1=inverse filter
    double dQdx_scale{1.0};
    double dQdx_offset{0.0};
    bool use_associate_points{false}; // use dpcloud("associate_points") + shower charge
    bool use_graph_vertices{false};   // dump PR::Vertex positions instead of PC
};
```

### `fill_bee_points(name, grouping)`

```
clus/src/MultiAlgBlobClustering.cxx line 380
```

Iterates all clusters in `grouping`. For each cluster calls
`fill_bee_points_from_cluster()` which:

1. Reads the named point cloud (`pcname`) from the cluster's local PC storage.
2. For each point, computes an encoded cluster ID:
   ```cpp
   encoded_id = cluster_id * 1000 + segment_graph_index
   ```
   This encoding lets the Bee viewer link each point back to a specific
   segment for selection/highlighting.
3. Applies `dQdx_scale` and `dQdx_offset` to the charge value before appending.
4. Appends `(x, y, z, charge, encoded_id)` to `Bee::Points`.

If `individual=true`, points are bucketed per APA+face combination so Bee can
display one readout plane at a time. If `individual=false`, all points go into
a single global `Bee::Points` object.

### `fill_bee_points_from_pr_graph(name, grouping)`

Reads the PR graph from the grouping instead of raw point clouds. For each
segment in the graph, uses the `"fit"` point cloud (fitted track positions).
Charge is set to the dQ value from the fit multiplied by the calibration
factors.

### `fill_bee_vertices_from_pr_graph(name, grouping)`

Dumps vertex positions from the PR graph:
- Primary vertex (`kNeutrinoVertex`): charge = 15000 (appears bright in Bee)
- Other vertices: charge = 0

### Trigger timing

Each `BeePointsConfig` has a `visitor` field. After each visitor runs,
`MultiAlgBlobClustering::operator()` checks whether any `BeePointsConfig`
names that visitor and, if so, calls `fill_bee_points` immediately. This
allows dumping intermediate states (e.g., after clustering, after Steiner,
after PR) into separate Bee files.

---

## Particle-Flow Tree

### Configuration (`BeePFConfig`)

```cpp
struct BeePFConfig {
    std::string name{"mc"};       // output file name in the ZIP
    std::string visitor;          // dump after this visitor (usually "" = at end)
    std::string grouping{"live"}; // grouping to read PR graph from
};
```

The default name `"mc"` matches what the Bee viewer expects for its
particle-flow panel.

### `fill_bee_pf_tree(cfg, grouping)`

```
clus/src/MultiAlgBlobClustering.cxx line 697
```

Builds a `Bee::ParticleTree` — a hierarchical JSON structure representing
the neutrino interaction — by combining the PR graph, shower information,
and kinematics.

#### Step 1 — Retrieve PR state

```cpp
auto pr_graph   = grouping.get_pr_graph();
auto tf         = grouping.get_track_fitting();
auto main_vertex = tf->get_main_vertex();
auto showers    = tf->get_showers();
auto pi0_showers = tf->get_pi0_showers();
```

#### Step 2 — BFS over the track skeleton

Starting from `main_vertex`, BFS traverses **track-only** (non-shower) segments.
Shower segments are pre-marked as `used_segs` so the BFS never crosses them.
The BFS builds four maps:

| Map | Key | Value |
|---|---|---|
| `seg_parent` | segment | parent segment (nullptr = root) |
| `seg_children` | segment | list of child segments |
| `seg_endpoints` | segment | `{near_vtx, far_vtx}` |
| `vtx_incoming_seg` | vertex | the segment that first arrived at this vertex |

#### Step 3 — Resolve shower attachment

Showers can attach at three types of parents:
1. **Root-level** (start vertex == main_vertex or in `root_reachable_vtxs`):
   the shower is a direct daughter of the neutrino interaction.
2. **Track-level** (start vertex in `vtx_incoming_seg`): the shower branches
   off a track segment.
3. **Nested** (start vertex belongs to another shower): shower-inside-shower
   (e.g., photon converting inside a different shower).

An iterative fixed-point loop (`while any_added`) propagates `vtx_incoming_seg`
into shower vertex sets until all showers are resolved.

Connection type 4 showers are excluded entirely from the output.

#### Step 4 — Build jsTree nodes

Each track segment becomes a node with:
- `id`: `cluster_id × 1000 + segment_graph_index`
- `text`: `"<particle_name> KE=<energy_MeV> MeV"`
- `data`: `{start: [x,y,z], end: [x,y,z]}` (from fitted endpoint positions)
- `children`: sub-track segments and attached showers

Each shower becomes a leaf node with:
- `id`: `shower_id` (mirrors start segment ID)
- `text`: `"<particle> KE=<kenergy_best_MeV> MeV"`
- For π0 pairs: a pseudo-gamma parent node is inserted with the π0 invariant
  mass, with the two photon showers as children.

#### Output format

The resulting `Bee::ParticleTree` is serialized as a JSON array:
```json
[
  {
    "id": 1042,
    "text": "mu- KE=312.4 MeV",
    "data": { "start": [10.3, -5.2, 44.1], "end": [82.1, -12.3, 55.6] },
    "children": [
      { "id": 1040, "text": "p KE=88.1 MeV", "data": {...}, "children": [] }
    ]
  },
  {
    "id": 2001,
    "text": "e- KE=245.0 MeV",
    "data": { "start": [10.5, -5.0, 44.3], "end": [55.2, 2.1, 71.8] },
    "children": []
  }
]
```

---

## Output ZIP File

All Bee output is written to the ZIP file configured by `bee_zip` in the
jsonnet:

```jsonnet
MultiAlgBlobClustering: {
    ...
    bee_zip: "bee-output.zip",
    bee_points: [ { name: "live", ... }, { name: "steiner", ... } ],
    bee_pf:     [ { name: "mc",   visitor: "" } ],
}
```

Inside the ZIP, each event produces files named:
```
<name>-<run>-<subrun>-<event>-*.json   (point cloud files)
mc-<run>-<subrun>-<event>.json         (particle flow tree)
```

The ZIP is flushed at the end of each event (`flush()`) and can be dragged
directly into the Bee web viewer.

---

## RSE (Run/Subrun/Event) Numbering

By default, the ident integer from the input tensor set is decoded as:
```cpp
run = (ident >> 16) & 0x7fff;
evt = ident & 0xffff;
```

If `use_config_rse=true` and `runNo`/`subRunNo`/`eventNo` are set in the
jsonnet configuration, those values override the ident decoding. This is
useful when the ident encodes something other than RSE.
