# Deep Dive: Data Structures — Facade Layer and PR Graph

**Source files:**
- `clus/inc/WireCellClus/Facade_Ensemble.h`
- `clus/inc/WireCellClus/Facade_Grouping.h`
- `clus/inc/WireCellClus/Facade_Cluster.h`
- `clus/inc/WireCellClus/Facade_Blob.h`
- `clus/inc/WireCellClus/Facade_Mixins.h`
- `clus/inc/WireCellClus/Facade_Util.h`
- `clus/inc/WireCellClus/PRVertex.h`
- `clus/inc/WireCellClus/PRSegment.h`
- `clus/inc/WireCellClus/PRGraph.h`
- `clus/inc/WireCellClus/PRGraphType.h`
- `clus/inc/WireCellClus/PRCommon.h`

---

## The Facade Pattern

The underlying data is stored in a **PC tree** (WireCell's point-cloud tree,
`WireCellUtil/PointTree.h`): a hierarchical structure where each node holds
a `PointCloud::Dataset` (a collection of named numeric arrays).

The Facade classes wrap these PC tree nodes and give them physics semantics
without copying the data. They provide convenient accessors, lazy-evaluated
cached properties, and named graph storage.

```
PC Tree node (raw)         Facade wrapper
────────────────────       ───────────────────
root node              →   Facade::Ensemble
  child node           →   Facade::Grouping  (e.g. "live", "dead")
    child node         →   Facade::Cluster
      child node       →   Facade::Blob
```

---

## Facade::Ensemble

```cpp
// clus/inc/WireCellClus/Facade_Ensemble.h
class Ensemble : public NaryTree::FacadeParent<Grouping, points_t> { ... };
```

The top-level container. Holds a set of `Grouping` objects. Accessed via:

```cpp
// Iterate groupings
for (auto& grouping : ensemble.children()) { ... }

// Find a named grouping
auto* live = ensemble.grouping("live");
```

---

## Facade::Grouping

```cpp
// clus/inc/WireCellClus/Facade_Grouping.h
class Grouping : public NaryTree::FacadeParent<Cluster, points_t> { ... };
```

A named collection of clusters. Standard groupings:

| Name | Contents |
|---|---|
| `"live"` | Clusters from active (non-dead) wires — main physics clusters |
| `"dead"` | Clusters reconstructed near dead wire regions |

The `FiducialUtils` object is attached here after `MakeFiducialUtils` runs.

```cpp
// Iterate clusters
for (auto& cluster : grouping.children()) { ... }

// Access the FiducialUtils attached to this grouping
auto& futil = grouping.fiducial_utils();
bool inside = futil.is_inside_fiducial(point);
```

---

## Facade::Cluster

```cpp
// clus/inc/WireCellClus/Facade_Cluster.h
class Cluster : public NaryTree::FacadeParent<Blob, points_t>
              , public Mixins::Cached<Cluster, ClusterCache>
              , public Mixins::Graphs
{ ... };
```

The central data object. A cluster wraps a PC tree node whose children are
blobs, and whose local dataset holds named arrays.

### Key accessors

```cpp
// Cluster identity
int  get_cluster_id() const;
void set_cluster_id(int id);

// T0 / timing
double get_cluster_t0() const;
void   set_cluster_t0(double t0);

// Coordinate scope management
void set_default_scope(const Tree::Scope& scope);
const Tree::Scope& get_default_scope() const;

// Named point clouds (local arrays on this cluster's node)
PointCloud::Dataset& local_pc(const std::string& name);
```

### Named point clouds on a cluster

| Name | Content | Set by |
|---|---|---|
| `"3d"` | Raw 3D blob positions | Upstream imaging |
| `"tagger_info"` | Numeric flag arrays from upstream taggers | Tagger stages |
| `"isolated"` | Points in "isolated" sub-cluster | Upstream clustering |
| `"perblob"` | One point per blob | Upstream clustering |
| `"associate_points"` | Points from shower/delta-ray clustering | `clustering_points()` |
| `"fit"` | Fitted 3D positions from TrackFitting | `TrackFitting` |

### Flags on a cluster

Flags are stored as a string-keyed map (base class from `NaryTree::Facade`):

```cpp
cluster.set_flag("beam_flash", true);
bool cosmic = cluster.get_flag("cosmic_tagger");
bool recovered = cluster.get_flag("associated_cluster");
```

### Cached properties (`ClusterCache`)

Expensive-to-compute properties (e.g., total length, bounding box, PCA
eigenvectors) are computed once and cached in a `ClusterCache` object.
The cache is invalidated whenever the cluster is mutated:

```cpp
void invalidate_cache() { clear_cache(); }
```

---

## Facade::Blob

```cpp
// clus/inc/WireCellClus/Facade_Blob.h
class Blob : public NaryTree::FacadeLeaf<points_t> { ... };
```

The leaf-level data object. A blob corresponds to one wire-cell volume: the
intersection of wire signals on the three planes during a time slice. It holds:

- The wire channel ranges for each plane that contributed.
- The time slice index.
- The charge measurement.
- A 3D position (centroid or sampling point).

---

## Mixins::Graphs

```cpp
// clus/inc/WireCellClus/Facade_Mixins.h
namespace Mixins {
    class Graphs {
    public:
        // Store and retrieve named PR::Graph objects
        void set_graph(const std::string& name, PR::Graph graph);
        PR::Graph& get_graph(const std::string& name);
        bool has_graph(const std::string& name) const;
        // ...
    };
}
```

The `Graphs` mixin is inherited by `Facade::Cluster`. It allows multiple
named graphs to coexist on a single cluster:

| Graph name | Type | Built by |
|---|---|---|
| `"relaxed_pid"` | `Graphs::Weighted` | `connect_graph_relaxed` visitor |
| `"steiner"` | `Graphs::Weighted` | `CreateSteinerGraph` visitor |
| `"pr"` (implicit) | `PR::Graph` | `TaggerCheckNeutrino` visitor |

---

## PR::Graph, PR::Vertex, PR::Segment

The final output of pattern recognition is stored as a `PR::Graph` — a
Boost `adjacency_list` where nodes hold `PR::Vertex` pointers and edges
hold `PR::Segment` pointers.

### PR::Graph type

```cpp
// clus/inc/WireCellClus/PRGraphType.h
using Graph = boost::adjacency_list<
    boost::vecS,        // edge container: vector
    boost::vecS,        // vertex container: vector
    boost::undirectedS, // undirected graph
    VertexBundle,       // node property: holds VertexPtr
    EdgeBundle          // edge property: holds SegmentPtr
>;
```

### PR::Vertex

```cpp
// clus/inc/WireCellClus/PRVertex.h
class Vertex
    : public Flagged<VertexFlags>        // flag bits
    , public Graphed<node_descriptor>    // knows its graph node
    , public HasCluster<Segment>         // back-pointer to Facade::Cluster*
{
    WCPoint m_wcpt;   // initial 3D position (from Steiner graph)
    Fit     m_fit;    // fitted 3D position (after TrackFitting)
    // ...
};

enum class VertexFlags {
    kUndefined      = 0,
    kNeutrinoVertex = 1<<1,   // this is the primary interaction vertex
};
```

**Key point**: `m_wcpt` is the initial estimate from the Steiner graph
topology; `m_fit` is the refined position after calling `fit_vertex()`.
The two can differ by a few mm.

### PR::Segment

```cpp
// clus/inc/WireCellClus/PRSegment.h
class Segment
    : public Flagged<SegmentFlags>    // kShowerTrajectory, kShowerTopology, etc.
    , public Graphed<edge_descriptor> // knows its graph edge
    , public HasCluster<Segment>      // back-pointer to Facade::Cluster*
    , public HasDPCs<Segment>         // holds DynamicPointClouds
{
    std::vector<WCPoint>  m_wcpts;       // original waypoints from Steiner path
    std::vector<Fit>      m_fits;        // fitted positions, one per waypoint
    int                   m_dirsign{0};  // +1/-1/0: direction relative to vertex
    std::shared_ptr<Aux::ParticleInfo> m_particle_info;  // PID result
    double                m_particle_score{100};          // PID score
    // ...
};

enum class SegmentFlags {
    kUndefined        = 0,
    kShowerTrajectory = 1<<1,   // EM shower (trajectory criterion)
    kShowerTopology   = 1<<2,   // EM shower (topology criterion)
    kAvoidMuonCheck   = 1<<3,   // exclude from muon test
    kFit              = 1<<4,   // fits are populated
};
```

### Graph construction functions

```cpp
// clus/inc/WireCellClus/PRGraph.h

// Create a vertex and add it to the graph
VertexPtr make_vertex(Graph& g, Args&&... args);

// Add a segment between two vertices
bool add_segment(Graph& g, SegmentPtr seg, VertexPtr v1, VertexPtr v2);

// Find the two endpoint vertices of a segment (ordered: first is closest
// to segment's initial wcpt)
std::pair<VertexPtr, VertexPtr> find_vertices(Graph& g, SegmentPtr seg);

// Given one endpoint vertex, find the vertex at the other end
VertexPtr find_other_vertex(Graph& g, SegmentPtr seg, VertexPtr vertex);
```

---

## PR::Fit (WCPoint and Fit)

The `WCPoint` type (defined in `PRCommon.h`) carries a 3D position plus
auxiliary wire-coordinate information used internally by `TrackFitting`.
The `Fit` struct carries the output of the fitting process:

```cpp
struct Fit {
    WireCell::Point point;   // fitted 3D position
    double range{0};         // arc length from segment start (cm)
    int    index{-1};        // index into the waypoint array
    bool   flag_fix{false};  // if true, held fixed during fit
    bool   flag_skip{false}; // if true, excluded from fit residual
    // ...
};
```

---

## PR::Shower

```cpp
// clus/inc/WireCellClus/PRShower.h
// clus/src/PRShower.cxx
```

A `PR::Shower` groups one or more `PR::Segment` objects into a coherent EM
shower. It stores:
- The shower's starting vertex (where it attaches to the main topology).
- The shower axis direction (from PCA of shower points).
- The total reconstructed energy (`kine_charge`).
- `ParticleInfo` holding the particle type (electron / photon) and mass.

Multiple showers can originate from the same vertex (e.g., the two photons
from a π0 decay). The π0 identification logic (in `NeutrinoTaggerPi0`) looks
for pairs of showers from the same vertex whose energies and opening angle
reconstruct the π0 mass.

---

## Summary: Object Lifetime

```
MultiAlgBlobClustering::operator()()
  │
  ├─ Deserialize → Facade::Ensemble (per-event, on stack)
  │                  └─ Facade::Grouping (per-event)
  │                       └─ Facade::Cluster (per-event)
  │                            ├─ Facade::Blob[] (per-event)
  │                            ├─ PointCloud::Dataset[] (per-event)
  │                            └─ PR::Graph (per-event, built during visit())
  │                                 ├─ PR::Vertex (shared_ptr, per-event)
  │                                 └─ PR::Segment (shared_ptr, per-event)
  │
  └─ Serialize → output TensorSet → destroyed at end of operator()
```

All PR objects (vertices, segments, showers) are owned by `shared_ptr` and
live for the duration of the event processing. They are serialized into the
output tensor set and then freed.
