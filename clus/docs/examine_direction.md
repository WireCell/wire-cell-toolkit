# Deep Dive: examine_direction — BFS Direction Propagation

**Source file:** `clus/src/NeutrinoVertexFinder.cxx` (function at line ~1018)

---

## Purpose

`examine_direction` is called after the primary vertex is chosen for a cluster.
Its job is to **re-propagate segment directions and PDG codes outward from the
primary vertex** using a BFS traversal of the PR graph.

This is necessary because `determine_direction()` (in `NeutrinoTrackShowerSep.cxx`)
assigns directions locally — each segment is processed independently using only
its own `dQ/dx` profile and its two vertex degrees. It does not know which vertex
is "upstream" (the primary). `examine_direction` corrects this by walking from the
known primary vertex outward, so every segment's `dirsign` is consistent with
particles flowing away from the vertex.

---

## Function Signature

```cpp
bool PatternAlgorithms::examine_direction(
    Graph& graph,
    VertexPtr vertex,           // the chosen primary vertex (or local vertex)
    VertexPtr main_vertex,      // the global primary vertex (may equal vertex)
    IndexedVertexSet& vertices_in_long_muon,
    IndexedSegmentSet& segments_in_long_muon,
    const ParticleDataSet::pointer& particle_data,
    const IRecombinationModel::pointer& recomb_model,
    bool flag_final);           // true = more aggressive direction reset
// Returns: true if any direction or PDG assignment was changed
```

---

## Pre-BFS Cluster Statistics

Before starting the BFS, the function collects cluster-level information that
influences the PDG assignment rules:

```cpp
double max_vtx_length;     // length of longest segment at the input vertex
double min_vtx_length;     // length of shortest segment at the input vertex
int    num_total_segments; // total number of segments in this cluster
bool   flag_only_showers;  // true if the cluster has no clear track segments
```

`flag_only_showers` is determined by scanning all vertices:
- If any vertex has at least one connected track segment and is NOT the endpoint
  of that segment alone → `flag_only_showers = false`.
- Additionally, if the input vertex has ≥3 segments, or 2 segments where both
  are longer than 15 cm → `flag_only_showers = false`.

This flag relaxes the electron-to-track reclassification guard (see below).

---

## BFS Mechanics

### Queue and state

```cpp
std::vector<pair<VertexPtr, SegmentPtr>> segments_to_be_examined;
std::set<VertexPtr>  used_vertices;
std::set<SegmentPtr> used_segments;
```

The queue is seeded with all segments connected to `vertex`. Each queue entry
is `(prev_vtx, current_seg)` — the vertex we came from, and the segment we are
about to process.

### Shower boundary rule

**BFS does not cross shower segments.** When a shower segment is encountered,
its direction is updated, but the BFS does not add the far-side vertex's
segments to the queue. This keeps shower sub-topologies independent; they are
handled later by `shower_clustering_with_nv`.

### Processing each segment

For each `(prev_vtx, current_seg)` dequeued:

1. **Check for incoming showers at `prev_vtx`**:
   Scan all already-processed (`used_segments`) siblings of `current_seg` at
   `prev_vtx`. If any is an incoming shower (pointing toward `prev_vtx`), set
   `flag_shower_in = true`. The "already processed" guard prevents
   order-dependent results: only segments from an earlier BFS level count.

2. **Check if direction needs updating**:
   Enter the update block if any of these is true:
   - `dirsign == 0` (undetermined)
   - `dir_weak == true` (weakly determined, set by the consistency fixers)
   - `is_shower == true` (showers always get direction re-set from context)
   - `flag_final == true` (forced reset pass)

3. **Set direction** from geometry:
   ```cpp
   bool flag_start = (dist(wcpts.front(), prev_vtx) < dist(wcpts.back(), prev_vtx));
   current_seg->dirsign(flag_start ? 1 : -1);
   ```
   `dirsign=+1` means the segment points away from `prev_vtx` (particle
   originates near `prev_vtx` and travels outward).

4. **Assign PDG** (see next section).

5. **Enqueue far-side vertex**: add all segments at `curr_vertex` (the vertex
   at the other end of `current_seg`) to the next BFS level, unless
   `curr_vertex` is already in `used_vertices`.

---

## PDG Assignment Rules

After direction is set, PDG is (re-)assigned. The rules apply in priority order:

### Rule 1 — `flag_shower_in` cases

If an incoming shower was found at `prev_vtx`:

| Sub-condition | Result |
|---|---|
| `dirsign == 0 && !is_shower` | PDG=11 (electron) |
| `length < 2 cm && !is_shower` | PDG=11 (short stub near shower = delta ray) |
| `abs(pdg)==13 OR pdg==0` | PDG=11 (muon/unknown adjacent to shower → electron) |
| Otherwise | No change |

Rationale: a segment starting just after an EM shower is very likely an
electron or part of the shower, not an independent muon.

### Rule 2 — No incoming shower, many daughter showers downstream

```cpp
if (num_daughter_showers >= 4 ||
    (length_daughter_showers > 50*cm && num_daughter_showers >= 2))
```
And additionally, the angle between `current_seg` and connected segments at
the far vertex is > 155° (near back-scatter), or > 135° and drift-parallel.

→ PDG=11

### Rule 3 — Reclassify electron back to track

```cpp
if (pdg == 11 && num_daughter_showers <= 2 && !flag_shower_in
    && !kShowerTopology && !kShowerTrajectory
    && length > 10*cm && !flag_only_showers
    && particle_score <= 100)
```
And the segment is geometrically straight (`direct_length / length > 0.93` OR
`direct_length >= 34 cm`):

→ PDG=0 (reclassify; the electron label from earlier was likely wrong)

### Rule 4 — Proton refinement

If `current_seg` already has PDG=2212 (proton) AND `flag_shower_in` AND
`num_daughter_showers == 0`:
- If the incoming shower's median `dQ/dx > 1.3 × MIP` → re-label the
  **incoming shower** as PDG=2212 (proton, not EM shower).
- Otherwise → re-label the incoming shower as PDG=211 (pion).
- Unset `kShowerTrajectory` and `kShowerTopology` flags on the relabelled
  shower.

### Rule 5 — Long muon detection

After BFS is complete, scan all segments connected to `vertex`. For each
segment with `dQ/dx < 1.3 × MIP`:
- Follow `find_cont_muon_segment()` to find a continuous chain of low-dQ/dx
  segments (collinear within an angle threshold).
- If the chain total length > 40 cm, add all segments and vertices of the
  chain to `segments_in_long_muon` and `vertices_in_long_muon`.

---

## `flag_final` Mode

When called with `flag_final=true` (during the final structure pass after the
global vertex is determined), the function:
- Enters the update block for **all** segments, including those with a strong
  non-weak direction.
- Applies more aggressive angle cuts for electron reclassification.
- Sets `dir_weak=true` on all newly assigned directions (they are treated as
  "soft" corrections rather than confident assignments).

This cleans up directions that were set before the global vertex was known
and are now inconsistent with the true event topology.

---

## Return Value

Returns `true` if **any** segment's direction or PDG was changed during this
call. The caller (`determine_main_vertex` and the DL path in
`TaggerCheckNeutrino`) may loop `examine_direction` until it returns `false`
to achieve convergence.

---

## Relationship to `determine_direction`

| | `determine_direction` | `examine_direction` |
|---|---|---|
| Context | No vertex known | Primary vertex known |
| Traversal | Per-segment, independent | BFS from primary vertex |
| Direction source | `dQ/dx` profile + vertex degree | Geometry from `prev_vtx` |
| PDG source | KS-test track PID | Contextual rules (shower_in, daughters) |
| Shower boundary | No (processes all segments) | Yes (stops at shower boundary) |
| Called from | `NeutrinoTrackShowerSep.cxx` | `NeutrinoVertexFinder.cxx` |
