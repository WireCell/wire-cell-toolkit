# Steiner Graph Port Review

**Reviewed:** 2026-04-10  
**Functions:** `create_steiner_graph`, `recover_steiner_graph`, `Create_steiner_tree`, `form_cell_points_map`, `establish_same_mcell_steiner_edges`, `remove_same_mcell_steiner_edges` and all callees.  
**Scope:** Functional equivalence, bugs, efficiency, determinism, multi-TPC/APA/face correctness.

---

## 1. Function Map (prototype → toolkit)

| Prototype (`prototype_base/pid/`) | Toolkit (`clus/src/`) | Status |
|---|---|---|
| `create_steiner_graph` — `PR3DCluster_steiner.h:10-75` | `Steiner::CreateSteinerGraph::visit` — `CreateSteinerGraph.cxx:62-279` | Ported |
| `recover_steiner_graph` — `PR3DCluster_steiner.h:77-180` | **MISSING** | Gap — see §4.1 |
| `Create_steiner_tree` — `PR3DCluster_steiner.h:182-709` | `Steiner::Grapher::create_steiner_tree` — `SteinerGrapher.cxx:20-124` + `Weighted::create_enhanced_steiner_graph` — `SteinerGrapher.cxx:898-1086` | Ported |
| `form_cell_points_map` — `PR3DCluster_steiner.h:1038-1056` | `Steiner::Grapher::form_cell_points_map` — `SteinerGrapher.cxx:545-587` | Ported |
| `establish_same_mcell_steiner_edges (flag=1)` — `PR3DCluster_graph.h:28-85` | `Steiner::Grapher::establish_same_blob_steiner_edges` — `SteinerGrapher.cxx:589-704` | Ported |
| `establish_same_mcell_steiner_edges (flag=2)` — `PR3DCluster_graph.h:86-130` | `Weighted::establish_same_blob_steiner_edges_steiner_graph` — `SteinerGrapher.cxx:843-894` | Ported |
| `remove_same_mcell_steiner_edges` — `PR3DCluster_graph.h:10-26` | `Steiner::Grapher::remove_same_blob_steiner_edges` — `SteinerGrapher.cxx:707-738` | Ported |
| `find_steiner_terminals` + `find_peak_point_indices` — `PR3DCluster_steiner.h:711-780` | `Grapher::find_steiner_terminals` `:516-541`, `Grapher::find_peak_point_indices` `:310-513` | Ported |

---

## 2. Bugs

### 2.1 Missing null guard in `establish_same_blob_steiner_edges_steiner_graph` (FIXED)

**File:** `SteinerGrapher.cxx:857`  
**Prototype ref:** not applicable (new function)

```cpp
// Before (no null check):
const auto* blob = nodes[blob_idx]->value.facade<Facade::Blob>();
cell_points_map[blob].insert(new_index);  // crash / nullptr key if facade fails
```

The member `form_cell_points_map` (`:573-575`) already guards against `blob == nullptr`. The free function `establish_same_blob_steiner_edges_steiner_graph` did not, which could crash or insert a `nullptr` key on bad blob nodes.

**Fix applied:** Added `if (!blob) continue;` after the facade cast.

---

## 3. Determinism Issues (FIXED)

### 3.1 Pointer-keyed `blob_vertex_map`

**Files:** `SteinerGrapher.h:67`, `SteinerGrapher.cxx:545-587, 589-704, 843-894, 516-541`

The original typedef was:
```cpp
using blob_vertex_map = std::map<const Facade::Blob*, vertex_set>;
```
Because `std::map` orders by key, and blob pointer values depend on heap layout (ASLR / allocator state), the iteration order over blobs varied run-to-run. This caused `boost::add_edge` calls to be inserted in different orders, changing edge descriptor identities and subtly affecting downstream Dijkstra/Voronoi tie-breaking.

**Fix applied:** Changed key to `size_t` (blob node index from `sv.nodes()`). The node vector is traversal-ordered and stable across `sv3d()` calls:
```cpp
using blob_vertex_map = std::map<size_t, vertex_set>;  // key = blob node index
```
All call sites updated accordingly. A local `build_blob_node_index_map()` helper was added to `find_peak_point_indices(vector<Blob*>)` to map a `Blob*` → node index for lookups.

### 3.2 Non-deterministic sort of `tree_edges` on Boost edge descriptors

**File:** `SteinerGrapher.cxx:990`

```cpp
std::sort(tree_edges.begin(), tree_edges.end());  // edge_type uses pointer comparison
```

For `adjacency_list<vecS, vecS, undirectedS, ...>`, `edge_type` (`edge_desc_impl`) embeds an internal property `void*`, so `operator<` compares pointer addresses. The resulting order is run-to-run non-deterministic, which affects the order edges are added to `result.graph` and therefore the edge numbering of the reduced Steiner graph.

**Fix applied:** `tree_edges` changed to `std::vector<std::pair<vertex_pair, edge_type>>`. Sort and unique use the `vertex_pair` (stable `size_t` indices) as the comparison key:
```cpp
std::sort(tree_edge_pairs.begin(), tree_edge_pairs.end(),
    [](const auto& a, const auto& b){ return a.first < b.first; });
std::unique(..., [](const auto& a, const auto& b){ return a.first == b.first; });
```
Edges are then added to `result.graph` in deterministic vertex-index order.

### 3.3 `edge_set` for tracked edges changed to `std::vector`

**Files:** `SteinerGrapher.h:229,235`, `SteinerGrapher.cxx:602,693-695`

`std::set<edge_type>` (alias `edge_set`) used `edge_type` comparison (pointer-based) as the ordering key. Since `remove_same_blob_steiner_edges` only iterates and removes (order-insensitive), using a `std::vector` is both cleaner and removes the pointer-ordering dependency.

**Fix applied:** `m_added_edges_by_graph` changed from `std::map<std::string, edge_set>` to `std::map<std::string, std::vector<edge_type>>`. Local `added_edges` in `establish_same_blob_steiner_edges` changed from `edge_set` to `std::vector<edge_type>`.

---

## 4. Functional Gaps

### 4.1 `recover_steiner_graph` is absent (documented — deferred)

The prototype's `recover_steiner_graph` (`PR3DCluster_steiner.h:77-180`) performs a second Voronoi+Dijkstra pass over the already-reduced `graph_steiner`, followed by `boost::kruskal_minimum_spanning_tree`, and populates two sets:

- `steiner_graph_terminal_indices` — all terminal indices from `flag_steiner_terminal`
- `steiner_graph_selected_terminal_indices` — endpoints of the Kruskal MST edges

It is called unconditionally in `wire-cell-graph.cxx` and `wire-cell-tracking.cxx` after `create_steiner_graph`, and is commented out in `wire-cell-stm.cxx`, `wire-cell-prod-stm.cxx`, and `wire-cell-prod-stm-port.cxx`. It is **not** called from `NeutrinoID.cxx`.

**Status in toolkit:** `m_steiner_graph_terminal_indices` is declared in `SteinerGrapher.h:296` but is never populated. No toolkit consumer currently reads `steiner_graph_terminal_indices` or `steiner_graph_selected_terminal_indices` (grep confirmed). All existing NeutrinoTagger/pr code uses only `get_flag_steiner_terminal()`, which IS populated correctly.

**Action:** No port needed at this time. If a future consumer needs the selected-terminal index sets (e.g., a port of `wire-cell-graph.cxx` output or `wire-cell-tracking.cxx` track-fitting), implement `Grapher::recover_steiner_graph(const std::string& steiner_graph_name)` mirroring `PR3DCluster_steiner.h:77-180` (multi-source Dijkstra via PAAL `make_nearest_recorder` → `TerminalGraph` → `kruskal_minimum_spanning_tree` → walk `vpred` back to terminals).

---

## 5. Efficiency Improvements (APPLIED)

### 5.1 Early exit for blobs with no terminals in `establish_same_blob_steiner_edges`

**File:** `SteinerGrapher.cxx:630-690`

The O(M²) inner pair loop was entered even for blobs with zero Steiner terminals, where every pair failed the `flag_index1 || flag_index2` test. Added a pre-check:
```cpp
bool blob_has_terminal = false;
for (auto idx : points_vec) {
    if (terminal_set.count(idx)) { blob_has_terminal = true; break; }
}
if (!blob_has_terminal) continue;
```
For typical events, most blobs are non-terminal, so this skip avoids the majority of O(M²) iterations.

### 5.2 Single-terminal blobs — iterate only M-1 pairs

When exactly one point in a blob is a terminal, the only edges that would be added are those between the one terminal and the M-1 non-terminals (since pairs of two non-terminals are skipped). Rather than run the full O(M²) loop, a special path builds the edges in O(M):
```cpp
// find the single terminal index
size_t term_k = ...;
for (size_t j = 0; j < M; ++j) {
    if (j == term_k) continue;
    // add edge between term_k and j at weight 0.9 * distance
}
```

### 5.3 Symmetric loop in `find_peak_point_indices` connected-component step

**File:** `SteinerGrapher.cxx:456-465`

The double loop `for j in N: for k in N` checked `boost::edge()` symmetrically and added edges twice (both `(j,k)` and `(k,j)` for undirected). Changed to `k = j+1` to visit each pair once, halving the number of `boost::edge()` calls and `boost::add_edge()` calls.

### 5.4 Redundant `find()` in `form_cell_points_map` removed

**File:** `SteinerGrapher.cxx:578`

Changed:
```cpp
if (cell_points.find(blob_node_idx) == cell_points.end()) {
    cell_points[blob_node_idx] = vertex_set();
}
cell_points[blob_node_idx].insert(point_idx);
```
to:
```cpp
cell_points[blob_node_idx].insert(point_idx);
```
`operator[]` default-constructs the `vertex_set` on first access, making the `find()` redundant.

---

## 6. Semantic Alignment Check

### 6.1 Wire-range convention (not a bug)

`Facade_Cluster.cxx:3239-3241` (`check_wire_ranges_match`) uses strict `<` for the upper wire bound:
```cpp
if (current_wire_u >= u_min && current_wire_u < u_max && ...)
```
The prototype uses inclusive `<=` with `back()->index()`. This is a **deliberate convention difference**: the toolkit treats wire ranges as half-open `[min, max)` while the prototype uses inclusive `[min, max]`. Do not change this to `<=`.

### 6.2 Path-based 2D/3D filter matches prototype

`filter_by_path_constraints` (`SteinerGrapher.cxx:159-249`) replicates the prototype's 0.6 cm step interpolation along `path_wcps`, the same 6.0 cm 3D threshold, and the same 1.8 cm 2D threshold. Thresholds verified against `PR3DCluster_steiner.h:351-371`. **Matches.**

### 6.3 Voronoi reduction — no Kruskal (correct)

Prototype `Create_steiner_tree` has Kruskal commented out (lines 489-490 of `PR3DCluster_steiner.h`); it uses all terminal-pair edges from `map_saved_edge` directly. Toolkit `create_enhanced_steiner_graph` matches this behavior. **Correct.**

### 6.4 Edge weighting in post-reduction `establish` — matches prototype

Prototype flag=2 path (`PR3DCluster_graph.h:86-126`): raw distance, no 0.8/0.9 factor. Toolkit `establish_same_blob_steiner_edges_steiner_graph` `:882-889`: raw distance. **Matches.**

### 6.5 `calc_charge_wcp` charge threshold

Both prototype and toolkit use `Q0=10000`, `factor1=0.8`, `factor2=0.4`, `charge_threshold=4000`. Verified at `SteinerGrapher.cxx:85-89`. **Matches.**

---

## 7. Multi-TPC / Multi-APA / Multi-Face Audit

The prototype is single-TPC only. The toolkit port must support multiple APAs and faces (e.g., PDHD with up to 4 TPC faces per cluster).

| Layer | Multi-face handling | Status |
|---|---|---|
| `form_cell_points_map` | Uses `sv3d()` — one 3D k-d tree spanning all faces. Blobs are face-local (each has one `wpid()`), so blob-grouping is face-safe. | OK |
| `establish_same_blob_steiner_edges` | All geometry via `m_cluster.point3d()`. Blob grouping from `sv3d()` node index — face-safe. | OK |
| `establish_same_blob_steiner_edges_steiner_graph` | Same: `sv3d()` node index + geometry from point cloud coordinates. | OK |
| `filter_by_path_constraints` | Correctly enumerates unique `WirePlaneId`s via `m_cluster.wpids_blob()`, builds per-(face,apa) `DynamicPointCloud` with per-face wire angles from `IDetectorVolumes`. Iterates over `std::set<WirePlaneId>` (value-ordered, deterministic). **Only place in Steiner code with explicit multi-face logic.** | OK |
| `is_point_spatially_related_to_time_blobs` | `ref_time_blob_map` is keyed `apa → face → time_slice`; point's own `wpid()` selects the right bucket. | OK |
| `CreateSteinerGraph::visit` | Iterates `grouping.children()` (vector, deterministic). Per-face handling is inside Grapher/Cluster helpers. | OK |
| `create_enhanced_steiner_graph` | Uses a single base graph `ctpc_ref_pid`. Whether this graph spans multiple faces correctly depends on `Cluster::find_graph("ctpc_ref_pid", ...)` — upstream of Steiner code, out of scope for this review. Should be smoke-tested on PDHD 4-face clusters. | Needs test |

**Summary:** The port correctly delegates per-face geometry to `filter_by_path_constraints` (using `IDetectorVolumes`) and to the `Cluster` facade (wire-index and time-slice checks). No per-face loops need to be added at the `Grapher` level. The blob-grouping approach is inherently face-safe.

---

## 8. Minor / Cosmetic Notes

- `steiner_flags_uint8` in `create_enhanced_steiner_graph` (`:1033`) uses `std::vector<int>` though the comment says uint8_t. Should be `std::vector<uint8_t>` to save 4× memory for large point clouds.
- `m_steiner_graph_terminal_indices` declared in `SteinerGrapher.h:296` is never written. Can be removed or repurposed as the backing store when `recover_steiner_graph` is ported.
- Commented-out debug blocks in `CreateSteinerGraph.cxx:280-348` can be cleaned up in a follow-up.

---

## 9. Verification Plan

1. **Build:** `ninja -C build clus` — must compile without errors or warnings in the changed files.
2. **Bit-reproducibility:** Run the Steiner stage twice on the same input; confirm `steiner_graph` vertex/edge counts and `flag_steiner_terminal` vector are identical across runs.
3. **Prototype parity (single-TPC):** On an ICARUS or uboone sample, compare `steiner_pc` point count, `flag_steiner_terminal` true-count, and `steiner_graph` edge count against the prototype's output for the same cluster. Allow for small float-rounding differences in edge weights; point/vertex counts should match exactly.
4. **Multi-face smoke test:** On a PDHD 4-face cluster, confirm no crashes, no empty `steiner_graph` where the prototype would produce one, and that `filter_by_path_constraints` branches through all unique `WirePlaneId`s.
