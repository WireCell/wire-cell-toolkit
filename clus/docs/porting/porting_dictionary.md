#  Between Clustering and Trajectoryand dQ/dx fitting (Current)

Missing info or questions are marked as: **bold** :warning:.  Some non-answers  are marked with :question: if I (bv) don't know or am unsure.

## [Opflash](./QLBundles/Opflash.md) (WCP) vs. optical arrays (WCT)

Covered in next item

## [FlashTPCBundle](./QLBundles/Bundle.md) (WCP) vs. optical arrays and flash association (WCT)

See link.

## [Examine_bundles](https://github.com/BNLIF/wire-cell-2dtoy/blob/master/docs/ExamineBundles.md) is a function that we need to port, [Original ExamineBundles Code](https://github.com/BNLIF/wire-cell-2dtoy/blob/a30305052fc54bbbbbd826b096066d6e8777b54d/src/ExamineBundles.cxx)
### `bundle->get_main_cluster()`  (WCP) 

WCT: see next.

### `bundle->get_other_clusters()` (WCP)

WCT: PC-tree / n-ary tree supports a 1-to-many "separation" and an many-to-1
"merge".  When used via the `Facade` classes from `clus/` a cluster may be
separated into a number of clusters.  Or, vice versa, a number of clusters may
be merged into a single cluster.

The user can represent "main" and "other" clusters in a few ways.

1. The user may simply keep "main" as a distinct cluster object and "other" in
   some collection of cluster objects.

2. The user may merge "main" and "other" clusters into a single cluster and have
   that merge recorded in the resulting (composite) cluster.  In that record,
   cluster "0" can be chosen to represent the "main" cluster.  This merge can be
   undone to get back to representation 1.

This and other related tree features are described in this presentation:

https://www.phy.bnl.gov/~bviren/talks/wire-cell/2024/12/10/narysepclus.pdf



### bundle->get_main_cluster()->get_mcells() (WCP) vs. [loop blob](https://github.com/WireCell/wire-cell-toolkit/blob/apply-pointcloud/clus/src/clustering_separate.cxx#L1377) (WCT)

See also.

https://www.phy.bnl.gov/~bviren/talks/wire-cell/2024/12/10/narysepclus.pdf



### Create a new temp cluster, PR3DCluster *new_cluster = new PR3DCluster(cluster_id) (WCP) vs. WCT

A WCT `Cluster` is a "facade" which wraps a n-ary tree node (with value type of
`Points`).  You can not create a temp `Cluster` per se but you can create a temp
node and get a `Cluster` facade from it.  This largely uses support in
`NaryTreeFacade` from which `Cluster` inherits.

A **truly temporary** node can be created **on the stack** and the cluster facade obtained with:

```c++
Points::node_t node(Points(...));
Cluster& cluster = node.facade<Cluster>()
```

This node **can not** be added to a tree.  

Make a temporary node that one day may be added to a tree, create it on the
heap, and preferably as a `unique_ptr`.

```c++
auto node = std::make_unique<Points::node_t>(Points(...));
// later you can add this to a parent ndoe
root->insert(node);
assert(node == nullptr);
```

If you know ahead of time the node is destined for the tree, you can make it in place

```c++
auto node = root->insert(Points(...)))
```

In both cases you can the `Cluster` facade regardless of heap vs stack:

```c++
auto& cluster = node->facade<Cluster>()
```

or if you have a `Grouping` facade object you can add a node it and get its facade without knowing anything about the underlying n-ary tree:

```c++
auto& cluster = grouping.make_facade();
```

More info on these items in WCT:
- https://github.com/WireCell/wire-cell-toolkit/blob/apply-pointcloud/util/test/doctest-pointtree-example.org  walks through many example using nodes.
- https://www.phy.bnl.gov/~bviren/talks/wire-cell/2024/12/10/narysepclus.pdf includes some info on facades.



### Add every blob to this cluster, `new_cluster->AddCell(blob, blob->GetTimeSlice())`  (WCP) vs. WCT

:question: This WCP example here implies to me you somehow have a set of "free
blob nodes".  I do not think that is usual in WCT.  Instead, all blob nodes
already reside in cluster nodes due to `PointTreeBuilding`.  And so the closest
WCT-equivalent is a "merge" as described above.


### std::vector<SMGCSelection> sep_mcells = new_cluster->Examine_graph(ct_point_cloud)  [code](https://github.com/BNLIF/wire-cell-data/blob/d5748d87c3113efcb44eed237bb48a10d60002d9/src/PR3DCluster.cxx#L2332) (WCP) vs. similar to [Create_graph](https://github.com/WireCell/wire-cell-toolkit/blob/apply-pointcloud/clus/src/Facade_Cluster.cxx#L1444) (WCT)

###  [Connect_graph_overclustering_protection](https://github.com/BNLIF/wire-cell-data/blob/d5748d87c3113efcb44eed237bb48a10d60002d9/src/PR3DCluster.cxx#L1853) (WCP)  vs. not existing, but example in [connect_graph](https://github.com/WireCell/wire-cell-toolkit/blob/apply-pointcloud/clus/src/Facade_Cluster.cxx#L1444) (WCT)

###  For each sep_mcells, we will form a new set of clusters, we need to find the cluster that overlapped the most with the original main cluster, and assign it as the new MAIN CLUSTER for the new bundle (achieved through compare blobs) (WCP) vs. WCT

Above we described how to interpret cluster "0" in a merged cluster as "main".
I think this same idea may apply here.

For example and algorithm may:

1. Consider a merged cluster.
2. Unmerge to get "main" cluster and the "other" clusters.
3. Example each of the "other" clusters for overlap with "main".
4. Remerge with the new main as cluster "0".


## [PRCluster->Create_Graph()](https://github.com/BNLIF/wire-cell/blob/master/uboone_nusel_app/apps/prod-wire-cell-matching-nusel.cxx#L817) (WCP) vs. [Create_graph](https://github.com/WireCell/wire-cell-toolkit/blob/apply-pointcloud/clus/src/Facade_Cluster.cxx#L1444) (WCT)

In addition, you can store a graph onto the PC-tree as a 2-array PC holding edges as a pair of (tail,head) array of indices.



## [cluster->get_highest_lowest_wcps()](https://github.com/BNLIF/wire-cell/blob/master/uboone_nusel_app/apps/prod-wire-cell-matching-nusel.cxx#L819C103-L819C130) (WCP) vs. [get_ghiehst_lowst_points](https://github.com/WireCell/wire-cell-toolkit/blob/apply-pointcloud/clus/src/Facade_Cluster.cxx#L1241) (WCT)

# WCT function only returns points, but we need the index to do the shortest path (following two functions), how to get point index from points?

A point has three indices:

- Its **absolute** index in an overall point cloud
- The **major** index of the node that provides the local PC with the point.
- The **minor** index of the the point in the local PC. 

It is possible to navigate between the absolute index and the pair of
major/minor indices.  But, perhaps there needs to be more helper methods in the
facades to do this in more ways.

Some facade methods that use this now:

- `Cluster::get_closest_blob`
- `Cluster::examine_graph`

More info:

- https://github.com/WireCell/wire-cell-toolkit/blob/apply-pointcloud/util/test/doctest-pointtree-example.org#scoped-k-d-tree examples using major/minor indices

- https://github.com/WireCell/wire-cell-toolkit/blob/apply-pointcloud/util/docs/nfkdvec-tree.org more examples.

## [cluster->->dijkstra_shortest_paths(wcps.first)](https://github.com/BNLIF/wire-cell/blob/master/uboone_nusel_app/apps/prod-wire-cell-matching-nusel.cxx#L822C25-L823C58) (WCP) vs. [dijkstra_shortest_paths using point index](https://github.com/WireCell/wire-cell-toolkit/blob/apply-pointcloud/clus/src/Facade_Cluster.cxx#L2631) (WCT)

## [cluster->cal_shortest_path(wcps.second)](https://github.com/BNLIF/wire-cell/blob/master/uboone_nusel_app/apps/prod-wire-cell-matching-nusel.cxx#L823) (WCP) vs. [cal_shortest_path using point index](https://github.com/WireCell/wire-cell-toolkit/blob/apply-pointcloud/clus/src/Facade_Cluster.cxx#L2664) (WCT)

## [Improve_PR3DCluster](https://github.com/BNLIF/wire-cell-2dtoy/blob/master/docs/Improve_PR3DCluster.md) is another function that we need to port. [Original Improve_PR3DCluster Code](https://github.com/BNLIF/wire-cell-2dtoy/blob/master/src/ImprovePR3DCluster.cxx)

### From existing clusters --> blobs --> time, fired channels --> get activities to be used to redo tiling (WCP) vs. WCT.

WCT: an algorithm that will "redo tiling" would make use of the `RayGrid::Tiling` class from `WireCellUtil/RayTiling.h` most usefully by calling 

https://github.com/WireCell/wire-cell-toolkit/blob/apply-pointcloud/util/inc/WireCellUtil/RayTiling.h#L189

Be sure that the first two "activities" span the horizontal/vertical "layers"
and then layers 2,3,4 are activities built from the wire bounds derived from
values in the PCs on the set of blobs in your clusters.

### Access shortest path std::list<WCPointCloud<double>::WCPoint>& wcps = cluster->get_path_wcps() (WCP) vs. [index of points along the path](https://github.com/WireCell/wire-cell-toolkit/blob/apply-pointcloud/clus/src/Facade_Cluster.cxx#L2689) (WCT)

### for a path point, use ct_point_cloud to convert the point into ch vs. time, and then judge if activities are available (WCP) vs. **xxx** :warning:  (WCT)

:question: I think this may be fore Haiwang.  I don't know enough about how the CT PC equivalent is stored/used.

### [update the activities according to the path_point's properties](https://github.com/BNLIF/wire-cell-2dtoy/blob/master/src/ImprovePR3DCluster.cxx#L136) (WCP) vs.  **xxx** :warning:  (WCT)

### [redo tiling using the newly created activities](https://github.com/BNLIF/wire-cell-2dtoy/blob/master/src/ImprovePR3DCluster.cxx#L203) (WCP) vs.  WCT

WCT: same answer as above.  Use `RayGrid::make_blobs()`.

### Compare the newly created blobs vs. the (old) existing blobs using the Overlap_fast function and judge if the blob should be kept or not, then create a new cluster from the remaining blobs (WCP)  [overlap_fast](https://github.com/WireCell/wire-cell-toolkit/blob/apply-pointcloud/clus/src/Facade_Cluster.cxx#L513)  (WCT)

WCT: as above, a cluster "merge" is likely used as "loose blobs" are not currently something I think would exist.

:question: Perhaps "redo the tiling" needs to be some larger chain which includes running `PointTreeBuilding` on the newly made blobs.


## [map from old to new cluster](https://github.com/BNLIF/wire-cell/blob/master/uboone_nusel_app/apps/prod-wire-cell-matching-nusel.cxx#L831) (WCP) ss. WCT

WCT: I can only think of some options so I add: :question:

- If the map is merely internal to an algorithm, use whatever you want (`std::map<Cluster*,Cluster*>`)

- If "old" and "new" cluster nodes are all held in the grouping, a cluster-to-cluster map can be made as a "edge" array and stored in a PC on the grouping node.

- An "old" and a "new" grouping (two trees) can be defined.  We may define a new root node (a "versions" node?) which has grouping nodes as children.  We can extend the facade to cover this new root.  If the clusters of the two groupings are 1-to-1 then their map is implicit.  Otherwise, you can store an "edge" array on the new root node.


# Trajectory Fitting and Pattern Recognition

The final portion of porting covers the transformation from clusters of blobs and their points to a "particle flow" description.  This section lists WCP data and algorithm summaries and design and links to development notes for the WCT equivalents (generally under [tjft/](./tjft/) directory).


## WCP Data

- [wcp-data-notes](./tjft/wcp-data-notes.org) notes and questions from reviewing the links to data summaries below.


### [ProtoSegment](https://github.com/BNLIF/wire-cell-pid/blob/537a3fd17f8a7b3cf5412594c4cc1775cb/docs/protosegment.md) (WCP) vs. **`WireCell::Clus::PR::Segment`** (WCT)

`PR::Segment` (`clus/inc/WireCellClus/PRSegment.h`, `clus/src/PRSegment.cxx`) is the WCT equivalent of `WCPPID::ProtoSegment`.  Key mapping:

| WCP (ProtoSegment) | WCT | Notes |
|---|---|---|
| `ProtoSegment(id, cluster, flag_shower)` | `PR::Segment` constructed via `PR::add_segment(graph, seg, vtx1, vtx2)` | Graph-owned; no raw constructor |
| `m_fit_pt_vec` | `PR::Segment::m_fits` (`std::vector<PR::Fit>`) | Same data; WCT adds `paf{apa,face}` per fit |
| `get_point_vec()` | `seg->fits()` | |
| `set_fit_associate_vec(pts, skip, idx)` | `seg->set_fit_associate_vec(fits, dv, cloud_name)` | Now takes by value + `IDetectorVolumes` |
| `get_closest_wcpt(point)` | `segment_get_closest_point(seg, point, cloud_name)` | In `PRSegmentFunctions.cxx` |
| `get_flag_shower()` | `seg->flags_any(SegmentFlags::kShowerTrajectory \| kShowerTopology) \|\| std::abs(pdg) == 11` | **THREE disjuncts, not two.** `ProtoSegment.cxx:1305` is `flag_shower_trajectory \|\| flag_shower_topology \|\| get_flag_shower_dQdx()`, and `get_flag_shower_dQdx()` (`:1309`) is *only* `fabs(particle_type)==11` — despite the name it tests nothing about dQ/dx. An earlier revision of this row omitted the third term, and doc pr/33 P6 traces two live toolkit sites to that omission |
| `get_flag_shower_dQdx()` | `std::abs(seg->particle_info()->pdg()) == 11` | Name is misleading in the prototype too; carries no dQ/dx test |
| `get_direct_length()` | `segment_track_direct_length(seg)` | In `PRSegmentFunctions.cxx` |
| `get_length()` | `segment_track_length(seg)` | |
| `search_kink(start, cloud, threshold)` | `segment_search_kink(seg, start, cloud, threshold)` | |
| `get_closest_2d_dis(point, apa, face)` | `segment_get_closest_2d_distances(seg, point, apa, face, cloud)` | |
| `update_point_cloud(map_seg_vtxs)` | `create_segment_point_cloud(seg, points, dv, cloud)` | Maps replaced by graph |

Full function-by-function table: `clus/docs/patternrecognition/prvertex_prsegment_prshower_review.md §1`.

PID / direction-finding / kinematics algorithms (`do_track_comp`, `eval_ks_ratio`, `do_track_pid`, `cal_kine_dQdx`, `cal_kine_range`, `cal_4mom`, `determine_dir_track`, `determine_dir_shower_trajectory`, `determine_shower_direction`, `break_segment_at_point`) — detailed review in `clus/docs/patternrecognition/pid_direction_kinematics_review.md`.

### [ProtoVertex](https://github.com/BNLIF/wire-cell-pid/blob/537a3fd17f8a7b3cf5412594c4cc1775cb/docs/protovertex.md) (WCP) vs. **`WireCell::Clus::PR::Vertex`** (WCT)

`PR::Vertex` (`clus/inc/WireCellClus/PRVertex.h`) is the WCT equivalent of `WCPPID::ProtoVertex`.

| WCP (ProtoVertex) | WCT | Notes |
|---|---|---|
| `ProtoVertex(id, wcpt, flag_vertex)` | `std::make_shared<PR::Vertex>(wcpt)` + `PR::add_vertex(graph, vtx)` | Graph-owned |
| `get_wcpt()` | `vtx->wcpt()` | |
| `get_fit_pt()` | `vtx->fit().point` | |
| `get_cluster()` | `vtx->cluster()` | Via `HasCluster<Vertex>` mixin |
| `add_segment(seg)` | `PR::add_segment(graph, seg, vtx1, vtx2)` | Graph edge; no per-vertex list |
| `get_all_segments()` | `ordered_out_edges(graph, vtx->get_descriptor())` (caller iterates) | |

### [WCShower](https://github.com/BNLIF/wire-cell-pid/blob/537a3fd17f8a7b3cf5412594c4cc1775cb/docs/wcshower.md) (WCP) vs. **`WireCell::Clus::PR::Shower`** (WCT)

`PR::Shower` (`clus/inc/WireCellClus/PRShower.h`, `clus/src/PRShower.cxx`) is the WCT equivalent of `WCPPID::WCShower`.

| WCP (WCShower) | WCT | Notes |
|---|---|---|
| `WCShower(seg, map_vtx_segs, map_seg_vtxs)` | `std::make_shared<PR::Shower>(graph)` then `shower->set_start_segment(seg)` | Maps eliminated |
| `get_start_vertex()` / `get_start_segment()` | `shower->start_vertex()` / `shower->start_segment()` | |
| `calculate_kinematics(particle_data, recomb)` | `shower->calculate_kinematics(particle_data, recomb)` | DI replaces singleton |
| `get_kine_range()` / `get_kine_dQdx()` / `get_kine_best()` | `shower->get_kine_range()` / `get_kine_dQdx()` / `get_kine_best()` | |
| `get_start_point()` / `get_end_point()` | `shower->data.start_point` / `shower->data.end_point` | Direct struct access |
| `fill_maps(map_vtx_shower, map_seg_shower)` | `shower->fill_maps()` (no out-params; callers iterate the view) | API change |
| `add_segment(seg, maps)` | `shower->add_segment(seg)` | Maps replaced by graph |
| `add_shower(other)` | `shower->add_shower(other)` | Batched DPC merge (improvement) |
| `complete_structure_with_start_segment(maps, used)` | `shower->complete_structure_with_start_segment(used_segs)` | |
| `get_total_length()` | `shower->get_total_length()` | Cached after S2 fixes |
| `get_num_segments()` | `shower->get_num_segments()` | O(1) via Boost `edges().size()` |

Full function-by-function table: `clus/docs/patternrecognition/prvertex_prsegment_prshower_review.md §1`.

### [Steiner Tree](https://github.com/BNLIF/wire-cell-pid/blob/537a3fd17f8a7b3cf5412594267c14c4cc1775cb/docs/PR3DCluster_steiner.md) (WCP) vs. **`WireCell::Clus::Steiner::Grapher`** (WCT)

`Steiner::Grapher` (`clus/src/SteinerGrapher.{h,cxx}`) plus the free functions in
`WireCell::Clus::Graphs::Weighted` are the WCT equivalent of
`WCPPID::PR3DCluster`'s Steiner half; the orchestration that WCP does inline in
`create_steiner_graph` lives in the `CreateSteinerGraph` visitor.

| WCP (`PR3DCluster_steiner.h` / `_graph.h`) | WCT | Notes |
|---|---|---|
| `create_steiner_graph(ct_point_cloud, gds, …)` | `Steiner::CreateSteinerGraph::visit` → `process_cluster_steiner` lambda | Component, not a member |
| `Improve_PR3DCluster_2(...)` + `calc_sampling_points` + `Create_point_cloud` | `m_grapher_config.retile->mutate(*src->node())` (`RetileCluster`) | The retiled cluster is a real grouping child, **not** deleted mid-algorithm |
| `Create_steiner_tree(pc_steiner, flags, gds, old_mcells, flag_path, disable_dead_mix_cell)` | `Grapher::create_steiner_tree(ref_cluster, path_indices, graph, steiner_graph, disable_dead_mix_cell, steiner_pc)` | Reference cluster passed explicitly instead of `old_mcells` |
| `find_steiner_terminals(gds, disable_dead_mix_cell)` | `Grapher::find_steiner_terminals(graph, disable_dead_mix_cell[, cell_points_map])` | Map may be passed in to avoid recomputing |
| `find_peak_point_indices(...)` | `Grapher::find_peak_point_indices(...)` (blob + live overloads) | |
| `calc_charge_wcp(wcp, gds, disable_dead_mix_cell)` | `Facade::Cluster::calc_charge_wcp(point_index, charge_cut, disable_dead_mix_cell)` | |
| `establish_same_mcell_steiner_edges(gds, ddmc)` — **flag=1** | `Grapher::establish_same_blob_steiner_edges(graph, ddmc)` + `form_cell_points_map()` | Both group by the **retiled** blob |
| `establish_same_mcell_steiner_edges(gds, true, 2)` — **flag=2** | `Graphs::Weighted::establish_same_blob_steiner_edges_steiner_graph(result, cluster)` | Grouping key differs — see K5 |
| `remove_same_mcell_steiner_edges()` | `Grapher::remove_same_blob_steiner_edges(graph)` | |
| `WCP::ToyPointCloud* point_cloud_steiner` (owns `WCPoint::mcell`) | a `"steiner_pc"` `PointCloud::Dataset` on the cluster | **No blob pointer per point** — see K5 |
| `std::vector<bool> flag_steiner_terminal` | `flag_steiner_terminal` column on `steiner_pc` | |
| `get_two_boundary_wcps(2)` / `get_extreme_wcps(2)` | `Facade::Cluster::get_two_boundary_steiner_graph_idx(...)` | Different construction — see K6 |
| `WCPPaal/steiner_tree_greedy.h` | `clus/src/PAAL.h` + `Graphs::Weighted::voronoi` | Character-for-character port |
| `MCUGraph` (`listS`/`setS`, `float` weights) | `Graphs::Weighted::graph_type` (`vecS`, `double` weights) | See K1, K2 |

Full audit, with per-line anchors on both sides:
`sbnd_xin/docs/pr/29_steiner-graph-build-port-audit.md` (wcp-porting-validation).

#### Known divergences — do NOT "correct" these back

**K1. Edge weights are `double` in WCT, `float` in WCP.** Deliberate. It blocks
bit-identicality with WCP output and is not a bug. (pr/29 D4)

**K2. Tree-edge dedup is by vertex pair, not by edge descriptor.** WCP sorts
edge descriptors whose `operator<` compares the edge-property **pointer**
(`PR3DCluster_steiner.h:505-506`), so its dedup result is allocation-order
dependent. WCT keys on `(source, target)`. **WCT is the deterministic side; keep
it.** The same applies to WCP's `std::map<SlimMergeGeomCell*, …>` in both
same-blob passes (`PR3DCluster_graph.h:60`, `:99`), which is *iterated* — the
edge set is order-independent but the insertion order is not. WCT's
`form_cell_points_map` keys on a `size_t` blob node index for exactly this
reason. (pr/29 D6, §6)

**K3. The Steiner terminal filter takes ±1 wire of slack; `get_extreme_wcps`
takes none.** WCP uses two different tolerances at two call sites
(`PR3DCluster_steiner.h:285-290` vs `PR3DCluster_path.h:111-119`); WCT serves
both from one helper, `Cluster::check_wire_ranges_match`, so the slack is a
**per-call-site parameter** (`terminal_wire_tol`, default 0), never an edit to
the shared body. **The arithmetic is the trap:** WCT blob wire ranges are
half-open `[min, max)` while WCP's are inclusive `[low, high]` with
`high == max - 1`, so WCP's `index <= high + 1` translates to `index < max + 1`,
**not** `<= max + 1` — the literal transcription is two wires loose on the high
side and one on the low. (pr/29 D1, and CLAUDE.md M7)

**K4. The adjacent-time-slice fallback must step by ticks-per-slice, not by 1.**
`time_blob_map` is keyed on `Blob::slice_index_min()`, which is in **ticks**
(`Facade_Blob.h:33`), so consecutive slice starts are `nticks_per_slice` apart
(SBND: 4). A literal `±1` names no real slice and the whole fallback is dead
code. WCT reads the stride from
`Grouping::get_nticks_per_slice().at(apa).at(face)` — never a hard-coded 4,
which would work on SBND and break elsewhere. (pr/29 D12)

**K5. The flag=2 same-blob pass groups by the RETILED blob in WCT and by the
ORIGINAL mcell in WCP — and WCT is correct.** WCP resolves each selected point
back to an original-cluster mcell by strict wire containment
(`PR3DCluster_steiner.h:542-563`) and skips points that match none
(`PR3DCluster_graph.h:90-91`). That is **not** a physics choice: `point_cloud_steiner`
is a member of the *original* cluster while the retiled cluster and its mcell
holder are `delete`d immediately after the tree is built (`:53-56`), so a
retiled `SlimMergeGeomCell*` stored on it would dangle, and a WCP `WCPoint` has
exactly one `mcell` field. WCP itself groups by the **retiled** blob one step
earlier, at flag=1. WCT's `steiner_pc` is a `Dataset` with no blob pointer, so
the constraint does not exist; the "skip points with no mcell" branch is that
lookup's failure path, not a rule. **Do not reintroduce it.** (pr/29 D3,
§10.2.9)

**K6. `get_two_boundary_wcps(2)`'s `mcell == 0` and
`Estimate_total_charge() < 1500` cuts are replaced, not dropped.** They depend
on the per-point original-mcell association K5 removes. WCT's
`get_two_boundary_steiner_graph_idx` scores boundaries on the **regular** point
cloud, where blob charge and dead-wire counts exist, then snaps each to the
nearest Steiner **terminal** (`Facade_Cluster.cxx:3423-3435`).

**K7. The path skeleton is resampled finer in WCT.** WCP uses
`num_steps = floor(dis/step)` and so *undershoots* its own 0.6 cm target by up
to a factor two; WCT uses `floor(dis/step) + 1` and never exceeds it
(`DynamicPointCloud.cxx:744`). A denser sample of the same curve is a more
accurate nearest-distance. WCP's duplicated endpoint is inert. **Keep WCT's.**
Note `make_points_cluster_skeleton` has a single body shared with
`clustering_deghost`, so any change there is not local. (pr/29 D7, §10.2.8)

**K8. `disable_dead_mix_cell` must reach the edge-weight charges.** WCP computes
the `Qs`/`Qt` in the edge weight with the same value `Create_steiner_tree` was
called with (`PR3DCluster_steiner.h:514`, `:521`; the chain passes `false`).
WCT's `create_enhanced_steiner_graph` defaults the parameter to `true`, so
dropping it at the call silently selects the other branch of `calc_charge_wcp`:
`true` sums all three planes then subtracts **dead** ones
(`charge_uncertainty > 1e10`), `false` sums only planes with a **nonzero**
charge value — independent predicates. Behind
`edge_charge_forward_dead_mix`. (pr/29 D2, §12)

**K9. Containment additionally requires the same APA and face in WCT.** Forced:
`time_blob_map` is `apa → face → tick → blobs`, and WCP is single-APA so it has
no counterpart. Correct as written. (pr/29 D8)

**K10. `recover_steiner_graph` is not ported.** It is WCP's only MST
(`PR3DCluster_steiner.h:77-180`) and no `wire-cell-pid` app calls it — verified
by an exhaustive `prototype_base/pid/apps/` grep. A gap on paper only. (pr/29
D11)

#### Config knobs that select WCP-faithful behaviour

All three default **OFF** in C++ (so no detector moves without opting in) and
are **ON in the SBND operating point** (`cfg/pgrapher/experiment/sbnd/wct-pr-perevt.jsonnet`):

| key on `CreateSteinerGraph` | restores | dictionary entry |
|---|---|---|
| `terminal_wire_tol` (int, 0) | WCP's ±1 wire slack in the terminal filter only | K3 |
| `terminal_adjacent_slice` (bool, false) | makes WCP's t±1 slice fallback actually resolve | K4 |
| `edge_charge_forward_dead_mix` (bool, false) | honours the caller's `disable_dead_mix_cell` in edge weights | K8 |

## WCP Algorithms

- [ ] :question: What is the overall "data flow graph" for these stages. 
    * Xin: For major components of the Pattern Recognition, the short answer is that it would require everything (3D points, saved in PCT, 2D measurements saved in CTPC, derived data 3D fitted trajectory and dQ/dx). These are needed to do trajectory/dQ/dx fitting, do PID, do energy reconstruction etc. The light information is not needed, but the matched time from flash is needed to position things at correct location. 
    * Xin: for the later feature extraction, it also need to access the above data to form various calculation. 
    * Xin: for the trajectory/dQ/dx fitting, the input consists of 2 parts: i) hypotheses, (derived from 3D points and expressed as ProtoSegment, ProtoVertex etc), ii) measurements (2D data).
- [ ] :question: What input data is required for each stage, what output is produced?  (can input/output be fully modeled as a graph or is there "extra" data that will not fit that model?)
    * Xin: There are two major outputs: 1. the particle flow (with a tree of particles, made by ProtoSegment, ProtoVertex, and WCShower in WCP), the different stages of the PAttern recognition is essentially figure out some information (e.g. PID, energy), as well as the order of things to form the particle flow tree. 
    * Xin, 2. in addition to the particle flow, we also extract various features for event selection. These features are based on various input (3D points, 2D measurements ...). These features in WCP were used to train BDTs for event selection.  
- [ ] :question: Is the data flow graph a linear pipeline or a more general DAG?
    * Xin: the WCP flow graph can be viewed as a linear pipeline. 
- [ ] :question: Will the PC-tree be required throughout all stages?
    * Xin: likely, since 3D points, which is needed, are stored in PC-tree  
- [ ] :question: Will the PC-tree be required to be output by the final algorithm that produces the particle flow data structure?
    * Xin:  I think the 3D points should be stored or associated with the particle flow data structure. For example, if I have a pion in the output, I do want to know which 3D points are associated with this pion. These would be useful in making features for latter usage. 
    * Xin: for this section, the performance is going to be the key, the computing should be less of a concern. This means that it would be great that the system can be expanded to incorporate more advanced algorithm (e.g. AI/ML).


### [multi dQ/dx fitting](https://github.com/BNLIF/wire-cell-pid/blob/537a3fd17f8a7b3cf5412594267c14c4cc1775cb/docs/PR3DCluster_multi_dQ_dx_fit.md) (WCP)

## Toolkit invariant: "perblob" node-local arrays are parallel to blob children order

Not a WCP concept (WCP keeps per-cell state on the cells themselves).  In the
toolkit, per-blob provenance (`isolated`, `assoc_cluster_id/main`,
`real_cluster_id/main`) lives in a cluster-level N-row `"perblob"` Dataset
where row *i* MUST describe blob child *i*.  The tree primitives
(`separate()`, `merge()`, `take_children()`) do NOT maintain this — any
separate-then-merge-back round trip silently permutes children while the
arrays keep their old row order.  Before writing or porting ANY code that
mutates a cluster's child set, read `clus/docs/perblob_invariant.md` (rules,
safe patterns, audit table; doc 52 §12-§13 in wcp-porting-validation has the
full bug history).

## `ProtoSegment::is_dir_weak()` (WCP) vs `PatternAlgorithms::seg_dir_weak()` (WCT)

WCP's raw `dir_weak` member is PRIVATE with no getter; the only public read
accessor is `is_dir_weak()`, which ORs the raw flag with particle-score
thresholds (mu: score>0.07 at >=5cm / >0.15 below; p: >0.13 / >0.27) — and the
score defaults to (and is reset to) the sentinel 100, so a typed-but-unscored
or PID-invalidated mu/p segment is ALWAYS weak in WCP.  The original port read
a raw `dir_weak()` getter at all ~84 sites (silent substitution, no doc trail;
even signed off as "equivalent" in a tagger review — corrected).  Since
2026-07-30 every read goes through `PatternAlgorithms::seg_dir_weak()`, gated
by the `TaggerCheckNeutrino` key `dir_weak_use_score` (C++ default false =
legacy raw reads, byte-identical; SBND module default true).  The faithful
predicate is `segment_is_dir_weak()` (`PRSegmentFunctions.cxx`).  Raw reads
remain ONLY in: that predicate's fall-through, `break_segment()` flag
propagation, and trace prints.  Full history: wcp-porting-validation
`sbnd_xin/docs/pr/6_is-dir-weak-port-divergence.md`.

## `pu/pv/pw/pt` published to Magnify-tracking: WCP writes DISPLAY coords, WCT writes INDEX coords

INTENTIONAL divergence — do not "restore" the missing `+ 0.5`.

WCP's `PR3DCluster::do_tracking()` publishes the projection of each fitted point
with a half added on every axis
(`prototype_base/pid/src/PR3DCluster_trajectory_fit.h:468-471`,
`pu.push_back(offset_u + 0.5 + ...)`, likewise `pv + 2400`, `pw + 4800`, `pt`).
That half is not geometry: it compensates the Magnify-tracking GUI, whose
projection histograms are `TH2F(n, 0, n)` filled with `SetBinContent(index+1,
...)`, so index *i* is painted centred on *i + 0.5*.  WCP therefore stores a
DISPLAY coordinate, not a wire index — `pu` there is half a wire off the wire it
names.

The toolkit writers (`UbooneMagnifyTrackingVisitor`, `SbndMagnifyTrackingVisitor`,
`SbndPrMagnifyTrackingVisitor`) publish `fit.pu` unshifted, i.e. an INDEX
coordinate in the same frame as `T_proj_data.channel`, consistent with
`Facade::point2wind` (`wind = round((y-center)/pitch - 0.5)`) and with
TrackFitting's own 2-D fit target.  Integer = wire centre, and `T_rec_charge.pu`
is directly comparable to `T_proj_data.channel`; same for `pt` and
`time_slice`.

The half bin lives in the VIEWER instead: `Magnify-tracking-SBND` bins its
projection pads with edges `lo-0.5 … hi-0.5` on both axes, so a bin is centred on
its index.  A port that "fixes" the writer to match WCP re-breaks that viewer by
exactly half a channel and half a slice.  Measurement and full history:
wcp-porting-validation `sbnd_xin/docs/pr/7_magnify-tracking-projection-alignment.md`.

## TrackFitting T frame: WCP fits in TIME-SLICE units, WCT fits in TICK units

INTENTIONAL divergence — do not "restore" the commented-out
`nticks_live_slice` factor in `clus/src/TrackFitting.cxx` (`time_slice_width =
md["tick_drift"]`).

WCP's dQ/dx fit works in time-slice units end to end: the data T coordinate is
`GetTimeSlice()` (one unit = nrebin = 4 ticks), the model uses `slope_xt =
1/time_slice_width` with `time_slice_width = mp.get_ts_width()` (= 4 ticks of
drift, 2.202 mm at uBooNE), and `sigma_L = hypot(diff_sigma_L,
add_sigma_L)/time_slice_width`
(`prototype_base/pid/src/PR3DCluster_dQ_dx_fit.h:303-305,435,612`).

WCT re-expresses the same fit in tick units end to end: the data T coordinate
is a tick index (`Grouping::convert_3Dpoint_time_ch`: `tind = round(time /
tick)`), the model uses `slope_t = 1/(xsign*drift_speed*tick)`
(`TrackFitting.cxx:517`), and sigma divides by the per-tick drift
(`tick_drift`).  Both coordinate and sigma scale by the same factor of 4, so
residual/sigma — the only thing the chi^2 sees — is identical.  The physical
add_sigma_L also matches: WCP `1.428249*ts_width/nrebin/0.5` = WCT preset
`1.428249*0.5505mm/0.5` = 1.5725 mm.

Verified at runtime (2026-07-30, instrumented build): model `central_T`
matches the data tick index to within rounding on both detectors (uBooNE
7052.28 vs 7052 with mm/T-unit = 0.5505 = tick_drift, not 2.202; SBND 1793.09
vs 1793 with 0.7815, not 3.126), and physical sigma_L is sane (2.22 mm uBooNE
/ 2.60 mm SBND).  Multiplying `time_slice_width` by `nticks_live_slice`
"for prototype parity" would inflate sigma 4x against tick-indexed data.
Full record: wcp-porting-validation
`sbnd_xin/docs/pr/2_uboone-chain-gap-analysis-and-validation-plan.md` sec 8.3a.

## Isochronous first-segment endpoints: `iso_endpoint` bypasses the boundary metric ON PURPOSE, knob-gated

WCP's `get_two_boundary_wcps` non-cosmic metric
(`prototype_base/pid/src/PR3DCluster_path.h:530-536`, ported verbatim to
`clus/src/Facade_Cluster.cxx calculate_boundary_metric`) deliberately zeroes
the wire-index-separation terms (`* 0`), scoring a candidate endpoint pair by
`|dx|/one_tick + live_wire_count(U+V+W)` only.  On an isochronous cluster
(charge in one narrow drift slab, imaged as a filled 2-D sheet) `|dx| ~ 0`, so
the winner is the widest wire-footprint pair — two corners on the sheet's
edges, not the tip of its axis — and the cheapest Steiner-graph paths between
such corners are the sheet's boundary edges (the two-edge-track fan, SBND evt
18255-271851, doc wcp-porting-validation `sbnd_xin/docs/pr/24` round 2).

The toolkit-only fix (`PatternAlgorithms::find_iso_first_segment_endpoints`,
`clus/src/NeutrinoPatternBase.cxx`, knob `iso_endpoint`, C++ default false)
gates on a quantile-trimmed corrected-frame drift-x extent and, when it
fires, takes the endpoints from the sheet's principal axis and SKIPS the
local-PCA endpoint refinement (itself a toolkit addition with no prototype
counterpart, commit 1eb097a9) on that branch.  Knob off = the prototype
metric byte-identically, `* 0.0` coefficients preserved.  Do NOT "fix" the
`* 0.0` coefficients in `calculate_boundary_metric` itself: the zeroing is
intentional prototype behavior and every non-isochronous cluster depends on
it.  Prototype precedents (neither reachable from NeutrinoID):
`adjust_wcpoints_parallel` (`prototype_base/data/src/PR3DCluster.cxx:428`,
cluster-separation only) and `search_for_connection_isochronous`
(`prototype_base/pid/src/PR3DCluster_graph.h:1445`, call site commented out).

Round 3 (doc pr/24 §15) additionally requires 2-D sheet-ness before the branch
may act (`iso_endpoint_min_aspect`, trimmed transverse/axial extent ratio,
default 0.12), so a 1-D track-like cluster is handed back to the prototype
metric rather than merely tolerated by it; and the endpoint is the UNTRIMMED
axial extreme (laterally centred only within a 3 cm end band), because a
trimmed extreme left a tip stub that `find_other_segments` then claimed as a
segment — a 0.9 deg "vertex" inside a straight track (SBND mcp1k evt 284794).
`iso_endpoint_tube_radius` is diagnostic only; a hard tube filter around the
straight axis line was measured to pull endpoints up to 28.6 cm INWARD on long
or curved clusters and was rejected.

## `examine_vertices_3` / `get_local_extension`: the "extension" recovery step has no outward check, in BOTH trees — `v3_extension_guard` is a prototype-limitation fix, not a port correction

`examine_vertices_3` (`clus/src/NeutrinoStructureExaminer.cxx`, a faithful
port of `NeutrinoID_proto_vertex.h:2412-2463`) revisits the main cluster's two
original `init_first_segment` endpoints and tries to push each one further
out via `get_local_extension` — a 10 cm-radius Hough-transform direction
estimate (`clus/src/NeutrinoStructureExaminer.cxx get_local_extension`, port
of `PR3DCluster_path.h:288-316`). Neither tree's caller ever checks that the
returned point actually moves the vertex FARTHER from the segment's other
(far) endpoint; the only checks are "different from both existing points"
(`wcp1.index == vtx / wcp2.index`, `PR3DCluster_path.h:2443`) and "the
rebuilt path isn't more than 2x longer" (`:2452`). At the axial extreme of an
isochronous (drift-perpendicular) sheet — exactly where `iso_endpoint` (round
2 above) picks its seed — the local 10 cm neighbourhood is dominated by the
sheet's transverse spread, so the Hough direction estimate is poorly
conditioned and can point back INTO the cluster. Measured on all three doc
pr/24 §18 (round 5) events (SBND 18259-42280, 18255-271851, 18255-350186),
the "extension" landed 7.5-8.9 cm closer to the far endpoint than the
original vertex, silently amputating the delivered trajectory by that much —
this is what produced round 4's 8.4-10.9 cm undershoot, NOT the endpoint pick
itself (round 4 had misattributed it to unidentified shared refinement code;
see doc pr/24 §17.5, retracted in §18).

Because the prototype has the identical gap (no distance-to-far-endpoint
check anywhere in `examine_vertices_3`), this is a **prototype limitation
(M15)**, not a port error — do NOT silently "fix" `get_local_extension` or
`examine_vertices_3` unconditionally. The fix, `v3_extension_guard`
(`clus/inc/WireCellClus/NeutrinoPatternBase.h`, `TaggerCheckNeutrino.h`), C++
default false, rejects a candidate unless it increases the vertex's distance
to the segment's far endpoint by more than `v3_extension_min_gain` (default
-1.0 cm, tolerating the legacy arm's own few-mm retreat from the same bug
while rejecting the multi-cm amputation). Knob off = both trees' unconditional
accept, byte-identical.

## `skip_trajectory_point`'s charge-consistency revert has no protection on an isochronous cluster, in BOTH trees — `skip_revert_iso_xext_cut` is a prototype-limitation fix, not a port correction

The multi-track trajectory fit's charge-consistency veto (toolkit
`clus/src/TrackFitting.cxx skip_trajectory_point`, a faithful port of
`prototype_base/pid/src/PR3DCluster_trajectory_fit.h:706-750`, called from
`PR3DCluster_multi_track_fitting.h:429`'s `examine_trajectory`, doc pr/28 §13
T1/T2) compares charge near a fitted point (`c1`) against charge near its
pre-fit position (`c2`) on each plane; if the ratio falls below threshold it
reverts the point (`p = ps_point`) with **no bound on how far that reverts
it** — the fit's own regularizer term toward the initial trajectory is
commented out in both trees (`PR3DCluster_trajectory_fit.h:302-304`), so
nothing else limits the jump either.

Measured (doc pr/28 §17, SBND 18255-271851 vs 18255-388): the individual
revert is *not* large (median 0.4 cm) — the failure is not one bad jump but
many small reverts applied consistently along a stretch whose pre-fit
positions already sit off the eventual smoothed line, which happens when the
cluster is **isochronous** (its charge sits in one narrow drift slab, same
geometry `iso_endpoint` above targets): `c1` and `c2` then integrate the same
overlapping 2-D charge blob rather than resolving two distinct samples of a
track, so the veto's premise — "the fitted point sits on less charge than a
resolvable alternative" — does not hold, and it blocks the fit from smoothing
those points onto the trunk. Two per-point candidate discriminators were
measured and both come back flat between the pathological cluster and a
normal one (§17.2): whether the fitted/reference projections share a time
tick, and the local displacement's angle to the drift axis. The measure that
DOES discriminate is a **cluster-level** one — blob-center drift-x extent,
the same measure `iso_band_like` above uses for its own, unrelated veto.

Because the prototype has the identical unbounded revert and the identical
missing regularizer (this is not a toolkit-introduced bug — `23bd6783`
correctly revived behavior the prototype has always had), this is a
**prototype limitation (M15)**, not a port error — do NOT "fix" the revert to
be bounded or conditional for every cluster; every non-isochronous cluster's
fit depends on the veto being unconditional (confirmed: SBND 18255-388, doc
pr/28 §13's own accepted case, is verified bit-identical on its own cluster
with the knob on — doc pr/28 §17.4 G3). The fix, `skip_revert_iso_xext_cut`
(`clus/inc/WireCellClus/TrackFitting.h`, reached through
`TaggerCheckNeutrino`/`TaggerCheckSTM`'s `trackfitting_config_file` JSON, NOT
an `IConfigurable`), C++ default -1 (off), abstains from the revert (keeps
the fitted point) only for a point whose segment's cluster has blob-center
drift-x extent below the cut. Knob off = both trees' unconditional revert,
byte-identical. A related dilution weakness in the same veto — a
non-informative plane contributes a free *agreement* vote rather than being
excluded — was found and is recorded, not fixed (doc pr/28 §17.5): the
naive fix (renormalize by informative-plane count) has the WRONG SIGN, since
it would remove the very term currently masking disagreement on the other
two planes.

## `NeutrinoTrackShowerSep.cxx`'s wholesale track-to-electron conversion sites and `segment_is_shower_topology` never consult the segment's own dQ/dx — `shower_reclass_dqdx_guard`/`shower_topo_dqdx_guard` are designed divergences, not port corrections

Doc pr/40 (owner report: 9 SBND events where a proton/pion/muon **track**
displays as an electron). Two mechanisms, both **prototype-faithful in both
trees**, not port bugs:

1. Three sites in `NeutrinoTrackShowerSep.cxx` —
   `examine_all_showers`'s `flag_change_showers` loop,
   `improve_maps_shower_in_track_out`'s two reclassify loops (out_tracks and
   no-direction segments), and `improve_maps_no_dir_tracks` Case E — convert
   a direction-weak, untyped, or vertex-topology-flagged segment to electron
   **unconditionally**, without re-checking the segment's own charge. Doc
   pr/9 §4 already established this is deliberate prototype design: *"in a
   shower-dominated cluster, only a strong-direction track survives as a
   track."* The prototype has the identical unconditional conversion at its
   own counterpart sites.
2. `segment_is_shower_topology` builds a per-point `vec_dQ_dx` array
   (normalized by `MIP_dQ_dx`) that is otherwise **dead** — read only as
   `.size()` (pr/31 GOTCHA 5, itself citing the prototype's identical
   dead-array shape). The whole shower/track call is decided by a 5-branch
   geometric spread test alone; the segment's charge never gets a vote,
   in either tree.

Measured (doc pr/40 Part 0 / G2): a segment with a genuinely hadron-like
median dQ/dx (proton ≥ 1.75× MIP, muon ≤ 1.2× MIP — the SAME thresholds
`segment_determine_dir_track`'s own short-track fallback already trusts) can
be swept into electron by either mechanism regardless of how confidently
hadron-like its own charge profile is. 8 of the 9 owner-reported cases
measured above 1.75× or matched the muon band; the population census found
129 segments across 42 events in the broader class.

Because both mechanisms are faithful to the prototype's own behaviour, the
fix is a **designed divergence (owner-approved this round), not a port
correction**: `shower_reclass_dqdx_guard` (all three sites in mechanism 1)
and `shower_topo_dqdx_guard` (mechanism 2) share one helper,
`segment_dqdx_spares_electron_reclass` (`PRSegmentFunctions.h/.cxx`) — spares
a segment from conversion / from having `kShowerTopology` set only when its
own median dQ/dx is decisively proton- or muon-like. C++ default `false` for
both = both trees' unconditional conversion, byte-identical. SBND
production default `true` since 2026-08-06 (doc pr/40 G1-G3 all pass:
48/48 byte-identical knob-off, 8/9 owner cases fixed, zero verdict
regression on the 48-event population).

**Do NOT make either guard unconditional or widen its threshold band.** The
gap between 1.2× and 1.75× MIP is deliberately left untouched — doc pr/40's
evt 256587 (median 1.26× MIP, still electron with the guard on) sits in that
band, and its own *intra-segment* charge is genuinely ambiguous, not a guard
failure. Widening the band to catch it would sweep a much larger,
less-certain population into an automatic override with no principled
stopping point.

**Superseded in part, doc pr/40 round 2**: evt 256587's own median dQ/dx really
is ambiguous, but its *topology* is not — segment 11079 starts exactly at
the neutrino vertex and its far end abuts a PID'd, charge-confirmed proton
(segment 11080, 3.72× MIP). An electron cannot father a proton. This is a
different evidence axis (topology, not the segment's own charge profile) and
does not widen the band above — see the new entry below for the
`shower_proton_daughter_pion` knob that resolves it.

## `set_default_shower_particle_info`'s electron default never consults graph topology — `shower_proton_daughter_pion` is a designed divergence, not a port correction

Doc pr/40 round 2 (owner: reviewing the pr/40 fix's own Bee display, "for
256587... in the end of the particle, there is a proton, which is high
dQ/dx... an electron cannot change to proton. So the fact that we identified
a proton should change it to pion instead of electron"). A third mechanism,
distinct from the two above: `set_default_shower_particle_info`
(`NeutrinoPatternBase.cxx`) is the single stage-4 choke point where any
`flag_shower` segment still missing `particle_info` gets defaulted to
electron (mirrors prototype `ProtoSegment::get_particle_type()`, which
**always** returns 11 for any shower segment, unconditionally — this default
itself is prototype-faithful and stays). Neither this function nor the
prototype it mirrors ever looks at the *graph* around the segment: not the
neutrino vertex, not the segment's neighbours' own PID. The prototype has
**no** proton-daughter veto anywhere — a designed divergence, not a port
correction.

Measured (48-event nueCC48 population, `work-pr40-on48` arm): of 2209
electron-labelled segments, a naive "has a >1.75× MIP neighbour anywhere"
rule fires on 348 — far too broad, since most electrons in a shower
legitimately sit next to a high-dQ/dx track. Requiring the segment (a) to
emanate from the neutrino vertex by graph identity (not a distance cut —
every measured case sits at exactly d=0.00 cm) **and** (b) its far end to be
a vertex whose out-edges include a segment already PID'd proton (2212) that
is *independently* charge-confirmed (its own median dQ/dx > 1.75× MIP)
narrows this to 5/2209. The neutrino-vertex requirement is what does the
work: it excludes the ordinary, correct nueCC topology where an electron and
a proton merely *share* the neutrino vertex as siblings, not parent/daughter.

`shower_proton_daughter_pion` (config key, threaded via
`m_shower_proton_daughter_pion`) relabels the candidate segment PION (211)
instead of electron when `segment_has_proton_daughter`
(`PRSegmentFunctions.h/.cxx`) fires. C++ default `false` = legacy
unconditional electron default, byte-identical.

**Round 2 found a fourth writer that reverted this end-to-end; round 3
closed it. SBND PRODUCTION DEFAULT ON.** The override fires correctly at
this choke point — traced, `pdg 11 -> 211` at `NeutrinoPatternBase.cxx` —
but `Shower::update_particle_type` (`PRShower.cxx`, called from 8 sites in
`NeutrinoShowerClustering.cxx`) unconditionally reasserted electron on a
shower's start segment whenever `shower_length > track_length`, with zero
awareness of PID or topology, reverting the override in the same pass (evt
256587 seg 11079: traced, `pdg 211 -> 11` at `PRShower.cxx:801`; population
census showed the override surviving end-to-end in only 1/2209 cases).
Round 3 (doc pr/40 round 3) threads `main_vertex` + `protect_proton_
daughter_pion` (both legacy-default) into `update_particle_type`, which
re-derives `segment_has_proton_daughter` on the SAME MIP scale
(`m_mip_dqdx_median`, not the function's own `mip_dqdx`, a different scale —
the two must not be conflated) and skips the reassignment when it fires.
Gate-clean: 48/48 byte-identical off; evt 256587 now reads pdg 211
end-to-end; population census shows exactly 2/2209 segments move, no
verdict regression. Flipped SBND ON, verified with a bare single-event
cfg-only run hash-matching the gated on-arm exactly — see doc pr/40 round 3.

**Do not widen this to "any high-dQ/dx neighbour."** The 348-vs-5 gap above
is the reason both the neutrino-vertex-emanation requirement and the
independent charge-confirmation requirement exist; dropping either
reintroduces the false-positive population this knob was designed to avoid.

## Relabelling a shower segment's PDG does not make it stop being a Shower — `shower_proton_daughter_pion_dissolve` (F7) is a designed divergence

Doc pr/40 round 4 (owner: reviewing the round-2/3 `shower_proton_daughter_pion`
fix's own Bee display, evt 256587: *"the particle flow do show the pion+, but
we do not see the proton after it in the particle flow... the EM shower were
modified as pion, but not on the individual tracks"*). F5 changes a segment's
*pdg* (11 → 211) but never touches `SegmentFlags::kShowerTrajectory` /
`kShowerTopology`. Those flags are what `shower_clustering_with_nv_in_main_
cluster` (`NeutrinoShowerClustering.cxx`) actually tests (`is_shower_seg =
flags_any(kShowerTrajectory) || flags_any(kShowerTopology) ||
|pdg|==11`) — with the flags still set, a pion-relabelled segment is still
rooted as a `Shower`. Two measured consequences on evt 256587 seg 11079:
`fill_bee_pf_tree` (`MultiAlgBlobClustering.cxx`) pre-claims every
shower-owned segment (`used_segs = shower_segs`), so the segment's own
charge-confirmed proton daughter (seg 11080, the very evidence F5 used to
relabel it) never gets its own particle-flow node; and the pi+ Bee node's
displayed extent is the *shower's* endpoint — a 0.35 cm fragment absorbed
from a different, non-main cluster (seg 81153) — not segment 11079's own end.

The prototype has **no** mechanism that reclassifies an already-formed shower
back into a track after a PDG override — designed divergence, not a port
correction. `shower_proton_daughter_pion_dissolve` (config key, threaded via
`m_shower_proton_daughter_pion_dissolve`) clears both shower flags in
`set_default_shower_particle_info` at the same site F5's override fires,
provided `shower_proton_daughter_pion` is also on. C++ default `false` =
legacy = byte-identical. Measured (48-event nueCC48 population,
`work-pr40r3-on48` arm): exactly 2/2209 electron-labelled segments carry
`pdg 211 && flag_shower` before this fix (evt 256587 seg 11079, evt 342199
seg 72098 — the latter on a non-main cluster, so it produces no `mc.json`
delta but is still a real behaviour change).

## A muon segment cannot terminate in a multi-proton hadronic vertex — `muon_multi_proton_pion` (F8) is a designed divergence

Doc pr/40 round 4 (owner, same Bee review round, evt 489330: *"there is one
muon → two protons. This is not physical, in this case, the muon should be
changed to pion"*). Measured: segment 4019 (mu-, 65.2 cm) has TWO
charge-confirmed proton daughters (segs 4018, 4044) at its far (non-neutrino-
vertex) end. The prototype has no proton-multiplicity veto on a track's PID
— designed divergence. `muon_multi_proton_pion` (config key, threaded via
`m_muon_multi_proton_pion`) relabels such a muon segment PION (211) via the
new pass `PatternAlgorithms::override_muon_multi_proton_pion`, called
immediately after `set_default_shower_particle_info` in `examine_direction`
(same per-cluster `main_vertex`, same last-word-before-shower-clustering
position). The topology test, `segment_at_multi_proton_vertex`
(`PRSegmentFunctions.h/.cxx`), generalizes F5's `segment_has_proton_daughter`:
same graph-identity main-vertex exclusion (a muon-plus-two-protons vertex AT
the neutrino vertex is the ordinary, correct numuCC topology and must not
fire), same 1.75× MIP independent charge-confirmation threshold, but
`min_protons=2` instead of F5's implicit 1, and it checks EITHER endpoint
other than `main_vertex` rather than only the segment's "far" end (a muon
segment has no owner-assigned direction the way an electron-daughter-of-the-
neutrino-vertex segment does). C++ default `false` = legacy = byte-identical.

**Owner decision: no propagation across a kink.** Evt 489330's muon segment
4019 sits behind a degree-2 vertex from a second muon segment (4043, 28.4 cm,
running to the neutrino vertex); only 4019 is relabelled, 4043 stays `mu-`.
Population (48-event nueCC48): exactly 1/N muon segments fires this rule
(evt 489330 seg 4019); 6 more muon segments have exactly one qualifying
proton at a non-main vertex and are deliberately left untouched — the
owner's "two protons" wording is read literally, `min_protons=2` is not a
placeholder.

## Segment-level straightness has no exemption from three shower-seeding tests — `track_pid_persist_dqdx_electron_guard`/`shower_connect_main_vertex_straight_guard`/`shower_traj_straight_guard` (F9/F10/F11) are designed divergences, **segment-level fix only, does not reach the displayed outcome, NOT flipped**

Doc pr/40 round 5 (owner, three new Bee cases: evts 18364-84229/18255-54341,
"electron → muon?", both read as stopping muon + Michel; evt 18255-55715,
"not electron → muon", read as an exiting muon behind a wrongly-labelled
pion). Three independent writer sites, none with a straightness exemption:
`segment_determine_dir_track`'s F1 persist-on-dQ/dx rescue fires
unconditionally once `pdg_code != 0` (F9); `shower_clustering_connecting_
to_main_vertex`'s three-branch skip has no branch for a long straight track
with no confident PID yet (F10); `segment_is_shower_trajectory` never got
the straightness exemption pr/40 F3 gave its topology sibling (F11). All
three key off `segment_is_straight_long_track` (new shared helper,
`PRSegmentFunctions.{h,cxx}`, `direct_length >= 34cm || direct_length >
0.93*length`, same threshold shape as the existing `NeutrinoVertexFinder.cxx:
1432-1447` demotion). None of the three sites has a prototype analog
(`prototype_base/pid/` checked directly) — designed divergences, not port
corrections. C++ default `false` on all three = legacy = byte-identical
(G1: 48/48 events, 96/96 archives, 0 mismatches).

**Measured: fixing the segment's own pdg does not change the Bee/mc.json
outcome, because the display is decided at the shower seeding/absorption
boundary, not at the segment.** `Shower::complete_structure_with_start_
segment` (`PRShower.cxx:337-408`) flood-fills a shower-seeded segment's
downstream sub-tree with no per-segment test, and a second seeding path,
`shower_clustering_with_nv_in_main_cluster`'s `is_shower_seg`
(`NeutrinoShowerClustering.cxx:116-119`), is untouched by F10 (which only
gates its sibling `shower_clustering_connecting_to_main_vertex`). Result on
all three owner cases with all three knobs on: evt 84229 seg 19038 becomes
correctly `pdg=13` but the `mc.json` node is unchanged (`19039 'e- 89 MeV'`,
still merged); evt 54341 seg 18005 does get excluded from the shower (split
shape achieved) but, unshielded from the shower path for the first time,
ordinary track PID calls it `proton` (2212) from its own elevated dQ/dx, not
`mu-` as intended — an open physics question, not a bug in these three
guards; evt 55715 seg 15007 becomes correctly `pdg=13`, but the shower
re-seeds one segment further up at seg 15005 (isolated to F11 alone via a
clean single-knob A/B: `pdg 211` all-off → `pdg 11` with F11 alone) — a
confirmed **regression against the owner's explicit round-5 planning answer**
that seg 15005 must stay untouched. **Not flipped.** All three knobs are
landed as gate-clean infrastructure (G1/G5/G6 pass) but G2 fails on every
owner-reported case; the fix implied by the root-cause finding is a change
to the seeding/absorption boundary itself (comparable in shape to round 4's
F7 shower-dissolve), scoped larger than this round and not attempted here.

## The shower flood-fill gains a per-segment track exclusion; the Michel rescue reaches confident-direction multi-prong stems — `shower_absorb_track_guard` (F12) + `michel_stem_muon_rescue` (F14) are the round-6 boundary-level fixes; `shower_connect_protected_pion_guard` (F13) is a measured-dead negative result, never flipped

doc pr/40 round 6 closes the round-5 negative result above. Three knobs, two
live, one dead:

- **F12 `shower_absorb_track_guard`** — `Shower::complete_structure_with_
  start_segment` (`PRShower.cxx`) gains a per-segment exclusion: a
  confidently PID'd non-electron (`pdg != 0 && |pdg| != 11`) that is long and
  straight (`segment_is_straight_long_track`) is not absorbed, the walk
  terminates there, the excluded segment is NOT claimed in `used_segments`,
  and long-muon pseudo-showers (`Shower::get_particle_type()==13`, set by the
  in_main_cluster seeder before completion) are exempt so broken-muon
  reassembly keeps working. One knob threads all 7 call sites in
  `NeutrinoShowerClustering.cxx`. The excluded segment automatically gets its
  own PF node: `fill_bee_pf_tree`'s suppression key is the shower's VIEW
  (`fill_sets()`), and the `pf_shower_vertex_barrier` orphan safety net (doc
  pr/38) covers the BFS-unreachable case. The prototype's counterpart
  flood-fill has no per-segment test either — designed divergence.
- **F14 `michel_stem_muon_rescue`** — new pass `override_michel_stem_muon`
  (`NeutrinoPatternBase.cxx`, called at F8's `examine_direction` call site).
  The toolkit's own Michel rescue ("a stopped proton cannot produce a Michel
  electron", `NeutrinoVertexFinder.cxx` weak-direction branch) is limited to
  `seg_dir_weak` segments with a degree-2 stopping vertex; the widened pass
  reaches a `pdg==2212`, straight-long, main-vertex-emanating stem with >=1
  shower-like sibling (`kShowerTrajectory || |pdg|==11`) at its stopping
  vertex, relabelling it mu-. The degree-2 restriction and the weak-dir gate
  are both prototype-faithful in the original — designed divergence.
- **F13 `shower_connect_protected_pion_guard`** — DEAD as shaped, kept as a
  documented negative result (doc pr/36 F2 precedent), excluded from the
  flip. The full-transition trace (`WCT_PID_WRITE_DEBUG=2`, a round-6
  widening of the round-1 probe) showed the motivating segment (55715/15005)
  is already `pdg 2212` at candidate-selection time; its baseline `211` was
  derivative of the very bug F11 fixes (wrong e- on 15007 → Michel rescue
  2212→13 → single-muon pion demotion 13→211). A `pdg==211` predicate
  cannot fire there, and F12 alone produces the intended display.

Measured on the three round-5/6 owner cases (all in doc pr/40 round 6's
Demonstration table): 84229 → `19038 'mu-'` own node + `19040 'e-'` child
(owner-accepted residual: parent stub 19039 reads `pi+` from the legacy
single-muon selection, its own honest call being a second muon); 54341 →
`18005 'mu-'` + children `e-`/`mu-` (exact owner request); 55715 → `15007
'mu-'` own node, `15006 'proton'` kept (owner-accepted: 15005 reads `proton`,
its own charge-based call, the baseline `pi+` having been derivative of the
fixed bug).
