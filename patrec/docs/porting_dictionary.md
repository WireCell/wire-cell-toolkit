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


### [ProtoSegment](https://github.com/BNLIF/wire-cell-pid/blob/537a3fd17f8a7b3cf5412594267c14c4cc1775cb/docs/protosegment.md) (WCP) vs. **xxx** :warning: (WCT)

### [ProtoVertex](https://github.com/BNLIF/wire-cell-pid/blob/537a3fd17f8a7b3cf5412594267c14c4cc1775cb/docs/protovertex.md) (WCP) vs. **xxx** :warning: (WCT)

### [WCShower](https://github.com/BNLIF/wire-cell-pid/blob/537a3fd17f8a7b3cf5412594267c14c4cc1775cb/docs/wcshower.md) (WCP) vs. **xxx** :warning: (WCT)

### [Steiner Tree](https://github.com/BNLIF/wire-cell-pid/blob/537a3fd17f8a7b3cf5412594267c14c4cc1775cb/docs/PR3DCluster_steiner.md) (WCP) vvs. **xxx** :warning: (WCT)

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
