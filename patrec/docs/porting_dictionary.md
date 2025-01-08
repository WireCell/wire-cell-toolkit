#  Between Clustering and Trajectoryand dQ/dx fitting (Current)

## [Opflash](./QLBundles/Opflash.md) (WCP) vs. <span style="color:red">xxx</span> (WCT)

## [FlashTPCBundle](./QLBundles/Bundle.md) (WCP) vs. <span style="color:red">xxx</span> (WCT)

## [Examine_bundles](https://github.com/BNLIF/wire-cell-2dtoy/blob/master/docs/ExamineBundles.md) is a function that we need to port, [Original ExamineBundles Code](https://github.com/BNLIF/wire-cell-2dtoy/blob/a30305052fc54bbbbbd826b096066d6e8777b54d/src/ExamineBundles.cxx)
- bundle->get_main_cluster()  (WCP) vs. <span style="color:red">xxx</span> (WCT)
- bundle->get_other_clusters() (WCP) vs. <span style="color:red">xxx</span> (WCT)
- bundle->get_main_cluster()->get_mcells() (WCP) vs. [loop blob](https://github.com/WireCell/wire-cell-toolkit/blob/apply-pointcloud/clus/src/clustering_separate.cxx#L1377) (WCT)
- Create a new temp cluster, PR3DCluster *new_cluster = new PR3DCluster(cluster_id) (WCP) vs.  <span style="color:red">xxx</span> (WCT)
- Add every blob to this cluster, new_cluster->AddCell(blob, blob->GetTimeSlice())  (WCP) vs. <span style="color:red">xxx</span> (WCT)
- std::vector<SMGCSelection> sep_mcells = new_cluster->Examine_graph(ct_point_cloud)  [code](https://github.com/BNLIF/wire-cell-data/blob/d5748d87c3113efcb44eed237bb48a10d60002d9/src/PR3DCluster.cxx#L2332) (WCP) vs. similar to [Create_graph](https://github.com/WireCell/wire-cell-toolkit/blob/apply-pointcloud/clus/src/Facade_Cluster.cxx#L1444) (WCT)
- [Connect_graph_overclustering_protection](https://github.com/BNLIF/wire-cell-data/blob/d5748d87c3113efcb44eed237bb48a10d60002d9/src/PR3DCluster.cxx#L1853) (WCP)  vs. not existing, but example in [connect_graph](https://github.com/WireCell/wire-cell-toolkit/blob/apply-pointcloud/clus/src/Facade_Cluster.cxx#L1444) (WCT)
- For each sep_mcells, we will form a new set of clusters, we need to find the cluster that overlapped the most with the original main cluster, and assign it as the new MAIN CLUSTER for the new bundle (achieved through compare blobs) (WCP) vs.  <span style="color:red">xxx</span> (WCT)

## [PRCluster->Create_Graph()](https://github.com/BNLIF/wire-cell/blob/master/uboone_nusel_app/apps/prod-wire-cell-matching-nusel.cxx#L817) (WCP) vs. [Create_graph](https://github.com/WireCell/wire-cell-toolkit/blob/apply-pointcloud/clus/src/Facade_Cluster.cxx#L1444) (WCT)

## [cluster->get_highest_lowest_wcps()](https://github.com/BNLIF/wire-cell/blob/master/uboone_nusel_app/apps/prod-wire-cell-matching-nusel.cxx#L819C103-L819C130) (WCP) vs. [get_ghiehst_lowst_points](https://github.com/WireCell/wire-cell-toolkit/blob/apply-pointcloud/clus/src/Facade_Cluster.cxx#L1241) (WCT)

<span style="color:red">WCT function only returns points, but we need the index to do the shortest path (following two functions), how to get point index from points? </span> 

## [cluster->->dijkstra_shortest_paths(wcps.first)](https://github.com/BNLIF/wire-cell/blob/master/uboone_nusel_app/apps/prod-wire-cell-matching-nusel.cxx#L822C25-L823C58) (WCP) vs. [dijkstra_shortest_paths using point index](https://github.com/WireCell/wire-cell-toolkit/blob/apply-pointcloud/clus/src/Facade_Cluster.cxx#L2631) (WCT)

## [cluster->cal_shortest_path(wcps.second)](https://github.com/BNLIF/wire-cell/blob/master/uboone_nusel_app/apps/prod-wire-cell-matching-nusel.cxx#L823) (WCP) vs. [cal_shortest_path using point index](https://github.com/WireCell/wire-cell-toolkit/blob/apply-pointcloud/clus/src/Facade_Cluster.cxx#L2664) (WCT)

## [Improve_PR3DCluster](https://github.com/BNLIF/wire-cell-2dtoy/blob/master/docs/Improve_PR3DCluster.md) is another function that we need to port. [Original Improve_PR3DCluster Code](https://github.com/BNLIF/wire-cell-2dtoy/blob/master/src/ImprovePR3DCluster.cxx)

- From existing clusters --> blobs --> time, fired channels --> get activities to be used to redo tiling (WCP) vs. <span style="color:red">xxx</span>  (WCT)
- Access shorest path std::list<WCPointCloud<double>::WCPoint>& wcps = cluster->get_path_wcps() (WCP) vs. [index of points along the path](https://github.com/WireCell/wire-cell-toolkit/blob/apply-pointcloud/clus/src/Facade_Cluster.cxx#L2689) (WCT)
- for a path point, use ct_point_cloud to convert the point into ch vs. time, and then judge if activities are available (WCP) vs. <span style="color:red">xxx</span>  (WCT)
- [update the activities according to the path_point's properties](https://github.com/BNLIF/wire-cell-2dtoy/blob/master/src/ImprovePR3DCluster.cxx#L136) (WCP) vs.  <span style="color:red">xxx</span>  (WCT)
- [redo tiling using the newly created activities](https://github.com/BNLIF/wire-cell-2dtoy/blob/master/src/ImprovePR3DCluster.cxx#L203) (WCP) vs.  <span style="color:red">xxx</span>  (WCT)
- Compare the newly created blobs vs. the (old) existing blobs using the Overlap_fast function and judge if the blob should be kept or not, then create a new cluster from the remaining blobs (WCP) vs. <span style="color:red">xxx</span>, [overlap_fast](https://github.com/WireCell/wire-cell-toolkit/blob/apply-pointcloud/clus/src/Facade_Cluster.cxx#L513)  (WCT)

## [map from old to new cluster](https://github.com/BNLIF/wire-cell/blob/master/uboone_nusel_app/apps/prod-wire-cell-matching-nusel.cxx#L831) (WCP) vs. <span style="color:red">xxx</span> (WCT)


# Include and after Trajectory Fitting (TBD)

## [ProtoSegment](https://github.com/BNLIF/wire-cell-pid/blob/537a3fd17f8a7b3cf5412594267c14c4cc1775cb/docs/protosegment.md) (WCP) vs. <span style="color:red">xxx</span> (WCT)

## [WCShower](https://github.com/BNLIF/wire-cell-pid/blob/537a3fd17f8a7b3cf5412594267c14c4cc1775cb/docs/wcshower.md) (WCP) vs. <span style="color:red">xxx</span> (WCT)

## [ProtoVertex](https://github.com/BNLIF/wire-cell-pid/blob/537a3fd17f8a7b3cf5412594267c14c4cc1775cb/docs/protovertex.md) (WCP) vs. <span style="color:red">xxx</span> (WCT)

## [Steiner Tree](https://github.com/BNLIF/wire-cell-pid/blob/537a3fd17f8a7b3cf5412594267c14c4cc1775cb/docs/PR3DCluster_steiner.md) (WCP) vvs. <span style="color:red">xxx</span> (WCT)
