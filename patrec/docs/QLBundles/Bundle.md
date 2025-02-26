## Describe FlashTPCBundle in WCP

See also the [OpFlash documentation](https://github.com/BNLIF/wire-cell-data/blob/master/docs/FlashTPCBundle.md) for specifics about that data product.  

## Example prototype jobs or files ...

A WCP rootfile can be found @ [this link](https://www.phy.bnl.gov/xqian/talks/wire-cell-porting/nuselEval_5384_137_6852.root)

Bundle information are saved in 
```cpp
  TTree *T_match1 = new TTree("T_match","T_match");
  T_match1->SetDirectory(file1);
  Int_t ncluster;
  T_match1->Branch("tpc_cluster_id",&ncluster,"tpc_cluster_id/I");  // TPC parent cluster id (see TC tree for more explainations)
  T_match1->Branch("flash_id",&flash_id,"flash_id/I"); // PMT Flash ID, see flash_id for more information
  T_match1->Branch("event_type",&event_type,"event_type/I"); // this is to save the event tagger information (saved as bits in WCP)
  Double_t flash_time;
  T_match1->Branch("flash_time",&flash_time,"flash_time/D"); // Flash time 
  cluster_length = 0;
  T_match1->Branch("cluster_length",&cluster_length,"cluster_length/D"); // cluster length for main cluster
```


TPC Blob information are saved in 
```cpp

  // load mcell
  TTree *TC = (TTree*)file->Get("TC");
  std::vector<int> *cluster_id_vec = new std::vector<int>;
  std::vector<int> *parent_cluster_id = new std::vector<int>;
  std::vector<int> *time_slice_vec = new std::vector<int>;
  std::vector<double> *q_vec = new std::vector<double>;
  std::vector<double> *uq_vec = new std::vector<double>;
  std::vector<double> *vq_vec = new std::vector<double>;
  std::vector<double> *wq_vec = new std::vector<double>;
  std::vector<double> *udq_vec = new std::vector<double>;
  std::vector<double> *vdq_vec = new std::vector<double>;
  std::vector<double> *wdq_vec = new std::vector<double>;

  std::vector<int> *nwire_u_vec = new  std::vector<int>;
  std::vector<int> *nwire_v_vec = new  std::vector<int>;
  std::vector<int> *nwire_w_vec = new  std::vector<int>;
  std::vector<int> *flag_u_vec = new  std::vector<int>;
  std::vector<int> *flag_v_vec = new  std::vector<int>;
  std::vector<int> *flag_w_vec = new  std::vector<int>;

  std::vector<std::vector<int>> *wire_index_u_vec = new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *wire_index_v_vec = new std::vector<std::vector<int>>;
  std::vector<std::vector<int>> *wire_index_w_vec = new std::vector<std::vector<int>>;
  std::vector<std::vector<double>> *wire_charge_u_vec = new std::vector<std::vector<double>>;
  std::vector<std::vector<double>> *wire_charge_v_vec = new std::vector<std::vector<double>>;
  std::vector<std::vector<double>> *wire_charge_w_vec = new std::vector<std::vector<double>>;
  std::vector<std::vector<double>> *wire_charge_err_u_vec = new std::vector<std::vector<double>>;
  std::vector<std::vector<double>> *wire_charge_err_v_vec = new std::vector<std::vector<double>>;
  std::vector<std::vector<double>> *wire_charge_err_w_vec = new std::vector<std::vector<double>>;
  
  TC->SetBranchAddress("cluster_id",&cluster_id_vec);  // actual cluster id where this blob belons
  TC->SetBranchAddress("parent_cluster_id",&parent_cluster_id);  // main cluster id that is used in T_match bundle
  TC->SetBranchAddress("time_slice",&time_slice_vec);
  TC->SetBranchAddress("q",&q_vec);
  TC->SetBranchAddress("uq",&uq_vec);
  TC->SetBranchAddress("vq",&vq_vec);
  TC->SetBranchAddress("wq",&wq_vec);
  TC->SetBranchAddress("udq",&udq_vec);
  TC->SetBranchAddress("vdq",&vdq_vec);
  TC->SetBranchAddress("wdq",&wdq_vec);
  TC->SetBranchAddress("nwire_u",&nwire_u_vec);
  TC->SetBranchAddress("nwire_v",&nwire_v_vec);
  TC->SetBranchAddress("nwire_w",&nwire_w_vec);
  TC->SetBranchAddress("flag_u",&flag_u_vec);
  TC->SetBranchAddress("flag_v",&flag_v_vec);
  TC->SetBranchAddress("flag_w",&flag_w_vec);
  TC->SetBranchAddress("wire_index_u",&wire_index_u_vec);
  TC->SetBranchAddress("wire_index_v",&wire_index_v_vec);
  TC->SetBranchAddress("wire_index_w",&wire_index_w_vec);
  TC->SetBranchAddress("wire_charge_u",&wire_charge_u_vec);
  TC->SetBranchAddress("wire_charge_v",&wire_charge_v_vec);
  TC->SetBranchAddress("wire_charge_w",&wire_charge_w_vec);
  TC->SetBranchAddress("wire_charge_err_u",&wire_charge_err_u_vec);
  TC->SetBranchAddress("wire_charge_err_v",&wire_charge_err_v_vec);
  TC->SetBranchAddress("wire_charge_err_w",&wire_charge_err_w_vec);
```


Opflash are saved in 
```cpp
  TTree *T_flash = (TTree*)file->Get("T_flash");
  Double_t time;
  Int_t type;
  Int_t flash_id;
  Int_t temp_run_no, temp_subrun_no, temp_event_no;
  T_flash->SetBranchAddress("runNo",&temp_run_no);
  T_flash->SetBranchAddress("subRunNo",&temp_subrun_no);
  T_flash->SetBranchAddress("eventNo",&temp_event_no);
  T_flash->SetBranchAddress("time",&time);
  T_flash->SetBranchAddress("type",&type);  // flash type, full waveform or Cosmic mode, two different types in MicroBooNE
  T_flash->SetBranchAddress("flash_id",&flash_id);   // this id is useful for matching with TPC object in bundle
  Double_t low_time, high_time, total_PE;
  Double_t temp_PE[32], temp_PE_err[32];
  std::vector<int> *fired_channels = new std::vector<int>;
  std::vector<double> *l1_fired_time = new std::vector<double>;
  std::vector<double> *l1_fired_pe = new std::vector<double>;
  T_flash->SetBranchAddress("low_time",&low_time);   // start time of flash
  T_flash->SetBranchAddress("high_time",&high_time);  // end time of flash
  T_flash->SetBranchAddress("total_PE",&total_PE);   // total PE
  T_flash->SetBranchAddress("PE",temp_PE);           // PE for each PMT
  T_flash->SetBranchAddress("PE_err",temp_PE_err);   // PE_err for each PMT
  T_flash->SetBranchAddress("fired_channels",&fired_channels);    // which channel are included in flash
  T_flash->SetBranchAddress("l1_fired_time",&l1_fired_time);      // advanced flash info 
  T_flash->SetBranchAddress("l1_fired_pe",&l1_fired_pe);          // advanced flash info
```

## Describe WCT version

### Review uboone blob loading

The `Trun`, `TC` and `TDC` ROOT trees may be loaded into WCT via [`UbooneBlobSource`](../../../root/src/UbooneBlobSource.cxx) to produce an `IBlobSet` of either "live" or "dead" blobs.  Two sources are needed to load both "live" and "dead" blobs simultaneously.

The [`uboone-blobs.smake` Snakemake workflow](../../../root/test/uboone-blobs.smake) runs `wire-cell` on [`uboone-blobs.jsonnet`](../../../root/test/uboone-blobs.jsonnet) in three kinds of modes: "live", "dead" and "clus".  

The "live" and "dead" modes each produce a *WCT cluster file* with a graph like:

```
UbooneBlobSources -> BlobClustering -> GlobalGeomClustering -> ClusterFileSink
```

The "live" mode has four sources `live-{uvw,uv,vw,wu}` and "dead" has three sources `dead-{uv,vw,wu}`.  

The "clus" job is like:

```
ClusterFileSources -> PointTreeBuilding -> MultiAlgBlobClustering -> TensorFileSink
```
There are two sources here loading each of "live" and "dead" cluster files.

The output of ["MABC"](../../../clus/src/MultiAlgBlobClustering.cxx) is a point-cloud tree in *WCT tensor data model* representation (`ITensorSet`).

:warning: This workflow uses uboone blobs but ignores uboone cluster info (`cluster_id`).

### Uboone cluster and flash loading

In order to consider `cluster_id` when loading `TC` and `T_Flash` information into WCT requires a complex procedure.  This is in large part due to the fact that WCT clusters are not "identified" per se but rather are emergent from the "connected components" graph operation.  After the operation we may identify a cluster by its ID in the "connected components" array which is effectively arbitrary.  Thus we must form WCT clusters in a context that still retains the association between blobs and their `cluster_id`.  But further, to represent the cluster-flash associations requires point-cloud level data tier.  And this requires blob sampling.  

We will call this new component `UbooneClusterSource` which supplies these operations:

- Job has the usual 4 `UbooneBlobSource`'s to get `IBlobSet`'s spanning uvw/uv/vw/wu view cases.
  - The `IBLobs` carry the ROOT `TTree::entry` number for their `TC` origin giving "blob order".
- A `BlobSetMerge` follows to provide 4-to-1 fan-in.

Then comer `UbooneClusterSource`.
- Each input `IBlobSet` is buffered (reminder: it is per-slice).
- Produces an output queue of `ITensorSet`.
- The output queue is expected to hold either zero or one `ITensorSet`.
- On EOS or change in frame ident, buffered `IBlobs` are converted and flushed to `ITensorSet`.
  - This is same behavior as `BlobClustering`

The flushing entails:

- Loads `TC` from file to get blob info including `cluster_id`
- Build CLUSTER ID ORDER as ordered  `cluster_id` set and map to index.
- Copy-paste parts of `PointTreeBuilding::sample_live()` and use to produce initial pc-tree.
  - We must define clusters based on `cluster_id` and not "connected components".
  - Cluster nodes are added to the grouping node in CLUSTER ID ORDER.
  - Map `cluster_id` to cluster `node_t*`.
  - Walk blobs, using their `tpc_cluster_id` to find the cluster node to receive them.
- Load `T_flash` ROOT TTree to produce **light**, **flash** and **flashlight** data as described in the tensor-data-model.
  - Store each of these three as arrays in a "local" PC on the "grouping" pc-tree node.
- Load `T_match` to get cluster-flash associations.
  - Form map from `cluster_id` to INDEX of **flash** array.
  - Store this index as a scalar `flash` array in a `cluster_scalar` local PC.  
  - Flash-less clusters get index of -1.
  - Also store `cluster_id` in this `cluster_scalar` cluster-local PC.
- Convert pc-tree to `ITensorSet` and output.
  - Reminder, all local PCs of the same name must be the same "shape".  Blob-nodes have a `scalar` local PC thus the name difference with `cluster_scaler`.
  
Overall graph configuration constraints:

- User **must** configure `Uboone*Source` to use the same stream of ROOT files.



## Example WCT jobs or files ... 

See [uboone.org](../../../root/docs/uboone.org).  In particular:

```
wire-cell -l stderr -L debug -A infiles=nuselEval_5384_137_6852.root root/test/uboone-clusters.jsonnet
```

And see `aux/src/ClusterFlashDump.cxx` as an example of a starting point for a real next-stage.  

Some comments on next stage:

- The `uboone-clusters.jsonnet` graph represents a single face pipeline (that's all of uboone).
- To do cluster-flash matching on other dets, we must merge across face+APAs.
  - This merge needs w.rt. optical `light/flash/flashlight` arrays.
  - On one hand, these are common across all face+APAs in which case we just take one.
  - OTOH, if they differ they must be appended **and** we must rewrite the `flash` scalar array in the `cluster_scalar` PCs.
  - So, the merge needs to know what previous people did to prepare the data.

- The `uboone-clusters.jsonnet` ends with `ClusterFlashDump`.
  - Obviously, replaced this with an `ITensorSet -> ITensorSet` filter is needed for cluster-flash making or other refinement.
  - This replacement node likely succeed in the graph the "big merge" node described above.
  - This new node type should probably look similar to `MABC` in that it allows a pc-tree to be passed through a pipeline of functions without requiring I/O through the TDM.
    - This pipeline could be hoisted up to the flow graph by making a new `IData` that passes a pc-tree w/out requiring a round trip through the TDM.

## WCP's requirements 

N/A

## WCT's questions to confirm functionality 

