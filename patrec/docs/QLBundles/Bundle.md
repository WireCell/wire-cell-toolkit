## Describe FlashTPCBundle in WCP

For more detailed information, please refer to the [OpFlash documentation](https://github.com/BNLIF/wire-cell-data/blob/master/docs/FlashTPCBundle.md).

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

## Example WCT jobs or files ... 


## WCP's requirements 

N/A

## WCT's questions to confirm functionality 

