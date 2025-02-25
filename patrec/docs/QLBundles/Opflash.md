## Describe OpFlash in WCP

For more detailed information, please refer to the [OpFlash documentation](https://github.com/BNLIF/wire-cell-data/blob/master/docs/OpFlash.md).

## Example prototype jobs or files ...

A WCP rootfile can be found @ [this link](https://www.phy.bnl.gov/xqian/talks/wire-cell-porting/nuselEval_5384_137_6852.root)


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

See "light" and "flash" sections of the [WCT tensor data model document](../../../aux/docs/tensor-data-model.org).

More discussion is found in the [`Bundle`](./Bundle.md) document.

## Example WCT jobs or files ... 

## WCP's requirements 

N/A

## WCT's questions to confirm functionality 

