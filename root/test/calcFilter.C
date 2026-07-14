

TH1F* get_power_density(TH2F* hraw, std::vector<int> channels,
                       std::string name,
		       int tick_offset=1700,int tick_offset2=1700, int nticks=200){

  TH1F* hfft= new TH1F(name.c_str(),"", nticks,0,1./0.512); // MHz, 512ns period
  hfft->GetXaxis()->SetTitle("Frequency (MHz)");
  hfft->GetYaxis()->SetTitle("|C(#omega)^{2}|");
    double k = double((tick_offset2-tick_offset))/double(channels.size());
    cout<<"k = "<<(tick_offset2-tick_offset)<<"  "<<channels.size()<<"  "<<k<<endl;
  for(auto channel: channels){

    int start_bin = k*(channel-channels[0])+tick_offset;
    cout<<"start_bin  "<<start_bin<<endl;
    TH1F *h2 = new TH1F(Form("Channel%d", channel),Form("Channel %d", channel),nticks,0,nticks);
    for (int i=0; i!=nticks; i++){
      int ibin = hraw->GetXaxis()->FindBin(channel);
      h2->SetBinContent(i+1,hraw->GetBinContent(ibin,i+1 + start_bin));
    }

    TH1F* hmag = (TH1F*)h2->FFT(0, "MAG");
    for(int i=1; i<nticks; i++){
      hfft->AddBinContent(i+1, hmag->GetBinContent(i+1) );
    }
    hmag->Delete();
  }
  hfft->Scale(1./channels.size());
  hfft->Scale(1400.0/(4096*4)); // convert ADC to mV

  for(int i=0; i<nticks; i++) {
    hfft->SetBinContent(i+1, std::pow(hfft->GetBinContent(i+1),2) );
  }

  return hfft;
}

void fill_range(vector<int>& v, int start, int end){
  v.resize(end-start+1);
  std::iota(v.begin(), v.end(), start);
}

TH1F* get_wiener_filter(TH1F* hsplusn_powden, TH1F* hnoise_powden) {
   const int nbins = hsplusn_powden->GetNbinsX();
   TH1F* hwiener = new TH1F("hwiener", "", nbins,
                                           hsplusn_powden->GetBinLowEdge(1),
                                           hsplusn_powden->GetBinLowEdge(nbins)
                                           + hsplusn_powden->GetBinWidth(nbins) );

   hwiener->GetXaxis()->SetTitle("Frequency (MHz)");
   for (int i=0; i<nbins; i++) {
       float M2 = hsplusn_powden->GetBinContent(i+1);
       float N2 = hnoise_powden->GetBinContent(i+1);
       hwiener->SetBinContent(i+1, M2==0? 0:(M2-N2)/M2);
   }
   return hwiener;
}

// main
void calcFilter(std::string magnify="/exp/dune/data/users/xning/data/protodunevd-sim-check-exa-for-wienerfilter.root"){

  gStyle->SetOptStat(0);

  std::vector<int> splusn_chans;
  std::vector<int> splusn_chans_data;
  std::vector<int> noise_chans;
  // selected channel number range
  // These numbers are manually selected according to the input file.(data/simulation)
  fill_range(splusn_chans, 1050, 1250);
  fill_range(noise_chans, 1050, 1250);


  TFile *file = TFile::Open(magnify.c_str());
   TH2F *hraw = (TH2F*)file->Get("hv_raw0");


  // lower edge of the track time: time_start=3550, time_end=3550; and time_window=200;
  // Usually for a track parallel to the plane and pependicular wires, the time_start and time_end is the same.
  // Some special cases will need to select tracks that are not parallel to the plane.
  // These numbers are manually selected according to the data/simulation

  // signal+noise
  TH1F* hsplusn_powden = (TH1F*)get_power_density(hraw, splusn_chans, "splusn", 3550,3550, 200);
  // noise only
  TH1F* hnoise_powden = (TH1F*)get_power_density(hraw, noise_chans, "signal", 1000,1000, 200);



  auto c1 = new TCanvas("c1","c1",600,600);
  hsplusn_powden->SetLineColor(kRed);
  hnoise_powden->SetLineColor(kBlue);
  // hnoise_powden->SetLineColor(kBlue);
  hsplusn_powden->GetXaxis()->SetRangeUser(0,1);
  // hsplusn_powden->GetYaxis()->SetRangeUser(0,500);
  hsplusn_powden->Draw("hist");
  hnoise_powden->Draw("hist same");
  // hnoise_powden->Draw("hist same");
    // cout<<hsplusn_powden->Integral()<<endl;
    // cout<<hnoise_powden->Integral()<<endl;
    // cout<<hsplusn_powden->Integral()/hnoise_powden->Integral()<<endl;
//
  auto c2 = new TCanvas("c2","c2",600,600);
  TH1F* hwiener = (TH1F*)get_wiener_filter(hsplusn_powden, hnoise_powden);
  hwiener->GetXaxis()->SetRangeUser(0,1/0.512/2);
  // fit wiener filter
//  float firstbin = hwiener->GetBinContent(2);
//  std::cout << "Enforce 1st bin to be 1, scaling with factor: " << 1.0/firstbin << std::endl;
  //hwiener->Scale(1./firstbin);
  float max = hwiener->GetMaximum();
  std::cout << "Enforce scaling with factor: " << 1.0/max << std::endl;
  hwiener->Scale(1./max); 
  hwiener->Draw("hist");
  auto f1 = new TF1("f1","exp(-0.5*(x/[0])**[1])", 0,1);
  // f1->SetParameters(0.3, 2);
  f1->SetParameters(0.1199, 3.8965);
  //f1->SetParameters(0.125448, 5.27080);
   hwiener->Fit("f1","R");
  f1->Draw("same");
}
