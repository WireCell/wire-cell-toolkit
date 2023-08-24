#include "WireCellImg/FrameQualityTagging.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellIface/IWaveformMap.h"
#include "WireCellAux/FrameTools.h"

#include "WireCellUtil/Units.h"
#include "WireCellUtil/NamedFactory.h"

#include <sstream>

WIRECELL_FACTORY(FrameQualityTagging, WireCell::Img::FrameQualityTagging, WireCell::IFrameFilter, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Img;

FrameQualityTagging::FrameQualityTagging()
    : Aux::Logger("FrameQualityTagging", "img")
{
}

FrameQualityTagging::~FrameQualityTagging() {}

WireCell::Configuration FrameQualityTagging::default_configuration() const
{
    Configuration cfg;

    //
    cfg["anode"] = "AnodePlane";
    //
    cfg["gauss_trace_tag"] = m_gauss_trace_tag;
    //
    cfg["wiener_trace_tag"] = m_wiener_trace_tag;
    //
    cfg["tick0"] = m_tick0;
    cfg["nticks"] = m_nticks;

    cfg["nrebin"] = m_nrebin;
    cfg["length_cut"] = m_length_cut;
    cfg["time_cut"] = m_time_cut;

    cfg["ch_threshold"] = m_ch_threshold;

    cfg["n_cover_cut1"] = m_n_cover_cut1;
    cfg["n_fire_cut1"] = m_n_fire_cut1;
    cfg["n_cover_cut2"] = m_n_cover_cut2;
    cfg["n_fire_cut2"] = m_n_fire_cut2;
    cfg["fire_threshold"] = m_fire_threshold;

    for(int i=0;i!=3;i++){
      cfg["n_cover_cut3"][i] = m_n_cover_cut3[i];
      cfg["percent_threshold"][i] = m_percent_threshold[i];
      cfg["threshold1"][i] = m_threshold1[i];
      cfg["threshold2"][i] = m_threshold2[i];
    }

    cfg["min_time"] = m_min_time;
    cfg["max_time"] = m_max_time;
    cfg["flag_corr"] = m_flag_corr;

    cfg["global_threshold"] = m_global_threshold;
    
    return cfg;
}

void FrameQualityTagging::configure(const WireCell::Configuration& cfg)
{
    auto anode_tn = cfg["anode"].asString();
    m_anode = Factory::find_tn<IAnodePlane>(anode_tn);
    m_gauss_trace_tag = get<std::string>(cfg, "gauss_trace_tag", m_gauss_trace_tag);
    m_wiener_trace_tag = get<std::string>(cfg, "wiener_trace_tag", m_wiener_trace_tag);
    m_cm_tag = get<std::string>(cfg, "cm_tag", m_cm_tag);
    m_tick0 = get<int>(cfg, "tick0", m_tick0);
    m_nticks = get<int>(cfg, "nticks", m_nticks);

    m_nrebin = get<int>(cfg, "nrebin", m_nrebin);
    m_length_cut = get<int>(cfg, "length_cut", m_length_cut);
    m_time_cut = get<int>(cfg, "time_cut", m_time_cut);
    m_ch_threshold = get<int>(cfg, "ch_threshold", m_ch_threshold);

    m_n_cover_cut1 = get<int>(cfg, "n_cover_cut1", m_n_cover_cut1);
    m_n_fire_cut1 = get<int>(cfg, "n_fire_cut1", m_n_fire_cut1);
    m_n_cover_cut2 = get<int>(cfg, "n_cover_cut2", m_n_cover_cut2);
    m_n_fire_cut2 = get<int>(cfg, "n_fire_cut2", m_n_cover_cut2);
    m_fire_threshold = get<float>(cfg,"fire_threshold", m_fire_threshold);

    if (cfg.isMember("n_cover_cut3")){
      m_n_cover_cut3.clear();
      for (auto jplane : cfg["n_cover_cut3"]){
	m_n_cover_cut3.push_back(jplane.asInt());
      }
    }
    if (cfg.isMember("percent_threshold")){
      m_percent_threshold.clear();
      for (auto jplane : cfg["percent_threshold"]){
	m_percent_threshold.push_back(jplane.asFloat());
      }
    }
    if (cfg.isMember("threshold1")){
      m_threshold1.clear();
      for (auto jplane : cfg["threshold1"]){
	m_threshold1.push_back(jplane.asFloat());
      }
    }
    if (cfg.isMember("threshold2")){
      m_threshold2.clear();
      for (auto jplane : cfg["threshold2"]){
	m_threshold2.push_back(jplane.asFloat());
      }
    }
    
    
    m_min_time = get<int>(cfg,"min_time", m_min_time);
    m_max_time = get<int>(cfg,"max_time", m_max_time);
    m_flag_corr = get<int>(cfg,"flag_corr",m_flag_corr);

    m_global_threshold = get<float>(cfg,"global_threshold", m_global_threshold);
    
}

bool FrameQualityTagging::operator()(const input_pointer& in, output_pointer& out)
{
    // frame quality
    // no mechanism to consume it for now
    int frame_quality = 0;

    out = nullptr;
    if (!in) {
        log->debug("see EOS");
        return true;  // eos
    }

    // Prepare a cmm for output
    auto cmm = in->masks();
    if (cmm.find(m_cm_tag)==cmm.end()) {
        THROW(RuntimeError()<< errmsg{"no ChannelMask with name "+m_cm_tag});
    }
    auto& cm = cmm[m_cm_tag];
    log->debug("input: {} size: {}", m_cm_tag, cm.size());

    // Sort out channel list for each plane
    std::unordered_map<int, Aux::channel_list> chbyplane;
    m_anode->channels();
    for (auto channel : m_anode->channels()) {
        auto iplane = (m_anode->resolve(channel)).index();
        if (iplane < 0 || iplane >= 3) {
            THROW(RuntimeError() << errmsg{"Illegal wpid"});
        }
        // TODO: asumming this sorted?
        chbyplane[iplane].push_back(channel);
    }

    // Prepare Eigen matrices for gauss, wiener
    std::vector<size_t> nchannels = {chbyplane[0].size(), chbyplane[1].size(), chbyplane[2].size()};
    log->debug("nchannels: {} {} {}", nchannels[0], nchannels[1], nchannels[2]);
    std::vector<Array::array_xxf> arr_gauss;
    std::vector<Array::array_xxf> arr_wiener;
    auto trace_gauss = Aux::tagged_traces(in, m_gauss_trace_tag);
    auto trace_wiener = Aux::tagged_traces(in, m_wiener_trace_tag);
    for (int iplane=0; iplane< 3; ++iplane) {
        arr_gauss.push_back(Array::array_xxf::Zero(nchannels[iplane], m_nticks));
        Aux::fill(arr_gauss.back(), trace_gauss, chbyplane[iplane].begin(), chbyplane[iplane].end(), m_tick0);
        arr_wiener.push_back(Array::array_xxf::Zero(nchannels[iplane], m_nticks));
        Aux::fill(arr_wiener.back(), trace_wiener, chbyplane[iplane].begin(), chbyplane[iplane].end(), m_tick0);
    }
    //    std::cout << arr_gauss[0] << std::endl;

    // Prepare threshold for wiener
    auto threshold = in->trace_summary(m_wiener_trace_tag);
    log->debug("threshold.size() {}", threshold.size());

    // begin for Xin
    //    log->debug("{}", m_nticks);

    std::vector<std::vector<std::tuple<int,int,float>> > summary_plane(nchannels.size());

    for (int i=0;i<std::round(m_nticks/m_nrebin);i++){
      std::vector<std::vector<std::pair<int,int>> > ROI_plane(nchannels.size());

      // loop over each plane ...
      for (size_t iplane = 0; iplane != nchannels.size(); iplane ++){
	bool flag_roi = false;
	int start = 0, end = 0;

	// loop channels in a plane ...
	for (size_t j=0;j!=chbyplane[iplane].size();j++){
	  float content = 0;
	  for (int k=0;k!=m_nrebin;k++){
	    int time = i * m_nrebin + k;
	    if (time < m_nticks){
	      content += arr_wiener[iplane](j,time);
	    }
	  }
	  
	  float rms = threshold.at(chbyplane[iplane].at(j));
	  //if (content > 0 && i <=1) log->debug("{}, {}, {}, {}, {}",iplane, i, j, content, rms);

	  if (!flag_roi){
	    if (content > rms){
	      start = j;
	      flag_roi = true;
	    }
	  }else{
	    if (content <= rms || j + 1 == chbyplane[iplane].size()){
	      end = j;
	      flag_roi = false;
	      if (end - start > m_length_cut)
		ROI_plane[iplane].push_back(std::make_pair(start,end));
	    }
	  }
	}
      }

      //if (i<=10)
      //	log->debug("{}, {}, {}, {}", i, ROI_plane[0].size(), ROI_plane[1].size(), ROI_plane[2].size());

      std::vector<float> nfired;
      std::vector<int> ncover;
      nfired.resize(chbyplane.size(),0);
      ncover.resize(chbyplane.size(),0);
      
      for (size_t iplane = 0; iplane != nchannels.size(); iplane ++){
	if (ROI_plane[iplane].size()>0){
	  for (size_t j=0;j!=ROI_plane[iplane].size();j++){
	    nfired[iplane] += ROI_plane[iplane].at(j).second - ROI_plane[iplane].at(j).first;
	  }
	  ncover[iplane] = ROI_plane[iplane].back().second - ROI_plane[iplane].front().first;
	  nfired[iplane] /= ncover[iplane];

	  if (ncover[iplane] * nfired[iplane] > m_ch_threshold)
	    summary_plane[iplane].push_back(std::make_tuple(i, ncover[iplane], nfired[iplane]));
	  
	}
      }
    }

    //    log->debug("{},{},{}",summary_plane[0].size(), summary_plane[1].size(), summary_plane[2].size());
    
    std::vector<bool> flag_plane(nchannels.size(), false);
    std::vector<std::vector<std::pair<int, int> > > noisy_plane(nchannels.size());
    
    for (size_t iplane = 0; iplane!=nchannels.size(); iplane ++){
      int prev_time = -1;
      int n_cover = 0;
      int n_fire = 0;
      int start_time =0;
      int end_time = 0;
      // int acc_cover = 0;
      int acc_total = 0;
      int acc_fire = 0;

      for (size_t i=0;i!=summary_plane[iplane].size();i++){
	int time = std::get<0>(summary_plane[iplane].at(i));
	int ncover = std::get<1>(summary_plane[iplane].at(i));
	float percentage = std::get<2>(summary_plane[iplane].at(i));

	//	log->debug("{}, {}, {}", time, ncover, percentage);
	if (time > prev_time + m_time_cut){
	  if ((n_cover >= m_n_cover_cut1 && n_fire >= m_n_fire_cut1) || (n_cover >= m_n_cover_cut2 && n_fire >= m_n_fire_cut2 && acc_fire > m_fire_threshold * acc_total)){
	    end_time = prev_time;
	    log->debug("{}: {}, {}, {}, {}",iplane, n_cover, n_fire, start_time, end_time);
	    noisy_plane[iplane].push_back(std::make_pair(start_time, end_time));
	    flag_plane[iplane] = true;
	  }
	  n_cover = 0;
	  n_fire = 0;
	  // acc_cover = 0;
	  acc_total = 0;
	  acc_fire = 0;
	  start_time = time;
	  end_time = time;
	}

	if (ncover > m_n_cover_cut3[iplane] && (percentage > m_percent_threshold[iplane] || ncover * percentage > m_threshold1[iplane]))
	  n_cover ++;
	if (ncover * percentage > m_threshold2[iplane]){
	  n_fire ++;
	  acc_total += nchannels[iplane];
	  // acc_cover += ncover;
	  acc_fire  += ncover*percentage;
	}
	prev_time = time;
      }
      if ((n_cover >= m_n_cover_cut1 && n_fire >= m_n_fire_cut1) || (n_cover >= m_n_cover_cut2 && n_fire >= m_n_fire_cut2 && acc_fire > m_fire_threshold * acc_total)){
	end_time = prev_time;
	log->debug("Plane {}: {}, {}, {}, {}",iplane, n_cover, n_fire, start_time, end_time);
	noisy_plane[iplane].push_back(std::make_pair(start_time, end_time));
	flag_plane[iplane] = true;
      }
    }

    // log->debug("{},{},{}",flag_plane[0], flag_plane[1],flag_plane[2]);
    //    log->debug("{},{},{},{},{},{},{},{}", m_length_cut, m_time_cut, m_ch_threshold, m_n_cover_cut1, m_n_fire_cut1, m_n_cover_cut2, m_n_fire_cut2, m_fire_threshold);
    // for (int i=0;i!=3;i++){
    //   log->debug("{}, {}, {}, {}",m_n_cover_cut3[i], m_percent_threshold[i], m_threshold1[i], m_threshold2[i]);
    // }

    // test
    // noisy_plane[0].push_back(std::make_pair(100,200));
    // noisy_plane[1].push_back(std::make_pair(100,200));
    // noisy_plane[2].push_back(std::make_pair(100,200));
    
    // extend size ...
    for (size_t iplane = 0; iplane != nchannels.size(); iplane++){
      for (auto it = noisy_plane[iplane].begin(); it != noisy_plane[iplane].end(); it++){
	int start_time = it->first; 
	int end_time = it->second;
	
	int min_start_time = start_time;
	int max_end_time = end_time;

	for (size_t i=0;i!=nchannels[iplane];i++){
	  int j = start_time;
	  float content = 0;

	  for (int k=0;k!=m_nrebin;k++){
	    int time = j * m_nrebin + k;
	    if (time < m_nticks && time >=0){
	      content += arr_wiener[iplane](i,time);
	    }
	  }

	  while(content >0){
	    j--;
	    content = 0;
	    for (int k=0;k!=m_nrebin;k++){
	      int time = j * m_nrebin + k;
	      if (time < m_nticks && time >=0){
		content += arr_wiener[iplane](i,time);
	      }
	    }
	  }
	  if (j < min_start_time) min_start_time = j;

	  j = end_time;
	  content = 0;
	  for (int k=0;k!=m_nrebin;k++){
	    int time = j * m_nrebin + k;
	    if (time < m_nticks && time >=0){
	      content += arr_wiener[iplane](i,time);
	    }
	  }

	  while(content>0){
	    j++;
	    content =0;
	    for (int k=0;k!=m_nrebin;k++){
	      int time = j * m_nrebin + k;
	      if (time < m_nticks && time >=0){
		content += arr_wiener[iplane](i,time);
	      }
	    }
	  }

	  if (j > max_end_time) max_end_time = j;
	}

	if (min_start_time<0) min_start_time = 0;
	if (max_end_time>=std::round(m_nticks/m_nrebin)) max_end_time = std::round(m_nticks/m_nrebin) - 1;

	it->first = min_start_time;
	it->second = max_end_time;

	//	log->debug("Plane {}: {},{},{},{}",iplane, start_time, end_time, min_start_time, max_end_time);
	
      }
    }


    bool flag_active = false;
    for (size_t iplane = 0; iplane != nchannels.size(); iplane++){
      for (auto it = noisy_plane[iplane].begin(); it!= noisy_plane[iplane].end(); it++){
	int start_time = std::max( it->first * m_nrebin ,0);
	int end_time = std::min(it->second * m_nrebin + 3, m_nticks);
	if (end_time > m_min_time && start_time < m_max_time)
	  flag_active = true;

	//	log->debug("{}, {}", start_time, end_time);
	
	if (m_flag_corr==1){
	  // log->debug("{}", cm.size());
	  for (size_t jplane = 0; jplane != nchannels.size(); jplane ++){
	    for (size_t i=0;i!=	chbyplane[jplane].size(); i++){ // channel ...
	      int chid = chbyplane[jplane][i];

	      if (cm.find(chid)==cm.end()){
		cm[chid].push_back(std::make_pair(start_time, end_time));
	      }else{
		int tmp_time = cm[chid].front().first;
		cm[chid].front().first = std::min(start_time, tmp_time);
		tmp_time = cm[chid].front().second;
		cm[chid].front().second = std::max(tmp_time, end_time);
	      }
	    }
	  }
	  //log->debug("{}", cm.size());
	}
      }
    }


    if (flag_plane[0] || flag_plane[1] || flag_plane[2]){
      if (flag_active) frame_quality = 2;
      else frame_quality = 1;
    }
    
    {
      int ntime = std::round(m_nticks/m_nrebin);
      std::vector<float> sum_total_plane(nchannels.size(),0);
      std::vector<float> sum_fired_plane(nchannels.size(),0);
      for (size_t iplane = 0; iplane != nchannels.size(); iplane++){

	for (size_t i=0;i!=nchannels[iplane];i++){ // channel ...

	  for (int j=0; j!= ntime; j++){

	    float content = 0;
	    for (int k=0;k!=m_nrebin;k++){
	      int time = j * m_nrebin + k;
	      if (time < m_nticks && time >=0){
		content += arr_wiener[iplane](i,time);
	      }
	    }
	    
	    if(content >0) sum_fired_plane[iplane] ++;
	    sum_total_plane[iplane] ++;
	  }
	}
      }
      if (sum_fired_plane[0]/sum_total_plane[0] +  sum_fired_plane[1]/sum_total_plane[1] + sum_fired_plane[2]/sum_total_plane[2] > m_global_threshold){
	log->debug("Too busy: {}, {}, {}",sum_fired_plane[0]/ sum_total_plane[0], sum_fired_plane[1]/sum_total_plane[1], sum_fired_plane[2]/ sum_total_plane[2]);
	frame_quality = 3;

	if (m_flag_corr){
	  for (size_t jplane = 0; jplane != nchannels.size(); jplane ++){
	    for (size_t i=0;i!=	chbyplane[jplane].size(); i++){ // channel ...
	      int chid = chbyplane[jplane][i];

	      if (cm.find(chid)==cm.end()){
		cm[chid].push_back(std::make_pair(0, m_nticks));
	      }else{
		cm[chid].front() = std::make_pair(0, m_nticks);
	      }
	      
	    }
	  }
	}      
      }
    }
    
    
    
    // end for Xin

    // good frame, pass trhough
    if (frame_quality == 0 || frame_quality == 1) {
      log->debug("good frame, pass through, flag {}", frame_quality);
        // make a copy of all input trace pointers
        ITrace::vector out_traces(*in->traces());
        // Basic frame stays the same.
        auto sfout = new Aux::SimpleFrame(in->ident(), in->time(), out_traces, in->tick(), cmm);
        // passing through other parts of the original frame
        for (auto ftag : in->frame_tags()) {
            sfout->tag_frame(ftag);
        }
        for (auto inttag : in->trace_tags()) {
            const auto& traces = in->tagged_traces(inttag);
            const auto& summary = in->trace_summary(inttag);
            sfout->tag_traces(inttag, traces, summary);
        };
        out = IFrame::pointer(sfout);
        log->debug("output: {} size: {}", m_cm_tag, out->masks()[m_cm_tag].size());
        return true;
    }

    // bad frame, out = nullptr for now
    return true;
}
