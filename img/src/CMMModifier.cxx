#include "WireCellImg/CMMModifier.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/FrameTools.h"

#include "WireCellIface/IWaveformMap.h"

#include "WireCellUtil/Units.h"
#include "WireCellUtil/NamedFactory.h"

#include <sstream>

WIRECELL_FACTORY(CMMModifier, WireCell::Img::CMMModifier, WireCell::IFrameFilter, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Img;

CMMModifier::CMMModifier()
    : Aux::Logger("CMMModifier", "img")
{
}

CMMModifier::~CMMModifier() {}

WireCell::Configuration CMMModifier::default_configuration() const
{
    Configuration cfg;

    //
    cfg["anode"] = "AnodePlane";
    //
    cfg["cm_tag"] = m_cm_tag;
    //
    cfg["trace_tag"] = m_trace_tag;

    cfg["m_start"] = m_start;
    cfg["m_end"] = m_end;
    
    // bad channels ...
    cfg["ncount_cont_ch"] = m_ncount_cont_ch;
    for(int i=0;i!=m_ncount_cont_ch;i++){
      cfg["cont_ch_llimit"][i] = m_cont_ch_llimit[i];
      cfg["cont_ch_hlimit"][i] = m_cont_ch_hlimit[i];
    }

    cfg["ncount_veto_ch"] = m_ncount_veto_ch;
    for(int i=0;i!=m_ncount_veto_ch;i++){
      cfg["veto_ch_llimit"][i] = m_veto_ch_llimit[i];
      cfg["veto_ch_hlimit"][i] = m_veto_ch_hlimit[i];
    }

    cfg["dead_ch_ncount"] = m_dead_ch_ncount;
    cfg["dead_ch_charge"] = m_dead_ch_charge;
    cfg["ncount_dead_ch"] = m_ncount_dead_ch;
    for(int i=0;i!=m_ncount_dead_ch;i++){
      cfg["dead_ch_llimit"][i] = m_dead_ch_llimit[i];
      cfg["dead_ch_hlimit"][i] = m_dead_ch_hlimit[i];
    }

    cfg["ncount_org"] = m_dead_org;
    for (int i=0;i!=m_dead_org;i++){
      cfg["org_llimit"] = m_dead_org_llimit[i];
      cfg["org_hlimit"] = m_dead_org_hlimit[i];
    }
    
    return cfg;
}

void CMMModifier::configure(const WireCell::Configuration& cfg)
{
    auto anode_tn = cfg["anode"].asString();
    m_anode = Factory::find_tn<IAnodePlane>(anode_tn);
    m_cm_tag = get<std::string>(cfg, "cm_tag", m_cm_tag);
    m_trace_tag = get<std::string>(cfg, "trace_tag", m_trace_tag);


    m_start = get<int>(cfg,"start",m_start);
    m_end = get<int>(cfg,"end",m_end);
    
    m_ncount_cont_ch = get<int>(cfg,"ncount_cont_ch",m_ncount_cont_ch);
    if (cfg.isMember("cont_ch_llimit")){
      m_cont_ch_llimit.clear();
      for (auto value: cfg["cont_ch_llimit"]){
	m_cont_ch_llimit.push_back(value.asInt());
      }
    }
    if (cfg.isMember("cont_ch_hlimit")){
      m_cont_ch_hlimit.clear();
      for (auto value: cfg["cont_ch_hlimit"]){
	m_cont_ch_hlimit.push_back(value.asInt());
      }
    }

    m_ncount_veto_ch = get<int>(cfg,"ncount_veto_ch",m_ncount_veto_ch);
    if (cfg.isMember("veto_ch_llimit")){
      m_veto_ch_llimit.clear();
      for (auto value: cfg["veto_ch_llimit"]){
	m_veto_ch_llimit.push_back(value.asInt());
      }
    }
    if (cfg.isMember("veto_ch_hlimit")){
      m_veto_ch_hlimit.clear();
      for (auto value: cfg["veto_ch_hlimit"]){
	m_veto_ch_hlimit.push_back(value.asInt());
      }
    }

    m_dead_ch_ncount = get<int>(cfg,"dead_ch_ncount",m_dead_ch_ncount);
    m_dead_ch_charge = get<float>(cfg,"dead_ch_charge",m_dead_ch_charge);
    m_ncount_dead_ch = get<int>(cfg,"ncount_dead_ch",m_ncount_dead_ch);
    if (cfg.isMember("dead_ch_llimit")){
      m_dead_ch_llimit.clear();
      for (auto value: cfg["dead_ch_llimit"]){
	m_dead_ch_llimit.push_back(value.asInt());
      }
    }
    if (cfg.isMember("dead_ch_hlimit")){
      m_dead_ch_hlimit.clear();
      for (auto value: cfg["dead_ch_hlimit"]){
	m_dead_ch_hlimit.push_back(value.asInt());
      }
    }

    m_dead_org = get<int>(cfg,"ncount_org",m_dead_org);
    if (cfg.isMember("org_llimit")){
      m_dead_org_llimit.clear();
      for (auto value: cfg["org_llimit"]){
	m_dead_org_llimit.push_back(value.asInt());
      }
    }
    if (cfg.isMember("org_hlimit")){
      m_dead_org_hlimit.clear();
      for (auto value: cfg["org_hlimit"]){
	m_dead_org_hlimit.push_back(value.asInt());
      }
    }
    
    // log->debug("Xin: {} , {} ", m_dead_ch_ncount, m_dead_ch_charge);
}

bool CMMModifier::operator()(const input_pointer& in, output_pointer& out)
{
    ++ m_count;

    out = nullptr;
    if (!in) {
        log->debug("see EOS at call={}", m_count-1);
        return true;  // eos
    }

    log->debug("call={} input frame: {}", m_count-1, Aux::taginfo(in));

    // copy a CMM from input frame
    auto cmm = in->masks();

    // ------------ Debug Ewerton 2024-02-29 -----------
    //auto mycmm = in->masks();
    for(auto x : cmm) //cmm[m_cm_tag])
    {
      std::cout << "\n [CMMModifier-Ew] cmm.first=" << x.first << "\n";// " (" << x.second.first << ", " << x.second.second << ")\n";
    }
    std::cout << "\n [CMMModifier-Ew] end of inmasks loop" ;
    // ------------ End debug Ewerton 2024-02-29 -----------
                                   

    if (cmm.find(m_cm_tag)==cmm.end()) {
        log->warn("no ChannelMask with name \""+m_cm_tag+"\", will create one");
    }
    auto& cm = cmm[m_cm_tag];
    log->debug("input: {} size: {}", m_cm_tag, cm.size());

    // begin for Xin
    // some dummy ops for now
    // cm.clear();
    //    int nq=0;
    // int tbin=0;
    //for (auto trace : Aux::tagged_traces(in, m_trace_tag)) {
    //  tbin = trace->tbin();
    //     const int chid = trace->channel();
    //  const auto& charge = trace->charge();
    //  nq = charge.size();
    //  break;
    //     cm[chid].push_back({tbin,tbin+nq});
    // }
    // for (auto it = cm.begin(); it!=cm.end(); it++){
    //   log->debug("{}, {}, {}", it->first, it->second.at(0).first, it->second.at(0).second);
    // }

    // shorted wires ...
    for (int i = 0;i != m_ncount_cont_ch; i++){
      for (int ch = m_cont_ch_llimit[i]; ch != m_cont_ch_hlimit[i]; ch++){
	//	log->debug("{}", ch);
	if (cm.find(ch) == cm.end()){
	  if ((cm.find(ch-1) != cm.end() && cm.find(ch+2) != cm.end())){
	    log->debug("(Shorted) {} added to the bad channel list", ch);
	    cm[ch] = cm[ch-1];
	  } else if ((cm.find(ch-2) != cm.end() && cm.find(ch+1) != cm.end())){
	    log->debug("(Shorted) {} added to the bad channel list", ch);
	    cm[ch] = cm[ch+1];
	  }
	}
      }
    }

    // veto these wires anyway
    for (int i=0;i!=m_ncount_veto_ch;i++){
      for (int ch = m_veto_ch_llimit[i]; ch != m_veto_ch_hlimit[i]; ch++){
	if (cm.find(ch) == cm.end()){
	  cm[ch].push_back({m_start,m_end});
	  log->debug("(Veto) {}, added to the bad channel list with range ({}, {}) ",ch, m_start,m_end);
	}
      }
    } 

    // dead channels ...
    std::vector<int> temp_counts;
    if (m_ncount_dead_ch>0){
      temp_counts.resize(m_ncount_dead_ch,0);
      
      for (auto trace : Aux::tagged_traces(in, m_trace_tag)) {
	//	const int tbin = trace->tbin();
	const int chid = trace->channel();
	const auto& charge = trace->charge();
	const int nq = charge.size();

	int saved_group = -1;
	for (int i=0;i!=m_ncount_dead_ch;i++){
	  if (chid >= m_dead_ch_llimit[i] && chid < m_dead_ch_hlimit[i]){
	    saved_group = i;
	    break;
	  }
	}

	if (saved_group >=0){
	  for (int i=0;i!=nq;i++){
	    if (charge.at(i) > m_dead_ch_charge)
	      temp_counts.at(saved_group) ++;
	  }
	}
      }

      //      log->debug("(Dead) group: {},{}", temp_counts.at(0), temp_counts.at(1));

      for(int i=0;i!=m_ncount_dead_ch;i++){
	if (temp_counts.at(i) < m_dead_ch_ncount){
	  for (int ch = m_dead_ch_llimit[i]; ch != m_dead_ch_hlimit[i]; ch++){
	    if (cm.find(ch) == cm.end()){
	      cm[ch].push_back({m_start, m_end});
	      log->debug("(Dead) {}, added to the bad channel list with range ({}, {}) ",ch, m_start, m_end);
	    }
	  }
	}
      }
    }

    // add a new algorithm to align the dead time range ...
    for (auto it = cm.begin(); it!= cm.end(); it++){
      std::map<int, int> record;
      for (auto it1 = it->second.begin(); it1!=it->second.end();it1++){
	//	log->debug("{}", it1->first, it1->second);
	if (m_dead_org > 0){
	  for (int i = m_dead_org-1; i>=0; i--){
	    if (it1->first >= m_dead_org_llimit[i] ){
	      it1->first = m_dead_org_llimit[i];
	      break;
	    }
	  }
	  for (int i= 0;i!= m_dead_org; i++){
	    if (it1->second <= m_dead_org_hlimit[i]){
	      it1->second = m_dead_org_hlimit[i];
	      break;
	    }
	  }
	  record[it1->first] = 0; // start
	  record[it1->second] = 1; // end
 	}	
      }

      // clear existing dead channels ...
      it->second.clear();
      bool flag_start = false;
      int start = 0;
      int end = -1;
      for (auto it1 = record.begin(); it1 != record.end(); it1++){
      	if (flag_start){
      	  if (it1->second == 1){
      	    end = it1->first;
      	    flag_start = false;
      	  }
	  
      	}else{
      	  if (it1->second == 0 ) {
      	    if (end >= start &&  end >=0 ){
      	      it->second.push_back({start,end});
      	      start = 0;
      	      end = -1;
      	    }
      	    start = it1->first;
      	    flag_start = true;
      	  }else if (it1->second == 1){
      	    end = it1->first;
      	  }
      	}
      }
      if (end >= start && end >=0)
      	it->second.push_back({start,end});
      
    }
    
    
    // end for Xin

    // make a copy of all input trace pointers
    ITrace::vector out_traces(*in->traces());

    // Basic frame stays the same.
    auto sfout = new Aux::SimpleFrame(in->ident(), in->time(), out_traces, in->tick(), cmm);

    // ------------ Debug Ewerton 2024-02-29 -----------
    // ChannelMaskMap = typedef std::map<std::string, ChannelMasks>
    auto mycmm = in->masks();  
    for(auto x : mycmm) //cmm[m_cm_tag])
    {
      std::cout << "\n [CMMModifier-Ew] cmm.first=" << x.first << "\n";// " (" << x.second.first << ", " << x.second.second << ")\n";
    }
    // ------------ End debug Ewerton 2024-02-29 -----------

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

    log->debug("call={} output frame: {}", m_count-1, Aux::taginfo(out));
    log->debug("output: {} size: {}", m_cm_tag, out->masks()[m_cm_tag].size());

    return true;
}
