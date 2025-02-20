#include "WireCellImg/ChargeErrorFrameEstimator.h"
#include "WireCellIface/IWaveformMap.h"
#include "WireCellAux/FrameTools.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleTrace.h"

#include "WireCellUtil/Units.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"

#include <sstream>

WIRECELL_FACTORY(ChargeErrorFrameEstimator, WireCell::Img::ChargeErrorFrameEstimator, WireCell::IFrameFilter, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Img;

using WireCell::Aux::SimpleTrace;
using WireCell::Aux::SimpleFrame;

// One place for the defaults
static Binning default_bins(1001, 0, (1+1000)*0.5*units::us);

ChargeErrorFrameEstimator::ChargeErrorFrameEstimator()
    : Aux::Logger("ChargeErrorFrameEstimator", "img")
{
}

ChargeErrorFrameEstimator::~ChargeErrorFrameEstimator() {}

WireCell::Configuration ChargeErrorFrameEstimator::default_configuration() const
{
    Configuration cfg;

    // The "service" providing us ROI error "waveforms".
    cfg["errors"] = "WaveformMap";

    // This maps plane index to name used in the waveform map.
    // User config must make sure it is correct.
    cfg["planes"][0] = "u";
    cfg["planes"][1] = "v";
    cfg["planes"][2] = "w";
    cfg["tick"] = default_bins.binsize();
    cfg["nbins"] = default_bins.nbins();

    // Frame trace tags
    cfg["intag"] = m_intag;
    cfg["outtag"] = m_outtag;

    // rebin choice, need to match with the input file
    cfg["rebin"] = m_rebin;

    // additional fudge factors to update error calculations for each plane
    cfg["fudge_factors"][0] = m_fudge_factors[0];
    cfg["fudge_factors"][1] = m_fudge_factors[1];
    cfg["fudge_factors"][2] = m_fudge_factors[2];

    // additional time limits in ticks
    cfg["time_limits"][0] = m_time_limits.first;
    cfg["time_limits"][1] = m_time_limits.second;
    
    cfg["anode"] = "AnodePlane";

    // add any other parameters you need.


    return cfg;
}

void ChargeErrorFrameEstimator::configure(const WireCell::Configuration& cfg)
{
    auto anode_tn = cfg["anode"].asString();
    m_anode = Factory::find_tn<IAnodePlane>(anode_tn);

    m_intag = get<std::string>(cfg, "intag", m_intag);
    m_outtag = get<std::string>(cfg, "outtag", m_outtag);

    m_rebin = get<int>(cfg,"rebin",m_rebin);
    if (cfg.isMember("fudge_factors")) {
        m_fudge_factors.clear();
        for (auto jplane : cfg["fudge_factors"]) {
	  m_fudge_factors.push_back(jplane.asFloat());
        }
    }
    
    if (cfg.isMember("time_limits")) {
      if (cfg["time_limits"].size()!=2)
	THROW(ValueError() << errmsg{"Formot of TIME LIMITS is Wrong!"});
      
      m_time_limits.first = cfg["time_limits"][0].asFloat();
      m_time_limits.second = cfg["time_limits"][1].asFloat();
    }
    
    //    log->debug("Rebin = {}, Fudge Factors = {}, Time Limits = {}",m_rebin, m_fudge_factors[0], m_time_limits.first);
    
    // Apply user binning
    int nbins = get(cfg, "nbins", default_bins.nbins());
    double tick = get(cfg, "tick", default_bins.binsize());
    m_bins.set(nbins + 1, 0, (1 + nbins) * tick);

    // Load ROI error data.
    auto tn = cfg["errors"].asString();
    m_errs.resize((size_t)cfg["planes"].size());
    auto wfm = Factory::find_tn<IWaveformMap>(tn);
    for (int iplane=0; iplane< (int)cfg["planes"].size(); ++iplane) {
        auto plane = cfg["planes"][iplane].asString();
        auto wf = wfm->waveform_lookup(plane);
        m_errs[iplane] = wf->waveform_samples(m_bins);

	// judge if the rebin is consistent with the error data
	if (m_rebin != std::round(wf->waveform_period()/tick))
	  log->error("Rebin {} does not match with Error Data period length {} with tick length {}", m_rebin, wf->waveform_period(), tick);
	// warn, inform, debug
    }
}

ITrace::vector ChargeErrorFrameEstimator::error_traces(const ITrace::vector& intraces) const {
    ITrace::vector ret;
    // Loop over tagged traces from input frame.
    // ---------- begin of the stub for Xin ----------
    for (auto itrace : intraces) {
        auto ich = m_anode->channel(itrace->channel());
        auto plane_index = ich->planeid().index();

	double fudge_factor = m_fudge_factors[plane_index];
	
        const auto& errs = m_errs[plane_index];
        auto& charge = itrace->charge();
        size_t nq = charge.size();

	
	//log->debug(nq);
	
        size_t roistart = 0;
        bool inroi = false;
        std::vector<float> error_charge(nq, 0);
        for (size_t iq = 0; iq < nq; iq++) {
            auto& q = charge[iq];
            // start an ROI
            if (!inroi && q != 0) {
                inroi = true;
                roistart = iq;
                continue;
            }
            // end an ROI
            if (inroi && (q == 0 )) {
	      int roilen = iq - roistart;
	      if (roilen < m_time_limits.first) roilen = m_time_limits.first;
	      else if (roilen > m_time_limits.second) roilen = m_time_limits.second;
          if (errs[roilen] <= 0) {
            log->warn("Error {} <= 0 for plane {} channel {} roilen {}",
                      errs[roilen], plane_index, ich->ident(), roilen);
            /// TODO: what to do here?
          }
	      std::generate(error_charge.begin() + roistart, error_charge.begin() + iq,
			    [=]() mutable { return errs[roilen] * fudge_factor / sqrt(m_rebin); });
	      inroi = false;
	      continue;
            }else if (inroi && iq == nq - 1){
	      int roilen = iq + 1 - roistart;
	      if (roilen < m_time_limits.first) roilen = m_time_limits.first;
	      else if (roilen > m_time_limits.second) roilen = m_time_limits.second;
          if (errs[roilen] <= 0) {
            log->warn("Error {} <= 0 for plane {} channel {} roilen {}",
                      errs[roilen], plane_index, ich->ident(), roilen);
            /// TODO: what to do here?
          }
	      std::generate(error_charge.begin() + roistart, error_charge.begin() + iq + 1,
			    [=]() mutable { return errs[roilen] * fudge_factor / sqrt(m_rebin); });
	      inroi = false;
	      continue;
	    }
        }

	// if (plane_index==1 && ich->ident() == 1758+2400){
	//   for (int i=1850; i!= 1950;i++){
	//     log->debug("{}, {}, {}",i, charge[i], error_charge[i]);
	//   }
	// }

        ret.push_back(
            std::make_shared<SimpleTrace>(ich->ident(), itrace->tbin(), error_charge));

    }
    // ---------- end of the stub for Xin ---------- 
    return ret;
}

bool ChargeErrorFrameEstimator::operator()(const input_pointer& in, output_pointer& out)
{
    out = nullptr;
    if (!in) {
        log->debug("see EOS");
        return true;  // eos
    }

    // calculate error traces
    ITrace::vector in_traces = Aux::tagged_traces(in, m_intag);
    ITrace::vector err_traces = error_traces(in_traces);

    // make a copy of all input trace pointers
    ITrace::vector out_traces(*in->traces());

    // list used for the new error trace tagging
    IFrame::trace_list_t out_trace_indices(err_traces.size());
    std::generate(out_trace_indices.begin(), out_trace_indices.end(), [n = out_traces.size()] () mutable {return n++;});

    // merge error traces to the output copy
    out_traces.insert(out_traces.end(), err_traces.begin(), err_traces.end());

    // Basic frame stays the same.
    auto sfout = new SimpleFrame(in->ident(), in->time(), out_traces, in->tick(), in->masks());

    // tag err traces
    sfout->tag_traces(m_outtag, out_trace_indices);

    // passing through other parts of the original frame
    for (auto ftag : in->frame_tags()) {
        sfout->tag_frame(ftag);
    }

    for (auto inttag : in->trace_tags()) {
        const auto& traces = in->tagged_traces(inttag);
        const auto& summary = in->trace_summary(inttag);
        if (traces.empty()) {
            log->warn("no traces for requested tag \"{}\"", inttag);
            continue;
        }
        sfout->tag_traces(inttag, traces, summary);
    };

    out = IFrame::pointer(sfout);

    log->debug("input : {}", Aux::taginfo(in));
    log->debug("output: {}", Aux::taginfo(out));

    return true;
}
