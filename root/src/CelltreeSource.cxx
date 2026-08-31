#include "WireCellRoot/CelltreeSource.h"
#include "WireCellAux/SimpleTrace.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/FrameTools.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/String.h"

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TH1F.h"

WIRECELL_FACTORY(CelltreeSource, WireCell::Root::CelltreeSource,
                 WireCell::INamed,
                 WireCell::IFrameSource, WireCell::IConfigurable)

using namespace WireCell;

Root::CelltreeSource::CelltreeSource()
    : Aux::Logger("CelltreeSource", "root")
    , m_calls(0)
{
}

Root::CelltreeSource::~CelltreeSource() {}

void Root::CelltreeSource::configure(const WireCell::Configuration& cfg)
{
    if (cfg["filename"].empty()) {
        THROW(ValueError() << errmsg{"CelltreeSource: must supply input \"filename\" configuration item."});
    }
    m_cfg = cfg;
}

WireCell::Configuration Root::CelltreeSource::default_configuration() const
{
    Configuration cfg;
    // Give a URL for the input file.
    cfg["filename"] = "";

    // which entry of the /Event/Sim tree to process
    cfg["entry"] = 0;

    // Give a list of frame/tree tags.

    // input branch base name
    // raw, calibGaussian, calibWiener
    cfg["in_branch_base_names"][0] = "raw";

    // output frame tags, should be corresponding to each input
    // branch name, raw: orig, calibGaussian: gauss, calibWiener, wiener
    cfg["out_trace_tags"][0] = "orig";

    // FIXME: is "channelThreshold" for both wiener and gauss?
    cfg["in_branch_thresholds"][0] = "";

    cfg["time_scale"] = 4;

    return cfg;
}

bool Root::CelltreeSource::read_traces(ITrace::vector& all_traces,
                                       std::unordered_map<IFrame::tag_t, IFrame::trace_list_t>& tagged_traces,
                                       std::unordered_map<IFrame::tag_t, IFrame::trace_summary_t>& tagged_threshold,
                                       TTree* tree, const std::string& br_name,
                                       const std::string& br_name_threshold, const std::string& frametag,
                                       const unsigned int entry, const int time_scale) const
{
    bool trace_has_threshold = false;
    if (tree->GetBranch(br_name_threshold.c_str()) != nullptr) trace_has_threshold = true;

    std::vector<int>* channelid = new std::vector<int>;
    std::vector<double>* channel_threshold = new std::vector<double>;
    TClonesArray* esignal = new TClonesArray;

    auto br_name_channelId = String::format("%s_channelId", br_name);
    auto br_name_wf = String::format("%s_wf",br_name);

    // tree->SetBranchStatus(br_name_channelId.c_str(), 1);
    tree->SetBranchAddress(br_name_channelId.c_str(), &channelid);
    // tree->SetBranchStatus(br_name_wf.c_str(), 1);
    tree->SetBranchAddress(br_name_wf.c_str(), &esignal);
    if (trace_has_threshold) tree->SetBranchAddress(br_name_threshold.c_str(), &channel_threshold);

    int size = tree->GetEntry(entry);
    if (!(size>0)) {
        THROW(RuntimeError() << errmsg{"Failed loading entry!"});
    }

    int nchannels = channelid->size();

    int channel_number = 0;
    

    log->debug("loading frame \"{}\", nchannels={}, thresholds={}, ticks={}",
               frametag, nchannels, channel_threshold->size(), esignal->GetSize());
    // std::cout << "CelltreeSource: loading " << frametag << " " << nchannels << " channels \n";
    // std::cout << "CelltreeSource: channel_threshold->GetSize(): " << channel_threshold->size() << "\n";
    // std::cout << "CelltreeSource: esignal->GetSize(): " << esignal->GetSize() << "\n";
    if (nchannels != esignal->GetSize() || (trace_has_threshold && nchannels != (int) channel_threshold->size())) {
        THROW(RuntimeError() << errmsg{"nchannels != channel_threshold->size() || nchannels != esignal->GetSize()"});
    }

    // fill waveform ...
    TH1* signal=nullptr;

    for (int ind = 0; ind < nchannels; ++ind) {
        signal = dynamic_cast<TH1*>(esignal->At(ind));
        if (!signal) continue;
        channel_number = channelid->at(ind);

        if(ind==0) {
            log->debug("branch {} channel number: {}, signal integral: {}",
                       br_name, channel_number, signal->Integral());
            // std::cout << "yuhw: " << br_name << " channel_number: " << channel_number << " signal->Integral(): " << signal->Integral() << std::endl;
        }

        ITrace::ChargeSequence charges;
        int nticks = signal->GetNbinsX();
        for (int itickbin = 0; itickbin != nticks; itickbin++) {
            for (int i=0; i<time_scale; ++i) {
                charges.push_back(signal->GetBinContent(itickbin + 1)/time_scale);
            }
        }
        const size_t index = all_traces.size();
        tagged_traces[frametag].push_back(index);
        all_traces.push_back(std::make_shared<Aux::SimpleTrace>(channel_number, 0, charges));
        if(trace_has_threshold) {
            // std::cout << String::format("channel_number %d threshold %f\n",channel_number, channel_threshold->at(ind)/time_scale);
            tagged_threshold[frametag].push_back(channel_threshold->at(ind)/time_scale);}
    }

    delete channelid;
    delete channel_threshold;
    delete esignal;
    return true;
}

bool Root::CelltreeSource::read_cmm(WireCell::Waveform::ChannelMaskMap &cmm, TTree* tree, const unsigned int entry) const
{
    std::vector<int>* badChannel = new std::vector<int>;
    std::vector<int>* badBegin = new std::vector<int>;
    std::vector<int>* badEnd = new std::vector<int>;

    tree->SetBranchAddress("badChannel", &badChannel);
    tree->SetBranchAddress("badBegin", &badBegin);
    tree->SetBranchAddress("badEnd", &badEnd);

    int size = tree->GetEntry(entry);
    if (!(size>0)) {
        THROW(RuntimeError() << errmsg{"Failed loading entry!"});
    }

    int nchannels = badChannel->size();
    if (nchannels != (int)badBegin->size() || nchannels != (int)badEnd->size()) {
        THROW(RuntimeError() << errmsg{"Length doesn't match"});
    }

    for (int ind = 0; ind < nchannels; ++ind) {
        auto chid = badChannel->at(ind);
        WireCell::Waveform::BinRange br;
        br.first = badBegin->at(ind);
        br.second = badEnd->at(ind);
        cmm["bad"][chid].push_back(br);
    }

    delete badChannel;
    delete badBegin;
    delete badEnd;
    return true;
}

bool Root::CelltreeSource::operator()(IFrame::pointer& out)
{
    out = nullptr;
    ++m_calls;

    if (m_calls > 2) {
        log->debug("past EOS at call={}", m_calls-1);
        return false;
    }
    if (m_calls > 1) {
        log->debug("EOS at call={}", m_calls-1);
        return true;  // this is to send out EOS
    }

    std::string url = m_cfg["filename"].asString();
    TFile* tfile = TFile::Open(url.c_str());
    TTree* tree = (TTree*) tfile->Get("/Event/Sim");
    if (!tree) {
        THROW(IOError() << errmsg{"No tree: /Event/Sim in input file"});
    }

    // tree->SetBranchStatus("*", 0);

    int run_no, subrun_no, event_no;
    // tree->SetBranchStatus("eventNo", 1);
    tree->SetBranchAddress("eventNo", &event_no);
    // tree->SetBranchStatus("runNo", 1);
    tree->SetBranchAddress("runNo", &run_no);
    // tree->SetBranchStatus("subRunNo", 1);
    tree->SetBranchAddress("subRunNo", &subrun_no);

    unsigned int ent = m_cfg["entry"].asUInt();
    unsigned int entries = tree->GetEntries();
    if (ent >= entries) {
        THROW(ValueError() << errmsg{"entry number out of range"});
    }
    int siz = tree->GetEntry(ent);

    // need run number and subrunnumber
    int frame_ident = event_no;
    double frame_time = 0;
    int time_scale = m_cfg["time_scale"].asInt();

    // some output using eventNo, runNo, subRunNO, ...
    log->debug("frame_ident:{} size:{} runNo:{}, subrunNo:{}, eventNo:{}, time_scale:{}, ent:{}",
               frame_ident, siz, run_no, subrun_no, event_no, time_scale, ent);

    ITrace::vector all_traces;
    std::unordered_map<IFrame::tag_t, IFrame::trace_list_t> tagged_traces;
    std::unordered_map<IFrame::tag_t, IFrame::trace_summary_t> tagged_threshold;

    if (m_cfg["in_branch_base_names"].size() != m_cfg["out_trace_tags"].size()) {
        THROW(ValueError() << errmsg{"#in_branch_base_names != #out_trace_tags"});
    }

    for (Json::ArrayIndex i = 0; i < m_cfg["in_branch_base_names"].size(); ++i) {
        auto branch_base_name = m_cfg["in_branch_base_names"][i].asString();
        auto frame_tag = m_cfg["out_trace_tags"][i].asString();
        auto branch_threshold = m_cfg["in_branch_thresholds"][i].asString();
        read_traces(all_traces, tagged_traces, tagged_threshold, tree, branch_base_name, branch_threshold, frame_tag, ent, time_scale);
    }

    WireCell::Waveform::ChannelMaskMap cmm;
    read_cmm(cmm, tree, ent);

    tfile->Close();
    delete tfile;

    auto sframe = new Aux::SimpleFrame(frame_ident, frame_time, all_traces, 0.5 * units::microsecond, cmm);
    for (auto const& it : tagged_traces) {
        sframe->tag_traces(it.first, it.second, tagged_threshold[it.first]);
    }

    out = IFrame::pointer(sframe);
    log->debug("output: {}", Aux::taginfo(out));

    return true;
}
