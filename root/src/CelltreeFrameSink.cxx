#include "WireCellRoot/CelltreeFrameSink.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/BoundingBox.h"
#include "WireCellUtil/Waveform.h"

#include "TFile.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TTree.h"
#include "TH1F.h"
#include "TClonesArray.h"

#include <iostream>
#include <algorithm>
#include <unordered_map>

WIRECELL_FACTORY(CelltreeFrameSink, WireCell::Root::CelltreeFrameSink, WireCell::IFrameFilter, WireCell::IConfigurable)

using namespace std;
using namespace WireCell;

Root::CelltreeFrameSink::CelltreeFrameSink()
  : m_nrebin(1)
{
}

Root::CelltreeFrameSink::~CelltreeFrameSink() {}

WireCell::Configuration Root::CelltreeFrameSink::default_configuration() const
{
    Configuration cfg;

    cfg["anode"] = "AnodePlane";
    cfg["output_filename"] = "";
    cfg["frames"] = Json::arrayValue;
    cfg["cmmtree"] = Json::arrayValue;
    cfg["root_file_mode"] = "RECREATE";
    cfg["nsamples"] = 0;
    cfg["nrebin"] = 1;

    // List tagged traces from which to save the "trace summary"
    // vector into a 1D vector which will be named after the tag.
    // See "summary_operator".
    cfg["summaries"] = Json::arrayValue;

    // An object mapping tags to operators for aggregating trace
    // summary values on the same channel.  Operator may be "sum" to
    // add up all values on the same channel or "set" to assign values
    // to the channel bin (last one wins).  If a tag is not found, the
    // default operator is "sum".
    cfg["summary_operator"] = Json::objectValue;

    return cfg;
}

void Root::CelltreeFrameSink::configure(const WireCell::Configuration &cfg)
{
    std::string fn;

    fn = cfg["output_filename"].asString();
    if (fn.empty()) {
        THROW(ValueError() << errmsg{"Must provide output filename to CelltreeFrameSink"});
    }

    auto anode_tn = get<std::string>(cfg, "anode", "AnodePlane");
    m_anode = Factory::lookup_tn<IAnodePlane>(anode_tn);
    if (!m_anode) {
        cerr << "Root::CelltreeFrameSink: failed to get anode: \"" << anode_tn << "\"\n";
        return;
    }

    m_nsamples = get<int>(cfg, "nsamples", 0);
    if (m_nsamples == 0) {
        THROW(ValueError() << errmsg{"nsamples has to be configured"});
    }

    m_cfg = cfg;

    m_nrebin = get<int>(cfg, "nrebin", m_nrebin);
}

typedef std::unordered_set<std::string> string_set_tc;
string_set_tc cgetset(const WireCell::Configuration &cfg)
{
    string_set_tc ret;
    for (auto jone : cfg) {
        ret.insert(jone.asString());
    }
    return ret;
}

ITrace::vector cget_tagged_traces(IFrame::pointer frame, IFrame::tag_t tag)
{
    ITrace::vector ret;
    auto const &all_traces = frame->traces();
    for (size_t index : frame->tagged_traces(tag)) {
        ret.push_back(all_traces->at(index));
    }
    if (!ret.empty()) {
        return ret;
    }
    auto ftags = frame->frame_tags();
    if (std::find(ftags.begin(), ftags.end(), tag) == ftags.end()) {
        return ret;
    }
    return *all_traces;  // must make copy
}

bool Root::CelltreeFrameSink::operator()(const IFrame::pointer &frame, IFrame::pointer &out_frame)
{
    out_frame = frame;
    if (!frame) {
        std::cerr << "CelltreeFrameSink: EOS\n";
        return true;
    }

    const std::string ofname = m_cfg["output_filename"].asString();
    const std::string mode = m_cfg["root_file_mode"].asString();
    std::cerr << "CelltreeFrameSink: opening for output: " << ofname << " with \"" << mode << "\"\n";
    TFile *output_tf = TFile::Open(ofname.c_str(), mode.c_str());
    if (!output_tf->GetDirectory("Event")) output_tf->mkdir("Event");
    output_tf->cd("Event");

    // [HY] Simulation output celltree
    // Sim Tree
    TTree *Sim;
    Sim = (TTree *) output_tf->Get("Event/Sim");
    if (!Sim) {
        Sim = new TTree("Sim", "Wire-cell toolkit simulation output");
        Int_t runNo = 0;
        Int_t subRunNo = 0;
        Int_t eventNo = 0;
        Sim->Branch("runNo", &runNo, "runNo/I");
        Sim->Branch("subRunNo", &subRunNo, "subRunNo/I");
        Sim->Branch("eventNo", &eventNo, "eventNo/I");
        Sim->Fill();
    }

    // trace frames
    for (auto tag : cgetset(m_cfg["frames"])) {
        // string tag = "gauss";
        ITrace::vector traces = cget_tagged_traces(frame, tag);
        if (traces.empty()) {
            std::cerr << "CelltreeFrameSink: no tagged traces for\"" << tag << "\"\n";
            continue;
        }

        std::cerr << "CelltreeFrameSink: tag: \"" << tag << "\" with " << traces.size() << " traces\n";

        std::vector<int> *raw_channelId = new std::vector<int>;
        std::string channelIdname;
        if (!tag.compare("gauss")) channelIdname = "calibGaussian_channelId";
        if (!tag.compare("wiener")) channelIdname = "calibWiener_channelId";
        if (!tag.compare("orig")) channelIdname = "raw_channelId";
        // const std::string channelIdname = Form("%s_channelId", tag.c_str());
        TBranch *bchannelId = Sim->Branch(channelIdname.c_str(), &raw_channelId);

        TClonesArray *sim_wf = new TClonesArray("TH1F");
        TH1::AddDirectory(kFALSE);
        std::string wfname;
        if (!tag.compare("gauss")) wfname = "calibGaussian_wf";
        if (!tag.compare("wiener")) wfname = "calibWiener_wf";
        if (!tag.compare("orig")) wfname = "raw_wf";
        // const std::string wfname = Form("%s_wf", tag.c_str());
        TBranch *bwf = Sim->Branch(wfname.c_str(), &sim_wf, 256000, 0);

        int nsamples = m_nsamples;
        // std::cout<<"nsamples: "<<m_nsamples<<std::endl;

        int sim_wf_ind = 0;

        for (auto trace : traces) {
            int ch = trace->channel();
            // std::cout<<"channel number: "<<ch<<std::endl;
            raw_channelId->push_back(ch);
            // fill raw_wf

            int temp_nbins = nsamples / m_nrebin;

            TH1F *htemp = new ((*sim_wf)[sim_wf_ind]) TH1F("", "", temp_nbins, 0, nsamples);
            // TH1I *htemp = new ( (*sim_wf)[sim_wf_ind] ) TH1I("", "",  nsamples, 0,  nsamples);
            auto const &wave = trace->charge();
            const int nbins = wave.size();
            // std::cout<<"waveform size: "<<nbins<<std::endl;
            const int tmin = trace->tbin();
            // std::cout<<"tmin: "<<tmin<<std::endl;
            for (Int_t i = 0; i < nbins; i++) {
                if (tmin + i + 1 <= nsamples) {
                    int ibin = (tmin + i) / m_nrebin;
                    htemp->SetBinContent(ibin + 1, wave[i] + htemp->GetBinContent(ibin + 1));
                }
            }
            sim_wf_ind++;
        }  // traces

        cout << "channelId size: " << raw_channelId->size() << "\n";
        bchannelId->Fill();
        bwf->Fill();
        delete raw_channelId;
        delete sim_wf;
    }  // frames

    // trace summaries
    // currently only one option "threshold"
    for (auto tag : cgetset(m_cfg["summaries"])) {
        auto traces = cget_tagged_traces(frame, tag);
        if (traces.empty()) {
            std::cerr << "CelltreeFrameSink: warning: no traces tagged with \"" << tag << "\", skipping summary\n";
            continue;
        }
        auto const &summary = frame->trace_summary(tag);
        if (summary.empty()) {
            std::cerr << "CelltreeFrameSink: warning: empty summary tagged with \"" << tag << "\", skipping summary\n";
            continue;
        }

        // default: "sum"
        // to save thresholds: "set"
        std::string oper = get<std::string>(m_cfg["summary_operator"], tag, "sum");

        std::cerr << "CelltreeFrameSink: summary tag: \"" << tag << "\" with " << traces.size() << " traces\n";

        // vector<double> channelThreshold index <--> channelId
        std::vector<double> *channelThreshold = new std::vector<double>;
        TBranch *bthreshold = Sim->Branch("channelThreshold", &channelThreshold);

        const int ntot = traces.size();
        // Find max channel ID to size the vector correctly
        int max_chid = 0;
        for (int ind = 0; ind < ntot; ++ind) {
            int chid = traces[ind]->channel();
            if (chid > max_chid) max_chid = chid;
        }
        channelThreshold->resize(max_chid + 1, 0);
        for (int ind = 0; ind < ntot; ++ind) {
            const int chid = traces[ind]->channel();
            const double val = summary[ind];
            if (oper == "set") {
                channelThreshold->at(chid) = val;
            }
            else {
                channelThreshold->at(chid) += val;
            }
        }

        // debug
        // for (int i=0; i<channelThreshold->size(); i++)
        //{
        //    cout<<" chid: "<< i << " threshold: "<< channelThreshold->at(i) << "\n";
        //}

        bthreshold->Fill();
        delete channelThreshold;
    }

    // bad channels
    // vector<int> badChannel
    // vector<int> badBegin
    // vector<int> badEnd
    for (auto tag : cgetset(m_cfg["cmmtree"])) {
        Waveform::ChannelMaskMap input_cmm = frame->masks();
        for (auto const &it : input_cmm) {
            auto cmmkey = it.first;

            if (tag.compare(cmmkey) != 0) continue;

            std::cerr << "CelltreeFrameSink: saving channel mask \"" << cmmkey << "\"\n";

            std::vector<int> *Channel = new std::vector<int>;
            std::vector<int> *Begin = new std::vector<int>;
            std::vector<int> *End = new std::vector<int>;
            const std::string Channelname = Form("%sChannel", tag.c_str());
            const std::string Beginname = Form("%sBegin", tag.c_str());
            const std::string Endname = Form("%sEnd", tag.c_str());
            TBranch *bch = Sim->Branch(Channelname.c_str(), &Channel);
            TBranch *bb = Sim->Branch(Beginname.c_str(), &Begin);
            TBranch *be = Sim->Branch(Endname.c_str(), &End);

            for (auto const &chmask : it.second) {
                auto mask = chmask.second;
                if (mask.size() != 1) {
                    std::cerr << "CelltreeFrameSink: Warning: channel mask: " << chmask.first
                              << " has >1 dead period [begin, end] \n";
                    continue;
                }
                Channel->push_back(chmask.first);
                for (size_t ind = 0; ind < mask.size(); ++ind) {
                    Begin->push_back(mask[ind].first);
                    End->push_back(mask[ind].second);
                }
            }

            bch->Fill();
            bb->Fill();
            be->Fill();
            delete Channel;
            delete Begin;
            delete End;
        }
    }

    std::cerr << "CelltreeFrameSink: closing output file " << ofname << std::endl;
    auto count = output_tf->Write();
    std::cerr << "\twrote " << count << " bytes." << std::endl;
    output_tf->Close();
    delete output_tf;
    output_tf = nullptr;
    return true;
}
