#include "WireCellRoot/MagnifySource.h"
#include "WireCellAux/SimpleTrace.h"
#include "WireCellAux/SimpleFrame.h"

#include "WireCellUtil/NamedFactory.h"

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"

WIRECELL_FACTORY(MagnifySource, WireCell::Root::MagnifySource, WireCell::IFrameSource, WireCell::IConfigurable)

using namespace WireCell;

Root::MagnifySource::MagnifySource()
  : m_calls(0)
{
}

Root::MagnifySource::~MagnifySource() {}

void Root::MagnifySource::configure(const WireCell::Configuration& cfg)
{
    if (cfg["filename"].empty()) {
        THROW(ValueError() << errmsg{"MagnifySource: must supply input \"filename\" configuration item."});
    }
    m_cfg = cfg;
}

WireCell::Configuration Root::MagnifySource::default_configuration() const
{
    Configuration cfg;
    // Give a URL for the input file.
    cfg["filename"] = "";

    // Give a list of frame tags.  These are translated to histogram
    // names h[uvw]_<tag> to be used for input.  Missing histograms
    // raise an exception.  Each frame found will be loaded as tagged
    // traces into the output frame.  The ensemble of traces will be
    // tagged by both <tag> and the per-plane subset as
    // <planeletter>plane.
    cfg["frames"][0] = "raw";

    // A list of pairs mapping a cmm key name to a ttree name.
    cfg["cmmtree"] = Json::arrayValue;

    return cfg;
}

bool Root::MagnifySource::operator()(IFrame::pointer& out)
{
    out = nullptr;
    ++m_calls;
    std::cerr << "MagnifySource: called " << m_calls << " times\n";
    if (m_calls > 2) {
        std::cerr << "MagnifySource: past EOS\n";
        return false;
    }
    if (m_calls > 1) {
        std::cerr << "MagnifySource: EOS\n";
        return true;  // this is to send out EOS
    }

    std::string url = m_cfg["filename"].asString();

    TFile* tfile = TFile::Open(url.c_str());
    if (!tfile || tfile->IsZombie()) {
        THROW(IOError() << errmsg{"MagnifySource: failed to open input file: " + url});
    }

    int frame_ident = 0;
    int nticks = 0;
    double frame_time = 0;
    {
        TTree* trun = (TTree*) tfile->Get("Trun");
        if (!trun) {
            std::cerr << "No tree: Trun in input file \n";
        }
        else {
            // runNo, subRunNo??
            trun->SetBranchAddress("eventNo", &frame_ident);
            trun->SetBranchAddress("total_time_bin", &nticks);
            // trun->SetBranchAddress("time_offset", &frame_time); use this??
            trun->GetEntry(0);
        }
    }
    std::cerr << "MagnifySource: frame ident=" << frame_ident << ", time=0, nticks=" << nticks << "\n";

    WireCell::Waveform::ChannelMaskMap cmm;
    for (auto jcmmtree : m_cfg["cmmtree"]) {
        auto cmmkey = jcmmtree[0].asString();
        auto treename = jcmmtree[1].asString();
        TTree* tree = dynamic_cast<TTree*>(tfile->Get(treename.c_str()));
        if (!tree) {
            std::cerr << "MagnifySource: failed to find tree \"" << treename << "\" in " << tfile->GetName()
                      << std::endl;
            THROW(IOError() << errmsg{"MagnifySource: failed to find tree."});
        }

        std::cerr << "MagnifySource: loading channel mask \"" << cmmkey << "\" from tree \"" << treename << "\"\n";

        int chid = 0, plane = 0, start_time = 0, end_time = 0;
        tree->SetBranchAddress("chid", &chid);
        tree->SetBranchAddress("plane", &plane);
        tree->SetBranchAddress("start_time", &start_time);
        tree->SetBranchAddress("end_time", &end_time);

        const int nentries = tree->GetEntries();
        for (int ientry = 0; ientry < nentries; ++ientry) {
            tree->GetEntry(ientry);
            WireCell::Waveform::BinRange br;
            br.first = start_time;
            br.second = end_time;
            cmm[cmmkey][chid].push_back(br);
        }
    }

    ITrace::vector all_traces;
    std::unordered_map<IFrame::tag_t, IFrame::trace_list_t> tagged_traces;

    {
        for (auto jtag : m_cfg["frames"]) {
            auto frametag = jtag.asString();
            int channel_number = 0;
            for (int iplane = 0; iplane < 3; ++iplane) {
                std::string plane_name = Form("%cplane", 'u' + iplane);
                std::string hist_name = Form("h%c_%s", 'u' + iplane, frametag.c_str());
                std::cerr << "MagnifySource: loading " << hist_name << std::endl;
                TH2* hist = (TH2*) tfile->Get(hist_name.c_str());
                if (!hist) {
                    THROW(IOError() << errmsg{"MagnifySource: failed to find histogram: " + hist_name});
                }

                ITrace::vector plane_traces;

                int nchannels = hist->GetNbinsX();
                int nticks = hist->GetNbinsY();
                double qtot = 0;

                // loop over channels in plane
                for (int chip = 0; chip < nchannels; ++chip) {
                    const int ichbin = chip + 1;

                    ITrace::ChargeSequence charges;
                    for (int itickbin = 1; itickbin <= nticks; ++itickbin) {
                        auto q = hist->GetBinContent(ichbin, itickbin);
                        charges.push_back(q);
                        qtot += q;
                    }
                    const size_t index = all_traces.size();
                    tagged_traces[frametag].push_back(index);
                    all_traces.push_back(std::make_shared<Aux::SimpleTrace>(channel_number, 0, charges));

                    ++channel_number;
                }
                std::cerr << "MagnifySource: plane " << iplane << ": " << nchannels << " X " << nticks
                          << " qtot=" << qtot << std::endl;
            }  // plane
        }
    }

    auto sframe = new Aux::SimpleFrame(frame_ident, frame_time, all_traces, 0.5 * units::microsecond, cmm);
    for (auto const& it : tagged_traces) {
        sframe->tag_traces(it.first, it.second);
        std::cerr << "MagnifySource: tag " << it.second.size() << " traces as: \"" << it.first << "\"\n";
    }

    out = IFrame::pointer(sframe);

    tfile->Close();
    delete tfile;

    return true;
}
