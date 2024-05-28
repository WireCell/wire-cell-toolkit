#include "WireCellRoot/MagnifySink.h"
#include "WireCellIface/ITrace.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TFile.h"
#include "TTree.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellAux/FrameTools.h"

#include <vector>
#include <string>

WIRECELL_FACTORY(MagnifySink, WireCell::Root::MagnifySink, WireCell::IFrameFilter, WireCell::IConfigurable)

using namespace WireCell;

Root::MagnifySink::MagnifySink()
  : m_nrebin(1)
  , log(Log::logger("magnify"))
{
}

Root::MagnifySink::~MagnifySink() {}

void Root::MagnifySink::configure(const WireCell::Configuration& cfg)
{
    std::string fn;

    fn = cfg["input_filename"].asString();
    if (fn.empty() and !cfg["shunt"].empty()) {
        log->error("MagnifySink: asked to shunt but not given input file name");
        THROW(ValueError() << errmsg{"MagnifySink: must provide input filename to shunt objects to output"});
    }

    fn = cfg["output_filename"].asString();
    if (fn.empty()) {
        THROW(ValueError() << errmsg{"Must provide output filename to MagnifySink"});
    }

    auto anode_tn = get<std::string>(cfg, "anode", "AnodePlane");
    m_anode = Factory::find_tn<IAnodePlane>(anode_tn);

    m_cfg = cfg;

    m_nrebin = get<int>(cfg, "nrebin", m_nrebin);

    create_file();
}

WireCell::Configuration Root::MagnifySink::default_configuration() const
{
    Configuration cfg;

    cfg["anode"] = "AnodePlane";

    // fixme: this TOTALLY violates the design of wire cell DFP.
    cfg["input_filename"] = "";

    // List of TObjects to copy from input file to output file.
    cfg["shunt"] = Json::arrayValue;

    // Name of ROOT file to write.
    cfg["output_filename"] = "";

    // A list of trace tags defining which waveforms are saved to Magnify histograms.
    cfg["frames"] = Json::arrayValue;

    // If no tags for traces, i.e. trace_has_tag=false in a frame,
    // set desired tag to ""
    // as FrameTool::tagged_traces(frame, "") calls untagged_traces(frame).
    cfg["trace_has_tag"] = true;

    // A list of pairs mapping a cmm key name to a ttree name.
    cfg["cmmtree"] = Json::arrayValue;

    // The ROOT file mode with which to open the file.  Use "RECREATE"
    // to overrite an existing file.  This might be useful for the
    // first MagnifySink in a chain.  Use "UPDATE" for subsequent
    // sinks that add to the file.
    cfg["root_file_mode"] = "RECREATE";

    // If runinfo is given it should be a JSON object and its values
    // will be copied into the Trun tree.  If instead it is null AND
    // an input file is given AND it contains a Trun tree, it will be
    // copied to output.
    cfg["runinfo"] = Json::nullValue;

    cfg["nrebin"] = 1;

    // List tagged traces from which to save the "trace summary"
    // vector into a 1D histogram which will be named after the tag.
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

typedef std::unordered_set<std::string> string_set_t;
string_set_t getset(const WireCell::Configuration& cfg)
{
    string_set_t ret;
    for (auto jone : cfg) {
        ret.insert(jone.asString());
    }
    return ret;
}

std::vector<WireCell::Binning> collate_byplane(const ITrace::vector& traces, const IAnodePlane::pointer anode,
                                               ITrace::vector byplane[])
{
    std::vector<int> uvwt[4];
    for (auto trace : traces) {
        const int chid = trace->channel();
        auto wpid = anode->resolve(chid);
        const int iplane = wpid.index();
        if (iplane < 0 || iplane >= 3) {
            THROW(RuntimeError() << errmsg{"Illegal wpid"});
        }
        uvwt[iplane].push_back(chid);
        byplane[iplane].push_back(trace);
        uvwt[3].push_back(trace->tbin());
        uvwt[3].push_back(trace->tbin() + trace->charge().size());
    }

    std::vector<Binning> binnings(4);
    for (int ind = 0; ind < 4; ++ind) {
        auto const& one = uvwt[ind];
        // std::cerr << "[wgu] get ind: " << ind << " size: " << one.size() << std::endl;
        if (one.empty()) {
            // THROW(ValueError() << errmsg{"MagnifySink: bogus bounds"});
            std::cerr << "[wgu] plane: " << ind << " has not traces. " << std::endl;
        }
        else {
            auto mme = std::minmax_element(one.begin(), one.end());
            const int vmin = *mme.first;
            const int vmax = *mme.second;
            if (ind == 3) {
                const int n = vmax - vmin;
                // binnings.push_back(Binning(n, vmin, vmax););
                binnings.at(ind) = Binning(n, vmin, vmax);
            }
            else {
                // Channel-centered binning
                const double diff = vmax - vmin;
                // binnings.push_back(Binning(diff+1, vmin-0.5, vmax+0.5));
                binnings.at(ind) = Binning(diff + 1, vmin - 0.5, vmax + 0.5);
            }
        }
    }
    return binnings;
}

void Root::MagnifySink::create_file()
{
    const std::string ofname = m_cfg["output_filename"].asString();
    const std::string mode = "RECREATE";
    TFile* output_tf = TFile::Open(ofname.c_str(), mode.c_str());
    output_tf->Close("R");
    delete output_tf;
    output_tf = nullptr;
}

void Root::MagnifySink::do_shunt(TFile* output_tf)
{
    std::stringstream ss;
    auto truncfg = m_cfg["runinfo"];
    if (!truncfg.empty()) {
        TTree* rtree = new TTree("Trun", "Trun");
        rtree->SetDirectory(output_tf);
        ss << "MagnifySink: making Tree RunInfo:\n";

        std::vector<int> ints;
        // issue to be fixed:
        // the first element seems to be misconnected
        // to a wrong address in the rtree->Branch when rtree->Fill
        // So kind of initilized to the vector could solve this issues
        ints.push_back(0);
        std::vector<float> floats;
        floats.push_back(0.0);

        int frame_number = 0;
        int celltree_input = 0;
        for (auto name : truncfg.getMemberNames()) {
            auto jval = truncfg[name];
            ss << "\t" << name << " = " << jval << "\n";
            if (name == "eventNo") frame_number = std::stoi(jval.asString());
            if (name == "Celltree") {
                celltree_input = jval.asInt();
                continue;
            }
            if (jval.isInt()) {
                ints.push_back(jval.asInt());
                rtree->Branch(name.c_str(), &ints.back(), (name + "/I").c_str());
                continue;
            }
            if (jval.isDouble()) {
                floats.push_back(jval.asFloat());
                rtree->Branch(name.c_str(), &floats.back(), (name + "/F").c_str());
                continue;
            }
            if (jval.isString()) {
                ints.push_back(std::stoi(jval.asString()));
                rtree->Branch(name.c_str(), &ints.back(), (name + "/I").c_str());
                continue;
            }
            log->warn("MagnifySink: warning: got unknown type for run info entry: \"{}\" = {}", name, jval);
        }
        if (ss.str().size() > 0) {
            log->debug(ss.str());
        }

        if (celltree_input) {
            // runNo and subRunNo, perhaps other info in the future
            TFile* input_runinfo = TFile::Open((m_cfg["input_filename"].asString()).c_str());
            TTree* run = (TTree*) input_runinfo->Get("/Event/Sim");
            if (!run) {
                log->warn("MagnifySink: runinfo: no tree: /Event/Sim in input file");
            }
            else {
                run->SetBranchStatus("*", 0);

                int run_no, subrun_no, event_no;
                run->SetBranchStatus("eventNo", 1);
                run->SetBranchAddress("eventNo", &event_no);
                run->SetBranchStatus("runNo", 1);
                run->SetBranchAddress("runNo", &run_no);
                rtree->Branch("runNo", &run_no, "runNo/I");
                run->SetBranchStatus("subRunNo", 1);
                run->SetBranchAddress("subRunNo", &subrun_no);
                rtree->Branch("subRunNo", &subrun_no, "subRunNo/I");

                unsigned int entries = run->GetEntries();
                bool legalevt = false;
                for (unsigned int ent = 0; ent < entries; ent++) {
                    int siz = run->GetEntry(ent);
                    if (siz > 0 && event_no == frame_number) {
                        legalevt = true;
                        break;
                    }
                }
                if (!legalevt) {
                    THROW(ValueError() << errmsg{"MagnifySink: event number out of range!"});
                }
                delete input_runinfo;
                input_runinfo = nullptr;
            }  // Event/Sim found
        }      // celltree input

        rtree->Fill();
    }

    // Now deal with "shunting" input Magnify data to output.
    std::string ifname = m_cfg["input_filename"].asString();
    if (ifname.empty()) {
        // good, we shouldn't be peeking into the input file anyways.
        return;
    }
    auto toshunt = getset(m_cfg["shunt"]);
    if (toshunt.empty()) {
        log->warn("MagnifySink: no objects to copy but given input: {}", ifname);
        return;
    }
    log->debug("MagnifySink: sneaking peaks into input file: {}", ifname);

    TFile* input_tf = TFile::Open(ifname.c_str());
    for (auto name : toshunt) {
        if (name == "Trun" && !truncfg.empty()) {  // no double dipping
            continue;
        }
        TObject* obj = input_tf->Get(name.c_str());

        if (!obj) {
            log->warn("MagnifySink: warning: failed to find input object: \"{}\" for copying to output", name);
        }

        TTree* tree = dynamic_cast<TTree*>(obj);
        if (tree) {
            log->debug("MagnifySink: copying tree: \"{}\"", name);
            tree = tree->CloneTree();
            tree->SetDirectory(output_tf);
            continue;
        }

        TH1* hist = dynamic_cast<TH1*>(obj);
        if (hist) {
            log->debug("MagnifySink: copying hist: \"{}\"", name);
            hist->SetDirectory(output_tf);
            continue;
        }

        log->warn("MagnifySink: warning: not copying object of unknown type: \"{}\"", name);
    }

    delete input_tf;
    input_tf = nullptr;
}

bool Root::MagnifySink::operator()(const IFrame::pointer& frame, IFrame::pointer& out_frame)
{
    out_frame = frame;
    if (!frame) {
        // eos
        log->debug("MagnifySink: EOS");
        return true;
    }
    if (frame->traces()->empty()) {
        log->debug("MagnifySink: passing through empty frame ID {}", frame->ident());
        return true;
    }

    const std::string ofname = m_cfg["output_filename"].asString();
    const std::string mode = m_cfg["root_file_mode"].asString();
    log->debug("MagnifySink: opening for output: {} with mode {}", ofname, mode);
    TFile* output_tf = TFile::Open(ofname.c_str(), mode.c_str());

    for (auto tag : getset(m_cfg["frames"])) {
        auto trace_tag = tag;
        auto trace_has_tag = m_cfg["trace_has_tag"].asBool();
        if (!trace_has_tag) {
            trace_tag = "";
            log->debug("MagnifySink: set desired trace tag to \"\" as cfg::trace_has_tag=false");
        }

        ITrace::vector traces_byplane[3], traces = Aux::tagged_traces(frame, trace_tag);
        if (traces.empty()) {
            log->warn("MagnifySink: no tagged traces for \"{}\"", tag);
            continue;
        }

        log->debug("MagnifySink: tag: \"{}\" with {} traces", tag, traces.size());

        auto binnings = collate_byplane(traces, m_anode, traces_byplane);
        // if (binnings.size() == 4) {
        Binning tbin = binnings[3];
        for (int iplane = 0; iplane < 3; ++iplane) {
            if (traces_byplane[iplane].empty()) continue;
            const std::string name = Form("h%c_%s", 'u' + iplane, tag.c_str());
            Binning cbin = binnings[iplane];
            std::stringstream ss;
            ss << "MagnifySink:"
               << " cbin:" << cbin.nbins() << "[" << cbin.min() << "," << cbin.max() << "]"
               << " tbin:" << tbin.nbins() << "[" << tbin.min() << "," << tbin.max() << "]";
            log->debug(ss.str());

            // consider to add nrebin ...
            int nbins = tbin.nbins() / m_nrebin;

            TH2F* hist = new TH2F(name.c_str(), name.c_str(), cbin.nbins(), cbin.min(), cbin.max(), nbins, tbin.min(),
                                  tbin.max());

            hist->SetDirectory(output_tf);

            for (auto trace : traces_byplane[iplane]) {
                const int tbin1 = trace->tbin();
                const int ch = trace->channel();
                auto const& charges = trace->charge();
                for (size_t itick = 0; itick < charges.size(); ++itick) {
                    // Using what should be an identical call to
                    // Fill() ends up with a file that is more than
                    // two times bigger:
                    // 826M Jul 27 18:22 orig-bl-nf-fill-tbin.root
                    // 342M Jul 27 18:28 orig-bl-nf-setbincontent-tplus1.root
                    // hist->Fill(ch, tbin+itick+0.5, charges[itick]);
                    // edit: it's due to saving errors.

                    int ibin = (tbin1 - tbin.min() + itick) / m_nrebin;

                    hist->SetBinContent(cbin.bin(ch) + 1, ibin + 1,
                                        charges[itick] + hist->GetBinContent(cbin.bin(ch) + 1, ibin + 1));
                }
            }
        }
        // }
    }

    // Handle any trace summaries
    for (auto tag : getset(m_cfg["summaries"])) {
        // auto traces = get_tagged_traces(frame, tag);
        auto traces = Aux::tagged_traces(frame, tag);
        if (traces.empty()) {
            log->warn("MagnifySink: no traces tagged with \"{}\", skipping summary", tag);
            continue;
        }
        auto const& summary = frame->trace_summary(tag);
        if (summary.empty()) {
            log->warn("MagnifySink: warning: empty summary tagged with \"{}\", skipping summary", tag);
            continue;
        }

        std::string oper = get<std::string>(m_cfg["summary_operator"], tag, "set");

        log->debug("MagnifySink: saving summaries tagged with \"{}\" into per-plane hists", tag);

        // warning: this is going to get ugly.  we jump through these
        // hoops to avoid hard coding microboone channel ranges and
        // allow for sparse summaries on the domain of channels;
        const int ntot = traces.size();
        std::vector<int> perplane_channels[3];
        std::vector<double> perplane_values[3];
        for (int ind = 0; ind < ntot; ++ind) {
            const int chid = traces[ind]->channel();
            const int iplane = m_anode->resolve(chid).index();
            perplane_channels[iplane].push_back(chid);
            perplane_values[iplane].push_back(summary[ind]);
        }
        for (int iplane = 0; iplane < 3; ++iplane) {
            std::vector<int>& chans = perplane_channels[iplane];
            std::vector<double>& vals = perplane_values[iplane];
            if (!chans.empty()) {
                auto mme = std::minmax_element(chans.begin(), chans.end());
                const int ch0 = *mme.first;
                const int chf = *mme.second;
                const std::string hname = Form("h%c_%s", 'u' + iplane, tag.c_str());
                TH1F* hist = new TH1F(hname.c_str(), hname.c_str(), chf - ch0 + 1, ch0, chf);
                for (size_t ind = 0; ind < chans.size(); ++ind) {
                    const int x = chans[ind] + 0.5;
                    const double val = vals[ind];
                    if (oper == "set") {
                        int bin = hist->FindBin(x);
                        hist->SetBinContent(bin, val);
                    }
                    else {
                        hist->Fill(x, val);
                    }
                }
                hist->SetDirectory(output_tf);
            }
        }
    }

    {
        std::unordered_map<std::string, std::string> cmmkey2treename;
        for (auto jcmmtree : m_cfg["cmmtree"]) {
            cmmkey2treename[jcmmtree[0].asString()] = jcmmtree[1].asString();
        }

        Waveform::ChannelMaskMap input_cmm = frame->masks();
        for (auto const& it : input_cmm) {
            auto cmmkey = it.first;
            auto ct = cmmkey2treename.find(cmmkey);
            if (ct == cmmkey2treename.end()) {
                log->warn("MagnifySink: warning: found channel mask \"{}\", but no tree configured to accept it",
                          cmmkey);
                continue;
            }

            auto treename = ct->second;

            log->debug("MagnifySink: saving channel mask \"{}\" to tree \"{}\"", cmmkey, treename);

            TTree* tree = new TTree(treename.c_str(), treename.c_str());
            int chid = 0, plane = 0, start_time = 0, end_time = 0;
            tree->Branch("chid", &chid, "chid/I");
            tree->Branch("plane", &plane, "plane/I");
            tree->Branch("start_time", &start_time, "start_time/I");
            tree->Branch("end_time", &end_time, "end_time/I");
            tree->SetDirectory(output_tf);

            for (auto const& chmask : it.second) {
                chid = chmask.first;
                plane = m_anode->resolve(chid).index();
                auto mask = chmask.second;
                for (size_t ind = 0; ind < mask.size(); ++ind) {
                    start_time = mask[ind].first;
                    end_time = mask[ind].second;
                    tree->Fill();
                }
            }
        }
    }

    do_shunt(output_tf);

    auto count = output_tf->Write();
    log->debug("MagnifySink: closing output file {}, wrote {} bytes", ofname, count);
    output_tf->Close();
    delete output_tf;
    output_tf = nullptr;
    return true;
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
