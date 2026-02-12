#include "WireCellRoot/UbooneMagnifyTrackingSink.h"

#include "TFile.h"
#include "TTree.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellAux/TensorDMpointtree.h"
#include "WireCellClus/Facade_Grouping.h"
// #include "WireCellIface/IAnodePlane.h"

WIRECELL_FACTORY(UbooneMagnifyTrackingSink, WireCell::Root::UbooneMagnifyTrackingSink,
                 WireCell::ITensorSetFilter, WireCell::IConfigurable)

using namespace WireCell;
using WireCell::Aux::TensorDM::as_pctree;

Root::UbooneMagnifyTrackingSink::UbooneMagnifyTrackingSink()
  : log(Log::logger("tracking"))
  , m_runNo(0)
  , m_subRunNo(0)
  , m_eventNo(0)
{
}

Root::UbooneMagnifyTrackingSink::~UbooneMagnifyTrackingSink() {}

void Root::UbooneMagnifyTrackingSink::configure(const WireCell::Configuration& cfg)
{
    std::string fn = cfg["output_filename"].asString();
    if (fn.empty()) {
        THROW(ValueError() << errmsg{"Must provide output filename to UbooneMagnifyTrackingSink"});
    }

    m_runNo = get<int>(cfg, "runNo", 0);
    m_subRunNo = get<int>(cfg, "subRunNo", 0);
    m_eventNo = get<int>(cfg, "eventNo", 0);

    // Configure anodes
    auto anode_tns = cfg["anodes"];
    for (auto anode_tn : anode_tns) {
        auto anode = Factory::find_tn<IAnodePlane>(anode_tn.asString());
        m_anodes.push_back(anode);
    }

    // Configure detector volumes
    m_dv = Factory::find_tn<IDetectorVolumes>(cfg["detector_volumes"].asString());

    m_cfg = cfg;

    create_file();
}

WireCell::Configuration Root::UbooneMagnifyTrackingSink::default_configuration() const
{
    Configuration cfg;

    // Name of ROOT file to write.
    cfg["output_filename"] = "track_com.root";

    // The ROOT file mode with which to open the file.
    cfg["root_file_mode"] = "RECREATE";

    // Datapath to find the grouping in tensor set
    cfg["datapath"] = "pointtrees/%d";

    // Which grouping to read (e.g., "live", "dead")
    cfg["grouping"] = "live";

    // Anodes (required)
    cfg["anodes"] = Json::arrayValue;

    // Detector volumes (required)
    cfg["detector_volumes"] = "";

    // Run, subrun, and event numbers
    cfg["runNo"] = 0;
    cfg["subRunNo"] = 0;
    cfg["eventNo"] = 0;

    return cfg;
}

void Root::UbooneMagnifyTrackingSink::create_file()
{
    const std::string ofname = m_cfg["output_filename"].asString();
    const std::string mode = "RECREATE";
    log->debug("UbooneMagnifyTrackingSink: creating output file: {}", ofname);
    TFile* output_tf = TFile::Open(ofname.c_str(), mode.c_str());
    output_tf->Close("R");
    delete output_tf;
    output_tf = nullptr;
}

void Root::UbooneMagnifyTrackingSink::write_bad_channels(TFile* output_tf, const ITensorSet::pointer& ts)
{
    if (!ts) {
        log->warn("UbooneMagnifyTrackingSink: null tensor set");
        return;
    }

    // Get the datapath from configuration
    std::string datapath = get<std::string>(m_cfg, "datapath", "pointtrees/%d");
    std::string grouping_name = get<std::string>(m_cfg, "grouping", "live");

    const int ident = ts->ident();
    if (datapath.find("%") != std::string::npos) {
        datapath = String::format(datapath, ident);
    }
     // Append the grouping name to the datapath
    datapath = datapath + "/" + grouping_name;

    // Get the point cloud tree from tensor set
    const auto& tens = *ts->tensors();


    std::unique_ptr<PointCloud::Tree::Points::node_t> root_ptr;
    try {
        root_ptr = as_pctree(tens, datapath);
    }
    catch (WireCell::KeyError& err) {
        log->warn("UbooneMagnifyTrackingSink: no pc-tree at tensor datapath {}", datapath);
        return;
    }

    if (!root_ptr) {
        log->warn("UbooneMagnifyTrackingSink: failed to get pc-tree root");
        return;
    }

    // Create Grouping facade from the node
    auto* grouping = root_ptr->value.facade<Clus::Facade::Grouping>();

    if (!grouping) {
        log->warn("UbooneMagnifyTrackingSink: failed to create Grouping facade");
        return;
    }

    // Set anodes and detector volumes on the grouping (required for accessing channel information)
    grouping->set_anodes(m_anodes);
    grouping->set_detector_volumes(m_dv);

    // Create the T_bad_ch tree
    TTree* tree = new TTree("T_bad_ch", "T_bad_ch");
    tree->SetDirectory(output_tf);

    // Branch variables
    int chid = 0;
    int plane = 0;
    int start_time = 0;
    int end_time = 0;
    int runNo = m_runNo;
    int subRunNo = m_subRunNo;
    int eventNo = m_eventNo;

    // Create branches
    tree->Branch("chid", &chid, "chid/I");
    tree->Branch("plane", &plane, "plane/I");
    tree->Branch("start_time", &start_time, "start_time/I");
    tree->Branch("end_time", &end_time, "end_time/I");
    tree->Branch("runNo", &runNo, "runNo/I");
    tree->Branch("subRunNo", &subRunNo, "subRunNo/I");
    tree->Branch("eventNo", &eventNo, "eventNo/I");

    // Get the wpids to determine which APAs/faces exist
    auto wpids = grouping->wpids();

    std::set<std::pair<int, int>> apa_face_set;
    for (const auto& wpid : wpids) {
        apa_face_set.insert({wpid.apa(), wpid.face()});
    }

    log->debug("UbooneMagnifyTrackingSink: found {} APA/face combinations", apa_face_set.size());

    // Iterate through all APAs, faces, and planes
    for (const auto& [apa, face] : apa_face_set) {
        for (int pind = 0; pind < 3; ++pind) {  // 3 planes: U, V, W
            try {
                // Get all dead channels for this plane
                auto dead_chs = grouping->get_all_dead_chs(apa, face, pind);

                log->debug("UbooneMagnifyTrackingSink: APA={}, face={}, plane={}: {} dead channels",
                          apa, face, pind, dead_chs.size());

                // Fill tree with dead channel information
                plane = pind;
                for (const auto& [ch, time_range] : dead_chs) {
                    chid = ch;
                    start_time = time_range.first;
                    end_time = time_range.second;
                    tree->Fill();
                }
            }
            catch (const std::exception& e) {
                log->warn("UbooneMagnifyTrackingSink: failed to get dead channels for APA={}, face={}, plane={}: {}",
                         apa, face, pind, e.what());
            }
        }
    }

    log->debug("UbooneMagnifyTrackingSink: wrote {} entries to T_bad_ch", tree->GetEntries());
}

bool Root::UbooneMagnifyTrackingSink::operator()(const ITensorSet::pointer& in, ITensorSet::pointer& out)
{
    out = in;
    if (!in) {
        // eos
        log->debug("UbooneMagnifyTrackingSink: EOS");
        return true;
    }

    // std::cout << "Test: UbooneMagnifyTrackingSink: got tensor set with ident: " << in->ident() << std::endl;

    const std::string ofname = m_cfg["output_filename"].asString();
    const std::string mode = m_cfg["root_file_mode"].asString();
    log->debug("UbooneMagnifyTrackingSink: opening for output: {} with mode {}", ofname, mode);
    TFile* output_tf = TFile::Open(ofname.c_str(), mode.c_str());

    write_bad_channels(output_tf, in);

    auto count = output_tf->Write();
    log->debug("UbooneMagnifyTrackingSink: closing output file {}, wrote {} bytes", ofname, count);
    output_tf->Close();
    delete output_tf;
    output_tf = nullptr;

    return true;
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End: