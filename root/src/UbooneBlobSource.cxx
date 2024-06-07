#include "WireCellRoot/UbooneBlobSource.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleSlice.h"
#include "WireCellAux/SimpleBlob.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"

#include <algorithm> // minmax


WIRECELL_FACTORY(UbooneBlobSource,
                 WireCell::Root::UbooneBlobSource,
                 WireCell::INamed,
                 WireCell::IBlobSetSource,
                 WireCell::IConfigurable)


using namespace WireCell;
using namespace WireCell::Aux;


Root::UbooneBlobSource::UbooneBlobSource()
    : Aux::Logger("UbooneBlobSource", "root")
    , m_calls(0)
{
}

Root::UbooneBlobSource::~UbooneBlobSource()
{
}

void Root::UbooneBlobSource::configure(const WireCell::Configuration& cfg)
{
    auto anode_tn = get<std::string>(cfg, "anode", "AnodePlane");
    m_anode = Factory::find_tn<IAnodePlane>(anode_tn);
    m_iface = nullptr;
    for (auto& iface : m_anode->faces()) {
        if (iface) {            // take first non-null
            m_iface = iface;
            break;
        }
    }
    if (!m_iface) {
        log->critical("anode %s has no face", anode_tn);
        raise<ValueError>("UbooneBlobSource: anode face is required");
    }        

    m_kind = get(cfg, "kind", m_kind);
    auto input = cfg["input"];
    if (input.isNull()) {
        log->critical("input is required");
        raise<ValueError>("UbooneBlobSource: input is required");
    }
    if (input.isString()) {
        m_input.push_back(input.asString());
    }
    else {
        m_input.clear();
        for (const auto& one: input) {
            m_input.push_back(one.asString());
        }
    }
    // so we can pop_back.
    std::reverse(m_input.begin(), m_input.end());
}

WireCell::Configuration Root::UbooneBlobSource::default_configuration() const
{
    Configuration cfg;
    cfg["input"] = Json::arrayValue;
    cfg["kind"] = m_kind;
    return cfg;
}


/**
   Bags to hold the branch addresses relevant to building an IBlobSet.
*/

namespace WireCell::Root {
    struct UbooneBlobSourceTrees {


        // "header"
        // double eventTime{0};    // 1.45767e9. Unix time?
        double triggerTime{0};  // 1.48822e10. Units?
        // int runNo{0}; // 5384
        // int subrunNo{0}; // 130
        int eventNo{0}; // 6501
        int nrebin{0};          // 4
        // int frame_length{0};    // 3200 ?
        // float unit_dis{0};      // 1.10306???

        // "activity"
        std::vector<int> *timesliceId{nullptr}; // [0,2399]
        std::vector<std::vector<int>> *timesliceChannel{nullptr}; // [0,8255]
        // Units are number of signal electrons.
        std::vector<std::vector<int>> *raw_charge{nullptr};
        std::vector<std::vector<int>> *raw_charge_err{nullptr};

        void set_trun(TTree& trun) {
            // "header"
            trun.SetBranchAddress("triggerTime",&triggerTime); 
            trun.SetBranchAddress("eventNo",&eventNo);
            // trun.SetBranchAddress("runNo",&runNo);
            // trun.SetBranchAddress("subRunNo",&subrunNo);
            // trun.SetBranchAddress("eventTime",&eventTime);
            // trun.SetBranchAddress("unit_dis",&unit_dis);
            trun.SetBranchAddress("nrebin",&nrebin);
            // trun.SetBranchAddress("frame_length",&frame_length);

            // "activity"
            trun.SetBranchAddress("timesliceId",&timesliceId);
            trun.SetBranchAddress("timesliceChannel",&timesliceChannel);
            trun.SetBranchAddress("raw_charge",&raw_charge);
            trun.SetBranchAddress("raw_charge_err",&raw_charge_err);
        }

        // "blob"
        // std::vector<int> *cluster_id_vec = new std::vector<int>;
        std::vector<int> *time_slice_vec = new std::vector<int>;
        // std::vector<int> *nwire_u_vec{nullptr};
        // std::vector<int> *nwire_v_vec{nullptr};
        // std::vector<int> *nwire_w_vec{nullptr};
        // std::vector<int> *flag_u_vec{nullptr};
        // std::vector<int> *flag_v_vec{nullptr};
        // std::vector<int> *flag_w_vec{nullptr};
        std::vector<std::vector<int>> *wire_index_u_vec{nullptr};
        std::vector<std::vector<int>> *wire_index_v_vec{nullptr};
        std::vector<std::vector<int>> *wire_index_w_vec{nullptr};

        // Only gets defined if "q" exists.  ie, for TC not TDC.
        std::vector<double> *q_vec{nullptr};

        // Can be called with TC or TDC tree
        void set_tc(TTree& tc) {
            // tc.SetBranchAddress("cluster_id",&cluster_id_vec);
            tc.SetBranchAddress("time_slice",&time_slice_vec);
            // tc.SetBranchAddress("nwire_u",&nwire_u_vec);
            // tc.SetBranchAddress("nwire_v",&nwire_v_vec);
            // tc.SetBranchAddress("nwire_w",&nwire_w_vec);
            // tc.SetBranchAddress("flag_u",&flag_u_vec);
            // tc.SetBranchAddress("flag_v",&flag_v_vec);
            // tc.SetBranchAddress("flag_w",&flag_w_vec);
            tc.SetBranchAddress("wire_index_u",&wire_index_u_vec);
            tc.SetBranchAddress("wire_index_v",&wire_index_v_vec);
            tc.SetBranchAddress("wire_index_w",&wire_index_w_vec);

            if (tc.GetBranch("q")) { 
                tc.SetBranchAddress("q",&q_vec);
            }
        }
    };
};




bool Root::UbooneBlobSource::next()
{
    if (m_tfile) {
        // We have an open file, try next entry.
        ++m_entry;
    }
    else {
        // We are starting a new file.
        m_entry=0;
        m_tc = nullptr;
        m_trun = nullptr;
        m_data = nullptr;

        if (m_input.empty()) {
            // We have exhausted the list of files.
            return false;
        }

        auto fname = m_input.back();
        m_input.pop_back();

        m_tfile.reset(TFile::Open(fname.c_str(), "READ"));
        if (!m_tfile) {
            log->error("failed to open {}, skipping", fname);
            return next();
        }
        m_trun = reinterpret_cast<TTree*>(m_tfile->Get("Trun"));
        if (!m_trun) {
            log->error("failed to get TTree Trun from {}, skipping", fname);
            m_tfile = nullptr;
            return next();
        }
        m_tc = reinterpret_cast<TTree*>(m_tfile->Get(m_kind.c_str()));
        if (!m_tc) {
            log->error("failed to get TTree {} from {}, skipping", m_kind, fname);
            m_tfile = nullptr;
            return next();
        }

        m_data = std::make_unique<UbooneBlobSourceTrees>();
        m_data->set_trun(*m_trun);
        m_data->set_tc(*m_tc);
    }

    // We have an open file with m_entry set to the one to try.
    
    if (m_entry < m_trun->GetEntries() && m_entry < m_tc->GetEntries()) {
        // We are still inside the current file.
        m_trun->GetEntry(m_entry);
        m_tc->GetEntry(m_entry);
        return true;
    }
        
    // We ran off the end of the file.
    m_tfile = nullptr;
    return next();
}


IFrame::pointer Root::UbooneBlobSource::gen_frame()
{
    // No way to make traces.
    return std::make_shared<SimpleFrame>(m_data->eventNo,
                                         m_data->triggerTime, 
                                         0.5*units::us);
}


IChannel::pointer Root::UbooneBlobSource::get_channel(int chanid)
{
    return m_anode->channel(chanid);
}

static std::pair<int,int> make_strip(const std::vector<int>& winds)
{
    const auto& [imin, imax] = std::minmax_element(winds.begin(), winds.end());
    return std::make_pair(*imin, *imax + 1);
}

IBlobSet::pointer Root::UbooneBlobSource::load()
{
    if (!m_data->timesliceId || !m_data->timesliceChannel) {
        log->error("no time slice information");
        return nullptr;
    }

    auto iframe = gen_frame();

    // Slices and their activity
    ISlice::vector islices;
    IChannel::vector ichannels;
    
    const auto& tids = m_data->timesliceId;
    const size_t nslices = tids->size();
    const double span = m_data->nrebin * 0.5 * units::us;
    for (size_t sind=0; sind<nslices; ++sind) {

        const auto& q = m_data->raw_charge->at(sind);
        const auto& dq = m_data->raw_charge_err->at(sind);
        const auto& chans = m_data->timesliceChannel->at(sind);;
        const size_t nchans = chans.size();
        
        ISlice::map_t activity;
        for (size_t cind=0; cind<nchans; ++cind) {
            auto ichan = get_channel(cind);
            activity[ichan] = ISlice::value_t(q[cind], dq[cind]);
        }
        int tsind = tids->at(sind);
        auto sslice = std::make_shared<SimpleSlice>(iframe, tsind, tsind*span, span, activity);
        islices.push_back(sslice);
    }            

    const RayGrid::Coordinates& coords = m_iface->raygrid();

    // Blobs
    const size_t nblobs = m_data->time_slice_vec->size();
    IBlob::vector iblobs;
    for (size_t bind=0; bind<nblobs; ++bind) {

        RayGrid::Blob blob;
        blob.add(coords, RayGrid::Strip{0, {0,1}});
        blob.add(coords, RayGrid::Strip{1, {0,1}});
        blob.add(coords, RayGrid::Strip{2, make_strip(m_data->wire_index_u_vec->at(bind))});
        blob.add(coords, RayGrid::Strip{2, make_strip(m_data->wire_index_v_vec->at(bind))});
        blob.add(coords, RayGrid::Strip{2, make_strip(m_data->wire_index_w_vec->at(bind))});
            
        int sind = m_data->time_slice_vec->at(bind);
        float blob_charge = 0;
        if (m_kind == "TC" && m_data->q_vec) {
            blob_charge = m_data->q_vec->at(bind);
        }

        auto iblob = std::make_shared<SimpleBlob>(bind, blob_charge, 0, blob, islices[sind], m_iface);
        iblobs.push_back(iblob);
    }
    return std::make_shared<SimpleBlobSet>(iframe->ident(), iblobs);
}

bool Root::UbooneBlobSource::operator()(IBlobSet::pointer& blobset)
{
    blobset = nullptr;

    if (! next()) {
        log->debug("EOS at call {}", ++m_calls);
        return true;
    }

    blobset = load();
    return true;
}

