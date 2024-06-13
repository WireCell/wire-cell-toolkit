#include "WireCellRoot/UbooneBlobSource.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleBlob.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"

#include "TInterpreter.h"

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

    // See in_views().
    m_views = 0;
    auto jviews = cfg["views"];
    if (jviews.isString() && jviews.size()) {
        std::string views = jviews.asString();
        for (char letter : views) {
            if (letter > 'Z') letter -= 32; // upper case
            if (letter == 'U') m_views |= 1;
            if (letter == 'V') m_views |= 2;
            if (letter == 'W') m_views |= 4;
        }
    }
    else {
        if (m_kind == "live") {
            m_views = 7;
        }
        else {
            m_views = 3;
        }
    }
    log->debug("loading {} blobs of views bits {}", m_kind, m_views);
}

WireCell::Configuration Root::UbooneBlobSource::default_configuration() const
{
    Configuration cfg;
    cfg["input"] = Json::arrayValue;
    cfg["kind"] = m_kind;
    cfg["views"] = Json::arrayValue;
    return cfg;
}



namespace WireCell::Root {

    /**
       Sole interface to data in ROOT TTrees.

       FIXME: consider T_bad_ch.
    */
    class UbooneBlobSourceTrees {

        std::unique_ptr<TFile> m_tfile;
        TTree* m_activity{nullptr}; // always load
        TTree* m_live{nullptr};     // always load
        TTree* m_dead{nullptr};     // load only for kind=="dead"

        Long64_t m_nentries{0}; // number of entries in the trees.
        Long64_t m_entry{-1};   // not yet loaded.  Only use next() to advance

    public:

        // The data is available after construction and refreshed on each call
        // to next().  The data is invalid when/if next() throws IndexError and
        // "this" instance should be dropped shortly after.

        struct Header {
            // double eventTime{0};    // 1.45767e9. Unix time?
            double triggerTime{0};  // 1.48822e10. Units?
            int eventNo{0};
            int nrebin{0};
        };
        Header header;

        struct Activity {

            // "activity"
            std::vector<int> *timesliceId{nullptr}; // [0,2399]
            std::vector<std::vector<int>> *timesliceChannel{nullptr}; // [0,8255]
            // Units are number of signal electrons.
            std::vector<std::vector<int>> *raw_charge{nullptr};
            std::vector<std::vector<int>> *raw_charge_err{nullptr};
        };
        Activity activity;

        // TC and TDC are nearly the same.
        struct Blob {
            std::vector<int> *cluster_id_vec{nullptr};

            // 0 value indicates the view for a blob is "dead"
            std::vector<int> *flag_u_vec{nullptr};
            std::vector<int> *flag_v_vec{nullptr};
            std::vector<int> *flag_w_vec{nullptr};

            std::vector<std::vector<int>> *wire_index_u_vec{nullptr};
            std::vector<std::vector<int>> *wire_index_v_vec{nullptr};
            std::vector<std::vector<int>> *wire_index_w_vec{nullptr};

            // Lord, forgive us our sins.
            std::vector<std::vector<int>*> flag_uvw;
        };
        struct TC : Blob { 
            std::vector<int> *time_slice_vec{nullptr};               // TC: singular
            std::vector<double> *q_vec{nullptr};
        };
        TC live;

        struct TDC : Blob {
            std::vector<std::vector<int>> *time_slices_vec{nullptr}; // TDC: plural
        };
        TDC dead;

    public:

        // Construct the interface to the trees.
        //
        // This does NOT load any entry.  Must call next() to load first entry, etc.
        UbooneBlobSourceTrees(const std::string& tfile_name, bool include_dead=false) {
            
            m_tfile.reset(TFile::Open(tfile_name.c_str(), "READ"));
            if (!m_tfile) {
                raise<IOError>("failed to open %s", tfile_name);
            }

            set_activty_address();

            m_live = open_ttree("TC");
            set_blob_address(m_live, live);
            m_live->SetBranchAddress("q", &live.q_vec);
            m_live->SetBranchAddress("time_slice", &live.time_slice_vec);

            if (include_dead) {
                m_dead = open_ttree("TDC");
                set_blob_address(m_dead, dead);
                m_dead->SetBranchAddress("time_slice", &dead.time_slices_vec);
            }
        }

        void next() {
            ++m_entry;
            if (m_entry < m_nentries) { 
                m_activity->GetEntry(m_entry);
                m_live->GetEntry(m_entry);
                if (m_dead) {
                    m_dead->GetEntry(m_entry);
                }
                return;
            }
            raise<IndexError>("attempt to get entry %d past end of TTree with %d",
                              m_entry, m_nentries);
        }

        // Get live or dead blob data
        const Blob& blob(const std::string& which) {
            if (which == "live") return live;
            return dead;
        }

        int nentries() const { return m_nentries; }
        int entry() const { return m_entry; }
        int nslices() const { return activity.timesliceId->size(); }

    private:


        void set_activty_address() {
            m_activity = open_ttree("Trun");
            if (!m_activity) {
                raise<IOError>("failed get activity tree");
            }
            m_nentries = m_activity->GetEntries();

            // "header"
            m_activity->SetBranchAddress("triggerTime", &header.triggerTime); 
            m_activity->SetBranchAddress("eventNo", &header.eventNo);
            m_activity->SetBranchAddress("nrebin", &header.nrebin);

            // "activity"
            m_activity->SetBranchAddress("timesliceId", &activity.timesliceId);
            m_activity->SetBranchAddress("timesliceChannel", &activity.timesliceChannel);
            m_activity->SetBranchAddress("raw_charge", &activity.raw_charge);
            m_activity->SetBranchAddress("raw_charge_err", &activity.raw_charge_err);
        }
        void set_blob_address(TTree* tree, Blob& blob) {
            if (!tree) {
                raise<IOError>("failed to open required blob tree (TC or TDC)");
            }
            tree->SetBranchAddress("cluster_id", &blob.cluster_id_vec); 
            tree->SetBranchAddress("flag_u", &blob.flag_u_vec);
            tree->SetBranchAddress("flag_v", &blob.flag_v_vec);
            tree->SetBranchAddress("flag_w", &blob.flag_w_vec);
            blob.flag_uvw = {blob.flag_u_vec, blob.flag_v_vec, blob.flag_w_vec};
            tree->SetBranchAddress("wire_index_u", &blob.wire_index_u_vec);
            tree->SetBranchAddress("wire_index_v", &blob.wire_index_v_vec);
            tree->SetBranchAddress("wire_index_w", &blob.wire_index_w_vec);
        }

        TTree* open_ttree(const std::string& tree_name) {
            auto ttree = reinterpret_cast<TTree*>(m_tfile->Get(tree_name.c_str()));
            if (!ttree) {
                raise<IOError>("failed to get TTree %s from %s", tree_name, m_tfile->GetName());
            }
            return ttree;
        }
    };
};


bool Root::UbooneBlobSource::next()
{
    if (!m_data) {
        // We are starting a new file.

        if (m_input.empty()) {
            // We have exhausted our input files.
            return false;
        }

        auto fname = m_input.back();
        m_input.pop_back();

        try {
            m_data = std::make_unique<UbooneBlobSourceTrees>(fname, m_kind == "dead");
        }
        catch (IOError& err) {
            log->warn("failed to open {}, skipping", fname);
            m_data = nullptr;
            return next();
        }
        catch (IndexError& err) {
            m_data = nullptr;
            return next();
        }

        // sanity check file
        const int nentries = m_data->nentries();
        if (! nentries) {
            log->warn("no entries {}, skipping", fname);
            m_data = nullptr;
            return next();
        }
        log->debug("starting {} with {} entries", fname, nentries);
    }

    // We have an open file with at least some entries, try to advance
    try {
        m_data->next();
    }
    catch (IndexError& err) {
        m_data = nullptr;
        return next();
    }

    // sanity check entry
    int nslices = m_data->nslices();
    if (!nslices) {
        log->warn("no slices in entry {}, skipping", m_data->entry());
        return next();
    }

    log->debug("read {} slices in entry {}", nslices, m_data->entry());
    return true;
}


IFrame::pointer Root::UbooneBlobSource::gen_frame()
{
    // We have no good way to make traces.  Best we might do is spread activity
    // over nrebin ticks.  Assuming the traces are not important for downstream
    // we make a frame that only holds metadata.
    return std::make_shared<SimpleFrame>(m_data->header.eventNo,
                                         m_data->header.triggerTime, 
                                         0.5*units::us);
}


IChannel::pointer Root::UbooneBlobSource::get_channel(int chanid)
{
    return m_anode->channel(chanid);
}


bool Root::UbooneBlobSource::in_views(int bind)
{
    // Inclusion test depends if we are loading live or dead blobs.
    const int want = m_kind == "live" ? 1 : 0;
    const auto& flag_uvw = m_data->blob(m_kind).flag_uvw;

    int cone = 0;
    for (int pln=0; pln<3; ++pln) {
        if (want == flag_uvw[pln]->at(bind)) {
            cone |= 1<<pln;
        }
    }

    return cone == m_views;
}


// If any view of the blob is "flagged" 0 then set its channels activity to special "bad" value. 
// FIXME: consider T_bad_ch
void Root::UbooneBlobSource::bodge_channels(SimpleSlicePtr slice, const RayGrid::Blob& blob, int bind)
{
    const auto& flag_uvw = m_data->blob(m_kind).flag_uvw;

    const std::vector<int> choff = {0,2400,4800};
    const auto& strips = blob.strips();
    for (int pln = 0; pln < 3; ++pln) {
        if (flag_uvw[pln]->at(bind) == 1) {
            continue;           // view is "good", keep whatever activity exists.
        }
        const auto& [beg,end] = strips[2+pln].bounds;
        for (int wire = beg; wire<end; ++wire) {
            const int ch = wire + choff[pln];
            auto ich = get_channel(ch);
            slice->activity()[ich] = {0, 1e12};
        }        
    }
}


/*

  WCP's wire-in plane (the "wire" column) follows WCT standard.

  Excerpt from ChannelWireGeometry_v2.txt.
  Edited to round and put into mm units.

  # ch   plane    wire    sx      sy    sz   ex      ey    ez
  0       0       0        0    1171     0    0    1175     5
  2399    0       2399     0   -1155 10365    0   -1152 10370
  2400    1       0       -3   -1152     0   -3   -1155     5
  4799    1       2399    -3    1174 10365   -3    1172 10370
  4800    2       0       -6   -1155     3   -6    1175     3
  8255    2       3455    -6   -1155 10368   -6    1175 10368

  Can also check with:

    wcwires -e 5e-2 -v microboone-celltree-wires-v2.1.json.bz2

  With a lower imprecision value, it will show pitch/direction precision errors,
  but no ordering errors.

  So, we may take the wire_index_{u,v,w}_vec values at face value!
*/
std::pair<int,int> Root::UbooneBlobSource::make_strip(const std::vector<int>& winds)
{
    const auto& [imin, imax] = std::minmax_element(winds.begin(), winds.end());
    return std::make_pair(*imin, *imax + 1);
}


// Map prototype's slice index to toolkit ISlice.
// Note, WCP's sole concept of "slice" is what we would call here just "live" slice.
// It is "nrebin*tick" long.  WCT also has "dead" slice that spans some number live slices.
// This method is used by both load_live() and load_dead().
Root::UbooneBlobSource::SliceMap Root::UbooneBlobSource::load_slices()
{
    SliceMap slices;

    auto iframe = gen_frame();

    const auto& act = m_data->activity; // abbrev

    const auto& tids = act.timesliceId;
    const size_t nslices = tids->size();
    const double span = m_data->header.nrebin * 0.5 * units::us;
    for (size_t sind=0; sind<nslices; ++sind) {

        const auto& q = act.raw_charge->at(sind);
        const auto& dq = act.raw_charge_err->at(sind);
        const auto& chans = act.timesliceChannel->at(sind);
        const size_t nchans = chans.size();
        
        ISlice::map_t activity;
        for (size_t cind=0; cind<nchans; ++cind) {
            auto ichan = get_channel(cind);
            activity[ichan] = ISlice::value_t(q[cind], dq[cind]);
        }
        int tsind = tids->at(sind);
        auto sslice = std::make_shared<SimpleSlice>(iframe, tsind, tsind*span, span, activity);
        slices[tsind] = sslice;
    }            
    return slices;
}


IBlobSet::pointer Root::UbooneBlobSource::load_live()
{
    auto iframe = gen_frame();

    SliceMap live_slices = load_slices();

    const auto& live_data = m_data->live;
    const RayGrid::Coordinates& coords = m_iface->raygrid();

    // Blobs
    const size_t nblobs = live_data.cluster_id_vec->size();
    IBlob::vector iblobs;
    for (size_t bind=0; bind<nblobs; ++bind) {

        if (! in_views(bind)) continue;

        RayGrid::Blob blob;
        blob.add(coords, RayGrid::Strip{0, {0,1}});
        blob.add(coords, RayGrid::Strip{1, {0,1}});
        blob.add(coords, RayGrid::Strip{2, make_strip(live_data.wire_index_u_vec->at(bind))});
        blob.add(coords, RayGrid::Strip{3, make_strip(live_data.wire_index_v_vec->at(bind))});
        blob.add(coords, RayGrid::Strip{4, make_strip(live_data.wire_index_w_vec->at(bind))});
            
        const int sid = live_data.time_slice_vec->at(bind);
        const float blob_charge = live_data.q_vec->at(bind);

        auto slice = live_slices[sid];
        bodge_channels(slice, blob, bind);

        auto iblob = std::make_shared<SimpleBlob>(bind, blob_charge, 0, blob, slice, m_iface);
        iblobs.push_back(iblob);
    }
    auto main_slice = live_slices.begin()->second;
    return std::make_shared<SimpleBlobSet>(iframe->ident(), main_slice, iblobs);
}


IBlobSet::pointer Root::UbooneBlobSource::load_dead()
{
    // There are some WCP/WCT concept mismatches that are critical here.  It is
    // assumed that all WCP dead blobs line up in time to form cleavage points
    // which we call "dead slices" in WCT.  That is, we assume any two WCP
    // "dead" blobs either totally overlap in time or do not overlap at all.

    SliceMap live_slices = load_slices();
    auto iframe = live_slices.begin()->second->frame();

    const RayGrid::Coordinates& coords = m_iface->raygrid();
    const auto& dead_data = m_data->dead;

    IBlob::vector iblobs;
    const size_t nblobs = dead_data.cluster_id_vec->size();
    if (!nblobs) {
        log->warn("no dead (TDC) blobs");
        return std::make_shared<SimpleBlobSet>(iframe->ident(), nullptr, iblobs);
    }

    SliceMap dead_slices;

    // pick up first slice to give to the blob set.
    ISlice::pointer main_slice;

    for (size_t bind=0; bind<nblobs; ++bind) {

        if (! in_views(bind)) continue;

        const auto& dead_tsvec = dead_data.time_slices_vec->at(bind);
        const int dead_tsind = dead_tsvec.front(); // assumes ordered....
        auto dead_sit = dead_slices.find(dead_tsind);

        // Either created or looked up
        SimpleSlicePtr dead_slice = nullptr;
        if (dead_sit == dead_slices.end()) {

            // First blob in this dead slice, create the ISlice
            const double dead_span = m_data->header.nrebin * 0.5 * units::us * dead_tsvec.size();
            const double dead_start = dead_tsind * dead_span;

            // Add in activity for all live slices spanned by this first dead blob.
            ISlice::map_t activity;
            for (int tsind : dead_tsvec) {
                const auto live_slice = live_slices[tsind];
                if (!live_slice) {
                    /// Don't warn.  The set of live slices may be sparse.
                    // log->warn("dead slice {}+{} spans a non-existent live slice {}",
                    //           dead_tsind, dead_tsvec.size(), tsind);
                    continue;
                }
                for (const auto& [ich, act] : live_slice->activity()) {
                    activity[ich] += act;
                }
            }
            dead_slice = std::make_shared<SimpleSlice>(iframe, dead_tsind, dead_start, dead_span, activity);
        }
        else {
            dead_slice = dead_sit->second;
        }
        if (! main_slice) {
            main_slice = dead_slice;
        }

        RayGrid::Blob blob;
        blob.add(coords, RayGrid::Strip{0, {0,1}});
        blob.add(coords, RayGrid::Strip{1, {0,1}});
        blob.add(coords, RayGrid::Strip{2, make_strip(dead_data.wire_index_u_vec->at(bind))});
        blob.add(coords, RayGrid::Strip{3, make_strip(dead_data.wire_index_v_vec->at(bind))});
        blob.add(coords, RayGrid::Strip{4, make_strip(dead_data.wire_index_w_vec->at(bind))});

        const float bval = 0;
        const float bunc = 0;

        bodge_channels(dead_slice, blob, bind);

        auto iblob = std::make_shared<SimpleBlob>(bind, bval, bunc, blob, dead_slice, m_iface);
        iblobs.push_back(iblob);
    }
    return std::make_shared<SimpleBlobSet>(iframe->ident(), main_slice, iblobs);
}


bool Root::UbooneBlobSource::operator()(IBlobSet::pointer& blobset)
{
    blobset = nullptr;

    if (m_done) {
        // log->debug("past EOS at call {}, stop calling me", ++m_calls);
        return false;
    }

    if (! next()) {
        log->debug("EOS at call {}", ++m_calls);
        m_done = true;
        return true;
    }

    if (m_kind == "live") {
        blobset = load_live();
    }
    else {
        blobset = load_dead();
    }

    return true;
}

