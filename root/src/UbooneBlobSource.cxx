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
            m_views = 7;        // uvw
        }
        else {
            m_views = 3;        // uv
        }
    }
    // Work out which planes get "bodged" (set as "bad"/"dead")
    // and maybe set as "dummy" for dead.  All this is only for 2-view.
    if (m_kind == "live") {
        if (m_views == 2+3) m_bodged = {0}; // U
        if (m_views == 3+1) m_bodged = {1}; // V
        if (m_views == 1+2) m_bodged = {2}; // W
    }
    else {
        auto iwplanes = m_iface->planes();

        if (m_views == 2+3) {
            m_bodged = {1,2}; // !U
            m_dummy = iwplanes[0];
        }
        if (m_views == 3+1) {
            m_bodged = {2,0}; // !V
            m_dummy = iwplanes[1];
        }
        if (m_views == 1+2) {
            m_bodged = {0,1}; // !W
            m_dummy = iwplanes[2];
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
    */
    class UbooneBlobSourceTrees {

        std::unique_ptr<TFile> m_tfile;
        TTree* m_activity{nullptr}; // always load
        TTree* m_live{nullptr};     // always load
        TTree* m_dead{nullptr};     // load only for kind=="dead"
        TTree* m_bad{nullptr};      // optional

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

            void set_addresses(TTree& tree) {
                tree.SetBranchAddress("triggerTime", &triggerTime);
                tree.SetBranchAddress("eventNo", &eventNo);
                tree.SetBranchAddress("nrebin", &nrebin);
            }
        };
        Header header;

        struct Activity {

            // "activity"
            std::vector<int> *timesliceId{nullptr}; // [0,2399]
            std::vector<std::vector<int>> *timesliceChannel{nullptr}; // [0,8255]
            // Units are number of signal electrons.
            std::vector<std::vector<int>> *raw_charge{nullptr};
            std::vector<std::vector<int>> *raw_charge_err{nullptr};

            void set_addresses(TTree& tree) {
                tree.SetBranchAddress("timesliceId", &timesliceId);
                tree.SetBranchAddress("timesliceChannel", &timesliceChannel);
                tree.SetBranchAddress("raw_charge", &raw_charge);
                tree.SetBranchAddress("raw_charge_err", &raw_charge_err);
            }
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

            void set_addresses(TTree& tree) {
                tree.SetBranchAddress("cluster_id", &cluster_id_vec); 
                tree.SetBranchAddress("flag_u", &flag_u_vec);
                tree.SetBranchAddress("flag_v", &flag_v_vec);
                tree.SetBranchAddress("flag_w", &flag_w_vec);
                tree.SetBranchAddress("wire_index_u", &wire_index_u_vec);
                tree.SetBranchAddress("wire_index_v", &wire_index_v_vec);
                tree.SetBranchAddress("wire_index_w", &wire_index_w_vec);

                // derived
                flag_uvw = {flag_u_vec, flag_v_vec, flag_w_vec};
            }
        };

        struct TC : Blob { 
            std::vector<int> *time_slice_vec{nullptr};               // TC: singular
            std::vector<double> *q_vec{nullptr};

            void set_addresses(TTree& tree) {
                Blob::set_addresses(tree);
                tree.SetBranchAddress("q", &q_vec);
                tree.SetBranchAddress("time_slice", &time_slice_vec);
            }
        };
        TC live;

        struct TDC : Blob {
            std::vector<std::vector<int>> *time_slices_vec{nullptr}; // TDC: plural

            void set_addresses(TTree& tree) {
                Blob::set_addresses(tree);
                tree.SetBranchAddress("time_slice", &time_slices_vec);
            };
        };
        TDC dead;

        // Bad channels may override an activity map entry or delete it.
        class Bad {
            int chid{0}, plane{0};
            int tick_beg{0}, tick_end{0}; // one past?

            Waveform::ChannelMasks _masks;

        public:
            void set_addresses(TTree& tree) {
                tree.SetBranchAddress("chid", &chid);
                tree.SetBranchAddress("plane", &plane);
                tree.SetBranchAddress("start_time", &tick_beg);
                tree.SetBranchAddress("end_time", &tick_end);

                // derived
                _masks.clear();  // this is a one-shot fill
                const int nentries = tree.GetEntries();
                for (int entry = 0; entry<nentries; ++entry) {
                    tree.GetEntry(entry);
                    auto& brl = _masks[entry];
                    brl.emplace_back(tick_beg, tick_end);
                }
            }

            const Waveform::ChannelMasks& masks() const { return _masks; }

        };
        Bad bad;

    public:

        // Construct the interface to the trees.
        //
        // This does NOT load any entry.  Must call next() to load first entry, etc.
        UbooneBlobSourceTrees(const std::string& tfile_name, bool include_dead=false) {
            
            m_tfile.reset(TFile::Open(tfile_name.c_str(), "READ"));
            if (!m_tfile) {
                raise<IOError>("failed to open %s", tfile_name);
            }

            m_activity = open_ttree("Trun");
            m_nentries = m_activity->GetEntries();
            activity.set_addresses(*m_activity);

            m_live = open_ttree("TC");
            live.set_addresses(*m_live);

            // Needed for live or dead 2-view but optional
            m_bad = open_ttree("T_bad_ch", false);
            if (m_bad) {
                bad.set_addresses(*m_bad);
            }

            if (include_dead) {
                // only needed for kind=="dead"
                m_dead = open_ttree("TDC");
                dead.set_addresses(*m_dead);

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

        TTree* open_ttree(const std::string& tree_name, bool required = true) {
            auto ttree = reinterpret_cast<TTree*>(m_tfile->Get(tree_name.c_str()));
            if (!ttree && required) {
                raise<IOError>("failed to get TTree %s from %s", tree_name, m_tfile->GetName());
            }
            return ttree;
        }
    };
}


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


// Mark all channels spanned by the blob in a bodged plane as bodged 
void Root::UbooneBlobSource::bodge_activity(ISlice::map_t& activity, const RayGrid::Blob& blob)
{
    // hard-wire uboone's wire-to-channel map.
    const std::vector<int> choff = {0,2400,4800};
    const auto& strips = blob.strips();

    for (const int plane_index : m_bodged) {
        const int layer = plane_index + 2;
        const auto& [beg,end] = strips[layer].bounds;
        for (int wire = beg; wire<end; ++wire) {
            const int ch = wire + choff[plane_index];
            auto ich = get_channel(ch);
            activity[ich] = m_bodge;
        }        
    }
}

// Mark all channels in dummy planes.  Only call for 2-view
void Root::UbooneBlobSource::dummy_activity(ISlice::map_t& activity)
{
    for (auto ich : m_dummy->channels()) {
        activity[ich] = m_bodge;
        // fixme: MaskSlices in principle can use different values for "dummy"
        // and "masked" aka "bad" activity.
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


IBlobSet::pointer Root::UbooneBlobSource::load_live()
{
    // For 3-view, we accept live activity as-is.
    // For 2 view, we ignore live activity in 3rd view and replace it with "bad" activity.

    SliceMap slices;
    ISlice::pointer main_slice;

    // activity

    auto iframe = gen_frame();
    const auto& act = m_data->activity; // abbrev
    const auto& tids = act.timesliceId;
    const size_t nslices = tids->size();
    const int nrebin = m_data->header.nrebin;
    const double span = nrebin * m_tick;

    for (size_t sind=0; sind<nslices; ++sind) {

        const auto& q = act.raw_charge->at(sind);
        const auto& dq = act.raw_charge_err->at(sind);
        const auto& chans = act.timesliceChannel->at(sind);
        const size_t nchans = chans.size();
        int tsind = tids->at(sind);
        
        ISlice::map_t activity;
        for (size_t cind=0; cind<nchans; ++cind) {
            auto ichan = get_channel(cind);
            if (m_views < 7 && ichan->planeid().index() == m_bodged[0]) {
                continue;       // ignore real activity, will use bad CMM
            }
            activity[ichan] = ISlice::value_t(q[cind], dq[cind]);
        }

        if (m_views < 7) {      // 2-view
            for (const auto& [ch, brl] : m_data->bad.masks()) {
                auto ichan = get_channel(ch);
                if (ichan->planeid().index() != m_bodged[0]) {
                    continue;                    
                }
                // If slice overlaps with any of the bin ranges.
                for (const auto& tt : brl) {
                    if (tt.second < tsind*nrebin || tt.first > (tsind+1)*nrebin) continue;
                    activity[ichan] = m_bodge;
                }
            }
        }

        auto sslice = std::make_shared<SimpleSlice>(iframe, tsind, tsind*span, span, activity);
        slices[tsind] = sslice;
        if (!main_slice) main_slice = sslice;
    }            

    // blobs

    const auto& live_data = m_data->live;
    const RayGrid::Coordinates& coords = m_iface->raygrid();
    const size_t nblobs = live_data.cluster_id_vec->size();
    IBlob::vector iblobs;
    for (size_t bind=0; bind<nblobs; ++bind) {

        // Only consider blobs that match our configured views.
        if (! in_views(bind)) {
            continue;
        }

        RayGrid::Blob blob;
        blob.add(coords, RayGrid::Strip{0, {0,1}});
        blob.add(coords, RayGrid::Strip{1, {0,1}});
        blob.add(coords, RayGrid::Strip{2, make_strip(live_data.wire_index_u_vec->at(bind))});
        blob.add(coords, RayGrid::Strip{3, make_strip(live_data.wire_index_v_vec->at(bind))});
        blob.add(coords, RayGrid::Strip{4, make_strip(live_data.wire_index_w_vec->at(bind))});
            
        const int sid = live_data.time_slice_vec->at(bind);
        const float blob_charge = live_data.q_vec->at(bind);

        auto slice = slices[sid];

        if (m_views < 7) {
            // In principle, this is redundant with considering the CMM above as WCP
            // 2-view live blobs should reflect the contents of the T_bad_ch.
            bodge_activity(slice->activity(), blob);
        }

        auto iblob = std::make_shared<SimpleBlob>(bind, blob_charge, 0, blob, slice, m_iface);
        iblobs.push_back(iblob);
    }

    return std::make_shared<SimpleBlobSet>(iframe->ident(), main_slice, iblobs);
}


IBlobSet::pointer Root::UbooneBlobSource::load_dead()
{
    auto iframe = gen_frame();

    const auto& dead_data = m_data->dead;

    IBlob::vector iblobs;
    const size_t nblobs = dead_data.cluster_id_vec->size();
    if (!nblobs) {
        log->warn("no dead (TDC) blobs, producing empty blob set with no slice");
        return std::make_shared<SimpleBlobSet>(iframe->ident(), nullptr, iblobs);
    }

    SliceMap slices;
    const RayGrid::Coordinates& coords = m_iface->raygrid();

    // pick up first slice to give to the blob set.
    ISlice::pointer main_slice;

    for (size_t bind=0; bind<nblobs; ++bind) {

        if (! in_views(bind)) continue;

        const auto& tsvec = dead_data.time_slices_vec->at(bind);
        const int tsind = tsvec.front(); // assumes ordered....
        auto sit = slices.find(tsind);

        // Either already created or we make it
        SimpleSlicePtr slice = nullptr;
        if (sit == slices.end()) {

            // First blob in this dead slice, create the ISlice
            const double span = m_data->header.nrebin * m_tick * tsvec.size();
            const double start = tsind * span;

            ISlice::map_t activity;
            dummy_activity(activity);
            slice = std::make_shared<SimpleSlice>(iframe, tsind, start, span, activity);
        }
        else {
            slice = sit->second;
        }
        if (! main_slice) {
            main_slice = slice;
        }

        RayGrid::Blob blob;
        blob.add(coords, RayGrid::Strip{0, {0,1}});
        blob.add(coords, RayGrid::Strip{1, {0,1}});
        blob.add(coords, RayGrid::Strip{2, make_strip(dead_data.wire_index_u_vec->at(bind))});
        blob.add(coords, RayGrid::Strip{3, make_strip(dead_data.wire_index_v_vec->at(bind))});
        blob.add(coords, RayGrid::Strip{4, make_strip(dead_data.wire_index_w_vec->at(bind))});

        const float bval = 0;
        const float bunc = 0;

        bodge_activity(slice->activity(), blob);

        auto iblob = std::make_shared<SimpleBlob>(bind, bval, bunc, blob, slice, m_iface);
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

    log->debug("blob set {} with {} blobs at call", blobset->ident(), blobset->blobs().size(), ++m_calls);

    return true;
}

