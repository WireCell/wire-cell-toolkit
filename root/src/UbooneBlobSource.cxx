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
    auto jviews = cfg["views"];
    if (jviews.isArray()) {
        m_views.clear();
        for (const auto& jone : jviews) {
            std::string view = jone.asString();
            char cone=0;
            for (char letter : view) {
                if (letter > 'Z') letter -= 32; // upper case
                if (letter == 'U') cone |= 1;
                if (letter == 'V') cone |= 2;
                if (letter == 'W') cone |= 4;
            }
            if (cone) {
                m_views.push_back(cone);
            }
        }
    }
    else {
        if (m_kind == "TC") {
            m_views = {3,5,6,7}; // 3x2-way + 3-way
        }
        else {
            m_views = {3,5,6};  // 3x2-way
        }
    }

    // ROOT does not always have this one predefined
    gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");
}

WireCell::Configuration Root::UbooneBlobSource::default_configuration() const
{
    Configuration cfg;
    cfg["input"] = Json::arrayValue;
    cfg["kind"] = m_kind;
    cfg["views"] = Json::arrayValue;
    return cfg;
}


/**
   Bags to hold the branch addresses relevant to building an IBlobSet.
*/

namespace WireCell::Root {

    // This hides TFile/TTree types.
    class UbooneBlobSourceTrees {

        std::unique_ptr<TFile> m_tfile;
        TTree* m_trun;
        TTree* m_tc;
        Long64_t m_entry{-1};   // not yet loaded.  Only use next()
        Long64_t m_nentries{0};


    public:

        // The data is available after construction and refreshed on each call
        // to next().  The data is invalid when/if next() throws IndexError.

        // "header"
        // double eventTime{0};    // 1.45767e9. Unix time?
        double triggerTime{0};  // 1.48822e10. Units?
        int eventNo{0};
        int nrebin{0};

        // "activity"
        std::vector<int> *timesliceId{nullptr}; // [0,2399]
        std::vector<std::vector<int>> *timesliceChannel{nullptr}; // [0,8255]
        // Units are number of signal electrons.
        std::vector<std::vector<int>> *raw_charge{nullptr};
        std::vector<std::vector<int>> *raw_charge_err{nullptr};

        // "slice + blob"
        std::vector<int> *cluster_id_vec{nullptr};
        std::vector<int> *time_slice_vec{nullptr};               // TC: singular
        std::vector<std::vector<int>> *time_slices_vec{nullptr}; // TDC: plural

        // 0 value indicates the view for a blob is "dead"
        std::vector<int> *flag_u_vec{nullptr};
        std::vector<int> *flag_v_vec{nullptr};
        std::vector<int> *flag_w_vec{nullptr};

        std::vector<std::vector<int>> *wire_index_u_vec{nullptr};
        std::vector<std::vector<int>> *wire_index_v_vec{nullptr};
        std::vector<std::vector<int>> *wire_index_w_vec{nullptr};

        // Lord, forgive us our sins.
        std::vector<std::vector<int>*> flag_uvw;

        // Only gets defined if "q" exists.  ie, for TC not TDC.
        std::vector<double> *q_vec{nullptr};

        bool is_tc{true};       // true, if our data is from TC.  else from TDC

    public:

        UbooneBlobSourceTrees(const std::string& tfile_name,
                              const std::string& tc_name = "TC",
                              const std::string& trun_name = "Trun") {
            
            m_tfile.reset(TFile::Open(tfile_name.c_str(), "READ"));
            if (!m_tfile) {
                raise<IOError>("failed to open %s", tfile_name);
            }
            m_tc = open_ttree(tc_name);
            m_trun = open_ttree(trun_name);
            m_nentries = m_trun->GetEntries();

            // "header"
            m_trun->SetBranchAddress("triggerTime",&triggerTime); 
            m_trun->SetBranchAddress("eventNo",&eventNo);
            m_trun->SetBranchAddress("nrebin",&nrebin);

            // "activity"
            m_trun->SetBranchAddress("timesliceId",&timesliceId);
            m_trun->SetBranchAddress("timesliceChannel",&timesliceChannel);
            m_trun->SetBranchAddress("raw_charge",&raw_charge);
            m_trun->SetBranchAddress("raw_charge_err",&raw_charge_err);

            // "blob"
            m_tc->SetBranchAddress("cluster_id",&cluster_id_vec); 
            if (m_tc->GetBranch("q")) { // this is TC
                is_tc = true;
                m_tc->SetBranchAddress("q",&q_vec);
                m_tc->SetBranchAddress("time_slice",&time_slice_vec);
            }
            else {
                is_tc = false;  // is TDC
                m_tc->SetBranchAddress("time_slice",&time_slices_vec);
            }
            m_tc->SetBranchAddress("flag_u",&flag_u_vec);
            m_tc->SetBranchAddress("flag_v",&flag_v_vec);
            m_tc->SetBranchAddress("flag_w",&flag_w_vec);
            flag_uvw = {flag_u_vec,flag_v_vec,flag_w_vec};

            m_tc->SetBranchAddress("wire_index_u",&wire_index_u_vec);
            m_tc->SetBranchAddress("wire_index_v",&wire_index_v_vec);
            m_tc->SetBranchAddress("wire_index_w",&wire_index_w_vec);
        }

        void next() {
            ++m_entry;
            if (m_entry < m_nentries) { 
                m_trun->GetEntry(m_entry);
                m_tc->GetEntry(m_entry);
                return;
            }
            raise<IndexError>("attempt to get entry past end of TTree");
        }

    private:

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
    if (m_data) {
        // We have an open file.
        try {
            m_data->next();
        }
        catch (IndexError& err) {
            m_data = nullptr;
            return next();
        }
        return true;
    }

    // We are starting a new file.

    if (m_input.empty()) {
        // We have exhausted our input files.
        return false;
    }

    auto fname = m_input.back();
    m_input.pop_back();

    try {
        m_data = std::make_unique<UbooneBlobSourceTrees>(fname, m_kind);
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

    log->debug("reading {}", fname);
    return true;
}


IFrame::pointer Root::UbooneBlobSource::gen_frame()
{
    // We have no way to make traces so frame holds just metadata.
    return std::make_shared<SimpleFrame>(m_data->eventNo,
                                         m_data->triggerTime, 
                                         0.5*units::us);
}


IChannel::pointer Root::UbooneBlobSource::get_channel(int chanid)
{
    return m_anode->channel(chanid);
}

bool Root::UbooneBlobSource::in_views(int bind)
{
    // The flag value for view inclusion.  Live wants flag=1 (good) and dead
    // wants flag=0 (bad).
    const int want = m_data->is_tc ? 1 : 0;

    int cone = 0;
    for (int pln=0; pln<3; ++pln) {
        if (want == m_data->flag_uvw[pln]->at(bind)) {
            cone |= 1<<pln;
        }
    }

    for (char one : m_views) {
        if (cone == one) return true;
    }
    return false;
}

// If any view of the blob is "flagged" 0 then set its channels activity to special "bad" value. 
void Root::UbooneBlobSource::bodge_channels(std::shared_ptr<SimpleSlice> slice, const RayGrid::Blob& blob, int bind)
{
    const std::vector<int> choff = {0,2400,4800};
    const auto& strips = blob.strips();
    for (int pln = 0; pln < 3; ++pln) {
        if (m_data->flag_uvw[pln]->at(bind) == 1) {
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

RayGrid::Blob Root::UbooneBlobSource::make_blob(int bind)
{
    RayGrid::Blob blob;
    const RayGrid::Coordinates& coords = m_iface->raygrid();
    blob.add(coords, RayGrid::Strip{0, {0,1}});
    blob.add(coords, RayGrid::Strip{1, {0,1}});
    blob.add(coords, RayGrid::Strip{2, make_strip(m_data->wire_index_u_vec->at(bind))});
    blob.add(coords, RayGrid::Strip{3, make_strip(m_data->wire_index_v_vec->at(bind))});
    blob.add(coords, RayGrid::Strip{4, make_strip(m_data->wire_index_w_vec->at(bind))});
    return blob;
}

IBlobSet::pointer Root::UbooneBlobSource::load_live()
{
    auto iframe = gen_frame();

    // Map prototype slice index to toolkit ISlice.
    std::unordered_map<int, std::shared_ptr<SimpleSlice>> slices;
    // A representative slice for the IBlobSet
    ISlice::pointer main_slice = nullptr;
    
    const auto& tids = m_data->timesliceId;
    const size_t nslices = tids->size();
    const double span = m_data->nrebin * 0.5 * units::us;
    for (size_t sind=0; sind<nslices; ++sind) {

        const auto& q = m_data->raw_charge->at(sind);
        const auto& dq = m_data->raw_charge_err->at(sind);
        const auto& chans = m_data->timesliceChannel->at(sind);
        const size_t nchans = chans.size();
        
        ISlice::map_t activity;
        for (size_t cind=0; cind<nchans; ++cind) {
            auto ichan = get_channel(cind);
            activity[ichan] = ISlice::value_t(q[cind], dq[cind]);
        }
        int tsind = tids->at(sind);
        auto sslice = std::make_shared<SimpleSlice>(iframe, tsind, tsind*span, span, activity);
        slices[tsind] = sslice;
        if (!main_slice) main_slice = sslice;
    }            

    // Blobs
    const size_t nblobs = m_data->cluster_id_vec->size();
    IBlob::vector iblobs;
    for (size_t bind=0; bind<nblobs; ++bind) {

        if (! in_views(bind)) continue;

        const RayGrid::Blob blob = make_blob(bind);
            
        const int sid = m_data->time_slice_vec->at(bind);
        const float blob_charge = m_data->q_vec->at(bind);

        auto slice = slices[sid];
        bodge_channels(slice, blob, bind);

        auto iblob = std::make_shared<SimpleBlob>(bind, blob_charge, 0, blob, slice, m_iface);
        iblobs.push_back(iblob);
    }
    return std::make_shared<SimpleBlobSet>(iframe->ident(), main_slice, iblobs);
}

IBlobSet::pointer Root::UbooneBlobSource::load_dead()
{
    log->warn("load_dead needs also load_live to provide valid charge in activity map");

    auto iframe = gen_frame();

    // Map initial prototype slice index to toolkit ISlice.
    std::unordered_map<int, std::shared_ptr<SimpleSlice>> slices;
    // A representative slice for the IBlobSet
    ISlice::pointer main_slice = nullptr;

    const size_t nblobs = m_data->cluster_id_vec->size();
    if (!nblobs) {
        log->warn("no dead (TDC) blobs");
        // fall through the following loop, will return valid but empty blobset.
    }

    IBlob::vector iblobs;
    for (size_t bind=0; bind<nblobs; ++bind) {

        if (! in_views(bind)) continue;

        const auto& tsvec = m_data->time_slices_vec->at(bind);
        const int tsind = tsvec.front(); // assumes ordered....
        auto sit = slices.find(tsind);

        // Either created or looked up
        std::shared_ptr<SimpleSlice> slice = nullptr;

        if (sit == slices.end()) {
            // first seen, create ISlice
            const double span = m_data->nrebin * 0.5 * units::us * tsvec.size();
            const double start = tsind * span;
            slice = std::make_shared<SimpleSlice>(iframe, tsind, start, span);
            // FIXME: we leave the activity map empty!  Does downstream care?
        }
        else {
            slice = sit->second;
        }

        if (! main_slice) {
            main_slice = slice;
        }

        const RayGrid::Blob blob = make_blob(bind);
        const float bval = 0;
        const float bunc = 0;

        bodge_channels(slice, blob, bind);

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

    if (m_data->is_tc) {
        blobset = load_live();
    }
    else {
        blobset = load_dead();
    }

    return true;
}

