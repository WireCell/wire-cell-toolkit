#include "WireCellRoot/UbooneBlobSource.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleBlob.h"
#include "WireCellAux/BlobTools.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Persist.h"
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

    auto check_file = [](const auto& jstr) -> std::string {
        auto got = Persist::resolve(jstr.asString());
        if (got.empty()) {
            raise<ValueError>("file not found: %s", jstr.asString());
        }
        return got;
    };

    auto input = cfg["input"];
    if (input.isNull()) {
        log->critical("input is required");
        raise<ValueError>("UbooneBlobSource: input is required");
    }
    std::vector<std::string> input_paths;
    if (input.isString()) {
        input_paths.push_back(check_file(input));
    }
    else {
        for (const auto& one: input) {
            input_paths.push_back(check_file(one));
        }
    }

    m_kind = get(cfg, "kind", m_kind);
    std::vector<std::string> kinds = {m_kind};
    m_files = std::make_unique<UbooneTFiles>(input_paths, kinds, log);


    // See in_views().
    m_views = 0;
    auto jviews = cfg["views"];
    if (jviews.isNull()) {
        if (m_kind == "live") {
            m_views = 7;        // uvw
        }
        else {
            m_views = 3;        // uv
        }
    }
    else {
        std::string views = jviews.asString();
        for (char letter : views) {
            if (letter > 'Z') letter -= 32; // upper case
            if (letter == 'U') m_views |= 1;
            if (letter == 'V') m_views |= 2;
            if (letter == 'W') m_views |= 4;
        }
    }

    // Work out which planes get "bodged" (set as "bad"/"dead")
    // and maybe set as "dummy" for dead.  All this is only for 2-view.
    if (m_kind == "live") {
        if (m_views == 2+4) m_bodged = {0}; // U
        if (m_views == 4+1) m_bodged = {1}; // V
        if (m_views == 1+2) m_bodged = {2}; // W
    }
    else {
        auto iwplanes = m_iface->planes();

        if (m_views == 2+4) {
            m_bodged = {1,2}; // !U
            m_dummy = iwplanes[0];
        }
        if (m_views == 4+1) {
            m_bodged = {2,0}; // !V
            m_dummy = iwplanes[1];
        }
        if (m_views == 1+2) {
            m_bodged = {0,1}; // !W
            m_dummy = iwplanes[2];
        }
    }

    m_frame_eos = get(cfg, "frame_eos", m_frame_eos);

    log->debug("loading {} blobs of views bits: {}", m_kind, m_views);
}

WireCell::Configuration Root::UbooneBlobSource::default_configuration() const
{
    Configuration cfg;
    cfg["input"] = Json::arrayValue;
    cfg["kind"] = m_kind;
    cfg["views"] = Json::arrayValue;
    return cfg;
}




IFrame::pointer Root::UbooneBlobSource::gen_frame()
{
    // We have no good way to make traces.  Best we might do is spread activity
    // over nrebin ticks.  Assuming the traces are not important for downstream
    // we make a frame that only holds metadata.
    int ident = ((m_files->trees->header.runNo & 0x7fff) << 16) | (m_files->trees->header.eventNo & 0xffff);
    return std::make_shared<SimpleFrame>(ident,
                                         m_files->trees->header.triggerTime, 
                                         0.5*units::us);
}


IChannel::pointer Root::UbooneBlobSource::get_channel(int chanid)
{
    auto ich = m_anode->channel(chanid);
    if (!ich) {
        log->error("No channel for ID {}, segfault to follow", chanid);
    }
    return ich;
}


bool Root::UbooneBlobSource::in_views(int bind)
{
    // Inclusion test depends if we are loading live or dead blobs.
    const int want = m_kind == "live" ? 1 : 0;
    const auto& flag_uvw = m_files->trees->blob(m_kind).flag_uvw;

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

// Need to keep the concrete type around as we build
using SimpleSlicePtr = std::shared_ptr<Aux::SimpleSlice>;
using SliceMap = std::map<int, SimpleSlicePtr>;
using SimpleBlobsetPtr = std::shared_ptr<Aux::SimpleBlobSet>;
using BlobsetMap = std::map<int, SimpleBlobsetPtr>;

void Root::UbooneBlobSource::load_live()
{
    auto iframe = gen_frame();
    const auto& act = m_files->trees->activity; // abbrev
    const auto& tids = act.timesliceId;
    const int n_slices_data = m_files->trees->nslices_data();
    const int n_slices_span = m_files->trees->nslices_span();
    const int nrebin = m_files->trees->header.nrebin;
    const double span = nrebin * m_tick;

    // premake
    BlobsetMap blobsets;
    SliceMap slices;
    for (int tsid = 0; tsid<n_slices_span; ++tsid) {
        // log->debug("SimpleSlice: {} start={} span={} nrebin={} tick={} frame={}",
        //            tsid, tsid*span, span, nrebin, m_tick, iframe->ident());
        auto sslice = std::make_shared<SimpleSlice>(iframe, tsid, tsid*span, span);
        slices[tsid] = sslice;
        blobsets[tsid] = std::make_shared<SimpleBlobSet>(tsid, sslice);
    }

    // activity
    std::map<int, std::set<IChannel::pointer> > dead_channels;
    for (int sind=0; sind<n_slices_data; ++sind) {

        int tsid = tids->at(sind); // time slice ID

        const auto& q = act.raw_charge->at(sind);
        const auto& dq = act.raw_charge_err->at(sind);
        const auto& chans = act.timesliceChannel->at(sind);
        const size_t nchans = chans.size();
        
        auto sslice = slices[tsid];
        auto& activity = sslice->activity();
        for (size_t cind=0; cind<nchans; ++cind) {
            const int chid = chans[cind];
            auto ichan = get_channel(chid);
            if (m_views < 7 && ichan->planeid().index() == m_bodged[0]) {
                continue;       // ignore real activity, will use bad CMM
            }
            activity[ichan] = ISlice::value_t(q[cind], dq[cind]);
        }

        if (m_views < 7) {      // 2-view
            for (const auto& [ch, brl] : m_files->trees->bad.masks()) {
                auto ichan = get_channel(ch);
                if (ichan->planeid().index() != m_bodged[0]) {
                    continue;                    
                }
                dead_channels[ichan->planeid().index()].insert(ichan);
                // If slice overlaps with any of the bin ranges.
                for (const auto& tt : brl) {
                    if (tt.second < tsid*nrebin || tt.first > (tsid+1)*nrebin) continue;
                    activity[ichan] = m_bodge;
                }
            }
        }
    }
    for (const auto& [pind, chans] : dead_channels) {
        log->debug("plane {} has {} dead channels", pind, chans.size());
    }

    // blobs
    const auto& live_data = m_files->trees->live;
    const RayGrid::Coordinates& coords = m_iface->raygrid();
    const size_t nblobs = live_data.cluster_id_vec->size();

    //std::cout << "Test: " << nblobs << std::endl;

    size_t n_blobs_loaded = 0;
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
            
        const float blob_charge = live_data.q_vec->at(bind);
        const int tsid = live_data.time_slice_vec->at(bind);

        auto bset = blobsets[tsid];
        auto sslice = slices[tsid];

        if (m_views < 7) {
            // In principle, this is redundant with considering the CMM above as WCP
            // 2-view live blobs should reflect the contents of the T_bad_ch.
            bodge_activity(sslice->activity(), blob);
        }

        bset->insert(std::make_shared<SimpleBlob>(bind, blob_charge, 0, blob, sslice, m_iface));
        ++n_blobs_loaded;
    }
    for (const auto& [_, bs] : blobsets) {
        m_queue.push_back(bs);
    }
    if (m_frame_eos) {
        m_queue.push_back(nullptr);
    }

    log->debug("live: loaded {} blobs in {} sets from entry {} of {}",
               n_blobs_loaded, blobsets.size(), m_files->trees->entry(), m_files->trees->nentries());
}


void Root::UbooneBlobSource::load_dead()
{
    auto iframe = gen_frame();

    const auto& dead_data = m_files->trees->dead;

    SliceMap slices;
    BlobsetMap blobsets;
    const RayGrid::Coordinates& coords = m_iface->raygrid();
    const size_t nblobs = dead_data.cluster_id_vec->size();
    size_t n_blobs_loaded = 0;

    // std::cout << "Test: " << nblobs << std::endl;

    for (size_t bind=0; bind<nblobs; ++bind) {

        if (! in_views(bind)) continue;

        const auto& tsvec = dead_data.time_slices_vec->at(bind);
        const int tsid = tsvec.front(); // assumes ordered....
        auto sit = slices.find(tsid);

        // Either already created or we make it
        SimpleSlicePtr slice = nullptr;
        SimpleBlobsetPtr bset = nullptr;
        if (sit == slices.end()) {

            // First blob in this dead slice, create the ISlice
            const double live_span = m_files->trees->header.nrebin * m_tick;
            const double start = live_span * tsid;
            const double span = live_span * tsvec.size();


            slice = std::make_shared<SimpleSlice>(iframe, tsid, start, span);
            dummy_activity(slice->activity());
            slices[tsid] = slice;
            blobsets[tsid] = bset = std::make_shared<SimpleBlobSet>(tsid, slice);
        }
        else {
            slice = sit->second;
            bset = blobsets[tsid];
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

        bset->insert(std::make_shared<SimpleBlob>(bind, bval, bunc, blob, slice, m_iface));
        ++n_blobs_loaded;
    }

    for (const auto& [_, bs] : blobsets) {
        m_queue.push_back(bs);
    }
    if (m_frame_eos) {
        m_queue.push_back(nullptr);
    }

    log->debug("dead: loaded {} blobs in {} sets from entry {} of {}",
               n_blobs_loaded, blobsets.size(), m_files->trees->entry(), m_files->trees->nentries());

}



// Try to fill the queue, or don't.
void Root::UbooneBlobSource::fill_queue()
{
    if (! m_files->next() ) { return; }

    if (m_kind == "live") {
        load_live();
    }
    else {
        load_dead();
    }
    if (m_frame_eos) {
        m_queue.push_back(nullptr);
    }
}

bool Root::UbooneBlobSource::operator()(IBlobSet::pointer& blobset)
{
    blobset = nullptr;

    if (m_done) {
        // log->debug("past EOS at call {}, stop calling me", m_calls++);
        return false;
    }

    if (! m_files) {
        return false;
    }

    if (m_queue.empty()) {
        fill_queue();
    }

    if (m_queue.empty()) {
        log->debug("EOS due to input exhaustion at call={}", m_calls++);
        m_done = true;          // next time we get angry.
        return true;
    }

    blobset = m_queue.front();
    m_queue.pop_front();

    if (! blobset) {
        log->debug("EOS due to frame end at call={}", m_calls++);
    }
    else {
        log->trace("blob set call={}: {}", m_calls++, dumps(blobset));
    }


    return true;
}

