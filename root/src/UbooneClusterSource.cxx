#include "WireCellRoot/UbooneClusterSource.h"
#include "WireCellAux/SamplingHelpers.h"
#include "WireCellAux/TensorDMpointtree.h"
#include "WireCellAux/SimpleBlob.h"
#include "WireCellAux/SimpleSlice.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/PointTree.h"

#include "TInterpreter.h"

#include <algorithm> // minmax


WIRECELL_FACTORY(UbooneClusterSource,
                 WireCell::Root::UbooneClusterSource,
                 WireCell::INamed,
                 WireCell::IBlobTensoring,
                 WireCell::IConfigurable)


using namespace WireCell;
using namespace WireCell::Aux;
using WireCell::PointCloud::Tree::Points;
using WireCell::PointCloud::Tree::named_pointclouds_t;
using WireCell::PointCloud::Dataset;
using WireCell::PointCloud::Array;
using WireCell::Aux::TensorDM::as_tensors;
using WireCell::Aux::TensorDM::as_tensorset;

Root::UbooneClusterSource::UbooneClusterSource()
    : Aux::Logger("UbooneClusterSource", "root")
    , m_calls(0)
{
}

Root::UbooneClusterSource::~UbooneClusterSource()
{
}

void Root::UbooneClusterSource::configure(const WireCell::Configuration& cfg)
{
    auto anode_tn = get<std::string>(cfg, "anode", "AnodePlane");
    m_anode = Factory::find_tn<IAnodePlane>(anode_tn);

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

    m_light_name = get<std::string>(cfg, "light", "");
    m_flash_name = get<std::string>(cfg, "flash", "");
    m_flashlight_name = get<std::string>(cfg, "flashlight", "");


    // A "live" kind is always implied in UboontTFiles so if user gives "dead",
    // this still does the right thing.
    std::vector<std::string> kinds = {get<std::string>(cfg, "kind", "live")};
    if (!m_light_name.empty() || !m_flash_name.empty() || m_flashlight_name.empty()) {
        kinds.push_back("light");
    }
    m_files = std::make_unique<UbooneTFiles>(input_paths, kinds, log);

    m_sampler.reset();
    auto sampler = get<std::string>(cfg, "sampler","");
    if (! sampler.empty()) {
        m_sampler = Factory::find_tn<IBlobSampler>(sampler);
    }
    // else {
    //     log->warn("no 'sampler' given, pc-tree will not have sampled points");
    // }
    m_datapath = get(cfg, "datapath", m_datapath);

    m_time_offset = get(cfg, "time_offset", m_time_offset);
    m_drift_speed = get(cfg, "drift_speed", m_drift_speed);

    
    m_angles = {
        get(cfg, "angle_u", m_angles[0]),
        get(cfg, "angle_v", m_angles[1]),
        get(cfg, "angle_w", m_angles[2])
    };
}

WireCell::Configuration Root::UbooneClusterSource::default_configuration() const
{
    Configuration cfg;
    cfg["input"] = Json::arrayValue; // required
    cfg["sampler"] = "";             // optional
    cfg["light"] = "";               // optional
    cfg["flash"] = "";               // optional
    cfg["flashlight"] = "";          // optional
    cfg["datapath"] = m_datapath;    // optional, default
    return cfg;
}


// dig out the frame ID
static int frame_ident(const IBlobSet::pointer& bs)
{
    return bs->slice()->frame()->ident();
}

bool Root::UbooneClusterSource::new_frame(const input_pointer& newbs) const
{
    if (m_cache.empty()) return false;
    return frame_ident(newbs) != frame_ident(m_cache[0]);
}


bool Root::UbooneClusterSource::operator()(const IBlobSet::pointer& in, output_queue& outq)
{
    // This flushes to the output queue on EOS or if the blobs' frame ID
    // changes.  A nullptr is appended to queue only on EOS.

    if (!in) {                  // eos
        bool ok = flush(outq);
        outq.push_back(nullptr); // forward eos
        log->debug("flush on eos at call {} okay:{}", m_calls, ok);
        ++m_calls;
        return ok;
    }

    if (new_frame(in)) {
        bool ok = flush(outq);
        log->debug("flush on new frame at call {} okay:{}", m_calls, ok);
        if (!ok) return ok;
    }

    m_cache.push_back(in);
    ++m_calls;
    return true;
}

IChannel::pointer Root::UbooneClusterSource::get_channel(int chanid) const
{
    auto ich = m_anode->channel(chanid);
    if (!ich) {
        log->error("No channel for ID {}, segfault to follow", chanid);
    }
    return ich;
}


// Fill slices and blobs from "dead" blob sets.  This simply extracts the blob
// and slice pointers from each blob set.  Output vectors are not cleared prior
// to filling.
void Root::UbooneClusterSource::extract_dead(
    const IBlobSet::vector& blobsets,
    ISlice::vector& out_slices,
    IBlob::vector& out_blobs) const
{
    for (const auto& ibs : blobsets) {
        const auto& fresh = ibs->blobs();
        out_blobs.insert(out_blobs.end(), fresh.begin(), fresh.end());
        out_slices.push_back(ibs->slice());
    }
}


// Map from time slice ID.  This MUST be used same as in UbooneBlobSource.
using SimpleSlicePtr = std::shared_ptr<Aux::SimpleSlice>;
using SliceMap = std::map<int, SimpleSlicePtr>;

// Fill slices and blobs from "live" blob sets.  This will produce new slices
// according to the ROOT data and new blobs same as old blobs that refer to
// corresponding new slices.  Output vectors are not cleared prior to filling.
void Root::UbooneClusterSource::extract_live(
    const IBlobSet::vector& blobsets,
    ISlice::vector& out_slices,
    IBlob::vector& out_blobs) const
{
    // Temporarily hold slices by their time slice id so we can look up later as
    // we patch in new slices to old blobs.
    SliceMap by_tid;

    // Get representative frame.  Note, this is likely a trace-free frame
    // manufactured by one of several upstream UbooneBlobSource's.  In
    // principle, right here we could replace it with a fully fleshed one built
    // from ROOT data.  But, it's a fully fleshed activity we currently seek to
    // keep we keep it as-is.
    auto iframe = blobsets[0]->slice()->frame(); 

    // WARNING: the code to fill the activity is ALMOST an exact copy-paste from
    // UbooneBlobSource.

    // Make all new slices, initially with no activity.
    const int n_slices_span = m_files->trees->nslices_span();
    const int nrebin = m_files->trees->header.nrebin;
    const double span = nrebin * m_tick;

    // Intermediate mutable version.
    std::vector<SimpleSlicePtr> slices(n_slices_span);
    for (int tsid = 0; tsid<n_slices_span; ++tsid) {
        auto sslice = std::make_shared<SimpleSlice>(iframe, tsid, tsid*span, span);
        slices[tsid] = sslice;
        // log->debug("simple slice: ident:{}, time:{}, span:{}, frame:{}",
        //            tsid, tsid*span, span, iframe->ident());
    }

    // Now fill in activity
    const auto& act = m_files->trees->activity; // abbrev
    const auto& tids = act.timesliceId;
    const int n_slices_data = m_files->trees->nslices_data();
    for (int sind=0; sind<n_slices_data; ++sind) {

        int tsid = tids->at(sind); // time slice ID

        const auto& q = act.raw_charge->at(sind);
        const auto& dq = act.raw_charge_err->at(sind);
        const auto& chans = act.timesliceChannel->at(sind);
        const size_t nchans = chans.size();
        
        auto& sslice = slices[tsid];
        if (!sslice) {
            raise<RuntimeError>("UbooneClusterSource given ROOT files inconsistent with input blobs");
        }
        auto& activity = sslice->activity();

        // Fill in known channel measures in this time slice.
        for (size_t cind=0; cind<nchans; ++cind) {
            const int chid = chans[cind];
            auto ichan = get_channel(chid);
            activity[ichan] = ISlice::value_t(q[cind], dq[cind]);
        }

        // Fill in known dead channels in this time slice
        for (const auto& [ch, brl] : m_files->trees->bad.masks()) {
            auto ichan = get_channel(ch);
            // If slice overlaps with any of the bin ranges.
            for (const auto& tt : brl) {
                if (tt.second < tsid*nrebin || tt.first > (tsid+1)*nrebin) continue;
                activity[ichan] = m_bodge;
            }
        }
    }
    out_slices.insert(out_slices.end(), slices.begin(), slices.end());

    // Now freshen the blobs with new slice and everything else copied from old.
    for (auto blobset : blobsets) {
        int tsid = blobset->ident(); // UbooneBlobSource must continue to follow this convention.
        for (auto old_blob : blobset->blobs()) {

            auto sslice = out_slices[tsid];
            if (!sslice) {
                raise<RuntimeError>("UbooneClusterSource given ROOT files inconsistent with input blobs");
            }

            out_blobs.push_back(
                std::make_shared<SimpleBlob>(
                    old_blob->ident(),
                    old_blob->value(),
                    old_blob->uncertainty(),
                    old_blob->shape(),
                    sslice,
                    old_blob->face()));
        }
    }
}


bool Root::UbooneClusterSource::flush(output_queue& outq)
{
    const bool load_ok = m_files->next();
    if (!load_ok) {
        log->error("failed to load uboone cluster event at call {}", m_calls);
        return false;
    }

    // The root node on which we grow the point tree.
    Points::node_t root;


    const auto& trees = *m_files->trees;

    // Create cluster nodes.  These start out empty/anonymous but we map them by
    // their uboone cluster ID for later personalizing.  This also puts them in
    // CLUSTER ID ORDER as defined by the UbooneTTrees.
    std::unordered_map<int, Points::node_t*> cnodes; 
    const auto& blob_cluster_ids = trees.blobs().cluster_ids();

    for (int cid : trees.cluster_ids) {
        auto* cnode = root.insert();
        cnodes[cid] = cnode;
        auto& spc = cnode->value.local_pcs()["cluster_scalar"];
        spc.add("flash", Array({(int)-1}));
        spc.add("ident", Array({cid}));

        apply_tagger_flags(cnode, cid);
    }

    // First, collect all the IBlobs and ISlices from all cached IBlobSets
    IBlob::vector iblobs;
    ISlice::vector islices;
    if (trees.is_live()) {
        extract_live(m_cache, islices, iblobs);
    }
    else {
        extract_dead(m_cache, islices, iblobs);
    }

    const int ident = islices[0]->frame()->ident();
    IAnodeFace::pointer iface = iblobs[0]->face();
    log->debug("is_live:{}, extracted {} slices, {} blobs from {} blob sets from frame ident {} face ident {} at call {}",
               trees.is_live(),
               
               islices.size(), iblobs.size(), m_cache.size(), ident, iface->ident(), m_calls);

    m_cache.clear();

    const size_t nublobs = blob_cluster_ids.size();
    const size_t niblobs = iblobs.size();

    log->debug("blobs: ub={} ib={} in {} clusters", nublobs, niblobs, cnodes.size());
    if (nublobs != niblobs) {
        raise<ValueError>("blob count mismatch, job is malformed input gives %d, root file gives %d",
                          niblobs, nublobs);
    }

    // Next dispatch the collected IBlobs and ISlices

    // First, dispatch the IBlobs
    const double tick = 500*units::ns;
    for (size_t bind=0; bind<niblobs; ++bind) {
        const IBlob::pointer iblob = iblobs[bind];

        // This relies on UbooneBlobSource to set the blob INDEX in the TTree
        // vectors to be the IBlob::ident().
        const int index = iblob->ident();

        const int cluster_id = blob_cluster_ids[index];
        auto cit = cnodes.find(cluster_id);
        if (cit == cnodes.end()) {
            raise<ValueError>("malformed job failed to find cluster node for cluster id %d", cluster_id);
        }
        auto* cnode = cit->second;

        // No sampling requested so we simply make an "empty" blob node.
        if (!m_sampler) {
            cnode->insert();
            continue;
        }

        if (trees.is_live()) {

            auto pcs = Aux::sample_live(m_sampler, iblob, m_angles, tick, bind);

            
            /// DO NOT EXTEND FURTHER! see #426, #430

            if (pcs.empty()) {
                SPDLOG_DEBUG("retile: skipping blob {} with no points", iblob->ident());
                continue;
            }

            cnode->insert(Points(std::move(pcs)));
        }
        else { // dead

            cnode->insert(Points(Aux::sample_dead(iblob, tick)));
            /// DO NOT EXTEND FURTHER! see #426, #430
        }
    }
    
    // Process individual cluster IDs for live data
    if (trees.is_live()) {
        const auto& individual_cluster_ids = trees.live.individual_cluster_ids();
        
        // Create isolated arrays for each cluster
        for (const auto& [cid, cnode] : cnodes) {
            auto& lpc = cnode->value.local_pcs();
            auto& pc = lpc["perblob"];
            
            std::vector<int> cluster_cc2;
            
            // Find blobs belonging to this cluster and assign group IDs
            for (size_t bind = 0; bind < niblobs; ++bind) {
                const int index = iblobs[bind]->ident();
                const int blob_cluster_id = blob_cluster_ids[index];
                
                if (blob_cluster_id == cid) {
                    const int individual_cluster_id = individual_cluster_ids[index];
                    
                    if (individual_cluster_id == cid) {
                        cluster_cc2.push_back(-1);  // Main cluster
                    } else {
                        cluster_cc2.push_back(individual_cluster_id);  // Sub-cluster
                    }
                }
            }
            
            // Add the isolated array to the cluster's PCTree
            if (!cluster_cc2.empty()) {
                PointCloud::Array::shape_t shape = {cluster_cc2.size()};
                pc.add("isolated", PointCloud::Array(cluster_cc2, shape, false));
            }
        }
    }

    // Add "flash"
    size_t nmatch=0;
    if (trees.is_live()) { 
        if (!m_light_name.empty()) {

            // This provides flash/light/flashlight arrays.
            root.value.local_pcs() = std::move(trees.optical);

            const auto& cf = trees.cluster_flash;
            nmatch = cf.size();

            for ( const auto& [cid, find] : cf) {
                auto* cnode = cnodes[cid];
                auto& spc = cnode->value.local_pcs()["cluster_scalar"];
                auto farr = spc.get("flash"); // initially set undefined/-1 above
                farr->element<int>(0) = find;

                //  std::cout << "Test: " << cid << " " << find << " " << farr->element<float>(1) << std::endl;
            }
        }
    }

    // Dispatch the ISlices
    if (trees.is_live()){
        Aux::add_ctpc(root, islices, iface, 0, m_time_offset, m_drift_speed);
        Aux::add_dead_winds(root, islices, iface, 0, m_time_offset, m_drift_speed);
    }


    for (const auto& [name, pc] : root.value.local_pcs()) {
        log->debug("contains point cloud {} size_major {}", name, pc.size_major());
    }

    // Serialize to tensors
    std::string datapath = m_datapath;
    if (datapath.find("%") != std::string::npos) {
        datapath = String::format(datapath, ident);
    }
    auto tens = as_tensors(root, datapath);
    log->debug("made pc-tree ncluster={} nblob={} nmatch={} in {} tensors at {} with ident {} in call {}",
               cnodes.size(), niblobs, nmatch, tens.size(), datapath, ident, m_calls);
    // for (auto& ten : tens) {
    //     log->debug("{}", ten->metadata()["datapath"]);
    // }

    auto out = as_tensorset(tens, ident);
    outq.push_back(out);
    return true;
}

void Root::UbooneClusterSource::apply_tagger_flags(WireCell::PointCloud::Tree::Points::node_t* cnode, int cluster_id) const
{
    if (!m_files || !m_files->trees) {
        return;
    }

    const auto& trees = *m_files->trees;
    
    if (!trees.is_live()) {
        return;
    }

    // Store tagger flag information in point cloud for later transfer to flags
    auto& tagger_pc = cnode->value.local_pcs()["tagger_info"];
    
    // Check beam flash coincidence
    bool beam_flash_coincident = trees.is_beam_flash_coincident(cluster_id);
    
    if (beam_flash_coincident) {
        // Store that this cluster should have beam_flash flag
        tagger_pc.add("has_beam_flash", PointCloud::Array({1}));
        
        // Get event type and extract individual flags
        int event_type = trees.get_event_type(cluster_id);
        // double cluster_length = trees.get_cluster_length(cluster_id);
        
        // Extract flags from event_type using bit operations
        int flag_tgm = (event_type >> 3) & 1U;
        int flag_low_energy = (event_type >> 4) & 1U;
        int flag_lm = (event_type >> 1) & 1U;
        int flag_fully_contained = (event_type >> 2) & 1U;
        int flag_stm = (event_type >> 5) & 1U;
        int flag_full_detector_dead = (event_type >> 6) & 1U;
        
        // // Store all flags in a single array to avoid potential point cloud issues
        // std::cout << "Xin: Adding flags array for cluster " << cluster_id << std::endl;
        // std::vector<int> flag_values = {flag_tgm, flag_low_energy, flag_lm, flag_fully_contained, flag_stm, flag_full_detector_dead};
        // tagger_pc.add("tagger_flags", PointCloud::Array(flag_values));
        
        // Also store individual flags for backward compatibility
        tagger_pc.add("has_tgm", PointCloud::Array({flag_tgm}));
        tagger_pc.add("has_low_energy", PointCloud::Array({flag_low_energy}));
        tagger_pc.add("has_light_mismatch", PointCloud::Array({flag_lm}));
        tagger_pc.add("has_fully_contained", PointCloud::Array({flag_fully_contained}));
        tagger_pc.add("has_short_track_muon", PointCloud::Array({flag_stm}));
        tagger_pc.add("has_full_detector_dead", PointCloud::Array({flag_full_detector_dead}));
        // std::cout << "Xin: Added all flags for cluster " << cluster_id << std::endl;
        
        // Store metadata
        // tagger_pc.add("event_type", PointCloud::Array({event_type}));
        // tagger_pc.add("cluster_length", PointCloud::Array({cluster_length}));
        
        // Debug: List all keys we just added
        // std::cout << "Xin: Final tagger_pc keys for cluster " << cluster_id << ": ";
        // for (const auto& key : tagger_pc.keys()) {
        //     std::cout << key << " ";
        // }
        // std::cout << std::endl;
        
        log->debug("Stored tagger flag data for cluster {}: beam_flash=1, event_type={}", 
                   cluster_id, event_type);

        // std::cout << "Xin: Cluster " << cluster_id 
        //           << " has beam flash coincident with event type: " << event_type 
        //           << ", cluster length: " << cluster_length << " " << flag_tgm << " " << flag_low_energy << " " << flag_lm << " " << flag_fully_contained << " " << flag_stm << " " << flag_full_detector_dead << std::endl;
    } else {
        // Store that this cluster should NOT have beam_flash flag
        tagger_pc.add("has_beam_flash", PointCloud::Array({0}));

        // Also store individual flags for backward compatibility
        tagger_pc.add("has_tgm", PointCloud::Array({0}));
        tagger_pc.add("has_low_energy", PointCloud::Array({0}));
        tagger_pc.add("has_light_mismatch", PointCloud::Array({0}));
        tagger_pc.add("has_fully_contained", PointCloud::Array({0}));
        tagger_pc.add("has_short_track_muon", PointCloud::Array({0}));
        tagger_pc.add("has_full_detector_dead", PointCloud::Array({0}));
        
        // Store metadata
        // tagger_pc.add("event_type", PointCloud::Array({0}));
        // tagger_pc.add("cluster_length", PointCloud::Array({0}));

        log->debug("Cluster {} not beam coincident", cluster_id);
    }
}
