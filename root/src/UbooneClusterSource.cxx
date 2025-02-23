#include "WireCellRoot/UbooneClusterSource.h"
#include "WireCellAux/SamplingHelpers.h"
#include "WireCellAux/TensorDMpointtree.h"

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

    m_angle_u = get(cfg, "angle_u", m_angle_u);
    m_angle_v = get(cfg, "angle_v", m_angle_v);
    m_angle_w = get(cfg, "angle_w", m_angle_w);
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


bool Root::UbooneClusterSource::flush(output_queue& outq)
{
    bool load_ok = m_files->next();
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
    }

    int ident = -1;


    // Collect all the IBlobs from all cached IBlobSets
    std::vector<IBlob::pointer> iblobs;
    IAnodeFace::pointer iface = nullptr;
    for (const auto& ibs : m_cache) {
        if (ident < 0) {
            ident = ibs->slice()->frame()->ident();
            log->debug("using ident {} from first blob set frame", ident);
        }
        // use the first iface we find
        if (!iface && ibs->blobs().size()) {
            iface = ibs->blobs().front()->face();
        }
        const auto& fresh = ibs->blobs();
        iblobs.insert(iblobs.end(), fresh.begin(), fresh.end());
    }

    size_t nublobs = blob_cluster_ids.size();
    size_t niblobs = iblobs.size();

    log->debug("blobs: ub={} ib={} in {} clusters", nublobs, niblobs, cnodes.size());
    if (nublobs != niblobs) {
        raise<ValueError>("blob count mismatch, job is malformed input gives %d, root file gives %d",
                          niblobs, nublobs);
    }

    //std::cout << "Test: " << niblobs << " " << nublobs << std::endl;

    const double tick = 500*units::ns;
    size_t n3dpoints_total = 0;
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

        named_pointclouds_t pcs;
        if (trees.is_live()) {
            auto [pc3d, aux] = m_sampler->sample_blob(iblob, bind);
            n3dpoints_total += pc3d.size_major();
            pcs.emplace("3d", pc3d);
            /// These seem unused and bring in horrible code
            pcs.emplace("2dp0", make2dds(pc3d, m_angle_u));
            pcs.emplace("2dp1", make2dds(pc3d, m_angle_v));
            pcs.emplace("2dp2", make2dds(pc3d, m_angle_w));
            const Point center = calc_blob_center(pcs["3d"]);
            auto scalar_ds = make_scalar_dataset(iblob, center, pcs["3d"].get("x")->size_major(), tick);
            int max_wire_interval = aux.get("max_wire_interval")->elements<int>()[0];
            int min_wire_interval = aux.get("min_wire_interval")->elements<int>()[0];
            int max_wire_type = aux.get("max_wire_type")->elements<int>()[0];
            int min_wire_type = aux.get("min_wire_type")->elements<int>()[0];
            scalar_ds.add("max_wire_interval", Array({(int)max_wire_interval}));
            scalar_ds.add("min_wire_interval", Array({(int)min_wire_interval}));
            scalar_ds.add("max_wire_type", Array({(int)max_wire_type}));
            scalar_ds.add("min_wire_type", Array({(int)min_wire_type}));
            pcs.emplace("scalar", std::move(scalar_ds));
        }
        else { // dead
            auto scalar_ds = make_scalar_dataset(iblob, {0,0,0}, 0, tick);
            scalar_ds.add("max_wire_interval", Array({(int)-1}));
            scalar_ds.add("min_wire_interval", Array({(int)-1}));
            scalar_ds.add("max_wire_type", Array({(int)-1}));
            scalar_ds.add("min_wire_type", Array({(int)-1}));
            pcs.emplace("scalar", scalar_ds);
            pcs.emplace("corner", make_corner_dataset(iblob));

        }
        cnode->insert(Points(std::move(pcs)));
    }
    log->debug("sampled {} points over {} blobs", n3dpoints_total, niblobs);
    

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

    std::string datapath = m_datapath;
    if (datapath.find("%") != std::string::npos) {
        datapath = String::format(datapath, ident);
    }

    Aux::add_ctpc(root, m_cache, iface, 0, m_time_offset, m_drift_speed);
    Aux::add_dead_winds(root, m_cache, iface, 0, m_time_offset, m_drift_speed);
    m_cache.clear();

    for (const auto& [name, pc] : root.value.local_pcs()) {
        log->debug("contains point cloud {} size_major {}", name, pc.size_major());
    }
    auto tens = as_tensors(root, datapath);

    log->debug("made pc-tree ncluster={} nblob={} nmatch={} in {} tensors at {} with ident {} in call {}",
               cnodes.size(), niblobs, nmatch, tens.size(), datapath, ident, m_calls);

    auto out = as_tensorset(tens, ident);
    outq.push_back(out);
    return true;
}

