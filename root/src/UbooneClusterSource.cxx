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

    std::vector<std::string> kinds = {"live"};
    if (!m_light_name.empty() || !m_flash_name.empty() || m_flashlight_name.empty()) {
        kinds.push_back("light");
    }
    
    m_files = std::make_unique<UbooneTFiles>(input_paths, kinds, log);

    m_sampler.reset();
    auto sampler = get<std::string>(cfg, "sampler","");
    if (! sampler.empty()) {
        m_sampler = Factory::find_tn<IBlobSampler>(sampler);
    }

    m_datapath = get(cfg, "datapath", m_datapath);

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


    // Create cluster nodes.  These start out empty/anonymous but we map them by
    // their uboone cluster ID for later personalizing.  This also puts them in
    // CLUSTER ID ORDER as defined by the UbooneTTrees.
    std::unordered_map<int, Points::node_t*> cnodes; 
    for (int cid : m_files->trees->cluster_ids) {
        auto* cnode = root.insert();
        cnodes[cid] = cnode;
        auto& spc = cnode->value.local_pcs()["cluster_scalar"];
        spc.add("flash", Array({(int)-1}));
        spc.add("cluster_id", Array({cid}));
    }

    int ident = -1;


    // Collect all the IBlobs from all cached IBlobSets
    std::vector<IBlob::pointer> iblobs;
    for (const auto& ibs : m_cache) {
        if (ident < 0) {
            ident = ibs->slice()->frame()->ident();
            log->debug("using ident {} from first blob set frame", ident);
        }
        const auto& fresh = ibs->blobs();
        iblobs.insert(iblobs.end(), fresh.begin(), fresh.end());
    }
    m_cache.clear();

    // From uboone TTrees
    const auto& tblob = m_files->trees->live;
    const auto& blob_cids = *tblob.cluster_id_vec; // spans blobs in "event"

    size_t nublobs = blob_cids.size();
    size_t niblobs = iblobs.size();

    log->debug("blobs: ub={} ib={} in {} clusters", nublobs, niblobs, cnodes.size());
    if (nublobs != niblobs) {
        raise<ValueError>("blob count mismatch, job is malformed input gives %d, root file gives %d",
                          niblobs, nublobs);
    }

    //std::cout << "Test: " << niblobs << " " << nublobs << std::endl;

    for (size_t bind=0; bind<niblobs; ++bind) {
        const IBlob::pointer iblob = iblobs[bind];
        // This MUST be the TTree entry number as set by UbooneBlobSource!
        const int entry = iblob->ident(); // HUGE TRUST HERE!!!
        const int cluster_id = blob_cids[entry];
        auto cit = cnodes.find(cluster_id);
        if (cit == cnodes.end()) {
            raise<ValueError>("malformed job");
        }
        auto* cnode = cit->second;

        if (m_sampler) {
            named_pointclouds_t pcs;
            auto [pc3d, aux] = m_sampler->sample_blob(iblob, bind);
            pcs.emplace("3d", pc3d);
            /// These seem unused and bring in horrible code
            // pcs.emplace("2dp0", make2dds(pc3d, angle_u));
            // pcs.emplace("2dp1", make2dds(pc3d, angle_v));
            // pcs.emplace("2dp2", make2dds(pc3d, angle_w));
            const Point center = calc_blob_center(pcs["3d"]);
            auto scaler_ds = make_scaler_dataset(iblob, center, pcs["3d"].get("x")->size_major(), 500*units::ns);
            int max_wire_interval = aux.get("max_wire_interval")->elements<int>()[0];
            int min_wire_interval = aux.get("min_wire_interval")->elements<int>()[0];
            int max_wire_type = aux.get("max_wire_type")->elements<int>()[0];
            int min_wire_type = aux.get("min_wire_type")->elements<int>()[0];
            scaler_ds.add("max_wire_interval", Array({(int)max_wire_interval}));
            scaler_ds.add("min_wire_interval", Array({(int)min_wire_interval}));
            scaler_ds.add("max_wire_type", Array({(int)max_wire_type}));
            scaler_ds.add("min_wire_type", Array({(int)min_wire_type}));
            pcs.emplace("scalar", std::move(scaler_ds));

            cnode->insert(Points(std::move(pcs)));
        }
    }
    

    size_t nmatch=0;
    if (!m_light_name.empty()) {
        root.value.local_pcs() = std::move(m_files->trees->optical);

        const auto& cf = m_files->trees->cluster_flash;
        nmatch = cf.size();

        for ( const auto& [cid, find] : cf) {
            auto* cnode = cnodes[cid];
            auto& spc = cnode->value.local_pcs()["cluster_scalar"];
            auto farr = spc.get("flash"); // initially set undefined/-1 above
            farr->element<int>(0) = find;

          //  std::cout << "Test: " << cid << " " << find << " " << farr->element<float>(1) << std::endl;
        }

    }

    std::string datapath = m_datapath;
    if (datapath.find("%") != std::string::npos) {
        datapath = String::format(datapath, ident);
    }
    auto tens = as_tensors(root, datapath);

    log->debug("made pc-tree ncluster={} nblob={} nmatch={} in {} tensors at {} with ident {} in call {}",
               cnodes.size(), niblobs, nmatch, tens.size(), datapath, ident, m_calls);
    auto out = as_tensorset(tens, ident);
    outq.push_back(out);
    return true;
}

