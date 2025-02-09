#include "WireCellRoot/UbooneClusterSource.h"
#include "WireCellAux/SamplingHelpers.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/PointTree.h"

#include "TInterpreter.h"

#include <algorithm> // minmax


WIRECELL_FACTORY(UbooneClusterSource,
                 WireCell::Root::UbooneClusterSource,
                 WireCell::INamed,
                 WireCell::IBlobSampling,
                 WireCell::IConfigurable)


using namespace WireCell;
using namespace WireCell::Aux;
using WireCell::PointCloud::Tree::Points;
using WireCell::PointCloud::Tree::named_pointclouds_t;
using WireCell::PointCloud::Dataset;
using WireCell::PointCloud::Array;

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


}

WireCell::Configuration Root::UbooneClusterSource::default_configuration() const
{
    Configuration cfg;
    cfg["input"] = Json::arrayValue; // required
    cfg["sampler"] = "";             // optional
    cfg["light"] = "";               // optional
    cfg["flash"] = "";               // optional
    cfg["flashlight"] = "";          // optional
    return cfg;
}



bool Root::UbooneClusterSource::operator()(const IBlobSet::pointer& in, ITensorSet::pointer& out)
{
    out = nullptr;
    if (!in) { return true; }   // eos

    bool load_ok = m_files->next();
    if (!load_ok) {
        log->error("failed to load uboone cluster event");
        return false;
    }

    const auto& tblob = m_files->trees->live;
    const auto& cluster_ids = *tblob.cluster_id_vec; // spans blobs in "event"

    // We make an initial pass to make a cluster node for each unique cluster ID.
    Points::node_ptr root = std::make_unique<Points::node_t>();
    // Will need to navigate from cluster_id -> pc tree cluster node below.
    std::unordered_map<int, Points::node_t*> cnodes; 
    for (int cluster_id : cluster_ids) {
        auto cit = cnodes.find(cluster_id);
        if (cit != cnodes.end()) { continue; }
        cnodes[cluster_id] = root->insert();
    }

    // Iterate on blobs, make their blob-node, sample, add to cluster node.
    size_t nblobs = cluster_ids.size();
    const auto& iblobs = in->blobs();
    for (size_t bind=0; bind<nblobs; ++bind) {
        const IBlob::pointer iblob = iblobs[bind];
        // This MUST be the TTree entry number as set by UbooneBlobSource!
        const int entry = iblob->ident(); 
        const int cluster_id = cluster_ids[entry];
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
    

    if (!m_light_name.empty()) {
        auto& optical = m_files->trees->optical;
        auto& cf = optical["clusterflash"];

        auto cf_cind = cf.get("cluster")->elements<int>();
        auto cf_find = cf.get("flash")->elements<int>();
        const int nc = cf_cind.size();
        auto rchildren = root->children();
        for (int ic = 0; ic < nc; ++ic) {
            const int iclus = cf_cind[ic];
            const int iflash = cf_find[ic];
            rchildren[iclus]->value.local_pcs()["scalar"].add("flash", Array({iflash}));
        }

        optical.erase("clusterflash");
        root->value.local_pcs() = std::move(m_files->trees->optical);
    }

    return true;
}

