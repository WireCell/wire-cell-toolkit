#include "WireCellClus/MultiAlgBlobClustering.h"
#include <WireCellClus/ClusteringFuncs.h>
#include "WireCellClus/Facade_Summary.h"


#include "WireCellAux/TensorDMpointtree.h"
#include "WireCellAux/TensorDMdataset.h"
#include "WireCellAux/TensorDMcommon.h"
#include "WireCellAux/SimpleTensorSet.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/ExecMon.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Graph.h"

#include <fstream>

WIRECELL_FACTORY(MultiAlgBlobClustering, WireCell::Clus::MultiAlgBlobClustering, WireCell::INamed,
                 WireCell::ITensorSetFilter, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Aux;
using namespace WireCell::Aux::TensorDM;
using namespace WireCell::PointCloud::Facade;
using namespace WireCell::PointCloud::Tree;

MultiAlgBlobClustering::MultiAlgBlobClustering()
  : Aux::Logger("MultiAlgBlobClustering", "clus")
//   , m_bee_dead("channel-deadarea", 1*units::mm, 3) // tolerance, minpts
{
}

void MultiAlgBlobClustering::configure(const WireCell::Configuration& cfg)
{
    m_inpath = get(cfg, "inpath", m_inpath);
    m_outpath = get(cfg, "outpath", m_outpath);

    if (cfg.isMember("bee_dir")) {
        log->warn("the 'bee_dir' option is no longer supported, instead use 'bee_zip' to name a .zip file");
    }
    std::string bee_zip = get<std::string>(cfg, "bee_zip", "mabc.zip");
    // Add new configuration option for initial index
    m_initial_index = get<int>(cfg, "initial_index", m_initial_index);

    //std::cout << "Xin: " << m_initial_index << " " << bee_zip << std::endl;
    m_sink.reset(bee_zip, m_initial_index);  // Use the new reset with initial index

    // Configure RSE numbers
    if (cfg.isMember("use_config_rse")) {
        m_use_config_rse = get(cfg, "use_config_rse", false);
        if (m_use_config_rse) {
            // Only read RSE if we're using configured values
            m_runNo = get(cfg, "runNo", m_runNo);
            m_subRunNo = get(cfg, "subRunNo", m_subRunNo);
            m_eventNo = get(cfg, "eventNo", m_eventNo);
            
             // Set RSE in sink during configuration
            m_sink.set_rse(m_runNo, m_subRunNo, m_eventNo);
        }
    }

    m_save_deadarea = get(cfg, "save_deadarea", m_save_deadarea);

    m_dead_live_overlap_offset = get(cfg, "dead_live_overlap_offset", m_dead_live_overlap_offset);

    m_func_cfgs = cfg["func_cfgs"];


    m_perf = get(cfg, "perf", m_perf);

    for (const auto& aname : cfg["anodes"]) {
        auto anode = Factory::find_tn<IAnodePlane>(aname.asString());
        m_anodes.push_back(anode);
    }

    m_dv = Factory::find_tn<IDetectorVolumes>(cfg["detector_volumes"].asString());

    m_dump_json = get<bool>(cfg, "dump_json", false);

    // Configure bee points sets
    if (cfg.isMember("bee_points_sets")) {
        auto bee_points_sets = cfg["bee_points_sets"];
        for (const auto& bps : bee_points_sets) {
            BeePointsConfig bpc;
            bpc.name = get<std::string>(bps, "name", "");
            bpc.detector = get<std::string>(bps, "detector", "uboone");
            bpc.algorithm = get<std::string>(bps, "algorithm", bpc.name);
            bpc.pcname = get<std::string>(bps, "pcname", "3d");
            
            // Get coordinates
            if (bps.isMember("coords")) {
                for (const auto& coord : bps["coords"]) {
                    bpc.coords.push_back(coord.asString());
                }
            } else {
                // Default coordinates
                bpc.coords = {"x", "y", "z"};
            }
            
            bpc.individual = get<bool>(bps, "individual", false);
            
            m_bee_points_configs.push_back(bpc);
            
            
            // If individual, also initialize bee points for each APA and face
            if (bpc.individual) {
                for (const auto& anode : m_anodes) {
                    int apa = anode->ident();
                    // Initialize the outer map if it doesn't exist
                    if (m_bee_points[bpc.name].by_apa_face.find(apa) == 
                        m_bee_points[bpc.name].by_apa_face.end()) {
                        m_bee_points[bpc.name].by_apa_face[apa] = std::map<int, Bee::Points>();
                    }
                    
                    // Initialize bee points for each face
                    for (size_t face_index = 0; face_index < anode->faces().size(); ++face_index) {
                        int face = anode->faces()[face_index]->which();
                        std::string algo_name = String::format("%s-apa%d-face%d", bpc.algorithm.c_str(), apa,  face);
                        // std::cout << "Test: Individual: " << algo_name << std::endl;
                        m_bee_points[bpc.name].by_apa_face[apa][face] =  Bee::Points(bpc.detector, algo_name);
                    }
                }
            }else{
                m_bee_points[bpc.name].global.detector(bpc.detector);
                m_bee_points[bpc.name].global.algorithm(String::format("%s-global", bpc.name));
                // std::cout << "Test: Global: " << m_bee_points[bpc.name].global.algorithm() << std::endl;
            }
            
            log->debug("Configured bee points set: {}, algorithm: {}, individual: {}", 
                        bpc.name, bpc.algorithm, bpc.individual ? "true" : "false");
        }
    } 

    // Initialize patches for each APA and face
    if (m_save_deadarea) {
        for (const auto& anode : m_anodes) {
            int apa = anode->ident();
            
            // Initialize the outer map if it doesn't exist
            if (m_bee_dead_patches.find(apa) == 
                m_bee_dead_patches.end()) {
                m_bee_dead_patches[apa] = std::map<int, Bee::Patches>();
            }
            
            // Initialize patches for each face
            for (size_t face_index = 0; face_index < anode->faces().size(); ++face_index) {
                int face = anode->faces()[face_index]->which();
                std::string name = String::format("channel-deadarea-apa%d-face%d", apa, face);
                m_bee_dead_patches[apa].insert({face,Bee::Patches(name, 1*units::mm, 3)}); // Same parameters as the global one
            }
        }
    }
}

WireCell::Configuration MultiAlgBlobClustering::default_configuration() const
{
    Configuration cfg;
    cfg["inpath"] = m_inpath;
    cfg["outpath"] = m_outpath;
    // cfg["bee_dir"] = m_bee_dir;
    cfg["bee_zip"] = "mabc.zip";
    cfg["save_deadarea"] = m_save_deadarea;

    // Add the new parameter to default configuration
    cfg["initial_index"] = m_initial_index;

    cfg["dead_live_overlap_offset"] = m_dead_live_overlap_offset;

    cfg["use_config_rse"] = false;  // By default, don't use configured RSE
    cfg["runNo"] = m_runNo;
    cfg["subRunNo"] = m_subRunNo;
    cfg["eventNo"] = m_eventNo;

    return cfg;
}

void MultiAlgBlobClustering::finalize()
{
    flush();
    m_sink.close();
}

static void reset_bee(int ident, WireCell::Bee::Points& bpts)
{
    int run=0, evt=0;
    if (ident > 0) {
        run = (ident >> 16) & 0x7fff;
        evt = (ident) & 0xffff;
    }
    bpts.reset(evt, 0, run);
}

void MultiAlgBlobClustering::flush(WireCell::Bee::Points& bpts, int ident)
{
    if (bpts.empty()) return;

    m_sink.write(bpts);
    reset_bee(ident, bpts);
}

void MultiAlgBlobClustering::flush(int ident)
{
    // flush(m_bee_img, ident);
    // flush(m_bee_ld,  ident);
     // Flush all bee points sets
     for (auto& [name, apa_bpts] : m_bee_points) {
        // Find the configuration for this name to check if it's individual
        auto it = std::find_if(m_bee_points_configs.begin(), m_bee_points_configs.end(),
                              [&name](const BeePointsConfig& cfg) { return cfg.name == name; });
        
        bool individual = (it != m_bee_points_configs.end()) ? it->individual : false;
        
        if (individual) {
            // Write individual bee points
            for (auto& [anode_id, face_map] : apa_bpts.by_apa_face) {
                for (auto& [face, bpts] : face_map) {
                    if (!bpts.empty()) {
                        m_sink.write(bpts);
                        // Clear after writing
                        int run = 0, evt = 0;
                        if (ident > 0) {
                            run = (ident >> 16) & 0x7fff;
                            evt = (ident) & 0xffff;
                        }
                        bpts.reset(evt, 0, run);
                    }
                }
            }
        } else {
            // Write global bee points
            if (!apa_bpts.global.empty()) {
                m_sink.write(apa_bpts.global);
                // Clear after writing
                int run = 0, evt = 0;
                if (ident > 0) {
                    run = (ident >> 16) & 0x7fff;
                    evt = (ident) & 0xffff;
                }
                apa_bpts.global.reset(evt, 0, run);
            }
        }
    }


    // if (m_save_deadarea && m_bee_dead.size()) {
    //     m_bee_dead.flush();
    //     m_sink.write(m_bee_dead);
    //     m_bee_dead.clear();
    // }
    if (m_save_deadarea) {
        
        // Flush individual patches
        for (auto& [apa, face_map] : m_bee_dead_patches) {
            for (auto& [face, patches] : face_map) {
                if (patches.size()) {
                    patches.flush();
                    m_sink.write(patches);
                    patches.clear();
                }
            }
        }
    }

    m_last_ident = ident;
}



// Helper function remains the same as in the previous response

void MultiAlgBlobClustering::fill_bee_points(const std::string& name, const Grouping& grouping)
{
    // std::cout << "Test: " << name << " " << grouping.wpids().size() << std::endl;
    
    if (m_bee_points.find(name) == m_bee_points.end()) {
        log->warn("Bee points set '{}' not found, skipping", name);
        return;
    }
    
    auto& apa_bpts = m_bee_points[name];
    
    // Find the configuration for this name
    auto it = std::find_if(m_bee_points_configs.begin(), m_bee_points_configs.end(),
                          [&name](const BeePointsConfig& cfg) { return cfg.name == name; });
    
    if (it == m_bee_points_configs.end()) {
        log->warn("Configuration for bee points set '{}' not found, skipping", name);
        return;
    }
    
    const auto& config = *it;
    
    // Reset RSE values for all points objects
    if (m_use_config_rse) {
        apa_bpts.global.rse(m_runNo, m_subRunNo, m_eventNo);
        for (auto& [apa, face_map] : apa_bpts.by_apa_face) {
            for (auto& [face, bpts] : face_map) {
                bpts.rse(m_runNo, m_subRunNo, m_eventNo);
            }
        }
    } else {
        // Use the default approach with ident
        int run = 0, evt = 0;
        if (m_last_ident > 0) {
            run = (m_last_ident >> 16) & 0x7fff;
            evt = (m_last_ident) & 0xffff;
        }
        apa_bpts.global.reset(evt, 0, run);
        for (auto& [anode_id, face_map] : apa_bpts.by_apa_face) {
            for (auto& [face, bpts] : face_map) {
                bpts.reset(evt, 0, run);
            }
        }
    }
    
    auto wpids = grouping.wpids();



    if (config.individual){ // fill in the individual APA
        for (auto wpid: wpids) {
            int apa = wpid.apa();
            int face = wpid.face();
            auto it = apa_bpts.by_apa_face.find(apa);
            if (it != apa_bpts.by_apa_face.end()) {
                auto it2 = it->second.find(face);
                if (it2 != it->second.end()) {
                    for (const auto* cluster : grouping.children()) {
                        fill_bee_points_from_cluster(it2->second, *cluster, config.pcname, config.coords);
                    }
                }
            }
        }
    }else{ // fill in the global
        // std::cout << "Test: " << name << " " << grouping.wpids().size() << " " << grouping.nchildren() << std::endl;

        for (const auto* cluster : grouping.children()) {
            fill_bee_points_from_cluster(apa_bpts.global, *cluster, config.pcname, config.coords);
        }
    }
}


// Helper function to fill bee points from a single cluster
void MultiAlgBlobClustering::fill_bee_points_from_cluster(
    Bee::Points& bpts, const Cluster& cluster, 
    const std::string& pcname, const std::vector<std::string>& coords)
{
    int clid = cluster.get_cluster_id(); //bpts.back_cluster_id() + 1;

    // std::cout << "Test: " << bpts.size() << " " << bpts.back_cluster_id() << " " <<  clid << std::endl;


    // Get the scope
    Scope scope = {pcname, coords};
    
    auto filter_scope = cluster.get_scope_filter(scope);

    // std::cout << "Test: " << cluster.get_cluster_id() << " " << clid << " " << filter_scope << std::endl;

    if(filter_scope){
        // Access the points through the cluster's scoped view
        const WireCell::PointCloud::Tree::ScopedView<double>& sv = cluster.sv<double>(scope);
        const auto& spcs = sv.pcs();
        const auto& nodes = sv.nodes(); // Get the nodes in the scoped view

        // Create a map to cache blob information to avoid recalculating for points in the same blob
        std::unordered_map<const WireCell::PointCloud::Facade::Blob*, std::pair<double, size_t>> blob_info;

        // std::cout << "Test: " << cluster.get_cluster_id() << " " << spcs.size() << std::endl;


        // For each scoped pointcloud (each corresponds to a blob)
        for (size_t spc_idx = 0; spc_idx < spcs.size(); ++spc_idx) {
            const auto& spc = spcs[spc_idx];
            auto x = spc.get().get(coords[0])->elements<double>();
            auto y = spc.get().get(coords[1])->elements<double>();
            auto z = spc.get().get(coords[2])->elements<double>();
            
            // Get the blob associated with this spc
            // The node_with_major() function gets the node for this major index (blob)
            const auto* node = nodes[spc_idx];
            const auto* blob = node->value.facade<WireCell::PointCloud::Facade::Blob>();
            
            // Calculate blob information if not already cached
            if (blob_info.find(blob) == blob_info.end()) {
                double blob_charge = blob->charge();
                size_t blob_npoints = blob->npoints();
                blob_info[blob] = {blob_charge, blob_npoints};
            }
            
            // Get cached blob info
            const auto& [blob_charge, blob_npoints] = blob_info[blob];
            
            // Calculate charge per point
            double point_charge = 0.0;
            if (blob_npoints > 0) {
                point_charge = blob_charge / blob_npoints;
            }
            
            const size_t size = x.size();
            for (size_t ind = 0; ind < size; ++ind) {
                // Use the calculated point_charge instead of the original charge
                bpts.append(Point(x[ind], y[ind], z[ind]), point_charge, clid, clid);
            }
        }

    }

}




void MultiAlgBlobClustering::fill_bee_patches_from_grouping(
    const WireCell::PointCloud::Facade::Grouping& grouping)
{
    // auto wpids = grouping.wpids();

    // For each cluster in the grouping
    for (const auto* cluster : grouping.children()) {
        // Get the wpids to determine which APA and face this cluster belongs to

        fill_bee_patches_from_cluster(*cluster);

        
        // if (!wpids.empty()) {
        //     // Store patches by APA and face
        //     for (auto wpid : wpids) {
        //         int apa = wpid.apa();
        //         int face = wpid.face();
        //        
        //     }
        // } 
    }
}


// Helper function to fill patches from a single cluster
void MultiAlgBlobClustering::fill_bee_patches_from_cluster(
    const WireCell::PointCloud::Facade::Cluster& cluster)
{
    int first_slice = -1;
    
    // Get the underlying node that contains this cluster
    const auto* cluster_node = cluster.node();
    if (!cluster_node) {
        log->warn("Cannot access node for cluster");
        return;
    }
    
    // Iterate through child nodes (blobs)
    for (const auto* bnode : cluster_node->children()) {
        auto wpid = bnode->value.facade<Blob>()->wpid();
        int apa = wpid.apa();
        int face = wpid.face();


        auto it_apa = m_bee_dead_patches.find(apa);
        if (it_apa != m_bee_dead_patches.end()) {
            auto it_face = it_apa->second.find(face);
            if (it_face != it_apa->second.end()) {
                auto & patches = it_face->second;

                // Access the local point clouds in the node
                const auto& lpcs = bnode->value.local_pcs();
                
                // Get the scalar PC to find the slice index
                if (lpcs.find("scalar") == lpcs.end()) {
                    continue;  // Skip if no scalar PC
                }
                const auto& pc_scalar = lpcs.at("scalar");
                
                // Get slice_index_min
                if (!pc_scalar.get("slice_index_min")) {
                    continue;  // Skip if no slice_index_min
                }
                int slice_index_min = pc_scalar.get("slice_index_min")->elements<int>()[0];
                
                // Set first_slice if not already set
                if (first_slice < 0) {
                    first_slice = slice_index_min;
                }
                
                // Skip blobs not on the first slice
                if (slice_index_min != first_slice) continue;
                
                // Access the corner point cloud
                if (lpcs.find("corner") == lpcs.end()) {
                    continue;  // Skip if no corner PC
                }
                const auto& pc_corner = lpcs.at("corner");
                
                // Get y and z coordinates
                if (!pc_corner.get("y") || !pc_corner.get("z")) {
                    continue;  // Skip if missing y or z
                }
                const auto& y = pc_corner.get("y")->elements<double>();
                const auto& z = pc_corner.get("z")->elements<double>();
                
                // Add to patches
                patches.append(y.begin(), y.end(), z.begin(), z.end());
            }
        }        
    }
}


struct Perf {
    bool enable;
    Log::logptr_t log;
    ExecMon em;

    Perf(bool e, Log::logptr_t l, const std::string& t = "starting MultiAlgBlobClustering")
      : enable(e)
      , log(l)
      , em(t)
    {
    }

    ~Perf()
    {
        if (!enable) return;
        log->debug("MultiAlgBlobClustering performance summary:\n{}", em.summary());
    }

    void operator()(const std::string& ctx)
    {
        if (!enable) return;
        em(ctx);
    }

    void dump(const std::string& ctx, const Grouping& grouping, bool shallow = true, bool mon = true)
    {
        if (!enable) return;
        if (mon) (*this)(ctx);
        log->debug("{} grouping {}", ctx, grouping);
        if (shallow) return;
        auto children = grouping.children();  // copy
        sort_clusters(children);
        size_t count = 0;
        for (const auto* cluster : children) {
            bool sane = cluster->sanity(log);
            log->debug("{} cluster {} {} sane:{}", ctx, count++, *cluster, sane);
        }
    }
};

bool MultiAlgBlobClustering::operator()(const input_pointer& ints, output_pointer& outts)
{
    outts = nullptr;
    if (!ints) {
        flush();
        log->debug("EOS at call {}", m_count++);
        return true;
    }

    Perf perf{m_perf, log};

    const int ident = ints->ident();
    log->debug("loading tensor set ident={} (last={})", ident, m_last_ident);
    if (m_last_ident < 0) {     // first time.
        if (m_use_config_rse) {
            // Set RSE in the sink
            m_sink.set_rse(m_runNo, m_subRunNo, m_eventNo);
        }
        // Use default behavior
        // reset_bee(ident, m_bee_img);
        // reset_bee(ident, m_bee_ld);
        m_last_ident = ident;
    }
    else if (m_last_ident != ident) {
        flush(ident);
        if (m_use_config_rse) {
            // Update event number for next event
            m_eventNo++;
            // Update RSE in sink
            m_sink.set_rse(m_runNo, m_subRunNo, m_eventNo);
        }
    }
    // else do nothing when ident is unchanged.

    std::string inpath = m_inpath;
    if (inpath.find("%") != std::string::npos) {
        inpath = String::format(inpath, ident);
    }

    const auto& intens = *ints->tensors();
    auto root_live = std::move(as_pctree(intens, inpath + "/live"));
    if (!root_live) {
        log->error("Failed to get dead point cloud tree from \"{}\"", inpath);
        raise<ValueError>("Failed to get live point cloud tree from \"%s\"", inpath);
    }
    auto grouping = root_live->value.facade<Grouping>();
    grouping->set_anodes(m_anodes);
    grouping->set_detector_volumes(m_dv);
    // grouping->set_params(m_geomhelper->get_params(m_anodes.front()->ident(), m_face));
    perf("loaded live clusters");
    {
        size_t npoints_total = 0;
        size_t nzero = 0;
        for (const auto* cluster : grouping->children()) {
            int n = cluster->npoints();
            if (n == 0) {
                ++nzero;
            }
            npoints_total += n;
        }
        log->debug("loaded live grouping with {} clusters, {} points, and {} clusters with no points",
                   grouping->nchildren(), npoints_total, nzero);
        // It is probably an error if nzero is not zero.
    }


    // log->debug("Got live pctree with {} children", root_live->nchildren());
    // log->debug(em("got live pctree"));
    log->debug("as_pctree from \"{}\"", inpath + "/dead");
    const std::string deadinpath = inpath + "/dead";
    Points::node_ptr root_dead;
    try {
        root_dead = as_pctree(intens, deadinpath);
        perf("loaded dead clusters");
    }
    catch (WireCell::KeyError& err) {
        log->warn("No pc-tree at datapath {}, assuming no 'dead' clusters", deadinpath);
        root_dead = std::make_unique<Points::node_t>();
    }


    
    // if (m_save_deadarea) {
    //     // fill_bee_patches(m_bee_dead, *root_dead.get());
    //     perf("loaded dump dead regions to bee");
    // }
    // log->debug("will {} {} dead patches", m_save_deadarea ? "save" : "not save", m_bee_dead.size());

    cluster_set_t cluster_connected_dead;

    // initialize clusters ...
    Grouping& live_grouping = *root_live->value.facade<Grouping>();

    Grouping& dead_grouping = *root_dead->value.facade<Grouping>();
    dead_grouping.set_anodes(m_anodes);
    dead_grouping.set_detector_volumes(m_dv);
    // dead_grouping.set_params(m_geomhelper->get_params(m_anodes.front()->ident(), m_face));
    
    perf("loaded dump live clusters to bee");
    if (m_save_deadarea) {
        // Fill patches from the dead grouping
        fill_bee_patches_from_grouping(dead_grouping); // true means use individual patches by APA/face
        perf("loaded dump dead regions to bee");
    }

    //perf.dump("original live clusters", live_grouping, false, false);
    //perf.dump("original dead clusters", dead_grouping, false, false);

    perf.dump("pre clustering", live_grouping);

    // set cluster id ... 
    int cluster_id = 1;
    for (auto* cluster : live_grouping.children()) {
        cluster->set_cluster_id(cluster_id++);
    }

    for (const auto& config : m_bee_points_configs) {
        if(config.name == "img")
            fill_bee_points(config.name, live_grouping);
    }

    for (const auto& func_cfg : m_func_cfgs) {
        // std::cout << "func_cfg: " << func_cfg << std::endl;
        auto func = getClusteringFunction(func_cfg);

        func(live_grouping, dead_grouping, cluster_connected_dead);

        perf.dump(func_cfg["name"].asString(), live_grouping);
    }

    // Fill all configured bee points sets
    for (const auto& config : m_bee_points_configs) {
        if(config.name != "img")
            fill_bee_points(config.name, live_grouping);
            // fill_bee_points(config.name, dead_grouping); // hack to check ClusteringRetile, shad_grouping is loaded as dead_grouping
    }

    perf("dump live clusters to bee");

    if (m_dump_json) {
        Persist::dump(String::format("live-summary-%d.json", ident), json_summary(live_grouping), true);
        Persist::dump(String::format("dead-summary-%d.json", ident), json_summary(dead_grouping), true);
    }


    std::string outpath = m_outpath;
    if (outpath.find("%") != std::string::npos) {
        outpath = String::format(outpath, ident);
    }
    auto outtens = as_tensors(*root_live.get(), outpath + "/live");
    perf("output live clusters to tensors");
    auto outtens_dead = as_tensors(*root_dead.get(), outpath + "/dead");
    perf("output dead clusters to tensors");
    // Merge
    outtens.insert(outtens.end(), outtens_dead.begin(), outtens_dead.end());
    outts = as_tensorset(outtens, ident);
    perf("combine tensors");

    root_live = nullptr;
    root_dead = nullptr;
    perf("clear pc tree memory");

    return true;
}
