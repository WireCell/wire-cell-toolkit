#include "WireCellClus/MultiAlgBlobClustering.h"
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
#include "WireCellUtil/NamedFactory.h"

#include <map>
#include <fstream>

WIRECELL_FACTORY(MultiAlgBlobClustering, WireCell::Clus::MultiAlgBlobClustering, WireCell::INamed,
                 WireCell::ITensorSetFilter, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Aux;
using namespace WireCell::Aux::TensorDM;
using namespace WireCell::Clus::Facade;
using namespace WireCell::PointCloud::Tree;

MultiAlgBlobClustering::MultiAlgBlobClustering()
  : Aux::Logger("MultiAlgBlobClustering", "clus")
//   , m_bee_dead("channel-deadarea", 1*units::mm, 3) // tolerance, minpts
{
}


static
std::string format_path(
    std::string path,
    const std::string& name,
    int ident,
    const std::map<std::string, std::string> subpaths)
{
    auto it = subpaths.find(name);
    if (it == subpaths.end()) {
        path += "/" + name;
    }
    else {
        path += it->second;
    }
    if (path.find("%") == std::string::npos) {
        return path;
    }
    return String::format(path, ident);
}

std::string MultiAlgBlobClustering::inpath(const std::string& name, int ident)
{
    return format_path(m_inpath, name, ident, m_insubpaths);
}
std::string MultiAlgBlobClustering::outpath(const std::string& name, int ident)
{
    return format_path(m_outpath, name, ident, m_outsubpaths);
}


void MultiAlgBlobClustering::configure(const WireCell::Configuration& cfg)
{
    m_groupings = convert(cfg["groupings"], m_groupings);

    m_inpath = get(cfg, "inpath", m_inpath);
    m_outpath = get(cfg, "outpath", m_outpath);

    for (const auto& jsp : cfg["insubpaths"]) {
        m_insubpaths[jsp["name"].asString()] = jsp["subpath"].asString();
    }
    for (const auto& jsp : cfg["outsubpaths"]) {
        m_outsubpaths[jsp["name"].asString()] = jsp["subpath"].asString();
    }

    {
        auto jcid = cfg["cluster_id_order"];
        if (jcid.isString()) {
            m_clusters_id_order = jcid.asString();
        }
    }

    if (cfg.isMember("bee_dir")) {
        log->debug("the 'bee_dir' option is no longer supported, instead use 'bee_zip' to name a .zip file");
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

    for (auto jtn : cfg["pipeline"]) {
        std::string tn = jtn.asString();
        log->debug("configuring clustering method: {}", tn);
        auto imeth = Factory::find_tn<IEnsembleVisitor>(tn);
        m_pipeline.emplace_back(EnsembleVisitor{tn, imeth});
    }

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
            bpc.grouping = get<std::string>(bps, "grouping", "live");
            bpc.visitor = get<std::string>(bps, "visitor", "");
            bpc.filter = get<int>(bps, "filter", 1); // 1 for on, 0 for off, -1 for inverse filter
            
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

    assign(cfg["groupings"], m_groupings);

    cfg["inpath"] = m_inpath;
    cfg["outpath"] = m_outpath;

    // repeat defaults as literals just incase some "clever" person tries to
    // call this method AFTER configure() as that method mutates m_inlive, etc.
    cfg["inlive"] = "/live";
    cfg["outlive"] = "/live";
    cfg["indead"] = "/dead";
    cfg["outdead"] = "/dead";

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
         // C++17 can not use structured bindings in lambda capture list.
         const std::string the_name = name;

        // Find the configuration for this name to check if it's individual
        auto it = std::find_if(m_bee_points_configs.begin(), m_bee_points_configs.end(),
                              [&the_name](const BeePointsConfig& cfg) { return cfg.name == the_name; });
        
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
                        fill_bee_points_from_cluster(it2->second, *cluster, config.pcname, config.coords, config.filter);
                    }
                }
            }
        }
    }else{ // fill in the global
        // std::cout << "Test: " << name << " " << grouping.wpids().size() << " " << grouping.nchildren() << std::endl;

        for (const auto* cluster : grouping.children()) {
            fill_bee_points_from_cluster(apa_bpts.global, *cluster, config.pcname, config.coords, config.filter);
        }
    }
}


// Helper function to fill bee points from a single cluster
void MultiAlgBlobClustering::fill_bee_points_from_cluster(
    Bee::Points& bpts, const Cluster& cluster, 
    const std::string& pcname, const std::vector<std::string>& coords, int filter)
{
    int clid = cluster.get_cluster_id(); //bpts.back_cluster_id() + 1;

    // std::cout << "Test: " << bpts.size() << " " << bpts.back_cluster_id() << " " <<  clid << std::endl;

    if (pcname == "steiner_pc"){
        // Export Steiner points ... 
        // std::cout << "Exporting Steiner points for cluster ID: " << clid << " " << cluster.nchildren() << std::endl;

        auto& steiner_pc = cluster.get_pc(pcname);
        if (steiner_pc.empty()) {
            return;
        }
        // Get coordinate arrays from the point cloud
        const auto& x_coords = steiner_pc.get(coords.at(0))->elements<double>();
        const auto& y_coords = steiner_pc.get(coords.at(1))->elements<double>(); 
        const auto& z_coords = steiner_pc.get(coords.at(2))->elements<double>();
        const auto& flag_steiner_terminal = steiner_pc.get("flag_steiner_terminal")->elements<int>();

         for (size_t i = 0; i < x_coords.size(); ++i) {
            // Create point from steiner point cloud
            Point vtx(x_coords[i], y_coords[i], z_coords[i]);

            // Get the point index from the default scope
            auto point_index = cluster.get_closest_point_index(vtx);
            
            auto charge_result = cluster.calc_charge_wcp(point_index, 4000, true);
            double point_charge = charge_result.second; // Extract the charge value from the pair

            if (flag_steiner_terminal[i]) {
                bpts.append(Point(x_coords[i], y_coords[i], z_coords[i]), point_charge, 1, 1);  // terminals  ... 
            }else{
                bpts.append(Point(x_coords[i], y_coords[i], z_coords[i]), point_charge, 0, 0); // non-terminals ...
            }
         }


    }else{
        // Get the scope
        Scope scope = {pcname, coords};
        
        auto filter_scope = cluster.get_scope_filter(scope);

        // std::cout << "Test: " << cluster.get_cluster_id() << " " << clid << " " << scope << " " << filter_scope << std::endl;

        bool use_scope = true;
        if (filter == 1) {
            use_scope = filter_scope;
        }
        else if (filter == 0) {
            use_scope = true; // ignore filter_scope, always true
        }
        else if (filter == -1) {
            use_scope = !filter_scope;
        }

        if (use_scope) {
            // Access the points through the cluster's scoped view
            const WireCell::PointCloud::Tree::ScopedView<double>& sv = cluster.sv<double>(scope);
            const auto& spcs = sv.pcs();
            const auto& nodes = sv.nodes(); // Get the nodes in the scoped view

            // Create a map to cache blob information to avoid recalculating for points in the same blob
            std::unordered_map<const WireCell::Clus::Facade::Blob*, std::pair<double, size_t>> blob_info;

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
                const auto* blob = node->value.facade<WireCell::Clus::Facade::Blob>();
                
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

}




void MultiAlgBlobClustering::fill_bee_patches_from_grouping(
    const WireCell::Clus::Facade::Grouping& grouping)
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
    const WireCell::Clus::Facade::Cluster& cluster)
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

    void dump(const std::string& ctx, const Ensemble& ensemble, bool shallow = true, bool mon = true)
    {
        if (!enable) return;
        if (mon) (*this)(ctx);

        log->debug("{} ensemble with {} groupings:", ctx, ensemble.nchildren());

        for (const auto* grouping : ensemble.children()) {

            {
                size_t npoints_total = 0;
                size_t nzero = 0;
                size_t count = 0;
                for (const auto* cluster : grouping->children()) {
                    int n = cluster->npoints();
                    if (n == 0) {
                        ++nzero;
                    }
                    npoints_total += n;
                    // log->debug("loaded cluster {} with {} points out of {}", count, n, npoints_total);
                    ++count;
                }

                auto name = grouping->get_name();

                log->debug("\tgrouping \"{}\": {}, {} points and {} clusters with no points",
                           name, *grouping, npoints_total, nzero);
            }

            if (shallow) continue;

            auto children = grouping->children();  // copy
            sort_clusters(children);
            size_t count = 0;
            for (const auto* cluster : children) {
                bool sane = cluster->sanity(log);
                log->debug("\t\tcluster {} {} sane:{}", count++, *cluster, sane);
            }
        }
    }
};


Grouping& MultiAlgBlobClustering::load_grouping(
    Ensemble& ensemble,
    const std::string& name,
    const std::string& path,
    const ITensorSet::pointer ints)
{
    const auto& tens = *ints->tensors();
    try {
        ensemble.add_grouping_node(name, as_pctree(tens, path));
    }
    catch (WireCell::KeyError& err) {
        log->warn("No pc-tree at tensor datapath {}, making empty", path);
        ensemble.make_grouping(name);
    }
        
    Grouping* grouping = ensemble.with_name(name).at(0);
    if (!grouping) {
        raise<KeyError>("failed to make grouping node %s at %s", name, path);
    }

    grouping->enumerate_idents();
    grouping->set_anodes(m_anodes);
    grouping->set_detector_volumes(m_dv);
    return *grouping;
}

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


    Points::node_t root;
    Ensemble& ensemble = *root.value.facade<Ensemble>();

    for (const auto& gname : m_groupings) {
        const auto datapath = inpath(gname, ident);
        load_grouping(ensemble, gname, datapath, ints);
        perf.dump("loaded " + gname, ensemble);
    }    

    if (m_save_deadarea) {
        auto gs = ensemble.with_name("live");
        if (gs.size()) {
            // Fill patches from the dead grouping
            fill_bee_patches_from_grouping(*gs[0]);
            perf("dump dead regions to bee");
        }
    }

    perf.dump("pre clustering", ensemble);

    for (const auto& config : m_bee_points_configs) {
        if (config.name != "img") {
            continue;
        }
        auto gs = ensemble.with_name("live");
        if (gs.empty()) {
            continue;
        }
        fill_bee_points(config.name, *gs[0]);
    }

    perf.dump("start clustering", ensemble);

    // THE MAIN LOOP
    for (const auto& cmeth : m_pipeline) {
        cmeth.meth->visit(ensemble);
        perf.dump(cmeth.name, ensemble);

        for (auto* grouping : ensemble.children()) {
            grouping->enumerate_idents(m_clusters_id_order);
        }

        // for (const auto& config : m_bee_points_configs) {
        //     if (config.name == "img") continue;
        //     if (config.visitor != cmeth.name) continue;
        //     auto gs = ensemble.with_name(config.grouping);
        //     if (gs.empty()) {
        //         continue;
        //     }
        //     fill_bee_points(config.name, *gs[0]);
        // }
    }

    //
    // At this point, the ensemble may have more or fewer groupings just "live"
    // and "dead" including no groupings at all.  But for now, we assume the
    // original "live" and "dead" still exist and with their original facades.
    // Famous last words....
    //
    

    // Fill all configured bee points sets
    for (const auto& config : m_bee_points_configs) {
        if(config.name == "img") continue;

        auto gs = ensemble.with_name(config.grouping);
        if (gs.empty()) {
            continue;
        }
        fill_bee_points(config.name, *gs[0]);

    }
    perf("dump live clusters to bee");

    auto grouping_names = ensemble.names();

    if (m_dump_json) {
        for (const auto& name : grouping_names) {
            auto gs = ensemble.with_name(name);
            Persist::dump(String::format("%s-summary-%d.json", name, ident),
                          json_summary(*gs[0]), true);
        }
    }

    log->debug("Produce pctrees with {} groupings", grouping_names.size());
    
    ITensor::vector outtens;
    for (const auto& name : grouping_names) {

        // This next bit may look a little weird and it is so some explanation
        // is warranted.  Originally, we had disembodied "root" grouping nodes,
        // live and dead.  To clean up the clustering api we added the
        // "ensemble" as root node with children consisting of grouping nodes.
        // At the time of writing, the as_tensors() does not like serializing
        // non-root nodes I do not want to debug right now.  And, I do not want
        // the "ensemble" concept to leak out from the MABC+clustering context.
        // So, I remove each grouping child node from the ensemble prior to
        // serializing.  The remove gives an auto_ptr so the node is destructed
        // as this loop progresses.
        auto gs = ensemble.with_name(name);
        auto& grouping = *gs[0];
        auto node = ensemble.remove_child(grouping);
        auto tens = as_tensors(*node, outpath(name, ident));
        outtens.insert(outtens.end(), tens.begin(), tens.end());
        log->debug("Produce {} tensors for grouping {}", tens.size(), name);
    }
    outts = as_tensorset(outtens, ident);

    perf("done");

    return true;
}
