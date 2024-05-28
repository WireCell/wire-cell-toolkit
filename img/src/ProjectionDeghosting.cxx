
#include "WireCellImg/ProjectionDeghosting.h"
#include "WireCellImg/CSGraph.h"
#include "WireCellImg/Projection2D.h"
#include "WireCellAux/ClusterShadow.h"
#include "WireCellAux/SimpleCluster.h"
#include "WireCellAux/ClusterHelpers.h"
#include "WireCellUtil/GraphTools.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Rectangles.h"

#include <iterator>
#include <chrono>
#include <fstream>

WIRECELL_FACTORY(ProjectionDeghosting, WireCell::Img::ProjectionDeghosting, WireCell::INamed, WireCell::IClusterFilter,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Img;
using namespace WireCell::Aux;

Img::ProjectionDeghosting::ProjectionDeghosting()
  : Aux::Logger("ProjectionDeghosting", "img")
{
}

Img::ProjectionDeghosting::~ProjectionDeghosting() {}

WireCell::Configuration Img::ProjectionDeghosting::default_configuration() const
{
    WireCell::Configuration cfg;
    cfg["verbose"] = m_verbose;
    cfg["nchan"] = (unsigned int) m_nchan;
    cfg["nslice"] = (unsigned int) (m_nslice);
    cfg["dryrun"] = m_dryrun;

    cfg["global_deghosting_cut_nparas"] = m_global_deghosting_cut_nparas;
    for (int i = 0; i != 3 * m_global_deghosting_cut_nparas; i++) {
        cfg["global_deghosting_cut_values"][i] = m_global_deghosting_cut_values[i];
    }
    cfg["uncer_cut"] = m_uncer_cut;
    cfg["dead_default_charge"] = m_dead_default_charge;
    for (int i = 0; i != 4; i++) {
        cfg["judge_alt_cut_values"][i] = m_judge_alt_cut_values[i];
    }
    return cfg;
}

void Img::ProjectionDeghosting::configure(const WireCell::Configuration& cfg)
{
    m_verbose = get<bool>(cfg, "verbose", m_verbose);
    m_nchan = get<size_t>(cfg, "nchan", m_nchan);
    m_nslice = get<size_t>(cfg, "nslice", m_nslice);
    m_dryrun = get<bool>(cfg, "dryrun", m_dryrun);

    m_global_deghosting_cut_nparas = get<int>(cfg, "global_deghosting_cut_nparas", m_global_deghosting_cut_nparas);
    if (cfg.isMember("global_deghosting_cut_values")) {
        m_global_deghosting_cut_values.clear();
        for (auto value : cfg["global_deghosting_cut_values"]) {
            m_global_deghosting_cut_values.push_back(value.asDouble());
        }
    }

    m_uncer_cut = get<double>(cfg, "uncer_cut", m_uncer_cut);
    if (cfg.isMember("judge_alt_cut_values")) {
        m_judge_alt_cut_values.clear();
        for (auto value : cfg["judge_alt_cut_values"]) {
            m_judge_alt_cut_values.push_back(value.asDouble());
        }
    }

    log->debug("{}", cfg);
}

namespace {
    // struct pair_hash
    // {
    //     template <class T1, class T2>
    //     std::size_t operator () (std::pair<T1, T2> const &pair) const
    //     {
    //         std::size_t h1 = std::hash<T1>()(pair.first);
    //         std::size_t h2 = std::hash<T2>()(pair.second);
    //         return h1 ^ h2;
    //     }
    // };

    using cluster_proj_t = std::unordered_map<ClusterShadow::vdesc_t, Projection2D::LayerProjection2DMap>;
    Projection2D::LayerProjection2DMap& get_projection(cluster_proj_t& id2lproj, const ClusterShadow::vdesc_t& id,
                                                       const WireCell::cluster_graph_t& cg,
                                                       const std::set<cluster_vertex_t>& group, const size_t nchan,
                                                       const size_t nslice, double uncer_cut,
                                                       double dead_default_charge)
    {
        if (id2lproj.find(id) != id2lproj.end()) {
            auto& lproj = id2lproj.at(id);
            return lproj;
        }
        else {
            auto lproj = Projection2D::get_projection(cg, group, nchan, nslice, uncer_cut, dead_default_charge);
            id2lproj[id] = lproj;
            return id2lproj[id];
        }
    }

    // run in keep mode if removal_mode is false
    cluster_graph_t remove_blobs(const cluster_graph_t& cg, std::unordered_set<cluster_vertex_t> blobs,
                                 bool removal_mode = true)
    {
        cluster_graph_t cg_out;

        size_t nblobs = 0;
        std::unordered_map<cluster_vertex_t, cluster_vertex_t> old2new;
        for (const auto& vtx : GraphTools::mir(boost::vertices(cg))) {
            const auto& node = cg[vtx];
            if (node.code() == 'b') {
                bool removal = blobs.find(vtx) != blobs.end();
                if (!removal_mode) removal = !removal;
                if (!removal) {
                    ++nblobs;
                    old2new[vtx] = boost::add_vertex(node, cg_out);
                }
            }
            else {
                old2new[vtx] = boost::add_vertex(node, cg_out);
            }
        }

        for (auto edge : GraphTools::mir(boost::edges(cg))) {
            auto old_tail = boost::source(edge, cg);
            auto old_head = boost::target(edge, cg);

            auto old_tit = old2new.find(old_tail);
            if (old_tit == old2new.end()) {
                continue;
            }
            auto old_hit = old2new.find(old_head);
            if (old_hit == old2new.end()) {
                continue;
            }
            boost::add_edge(old_tit->second, old_hit->second, cg_out);
        }
        // std::cout << "boost::num_vertices(cg_out): " << boost::num_vertices(cg_out) << std::endl;
        return cg_out;
    }
}  // namespace

bool Img::ProjectionDeghosting::operator()(const input_pointer& in, output_pointer& out)
{
    out = nullptr;
    if (!in) {
        log->debug("EOS");
        return true;
    }

    const auto in_graph = in->graph();
    log->debug("in_graph: {}", dumps(in_graph));
    // blob shadow ... vertex = blob, edge --> wire/channel, plane
    BlobShadow::graph_t bsgraph = BlobShadow::shadow(in_graph, 'w');  // or 'c'

    // cluster shadow map, blob to cluster map ???
    ClusterShadow::blob_cluster_map_t clusters;
    auto cs_graph = ClusterShadow::shadow(in_graph, bsgraph, clusters);

    // make a cluster -> blob map
    std::unordered_map<ClusterShadow::vdesc_t, std::set<cluster_vertex_t>> c2b;
    for (auto p : clusters) {
        c2b[p.second].insert(p.first);
    }
    log->debug("c2b.size(): {} ", c2b.size());

    // Count edges
    std::map<std::string, int> counters = {{"cs_edges", 0},      {"cs_vertices", 0},         {"nblobs", 0},
                                           {"nblobs_filter", 0}, {"nblobs_clusters_tmp", 0}, {"nblobs_clusters", 0}};

    //    std::ofstream fout("ProjectionDeghosting.log");
    // using pairset_t = std::unordered_set<std::pair<ClusterShadow::vdesc_t, ClusterShadow::vdesc_t>, pair_hash>;
    // using layer2pairset_t = std::unordered_map<WireCell::WirePlaneLayer_t, pairset_t>;
    // layer2pairset_t layer2compared;

    // cluster desc -> LayerProjection2DMap
    cluster_proj_t id2lproj;

    // TODO: figure out how to do this right
    std::unordered_set<cluster_vertex_t> tagged_bs;

    // how to build a 2D cluster --> 3D cluster vector map (for each view)???
    // tail and head --> key of map -->
    // LayerProjection2DMap& aaa = tail->get_projection();
    // modify --> ClusterShadow ???
    // 2D projection --> (3D cluster, view)
    // also need clusters not in this cluster shadow things ...

    // figure out the connected components --> get projection for every vertex
    // store group number --> set of vertices & set of vertices and group number ...
    std::unordered_map<ClusterShadow::vdesc_t, int> cluster2id;
    // auto nclust =
    boost::connected_components(cs_graph, boost::make_assoc_property_map(cluster2id));

    // store the three maps in terms of group number?
    std::unordered_map<WireCell::WirePlaneLayer_t,
                       std::unordered_map<ClusterShadow::vdesc_t, std::vector<ClusterShadow::vdesc_t>>>
        wp_2D_3D_clus_map;
    {
        std::unordered_map<ClusterShadow::vdesc_t, std::vector<ClusterShadow::vdesc_t>> u_2D_3D_clus_map;
        std::unordered_map<ClusterShadow::vdesc_t, std::vector<ClusterShadow::vdesc_t>> v_2D_3D_clus_map;
        std::unordered_map<ClusterShadow::vdesc_t, std::vector<ClusterShadow::vdesc_t>> w_2D_3D_clus_map;
        wp_2D_3D_clus_map[kUlayer] = u_2D_3D_clus_map;
        wp_2D_3D_clus_map[kVlayer] = v_2D_3D_clus_map;
        wp_2D_3D_clus_map[kWlayer] = w_2D_3D_clus_map;
    }

    // run algorithm ...
    for (auto cs_cluster : GraphTools::mir(boost::vertices(cs_graph))) {
        auto b_cluster = c2b[cs_cluster];          // cluster_vertex_t
        int gid_cluster = cluster2id[cs_cluster];  // cluster id ...
        Projection2D::LayerProjection2DMap& proj_cluster =
            get_projection(id2lproj, cs_cluster, in_graph, b_cluster, m_nchan, m_nslice, m_uncer_cut,
                           m_dead_default_charge);  // 3 views are all here ...

        // empty ...
        if (proj_cluster.m_number_slices == 0) continue;

        // loop over each plane ...
        for (auto it = wp_2D_3D_clus_map.begin(); it != wp_2D_3D_clus_map.end(); it++) {
            WireCell::WirePlaneLayer_t layer_cluster = it->first;

            //
            if (proj_cluster.m_number_layer_slices[layer_cluster] == 0) continue;

            std::unordered_map<ClusterShadow::vdesc_t, std::vector<ClusterShadow::vdesc_t>>& clus_2D_3D_map =
                it->second;
            Projection2D::Projection2D& proj2D_cluster = proj_cluster.m_layer_proj[layer_cluster];

            // flag about saving or not ...
            bool flag_save = true;
            std::vector<ClusterShadow::vdesc_t> to_be_removed;

            for (auto it1 = clus_2D_3D_map.begin(); it1 != clus_2D_3D_map.end(); it1++) {
                ClusterShadow::vdesc_t clust3D = it1->first;
                std::vector<ClusterShadow::vdesc_t>& vec_clust3D = it1->second;
                int gid_clust3D = cluster2id[clust3D];
                // not in the same connected component, remove ...
                if (gid_clust3D != gid_cluster) continue;

                // judge if an edge (with the proper plan number) existed between cs_cluster and clust3D ...
                //	  auto edge_ranges = boost::edge_range(cs_cluster, clust3D, cs_graph);
                //	  for (auto it2 = edge_ranges.first; it2 != edge_ranges.second; it2++){
                for (const auto& cedge : GraphTools::mir(boost::edge_range(cs_cluster, clust3D, cs_graph))) {
                    // if not the same plane, skip ...
                    const auto& eobj = cs_graph[cedge];
                    if (layer_cluster != eobj.wpid.layer()) continue;

                    auto b_clust3D = c2b[clust3D];
                    Projection2D::LayerProjection2DMap& proj_clust3D = get_projection(
                        id2lproj, clust3D, in_graph, b_clust3D, m_nchan, m_nslice, m_uncer_cut, m_dead_default_charge);

                    Projection2D::Projection2D& proj2D_clust3D = proj_clust3D.m_layer_proj[layer_cluster];

                    int coverage = Projection2D::judge_coverage(proj2D_clust3D, proj2D_cluster,
                                                                m_uncer_cut);  // ref is eixsting, tar is new clust ...
                    // if (coverage == 1){ // tar is part of ref
                    if (coverage == Projection2D::REF_COVERS_TAR) {
                        flag_save = false;
                        break;
                        //	    }else if (coverage == 2){ // tar is same as ref
                    }
                    else if (coverage == Projection2D::REF_EQ_TAR) {
                        flag_save = false;
                        vec_clust3D.push_back(cs_cluster);
                        break;
                        //	    }else if (coverage == -1){ // ref is part of tar
                    }
                    else if (coverage == Projection2D::TAR_COVERS_REF) {
                        to_be_removed.push_back(clust3D);
                    }

                }  // loop over edges ..

                if (!flag_save) break;

            }  // loop over existing ...

            for (auto it1 = to_be_removed.begin(); it1 != to_be_removed.end(); it1++) {
                clus_2D_3D_map.erase(*it1);
            }
            if (flag_save) {
                std::vector<ClusterShadow::vdesc_t> vec_3Dclus;
                vec_3Dclus.push_back(cs_cluster);
                clus_2D_3D_map[cs_cluster] = vec_3Dclus;
            }

        }  // loop over plane
    }

    log->debug("2D --> 3D size: {} {} {}", wp_2D_3D_clus_map[kUlayer].size(), wp_2D_3D_clus_map[kVlayer].size(),
               wp_2D_3D_clus_map[kWlayer].size());

    // secondary loop over plane
    for (auto it = wp_2D_3D_clus_map.begin(); it != wp_2D_3D_clus_map.end(); it++) {
        WireCell::WirePlaneLayer_t layer_cluster = it->first;
        std::unordered_map<ClusterShadow::vdesc_t, std::vector<ClusterShadow::vdesc_t>>& clus_2D_3D_map = it->second;

        // loop over cluster
        std::vector<ClusterShadow::vdesc_t> to_be_removed;
        for (auto it1 = clus_2D_3D_map.begin(); it1 != clus_2D_3D_map.end(); it1++) {
            ClusterShadow::vdesc_t clust3D = it1->first;
            if (find(to_be_removed.begin(), to_be_removed.end(), clust3D) == to_be_removed.end()) {
                int gid_clust3D = cluster2id[clust3D];
                auto b_clust3D = c2b[clust3D];
                Projection2D::LayerProjection2DMap& proj_clust3D = get_projection(
                    id2lproj, clust3D, in_graph, b_clust3D, m_nchan, m_nslice, m_uncer_cut, m_dead_default_charge);
                Projection2D::Projection2D& proj2D_clust3D = proj_clust3D.m_layer_proj[layer_cluster];
                auto it2 = it1;
                it2++;
                for (auto it3 = it2; it3 != clus_2D_3D_map.end(); it3++) {
                    ClusterShadow::vdesc_t comp_3Dclus = it3->first;
                    int gid_comp_3Dclus = cluster2id[comp_3Dclus];
                    if (gid_clust3D != gid_comp_3Dclus) continue;

                    // judge if an edge (with the proper plan number) existed between cs_cluster and clust3D ...
                    //	  auto edge_ranges = boost::edge_range(cs_cluster, clust3D, cs_graph);
                    //	  for (auto it2 = edge_ranges.first; it2 != edge_ranges.second; it2++){
                    for (const auto& cedge : GraphTools::mir(boost::edge_range(clust3D, comp_3Dclus, cs_graph))) {
                        // if not the same plane, skip ...
                        const auto& eobj = cs_graph[cedge];
                        if (layer_cluster != eobj.wpid.layer()) continue;

                        auto b_comp_3Dclus = c2b[comp_3Dclus];
                        Projection2D::LayerProjection2DMap& proj_comp_3Dclus =
                            get_projection(id2lproj, comp_3Dclus, in_graph, b_comp_3Dclus, m_nchan, m_nslice,
                                           m_uncer_cut, m_dead_default_charge);
                        Projection2D::Projection2D& proj2D_comp_3Dclus = proj_comp_3Dclus.m_layer_proj[layer_cluster];

                        // std::vector<double> cut_values{0.05, 0.33, 0.15, 0.33};
                        int coverage_alt = Projection2D::judge_coverage_alt(proj2D_comp_3Dclus, proj2D_clust3D,
                                                                            m_judge_alt_cut_values, m_uncer_cut);
                        // if (coverage_alt == 1){
                        if (coverage_alt == Projection2D::REF_COVERS_TAR) {
                            to_be_removed.push_back(clust3D);
                            //	      }else if (coverage_alt==-1){
                        }
                        else if (coverage_alt == Projection2D::TAR_COVERS_REF) {
                            to_be_removed.push_back(comp_3Dclus);
                        }

                    }  // loop over edges
                }
            }
        }  // loop clusters

        for (auto it1 = to_be_removed.begin(); it1 != to_be_removed.end(); it1++) {
            clus_2D_3D_map.erase(*it1);
        }

        for (auto it1 = clus_2D_3D_map.begin(); it1 != clus_2D_3D_map.end(); it1++) {
            for (auto it2 = it1->second.begin(); it2 != it1->second.end(); it2++) {
                ClusterShadow::vdesc_t clust3D = *it2;
                auto b_clust3D = c2b[clust3D];
                Projection2D::LayerProjection2DMap& proj_clust3D = get_projection(
                    id2lproj, clust3D, in_graph, b_clust3D, m_nchan, m_nslice, m_uncer_cut, m_dead_default_charge);
                proj_clust3D.m_saved_flag++;
            }
        }

    }  // loop over plane

    log->debug("2D --> 3D size: {} {} {}", wp_2D_3D_clus_map[kUlayer].size(), wp_2D_3D_clus_map[kVlayer].size(),
               wp_2D_3D_clus_map[kWlayer].size());

    // summarize the results ...
    for (auto it = wp_2D_3D_clus_map.begin(); it != wp_2D_3D_clus_map.end(); it++) {
        //      WireCell::WirePlaneLayer_t layer_cluster = it->first;
        std::unordered_map<ClusterShadow::vdesc_t, std::vector<ClusterShadow::vdesc_t>>& clus_2D_3D_map = it->second;
        for (auto it1 = clus_2D_3D_map.begin(); it1 != clus_2D_3D_map.end(); it1++) {
            int max_flag_saved = -1;
            for (auto it2 = it1->second.begin(); it2 != it1->second.end(); it2++) {
                ClusterShadow::vdesc_t clust3D = *it2;
                auto b_clust3D = c2b[clust3D];
                Projection2D::LayerProjection2DMap& proj_clust3D = get_projection(
                    id2lproj, clust3D, in_graph, b_clust3D, m_nchan, m_nslice, m_uncer_cut, m_dead_default_charge);
                if (proj_clust3D.m_saved_flag > max_flag_saved) max_flag_saved = proj_clust3D.m_saved_flag;
            }

            for (auto it2 = it1->second.begin(); it2 != it1->second.end(); it2++) {
                ClusterShadow::vdesc_t clust3D = *it2;
                auto b_clust3D = c2b[clust3D];
                Projection2D::LayerProjection2DMap& proj_clust3D = get_projection(
                    id2lproj, clust3D, in_graph, b_clust3D, m_nchan, m_nslice, m_uncer_cut, m_dead_default_charge);
                if (proj_clust3D.m_saved_flag != max_flag_saved) proj_clust3D.m_saved_flag_1++;
            }
        }
    }

    // delete blobs ...
    for (auto cs_cluster : GraphTools::mir(boost::vertices(cs_graph))) {
        auto b_cluster = c2b[cs_cluster];  // cluster_vertex_t
        //       int gid_cluster = cluster2id[cs_cluster];
        Projection2D::LayerProjection2DMap& proj_cluster = get_projection(
            id2lproj, cs_cluster, in_graph, b_cluster, m_nchan, m_nslice, m_uncer_cut, m_dead_default_charge);
        //       float total_charge = proj_cluster.m_estimated_total_charge;
        double min_charge = proj_cluster.m_estimated_minimum_charge;
        int flag_saved = proj_cluster.m_saved_flag;
        int flag_saved_1 = proj_cluster.m_saved_flag_1;
        int n_blobs = proj_cluster.m_number_blobs;
        int n_timeslices = proj_cluster.m_number_slices;

        int saved = 0;
        if (flag_saved - flag_saved_1 == 3) {
            // look at each cell level ...
            if (sqrt(pow(n_timeslices / m_global_deghosting_cut_values.at(0), 2) +
                     pow(min_charge / n_blobs / m_global_deghosting_cut_values.at(1), 2)) < 1 ||
                min_charge / n_blobs / m_global_deghosting_cut_values.at(2) < 1.) {
                saved = 0;
            }
            else {
                saved = 1;
            }
        }
        else if (flag_saved - flag_saved_1 == 2) {
            if (sqrt(pow(n_timeslices / m_global_deghosting_cut_values.at(m_global_deghosting_cut_nparas), 2) +
                     pow(min_charge / n_blobs / m_global_deghosting_cut_values.at(m_global_deghosting_cut_nparas + 1),
                         2)) < 1 ||
                min_charge / n_blobs / m_global_deghosting_cut_values.at(m_global_deghosting_cut_nparas + 2) < 1.) {
                saved = 0;
            }
            else {
                saved = 1;
            }
        }
        else if (flag_saved - flag_saved_1 == 1) {
            if (sqrt(pow(n_timeslices / m_global_deghosting_cut_values.at(2 * m_global_deghosting_cut_nparas), 2) +
                     pow(min_charge / n_blobs /
                             m_global_deghosting_cut_values.at(2 * m_global_deghosting_cut_nparas + 1),
                         2)) < 1 ||
                min_charge / n_blobs / m_global_deghosting_cut_values.at(2 * m_global_deghosting_cut_nparas + 2) < 1.) {
                saved = 0;
            }
            else {
                saved = 1;
            }
        }

        if (saved == 0) {
            for (auto& b : b_cluster) {
                tagged_bs.insert(b);
            }
        }
    }

    //    fout.close();
    // for (auto c : counters) {
    //    log->debug("{} : {} ", c.first, c.second);
    // }

    // true: remove tagged_bs; false: only keep tagged_bs
    auto out_graph = remove_blobs(in_graph, tagged_bs, true);

    // debug info
    log->debug("tagged_bs.size(): {}", tagged_bs.size());
    log->debug("in_graph: {}", dumps(in_graph));
    log->debug("out_graph: {}", dumps(out_graph));

    out = std::make_shared<Aux::SimpleCluster>(out_graph, in->ident());
    if (m_dryrun) {
        out = std::make_shared<Aux::SimpleCluster>(in_graph, in->ident());
    }
    return true;
}
