#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/ClusteringFuncsMixins.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"

class ClusteringSeparate;
WIRECELL_FACTORY(ClusteringSeparate, ClusteringSeparate,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)


using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Graphs;
using namespace WireCell::Clus::Facade;
using namespace WireCell::PointCloud::Tree;



static void clustering_separate(Grouping& live_grouping,
                                IDetectorVolumes::pointer dv,
                                IPCTransformSet::pointer pcts,
                                const Tree::Scope& scope, 
                                const bool use_ctpc);

class ClusteringSeparate : public IConfigurable, public Clus::IEnsembleVisitor, private NeedDV, private NeedPCTS, private NeedScope {
public:
    ClusteringSeparate() {}
    virtual ~ClusteringSeparate() {}

    void configure(const WireCell::Configuration& config) {
        NeedDV::configure(config);
        NeedPCTS::configure(config);
        NeedScope::configure(config);
        
        use_ctpc_ = get(config, "use_ctpc", true);
    }

    void visit(Ensemble& ensemble) const {
        auto& live = *ensemble.with_name("live").at(0);
        clustering_separate(live, m_dv, m_pcts, m_scope, use_ctpc_);
    }

private:
    double use_ctpc_{true};
};


// The original developers do not care.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"

// this algorithm should be able to handle multiple APA/face now ..
static void clustering_separate(
    Grouping& live_grouping,
    const IDetectorVolumes::pointer dv,                // detector volumes
    const IPCTransformSet::pointer pcts,
    const Tree::Scope& scope,
    const bool use_ctpc)
{
    // Check that live_grouping has exactly one wpid
	// if (live_grouping.wpids().size() != 1 ) {
	// 	throw std::runtime_error("Live or Dead grouping must have exactly one wpid");
	// }
    geo_point_t drift_dir_abs(1,0,0);

    std::vector<Cluster *> live_clusters = live_grouping.children();  // copy
    // sort the clusters by length using a lambda function
    std::sort(live_clusters.begin(), live_clusters.end(), [](const Cluster *cluster1, const Cluster *cluster2) {
        return cluster1->get_length() > cluster2->get_length();
    });

    auto wpids = live_grouping.wpids();

    WirePlaneId wpid_all(0);
    // double det_FV_ymin = dv->metadata(wpid_all)["FV_ymin"].asDouble();
    double det_FV_ymax = dv->metadata(wpid_all)["FV_ymax"].asDouble();

    geo_point_t beam_dir(0, 0, 1);
    geo_point_t vertical_dir(0, 1, 0);
    
    // Get vertical_dir from metadata
    Json::Value vertical_dir_json = dv->metadata(wpid_all)["vertical_dir"];
    Json::Value beam_dir_json = dv->metadata(wpid_all)["beam_dir"];

    if (!vertical_dir_json.isNull() && vertical_dir_json.isArray() && vertical_dir_json.size() >= 3) {
        vertical_dir = geo_point_t(
            vertical_dir_json[0].asDouble(),
            vertical_dir_json[1].asDouble(),
            vertical_dir_json[2].asDouble()
        );
    } 
    if (!beam_dir_json.isNull() && beam_dir_json.isArray() && beam_dir_json.size() >= 3) {
        beam_dir = geo_point_t(
            beam_dir_json[0].asDouble(),
            beam_dir_json[1].asDouble(),
            beam_dir_json[2].asDouble()
        );
    } 

    std::vector<Cluster *> new_clusters;
    std::vector<Cluster *> del_clusters;

    for (size_t i = 0; i != live_clusters.size(); i++) {
        Cluster *cluster = live_clusters.at(i);
        if (!cluster->get_scope_filter(scope)) continue;
        if (cluster->get_default_scope().hash() != scope.hash()) {
            cluster->set_default_scope(scope);
            // std::cout << "Test: Set default scope: " << pc_name << " " << coords[0] << " " << coords[1] << " " << coords[2] << " " << cluster->get_default_scope().hash() << " " << scope.hash() << std::endl;
        }

        if (cluster->get_length() > 100 * units::cm) {
            std::vector<geo_point_t> boundary_points;
            std::vector<geo_point_t> independent_points;

            bool flag_proceed =
                JudgeSeparateDec_2(cluster, dv, drift_dir_abs, boundary_points, independent_points, cluster->get_length());
        
            if (!flag_proceed && cluster->get_length() > 100 * units::cm &&
                JudgeSeparateDec_1(cluster, drift_dir_abs, cluster->get_length()) &&
                independent_points.size() > 0) {
                bool flag_top = false;
                for (size_t j = 0; j != independent_points.size(); j++) {
                    if (independent_points.at(j).y() > det_FV_ymax) {
                        flag_top = true;
                        break;
                    }
                }

                geo_point_t main_dir(cluster->get_pca().axis.at(0).x(), cluster->get_pca().axis.at(0).y(),
                                     cluster->get_pca().axis.at(0).z());

                if (flag_top) {
                    if (fabs(main_dir.angle(beam_dir) - 3.1415926 / 2.) / 3.1415926 * 180. < 16 ||
                        fabs(main_dir.angle(beam_dir) - 3.1415926 / 2.) / 3.1415926 * 180. < 33 &&
                            cluster->get_length() > 160 * units::cm ||
                        fabs(main_dir.angle(beam_dir) - 3.1415926 / 2.) / 3.1415926 * 180. < 40 &&
                            cluster->get_length() > 260 * units::cm ||
                        fabs(main_dir.angle(beam_dir) - 3.1415926 / 2.) / 3.1415926 * 180. < 65 &&
                            cluster->get_length() > 360 * units::cm ||
                        fabs(main_dir.angle(beam_dir) - 3.1415926 / 2.) / 3.1415926 * 180. < 45 &&
                            cluster->get_pca().values.at(1) > 0.75 * cluster->get_pca().values.at(0) ||
                        fabs(main_dir.angle(beam_dir) - 3.1415926 / 2.) / 3.1415926 * 180. < 40 &&
                            cluster->get_pca().values.at(1) > 0.55 * cluster->get_pca().values.at(0)) {
                        flag_proceed = true;
                    }
                    else {
                        if (fabs(main_dir.angle(beam_dir) - 3.1415926 / 2.) / 3.1415926 * 180. < 40 &&
                            cluster->get_pca().values.at(1) > 0.2 * cluster->get_pca().values.at(0)) {
                            // std::vector<Cluster *> temp_sep_clusters = Separate_2(cluster, 10 * units::cm);
                            const auto b2id = Separate_2(cluster, scope, 10 * units::cm);
                            std::set<int> ids;
                            for (const auto& id : b2id) {
                                ids.insert(id);
                            }
                            int num_clusters = 0;
                    
                            for (const auto id : ids) {
                                double length_1 = get_length(cluster, b2id, id);
                                if (length_1 > 60 * units::cm) num_clusters++;
                            }
                            if (num_clusters > 1) flag_proceed = true;
                        }
                    }
                }
                else {
                    if (fabs(main_dir.angle(beam_dir) - 3.1415926 / 2.) / 3.1415926 * 180. < 4 &&
                            cluster->get_length() > 170 * units::cm ||
                        fabs(main_dir.angle(beam_dir) - 3.1415926 / 2.) / 3.1415926 * 180. < 25 &&
                            cluster->get_length() > 210 * units::cm ||
                        fabs(main_dir.angle(beam_dir) - 3.1415926 / 2.) / 3.1415926 * 180. < 28 &&
                            cluster->get_length() > 270 * units::cm ||
                        fabs(main_dir.angle(beam_dir) - 3.1415926 / 2.) / 3.1415926 * 180. < 35 &&
                            cluster->get_length() > 330 * units::cm ||
                        fabs(main_dir.angle(beam_dir) - 3.1415926 / 2.) / 3.1415926 * 180. < 30 &&
                            cluster->get_pca().values.at(1) > 0.55 * cluster->get_pca().values.at(0)) {
                        flag_proceed = true;
                    }
                }
                //	std::cout << flag_top << " " << flag_proceed << std::endl;
            }
     

            if (flag_proceed) {
                if (JudgeSeparateDec_1(cluster, drift_dir_abs, cluster->get_length())) {
                    //	  std::cerr << em("sep prepare sep") << std::endl;

                    const size_t orig_nchildren = cluster->nchildren();
                    //std::cout << "Separate Cluster with " << orig_nchildren << " blobs (ctpc) length " << cluster->get_length() << std::endl;
                    std::vector<Cluster *> sep_clusters =
                        Separate_1(use_ctpc, cluster, boundary_points, independent_points,  cluster->get_length(), vertical_dir, beam_dir, dv, pcts, scope);
                    
                    std::cout << "Separate Separate_1 for " << orig_nchildren << " " << " returned " << sep_clusters.size() << " clusters" << std::endl;

                    Cluster *cluster1 = sep_clusters.at(0);
                    new_clusters.push_back(cluster1);

                    // del_clusters.push_back(cluster);
                    // delete cluster;

                    if (sep_clusters.size() >= 2) {  // 1
                        for (size_t k = 2; k < sep_clusters.size(); k++) {
                            new_clusters.push_back(sep_clusters.at(k));
                        }
                        //	  std::cerr << em("sep sep1") << std::endl;

                        std::vector<Cluster *> temp_del_clusters;
                        Cluster *cluster2 = sep_clusters.at(1);
                        double length_1 = cluster2->get_length();

                        Cluster *final_sep_cluster = cluster2;

                        if (length_1 > 100 * units::cm) {
                            boundary_points.clear();
                            independent_points.clear();

                            if (JudgeSeparateDec_1(cluster2, drift_dir_abs, length_1) &&
                                JudgeSeparateDec_2(cluster2, dv, drift_dir_abs, boundary_points, independent_points,
                                                   length_1)) {
                                std::vector<Cluster *> sep_clusters =
                                    Separate_1(use_ctpc, cluster2, boundary_points, independent_points, length_1, vertical_dir, beam_dir, dv, pcts, scope);

                                // std::cout << "Separate Separate_1 1 for " << orig_nchildren << " " << " returned " << sep_clusters.size() << " clusters" << std::endl;

                                Cluster *cluster3 = sep_clusters.at(0);
                                new_clusters.push_back(cluster3);

                                temp_del_clusters.push_back(cluster2);
                                // delete cluster2;

                                if (sep_clusters.size() >= 2) {  // 2
                                    for (size_t k = 2; k < sep_clusters.size(); k++) {
                                        new_clusters.push_back(sep_clusters.at(k));
                                    }

                                    Cluster *cluster4 = sep_clusters.at(1);
                                    final_sep_cluster = cluster4;
                                    length_1 = cluster4->get_length();

                                    if (length_1 > 100 * units::cm) {
                                        boundary_points.clear();
                                        independent_points.clear();
                                        if (JudgeSeparateDec_1(cluster4, drift_dir_abs, length_1) &&
                                            JudgeSeparateDec_2(cluster4, dv, drift_dir_abs, boundary_points, independent_points,
                                                               length_1)) {
                                            
                                                                // std::cout << "Separate 3rd level" << std::endl;

                                            std::vector<Cluster *> sep_clusters = Separate_1(
                                                use_ctpc, cluster4, boundary_points, independent_points, length_1, vertical_dir, beam_dir, dv, pcts, scope);

                                            //		  std::cerr << em("sep sep3") << std::endl;

                                            Cluster *cluster5 = sep_clusters.at(0);
                                            new_clusters.push_back(cluster5);

                                            temp_del_clusters.push_back(cluster4);
                                            //		  delete cluster4;
                                            if (sep_clusters.size() >= 2) {  // 3
                                                for (size_t k = 2; k < sep_clusters.size(); k++) {
                                                    new_clusters.push_back(sep_clusters.at(k));
                                                }
                                                Cluster *cluster6 = sep_clusters.at(1);
                                                final_sep_cluster = cluster6;
                                            }
                                            else {
                                                final_sep_cluster = 0;
                                            }
                                        }
                                    }
                                }
                                else {
                                    final_sep_cluster = 0;
                                }
                            }
                        }

                        if (final_sep_cluster != 0) {  // 1
                            // temporary
                            length_1 = final_sep_cluster->get_length();

                            //	std::cout << length_1/units::cm << std::endl;

                            if (length_1 > 60 * units::cm) {
                                boundary_points.clear();
                                independent_points.clear();
                                JudgeSeparateDec_1(final_sep_cluster, drift_dir_abs, length_1);
                                JudgeSeparateDec_2(final_sep_cluster, dv, drift_dir_abs, boundary_points, independent_points,
                                                   length_1);
                                if (independent_points.size() > 0) {
                                    
                                    // std::cout << "Separate final one" << std::endl;

                                    std::vector<Cluster *> sep_clusters = Separate_1(
                                        use_ctpc, final_sep_cluster, boundary_points, independent_points, length_1, vertical_dir, beam_dir, dv, pcts, scope);
                                    //std::cout << "Separate Separate_1 2 for " << orig_nchildren << " " << " returned " << sep_clusters.size() << " clusters" << std::endl;
                                    //	      std::cerr << em("sep sep4") << std::endl;

                                    Cluster *cluster5 = sep_clusters.at(0);
                                    new_clusters.push_back(cluster5);

                                    temp_del_clusters.push_back(final_sep_cluster);
                                    // delete final_sep_cluster;

                                    if (sep_clusters.size() >= 2) {  // 4
                                        for (size_t k = 2; k < sep_clusters.size(); k++) {
                                            new_clusters.push_back(sep_clusters.at(k));
                                        }

                                        Cluster *cluster6 = sep_clusters.at(1);
                                        final_sep_cluster = cluster6;
                                    }
                                    else {
                                        final_sep_cluster = 0;
                                    }
                                }
                            }

                            if (final_sep_cluster != 0) {  // 2
                                // std::vector<Cluster *> final_sep_clusters = Separate_2(final_sep_cluster);
                                const auto b2id = Separate_2(final_sep_cluster, scope);
                                auto scope_transform = final_sep_cluster->get_scope_transform(scope);
                                auto final_sep_clusters = live_grouping.separate(final_sep_cluster,b2id,true); 
                                 
                                // Apply the scope filter settings to all new clusters
                                for (auto& [id, new_cluster] : final_sep_clusters) {
                                    new_cluster->set_scope_filter(scope, true);
                                    new_cluster->set_default_scope(scope);
                                    new_cluster->set_scope_transform(scope,scope_transform);
                                }

                                assert(final_sep_cluster == nullptr);
                     
                            }
                        }

                    }
                }
                else if (cluster->get_length() < 6 * units::m) {
                    // const size_t orig_nchildren = cluster->nchildren();
                    //std::cout << "Stripping Cluster with " << orig_nchildren << " blobs (ctpc) length " << cluster->get_length() << std::endl;
                    // std::cout << boundary_points.size() << " " << independent_points.size() << std::endl;
                    std::vector<Cluster *> sep_clusters =
                        Separate_1(use_ctpc, cluster, boundary_points, independent_points, cluster->get_length(), vertical_dir, beam_dir, dv, pcts, scope);
                    
                        // std::cout << "Stripping Separate_1 for " << orig_nchildren << " returned " << sep_clusters.size() << " clusters" << std::endl;

                    Cluster *cluster1 = sep_clusters.at(0);
                    new_clusters.push_back(cluster1);

                    del_clusters.push_back(cluster);
                    // delete cluster;

                    if (sep_clusters.size() >= 2) {
                        for (size_t k = 2; k < sep_clusters.size(); k++) {
                            new_clusters.push_back(sep_clusters.at(k));
                        }

                        //std::vector<Cluster *> temp_del_clusters;
                        Cluster *cluster2 = sep_clusters.at(1);
                        // double length_1 = cluster2->get_length();
                        Cluster *final_sep_cluster = cluster2;

                        // std::vector<Cluster *> final_sep_clusters = Separate_2(final_sep_cluster);
                        const auto b2id = Separate_2(final_sep_cluster, scope);
                        auto scope_transform = final_sep_cluster->get_scope_transform(scope);
                        auto final_sep_clusters = live_grouping.separate(final_sep_cluster, b2id, true);
                        // Apply the scope filter settings to all new clusters
                        for (auto& [id, new_cluster] : final_sep_clusters) {
                            new_cluster->set_scope_filter(scope, true);
                            new_cluster->set_default_scope(scope);
                            new_cluster->set_scope_transform(scope, scope_transform);
                        }
                        assert(final_sep_cluster == nullptr);
                        cluster2 = final_sep_cluster = sep_clusters[1] = nullptr;

                    }
                }  // else ...
            }
        }
    }


//  {
//    auto live_clusters = live_grouping.children(); // copy
//     // Process each cluster
//     for (size_t iclus = 0; iclus < live_clusters.size(); ++iclus) {
//         Cluster* cluster = live_clusters.at(iclus);
//         auto& scope = cluster->get_default_scope();
//         std::cout << "Test: " << iclus << " " << cluster->nchildren() << " " << scope.pcname << " " << scope.coords[0] << " " << scope.coords[1] << " " << scope.coords[2] << " " << cluster->get_scope_filter(scope)<< " " << cluster->get_pca().center << std::endl;
//     }
//   }







}

/// @brief PCA based, drift_dir +x, -x the same ...
bool WireCell::Clus::Facade::JudgeSeparateDec_1(const Cluster* cluster, const geo_point_t& drift_dir_abs, const double length)
{
    // get the main axis
    geo_point_t dir1(cluster->get_pca().axis.at(0).x(), cluster->get_pca().axis.at(0).y(), cluster->get_pca().axis.at(0).z());
    geo_point_t dir2(cluster->get_pca().axis.at(1).x(), cluster->get_pca().axis.at(1).y(), cluster->get_pca().axis.at(1).z());
    geo_point_t dir3(cluster->get_pca().axis.at(2).x(), cluster->get_pca().axis.at(2).y(), cluster->get_pca().axis.at(2).z());

    double angle1 = fabs(dir2.angle(drift_dir_abs) - 3.1415926 / 2.) / 3.1415926 * 180.;

    /// CHECKME: is "time_slice_length" drift_speed * tick?
    auto points = cluster->get_earliest_latest_points();
    double temp_angle1 = asin(fabs(points.first.x() - points.second.x()) / length) / 3.1415926 * 180.;

    double angle2 = fabs(dir3.angle(drift_dir_abs) - 3.1415926 / 2.) / 3.1415926 * 180.;
    double ratio1 = cluster->get_pca().values.at(1) / cluster->get_pca().values.at(0);
    double ratio2 = cluster->get_pca().values.at(2) / cluster->get_pca().values.at(0);

    // std::cout << ratio1 << " " <<  pow(10, exp(1.38115 - 1.19312 * pow(angle1, 1. / 3.)) - 2.2)  << " "  << ratio1 << " " << pow(10, exp(1.38115 - 1.19312 * pow(temp_angle1, 1. / 3.)) - 2.2) << " " << ratio2 << " " << pow(10, exp(1.38115 - 1.19312 * pow(angle2, 1. / 3.)) - 2.2) << " " << ratio1 << " " << angle1 << " " << angle2 << std::endl;
    
    // add some additional check ... PDHD
    if (angle1 < 5 && ratio2 < 0.05 && cluster->get_length() > 300*units::cm) return false;

    if (ratio1 > pow(10, exp(1.38115 - 1.19312 * pow(angle1, 1. / 3.)) - 2.2) ||
        ratio1 > pow(10, exp(1.38115 - 1.19312 * pow(temp_angle1, 1. / 3.)) - 2.2) ||
        ratio2 > pow(10, exp(1.38115 - 1.19312 * pow(angle2, 1. / 3.)) - 2.2) || 
        ratio1 > 0.75)
        return true;
    return false;
}

bool WireCell::Clus::Facade::JudgeSeparateDec_2(const Cluster* cluster, const IDetectorVolumes::pointer dv, const geo_point_t& drift_dir_abs,
                               std::vector<geo_point_t>& boundary_points, std::vector<geo_point_t>& independent_points,
                               const double cluster_length)
{
    // const auto &mp = cluster->grouping()->get_params();

    auto wpids = cluster->grouping()->wpids();
    std::map<WirePlaneId, double> map_FV_xmin;
    std::map<WirePlaneId, double> map_FV_xmax;
    std::map<WirePlaneId, double> map_FV_xmin_margin;
    std::map<WirePlaneId, double> map_FV_xmax_margin;
    for (const auto& wpid : wpids) {
        map_FV_xmin[wpid] = dv->metadata(wpid)["FV_xmin"].asDouble() ;
        map_FV_xmax[wpid] = dv->metadata(wpid)["FV_xmax"].asDouble() ;
        map_FV_xmin_margin[wpid] = dv->metadata(wpid)["FV_xmin_margin"].asDouble() ;
        map_FV_xmax_margin[wpid] = dv->metadata(wpid)["FV_xmax_margin"].asDouble() ;
    }
    WirePlaneId wpid_all(0);
    double det_FV_xmin = dv->metadata(wpid_all)["FV_xmin"].asDouble();
    double det_FV_xmax = dv->metadata(wpid_all)["FV_xmax"].asDouble();
    double det_FV_ymin = dv->metadata(wpid_all)["FV_ymin"].asDouble();
    double det_FV_ymax = dv->metadata(wpid_all)["FV_ymax"].asDouble();
    double det_FV_zmin = dv->metadata(wpid_all)["FV_zmin"].asDouble();
    double det_FV_zmax = dv->metadata(wpid_all)["FV_zmax"].asDouble();
    double det_FV_xmin_margin = dv->metadata(wpid_all)["FV_xmin_margin"].asDouble();
    double det_FV_xmax_margin = dv->metadata(wpid_all)["FV_xmax_margin"].asDouble();
    // double det_FV_ymin_margin = dv->metadata(wpid_all)["FV_ymin_margin"].asDouble();
    double det_FV_ymax_margin = dv->metadata(wpid_all)["FV_ymax_margin"].asDouble();
    double det_FV_zmin_margin = dv->metadata(wpid_all)["FV_zmin_margin"].asDouble();
    double det_FV_zmax_margin = dv->metadata(wpid_all)["FV_zmax_margin"].asDouble();

    boundary_points = cluster->get_hull();
    std::vector<geo_point_t> hy_points;
    std::vector<geo_point_t> ly_points;
    std::vector<geo_point_t> hz_points;
    std::vector<geo_point_t> lz_points;
    std::vector<geo_point_t> hx_points;
    std::vector<geo_point_t> lx_points;

    std::set<int> independent_surfaces;

    for (size_t j = 0; j != boundary_points.size(); j++) {
        if (j == 0) {
            hy_points.push_back(boundary_points.at(j));
            ly_points.push_back(boundary_points.at(j));
            hz_points.push_back(boundary_points.at(j));
            lz_points.push_back(boundary_points.at(j));
            hx_points.push_back(boundary_points.at(j));
            lx_points.push_back(boundary_points.at(j));
        }
        else {
            geo_point_t test_p(boundary_points.at(j).x(), boundary_points.at(j).y(), boundary_points.at(j).z());
            if (cluster->nnearby(test_p, 15 * units::cm) > 75) {
                if (boundary_points.at(j).y() > hy_points.at(0).y()) hy_points.at(0) = boundary_points.at(j);
                if (boundary_points.at(j).y() < ly_points.at(0).y()) ly_points.at(0) = boundary_points.at(j);
                if (boundary_points.at(j).x() > hx_points.at(0).x()) hx_points.at(0) = boundary_points.at(j);
                if (boundary_points.at(j).x() < lx_points.at(0).x()) lx_points.at(0) = boundary_points.at(j);
                if (boundary_points.at(j).z() > hz_points.at(0).z()) hz_points.at(0) = boundary_points.at(j);
                if (boundary_points.at(j).z() < lz_points.at(0).z()) lz_points.at(0) = boundary_points.at(j);
            }
        }
    }

    bool flag_outx = false;
    if (hx_points.at(0).x() > det_FV_xmax + det_FV_xmax_margin || lx_points.at(0).x() < det_FV_xmin - det_FV_xmin_margin) flag_outx = true;

    if (hy_points.at(0).y() > det_FV_ymax) {
        for (size_t j = 0; j != boundary_points.size(); j++) {
            if (boundary_points.at(j).y() > det_FV_ymax) {
                bool flag_save = true;
                for (size_t k = 0; k != hy_points.size(); k++) {
                    double dis = sqrt(pow(hy_points.at(k).x() - boundary_points.at(j).x(), 2) +
                                      pow(hy_points.at(k).y() - boundary_points.at(j).y(), 2) +
                                      pow(hy_points.at(k).z() - boundary_points.at(j).z(), 2));
                    if (dis < 25 * units::cm) {
                        if (boundary_points.at(j).y() > hy_points.at(k).y()) hy_points.at(k) = boundary_points.at(j);
                        flag_save = false;
                    }
                }
                if (flag_save) hy_points.push_back(boundary_points.at(j));
            }
        }
    }

    if (ly_points.at(0).y() < det_FV_ymin) {
        for (size_t j = 0; j != boundary_points.size(); j++) {
            if (boundary_points.at(j).y() < det_FV_ymin) {
                bool flag_save = true;
                for (size_t k = 0; k != ly_points.size(); k++) {
                    double dis = sqrt(pow(ly_points.at(k).x() - boundary_points.at(j).x(), 2) +
                                      pow(ly_points.at(k).y() - boundary_points.at(j).y(), 2) +
                                      pow(ly_points.at(k).z() - boundary_points.at(j).z(), 2));
                    if (dis < 25 * units::cm) {
                        if (boundary_points.at(j).y() < ly_points.at(k).y()) ly_points.at(k) = boundary_points.at(j);
                        flag_save = false;
                    }
                }
                if (flag_save) ly_points.push_back(boundary_points.at(j));
            }
        }
    }
    if (hz_points.at(0).z() > det_FV_zmax) {
        for (size_t j = 0; j != boundary_points.size(); j++) {
            if (boundary_points.at(j).z() > det_FV_zmax) {
                bool flag_save = true;
                for (size_t k = 0; k != hz_points.size(); k++) {
                    double dis = sqrt(pow(hz_points.at(k).x() - boundary_points.at(j).x(), 2) +
                                      pow(hz_points.at(k).y() - boundary_points.at(j).y(), 2) +
                                      pow(hz_points.at(k).z() - boundary_points.at(j).z(), 2));
                    if (dis < 25 * units::cm) {
                        if (boundary_points.at(j).z() > hz_points.at(k).z()) hz_points.at(k) = boundary_points.at(j);
                        flag_save = false;
                    }
                }
                if (flag_save) hz_points.push_back(boundary_points.at(j));
            }
        }
    }
    if (lz_points.at(0).z() < det_FV_zmin) {
        for (size_t j = 0; j != boundary_points.size(); j++) {
            if (boundary_points.at(j).z() < det_FV_zmin) {
                bool flag_save = true;
                for (size_t k = 0; k != lz_points.size(); k++) {
                    double dis = sqrt(pow(lz_points.at(k).x() - boundary_points.at(j).x(), 2) +
                                      pow(lz_points.at(k).y() - boundary_points.at(j).y(), 2) +
                                      pow(lz_points.at(k).z() - boundary_points.at(j).z(), 2));
                    if (dis < 25 * units::cm) {
                        if (boundary_points.at(j).z() < lz_points.at(k).z()) lz_points.at(k) = boundary_points.at(j);
                        flag_save = false;
                    }
                }
                if (flag_save) lz_points.push_back(boundary_points.at(j));
            }
        }
    }

    int num_outside_points = 0;
    int num_outx_points = 0;

    for (size_t j = 0; j != hy_points.size(); j++) {
        if (hy_points.at(j).x() >= det_FV_xmin && hy_points.at(j).x() <= det_FV_xmax &&
            hy_points.at(j).y() >= det_FV_ymin && hy_points.at(j).y() <= det_FV_ymax &&
            hy_points.at(j).z() >= det_FV_zmin && hy_points.at(j).z() <= det_FV_zmax && (!flag_outx))
            continue;

        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis = sqrt(pow(hy_points.at(j).x() - independent_points.at(k).x(), 2) +
                              pow(hy_points.at(j).y() - independent_points.at(k).y(), 2) +
                              pow(hy_points.at(j).z() - independent_points.at(k).z(), 2));
            if (dis < 15 * units::cm) flag_save = false;
        }
        if (flag_save) {
            independent_points.push_back(hy_points.at(j));
            if (hy_points.at(j).y() > det_FV_ymax + det_FV_ymax_margin) {
                independent_surfaces.insert(0);
            }
            else if (hy_points.at(j).y() < det_FV_ymin) {
                independent_surfaces.insert(1);
            }
            else if (hy_points.at(j).z() > det_FV_zmax + det_FV_zmax_margin) {
                independent_surfaces.insert(2);
            }
            else if (hy_points.at(j).z() < det_FV_zmin - det_FV_zmin_margin) {
                independent_surfaces.insert(3);
            }
            else if (hy_points.at(j).x() > det_FV_xmax) {
                independent_surfaces.insert(4);
            }
            else if (hy_points.at(j).x() < det_FV_xmin) {
                independent_surfaces.insert(5);
            }

            if (hy_points.at(j).y() > det_FV_ymax + det_FV_ymax_margin || hy_points.at(j).y() < det_FV_ymin ||
                hy_points.at(j).z() < det_FV_zmin - det_FV_zmin_margin || hy_points.at(j).z() > det_FV_zmax + det_FV_zmax_margin ||
                hy_points.at(j).x() < det_FV_xmin || hy_points.at(j).x() > det_FV_xmax)
                num_outside_points++;
            if (hy_points.at(j).x() < det_FV_xmin - det_FV_xmin_margin || hy_points.at(j).x() > det_FV_xmax - det_FV_xmax_margin) num_outx_points++;
        }
    }
    for (size_t j = 0; j != ly_points.size(); j++) {
        if (ly_points.at(j).x() >= det_FV_xmin && ly_points.at(j).x() <= det_FV_xmax &&
            ly_points.at(j).y() >= det_FV_ymin && ly_points.at(j).y() <= det_FV_ymax &&
            ly_points.at(j).z() >= det_FV_zmin && ly_points.at(j).z() <= det_FV_zmax && (!flag_outx))
            continue;

        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis = sqrt(pow(ly_points.at(j).x() - independent_points.at(k).x(), 2) +
                              pow(ly_points.at(j).y() - independent_points.at(k).y(), 2) +
                              pow(ly_points.at(j).z() - independent_points.at(k).z(), 2));
            if (dis < 15 * units::cm) flag_save = false;
        }
        if (flag_save) {
            independent_points.push_back(ly_points.at(j));

            if (ly_points.at(j).y() < det_FV_ymin) {
                independent_surfaces.insert(1);
            }
            else if (ly_points.at(j).y() > det_FV_ymax + det_FV_ymax_margin) {
                independent_surfaces.insert(0);
            }
            else if (ly_points.at(j).z() > det_FV_zmax + det_FV_zmax_margin) {
                independent_surfaces.insert(2);
            }
            else if (ly_points.at(j).z() < det_FV_zmin - det_FV_zmin_margin) {
                independent_surfaces.insert(3);
            }
            else if (ly_points.at(j).x() > det_FV_xmax) {
                independent_surfaces.insert(4);
            }
            else if (ly_points.at(j).x() < det_FV_xmin) {
                independent_surfaces.insert(5);
            }

            if (ly_points.at(j).y() > det_FV_ymax + det_FV_ymax_margin || ly_points.at(j).y() < det_FV_ymin ||
                ly_points.at(j).z() < det_FV_zmin - det_FV_zmin_margin || ly_points.at(j).z() > det_FV_zmax + det_FV_zmax_margin ||
                ly_points.at(j).x() < det_FV_xmin || ly_points.at(j).x() > det_FV_xmax)
                num_outside_points++;
            if (ly_points.at(j).x() < det_FV_xmin - det_FV_xmin_margin || ly_points.at(j).x() > det_FV_xmax - det_FV_xmax_margin) num_outx_points++;
        }
    }
    for (size_t j = 0; j != hz_points.size(); j++) {
        if (hz_points.at(j).x() >= det_FV_xmin && hz_points.at(j).x() <= det_FV_xmax &&
            hz_points.at(j).y() >= det_FV_ymin && hz_points.at(j).y() <= det_FV_ymax &&
            hz_points.at(j).z() >= det_FV_zmin && hz_points.at(j).z() <= det_FV_zmax && (!flag_outx))
            continue;

        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis = sqrt(pow(hz_points.at(j).x() - independent_points.at(k).x(), 2) +
                              pow(hz_points.at(j).y() - independent_points.at(k).y(), 2) +
                              pow(hz_points.at(j).z() - independent_points.at(k).z(), 2));
            if (dis < 15 * units::cm) flag_save = false;
        }
        if (flag_save) {
            independent_points.push_back(hz_points.at(j));

            if (hz_points.at(j).z() > det_FV_zmax + det_FV_zmax_margin) {
                independent_surfaces.insert(2);
            }
            else if (hz_points.at(j).z() < det_FV_zmin - det_FV_zmin_margin) {
                independent_surfaces.insert(3);
            }
            else if (hz_points.at(j).y() > det_FV_ymax + det_FV_ymax_margin) {
                independent_surfaces.insert(0);
            }
            else if (hz_points.at(j).y() < det_FV_ymin) {
                independent_surfaces.insert(1);
            }
            else if (hz_points.at(j).x() > det_FV_xmax) {
                independent_surfaces.insert(4);
            }
            else if (hz_points.at(j).x() < det_FV_xmin) {
                independent_surfaces.insert(5);
            }

            if (hz_points.at(j).y() > det_FV_ymax + det_FV_ymax_margin || hz_points.at(j).y() < det_FV_ymin ||
                hz_points.at(j).z() < det_FV_zmin - det_FV_zmin_margin || hz_points.at(j).z() > det_FV_zmax + det_FV_zmax_margin ||
                hz_points.at(j).x() < det_FV_xmin || hz_points.at(j).x() > det_FV_xmax)
                num_outside_points++;
            if (hz_points.at(j).x() < det_FV_xmin - det_FV_xmin_margin || hz_points.at(j).x() > det_FV_xmax - det_FV_xmax_margin) num_outx_points++;
        }
    }
    for (size_t j = 0; j != lz_points.size(); j++) {
        if (lz_points.at(j).x() >= det_FV_xmin && lz_points.at(j).x() <= det_FV_xmax &&
            lz_points.at(j).y() >= det_FV_ymin && lz_points.at(j).y() <= det_FV_ymax &&
            lz_points.at(j).z() >= det_FV_zmin && lz_points.at(j).z() <= det_FV_zmax && (!flag_outx))
            continue;

        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis = sqrt(pow(lz_points.at(j).x() - independent_points.at(k).x(), 2) +
                              pow(lz_points.at(j).y() - independent_points.at(k).y(), 2) +
                              pow(lz_points.at(j).z() - independent_points.at(k).z(), 2));
            if (dis < 15 * units::cm) flag_save = false;
        }
        if (flag_save) {
            independent_points.push_back(lz_points.at(j));

            if (lz_points.at(j).z() < det_FV_zmin - det_FV_zmin_margin) {
                independent_surfaces.insert(3);
            }
            else if (lz_points.at(j).z() > det_FV_zmax + det_FV_zmax_margin) {
                independent_surfaces.insert(2);
            }
            else if (lz_points.at(j).y() > det_FV_ymax + det_FV_ymax_margin) {
                independent_surfaces.insert(0);
            }
            else if (lz_points.at(j).y() < det_FV_ymin) {
                independent_surfaces.insert(1);
            }
            else if (lz_points.at(j).x() > det_FV_xmax) {
                independent_surfaces.insert(4);
            }
            else if (lz_points.at(j).x() < det_FV_xmin) {
                independent_surfaces.insert(5);
            }

            if (lz_points.at(j).y() > det_FV_ymax + det_FV_ymax_margin || lz_points.at(j).y() < det_FV_ymin ||
                lz_points.at(j).z() < det_FV_zmin - det_FV_zmin_margin || lz_points.at(j).z() > det_FV_zmax + det_FV_zmax_margin ||
                lz_points.at(j).x() < det_FV_xmin || lz_points.at(j).x() > det_FV_xmax)
                num_outside_points++;
            if (lz_points.at(j).x() < det_FV_xmin - det_FV_xmin_margin || lz_points.at(j).x() > det_FV_xmax - det_FV_xmax_margin) num_outx_points++;
        }
    }
    for (size_t j = 0; j != hx_points.size(); j++) {
        if (hx_points.at(j).x() >= det_FV_xmin && hx_points.at(j).x() <= det_FV_xmax &&
            hx_points.at(j).y() >= det_FV_ymin && hx_points.at(j).y() <= det_FV_ymax &&
            hx_points.at(j).z() >= det_FV_zmin && hx_points.at(j).z() <= det_FV_zmax && (!flag_outx))
            continue;

        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis = sqrt(pow(hx_points.at(j).x() - independent_points.at(k).x(), 2) +
                              pow(hx_points.at(j).y() - independent_points.at(k).y(), 2) +
                              pow(hx_points.at(j).z() - independent_points.at(k).z(), 2));
            if (dis < 15 * units::cm) flag_save = false;
        }
        if (flag_save) {
            independent_points.push_back(hx_points.at(j));

            if (hx_points.at(j).y() > det_FV_ymax + det_FV_ymax_margin || hx_points.at(j).y() < det_FV_ymin ||
                hx_points.at(j).z() < det_FV_zmin - det_FV_zmin_margin || hx_points.at(j).z() > det_FV_zmax + det_FV_zmax_margin ||
                hx_points.at(j).x() < det_FV_xmin || hx_points.at(j).x() > det_FV_xmax)
                num_outside_points++;
            if (hx_points.at(j).x() < det_FV_xmin - det_FV_xmin_margin || hx_points.at(j).x() > det_FV_xmax - det_FV_xmax_margin) {
                num_outx_points++;
            }

            if (lx_points.at(j).x() > det_FV_xmax) {
                independent_surfaces.insert(4);
            }
            else if (lx_points.at(j).x() < det_FV_xmin) {
                independent_surfaces.insert(5);
            }
            else if (lx_points.at(j).y() > det_FV_ymax + det_FV_ymax_margin) {
                independent_surfaces.insert(0);
            }
            else if (lx_points.at(j).y() < det_FV_ymin) {
                independent_surfaces.insert(1);
            }
            else if (lx_points.at(j).z() > det_FV_zmax + det_FV_zmax_margin) {
                independent_surfaces.insert(2);
            }
            else if (lx_points.at(j).z() < det_FV_zmin - det_FV_zmin_margin) {
                independent_surfaces.insert(3);
            }
        }
    }
    for (size_t j = 0; j != lx_points.size(); j++) {
        if (lx_points.at(j).x() >= det_FV_xmin && lx_points.at(j).x() <= det_FV_xmax &&
            lx_points.at(j).y() >= det_FV_ymin && lx_points.at(j).y() <= det_FV_ymax &&
            lx_points.at(j).z() >= det_FV_zmin && lx_points.at(j).z() <= det_FV_zmax && (!flag_outx))
            continue;

        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis = sqrt(pow(lx_points.at(j).x() - independent_points.at(k).x(), 2) +
                              pow(lx_points.at(j).y() - independent_points.at(k).y(), 2) +
                              pow(lx_points.at(j).z() - independent_points.at(k).z(), 2));
            if (dis < 15 * units::cm) flag_save = false;
        }
        if (flag_save) {
            independent_points.push_back(lx_points.at(j));

            if (lx_points.at(j).y() > det_FV_ymax + det_FV_ymax_margin || lx_points.at(j).y() < det_FV_ymin ||
                lx_points.at(j).z() < det_FV_zmin - det_FV_zmin_margin || lx_points.at(j).z() > det_FV_zmax + det_FV_zmax_margin ||
                lx_points.at(j).x() < det_FV_xmin || lx_points.at(j).x() > det_FV_xmax)
                num_outside_points++;
            if (lx_points.at(j).x() < det_FV_xmin - det_FV_xmin_margin || lx_points.at(j).x() > det_FV_xmax - det_FV_xmax_margin) {
                num_outx_points++;
            }

            if (lx_points.at(j).x() < det_FV_xmin) {
                independent_surfaces.insert(5);
            }
            else if (lx_points.at(j).x() > det_FV_xmax) {
                independent_surfaces.insert(4);
            }
            else if (lx_points.at(j).y() > det_FV_ymax + det_FV_ymax_margin) {
                independent_surfaces.insert(0);
            }
            else if (lx_points.at(j).y() < det_FV_ymin) {
                independent_surfaces.insert(1);
            }
            else if (lx_points.at(j).z() > det_FV_zmax + det_FV_zmax_margin) {
                independent_surfaces.insert(2);
            }
            else if (lx_points.at(j).z() < det_FV_zmin - det_FV_zmin_margin) {
                independent_surfaces.insert(3);
            }
        }
    }

    int num_far_points = 0;

    if (independent_points.size() == 2 && (independent_surfaces.size() > 1 || flag_outx)) {
        geo_vector_t dir_1(independent_points.at(1).x() - independent_points.at(0).x(),
                           independent_points.at(1).y() - independent_points.at(0).y(),
                           independent_points.at(1).z() - independent_points.at(0).z());
        dir_1 = dir_1.norm();
        for (size_t j = 0; j != boundary_points.size(); j++) {
            geo_vector_t dir_2(boundary_points.at(j).x() - independent_points.at(0).x(),
                               boundary_points.at(j).y() - independent_points.at(0).y(),
                               boundary_points.at(j).z() - independent_points.at(0).z());
            double angle_12 = dir_1.angle(dir_2);
            geo_vector_t dir_3 = dir_2 - dir_1 * dir_2.magnitude() * cos(angle_12);
            double angle_3 = dir_3.angle(drift_dir_abs);
            // std::cout << dir_3.Mag()/units::cm << " " << fabs(angle_3-3.1415926/2.)/3.1415926*180. << " " <<
            // fabs(dir_3.X()/units::cm) << std::endl;
            if (fabs(angle_3 - 3.1415926 / 2.) / 3.1415926 * 180. < 7.5) {
                if (fabs(dir_3.x() / units::cm) > 14 * units::cm) num_far_points++;
                if (fabs(dir_1.angle(drift_dir_abs) - 3.1415926 / 2.) / 3.1415926 * 180. > 15) {
                    if (dir_3.magnitude() > 20 * units::cm) num_far_points++;
                }
            }
            else {
                if (dir_3.magnitude() > 20 * units::cm) num_far_points++;
            }
        }

        // find the middle points and close distance ...
        geo_point_t middle_point((independent_points.at(1).x() + independent_points.at(0).x()) / 2.,
                                 (independent_points.at(1).y() + independent_points.at(0).y()) / 2.,
                                 (independent_points.at(1).z() + independent_points.at(0).z()) / 2.);
        double middle_dis = cluster->get_closest_dis(middle_point);
        // std::cout << middle_dis/units::cm << " " << num_far_points << std::endl;
        if (middle_dis > 25 * units::cm) {
            num_far_points = 0;
        }
    }

    double max_x = -1e9, min_x = 1e9;
    double max_y = -1e9, min_y = 1e9;
    double max_z = -1e9, min_z = 1e9;
    for (auto it = independent_points.begin(); it != independent_points.end(); it++) {
        if ((*it).x() > max_x) max_x = (*it).x();
        if ((*it).x() < min_x) min_x = (*it).x();
        if ((*it).y() > max_y) max_y = (*it).y();
        if ((*it).y() < min_y) min_y = (*it).y();
        if ((*it).z() > max_z) max_z = (*it).z();
        if ((*it).z() < min_z) min_z = (*it).z();
        // std::cout << (*it).x()/units::cm << " " << (*it).y()/units::cm << " " << (*it).z()/units::cm << std::endl;
    }
    if (hx_points.size() > 0) {
        if (hx_points.at(0).x() > max_x + 10 * units::cm) max_x = hx_points.at(0).x();
    }
    if (lx_points.size() > 0) {
        if (lx_points.at(0).x() < min_x - 10 * units::cm) min_x = lx_points.at(0).x();
    }

    if (max_x - min_x < 2.5 * units::cm &&
        sqrt(pow(max_y - min_y, 2) + pow(max_z - min_z, 2) + pow(max_x - min_x, 2)) > 150 * units::cm) {
        independent_points.clear();
        return false;
    }
    if (max_x - min_x < 2.5 * units::cm && independent_points.size() == 2 && num_outx_points == 0) {
        independent_points.clear();
        return false;
    }

    if ((num_outside_points > 1 && independent_surfaces.size() > 1 ||
         num_outside_points > 2 && cluster_length > 250 * units::cm || num_outx_points > 0) &&
        (independent_points.size() > 2 || independent_points.size() == 2 && num_far_points > 0))
        return true;

    // about to return false ...
    independent_points.clear();

    for (size_t j = 0; j != hy_points.size(); j++) {
        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis = sqrt(pow(hy_points.at(j).x() - independent_points.at(k).x(), 2) +
                              pow(hy_points.at(j).y() - independent_points.at(k).y(), 2) +
                              pow(hy_points.at(j).z() - independent_points.at(k).z(), 2));
            if (dis < 15 * units::cm) flag_save = false;
        }
        if (flag_save) independent_points.push_back(hy_points.at(j));
    }

    for (size_t j = 0; j != ly_points.size(); j++) {
        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis = sqrt(pow(ly_points.at(j).x() - independent_points.at(k).x(), 2) +
                              pow(ly_points.at(j).y() - independent_points.at(k).y(), 2) +
                              pow(ly_points.at(j).z() - independent_points.at(k).z(), 2));
            if (dis < 15 * units::cm) flag_save = false;
        }
        if (flag_save) independent_points.push_back(ly_points.at(j));
    }

    for (size_t j = 0; j != hx_points.size(); j++) {
        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis = sqrt(pow(hx_points.at(j).x() - independent_points.at(k).x(), 2) +
                              pow(hx_points.at(j).y() - independent_points.at(k).y(), 2) +
                              pow(hx_points.at(j).z() - independent_points.at(k).z(), 2));
            if (dis < 15 * units::cm) flag_save = false;
        }
        if (flag_save) independent_points.push_back(hx_points.at(j));
    }

    for (size_t j = 0; j != lx_points.size(); j++) {
        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis = sqrt(pow(lx_points.at(j).x() - independent_points.at(k).x(), 2) +
                              pow(lx_points.at(j).y() - independent_points.at(k).y(), 2) +
                              pow(lx_points.at(j).z() - independent_points.at(k).z(), 2));
            if (dis < 15 * units::cm) flag_save = false;
        }
        if (flag_save) independent_points.push_back(lx_points.at(j));
    }

    for (size_t j = 0; j != hz_points.size(); j++) {
        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis = sqrt(pow(hz_points.at(j).x() - independent_points.at(k).x(), 2) +
                              pow(hz_points.at(j).y() - independent_points.at(k).y(), 2) +
                              pow(hz_points.at(j).z() - independent_points.at(k).z(), 2));
            if (dis < 15 * units::cm) flag_save = false;
        }
        if (flag_save) independent_points.push_back(hz_points.at(j));
    }

    for (size_t j = 0; j != lz_points.size(); j++) {
        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis = sqrt(pow(lz_points.at(j).x() - independent_points.at(k).x(), 2) +
                              pow(lz_points.at(j).y() - independent_points.at(k).y(), 2) +
                              pow(lz_points.at(j).z() - independent_points.at(k).z(), 2));
            if (dis < 15 * units::cm) flag_save = false;
        }
        if (flag_save) independent_points.push_back(lz_points.at(j));
    }

    return false;
}

#define _INDEV_
#ifdef _INDEV_

std::vector<Cluster *> WireCell::Clus::Facade::Separate_1(const bool use_ctpc, Cluster *cluster,
                                                     std::vector<geo_point_t> &boundary_points,
                                                     std::vector<geo_point_t> &independent_points,
                                                     double length, geo_point_t dir_cosmic, geo_point_t dir_beam, const IDetectorVolumes::pointer dv, const IPCTransformSet::pointer pcts, const Tree::Scope& scope)
{
    
    auto* grouping = cluster->grouping();

    auto gwpids = grouping->wpids();

    std::map<int, std::map<int, std::map<int, std::pair<double, double>>>> af_dead_u_index ;
    std::map<int, std::map<int, std::map<int, std::pair<double, double>>>> af_dead_v_index ;
    std::map<int, std::map<int, std::map<int, std::pair<double, double>>>> af_dead_w_index ;
    std::map<int, std::map<int, std::shared_ptr<Multi2DPointCloud>>> af_temp_cloud;
    for (auto wpid : gwpids) {
        int apa = wpid.apa();
        int face = wpid.face();
        af_dead_u_index[apa][face] = grouping->get_dead_winds(apa, face, 0); // raw 
        af_dead_v_index[apa][face] = grouping->get_dead_winds(apa, face, 1); // raw
        af_dead_w_index[apa][face] = grouping->get_dead_winds(apa, face, 2); // raw 

        // Create wpids for all three planes with this APA and face
        WirePlaneId wpid_u(kUlayer, face, apa);
        WirePlaneId wpid_v(kVlayer, face, apa);
        WirePlaneId wpid_w(kWlayer, face, apa);
     
        // Get wire directions for all planes
        Vector wire_dir_u = dv->wire_direction(wpid_u);
        Vector wire_dir_v = dv->wire_direction(wpid_v);
        Vector wire_dir_w = dv->wire_direction(wpid_w);

        // Calculate angles
        double angle_u = std::atan2(wire_dir_u.z(), wire_dir_u.y());
        double angle_v = std::atan2(wire_dir_v.z(), wire_dir_v.y());
        double angle_w = std::atan2(wire_dir_w.z(), wire_dir_w.y());

        af_temp_cloud[apa][face] = std::make_shared<Multi2DPointCloud>(angle_u, angle_v, angle_w); // 2D Dynamic Point Cloud
    }

    // std::cout << "Test: " << pc_name << " " << coords[0] << " " << coords[1] << " " << coords[2] << std::endl;

    geo_point_t cluster_center = cluster->get_pca().center;
    geo_point_t main_dir = cluster->get_pca().axis.at(0);
    geo_point_t second_dir = cluster->get_pca().axis.at(1);
    

    // special case, if one of the cosmic is very close to the beam direction
    if (cluster->get_pca().values.at(1) > 0.08 * cluster->get_pca().values.at(0) &&
        fabs(main_dir.angle(dir_beam) - 3.1415926 / 2.) > 75 / 180. * 3.1415926 &&
        fabs(second_dir.angle(dir_cosmic) - 3.1415926 / 2.) > 60 / 180. * 3.1415926) {
        main_dir = second_dir;
    }

    main_dir = main_dir.norm();
    if (main_dir.y() > 0)
        main_dir = main_dir * -1;  // make sure it is pointing down????

    geo_point_t start_wcpoint;
    geo_point_t end_wcpoint;
    geo_point_t drift_dir_abs(1, 0, 0);
    geo_point_t dir;

    double min_dis = 1e9;
    // double max_pca_dis;
    int min_index = 0;
    double max_dis = -1e9;
    // double min_pca_dis;
    int max_index = 0;
    // if (flag_debug_porting) {
    //     std::cout << "independent_points.size(): " << independent_points.size() << std::endl;
    // }
    for (size_t j = 0; j != independent_points.size(); j++) {
        geo_point_t dir(independent_points.at(j).x() - cluster_center.x(),
                        independent_points.at(j).y() - cluster_center.y(),
                        independent_points.at(j).z() - cluster_center.z());
        geo_point_t temp_p(independent_points.at(j).x(), independent_points.at(j).y(), independent_points.at(j).z());
        double dis = dir.dot(main_dir);
        // double dis_to_pca = dir.cross(main_dir).magnitude();
        // std::cout << j << " " << dis << " " << dir.Mag() << " " << sqrt(dir.Mag()*dir.Mag() - dis*dis) << std::endl;
        bool flag_connect = false;
        int num_points = cluster->nnearby(temp_p, 15 * units::cm);
        if (num_points > 100) {
            flag_connect = true;
        }
        else if (num_points > 75) {
            num_points = cluster->nnearby(temp_p, 30 * units::cm);
            if (num_points > 160) flag_connect = true;
        }

        // std::cout << dis / units::cm << " A " << cluster->get_num_points(temp_p,15*units::cm) << " " <<
        // cluster->get_num_points(temp_p,30*units::cm)  << std::endl;
        if (dis < min_dis && flag_connect) {
            min_dis = dis;
            min_index = j;
            // min_pca_dis = dis_to_pca;
        }
        if (dis > max_dis && flag_connect) {
            max_dis = dis;
            max_index = j;
            // max_pca_dis = dis_to_pca;
        }
    }
 


    size_t start_wcpoint_idx = 0;
    size_t end_wcpoint_idx = 0;
    {
        start_wcpoint = independent_points.at(min_index);

        // change direction if certain thing happened ...

        {
            geo_point_t p1(independent_points.at(max_index).x(), independent_points.at(max_index).y(),
                           independent_points.at(max_index).z());
            geo_point_t p2(independent_points.at(min_index).x(), independent_points.at(min_index).y(),
                           independent_points.at(min_index).z());
            geo_point_t temp_dir1 = cluster->vhough_transform(p1, 15 * units::cm);
            geo_point_t temp_dir2 = cluster->vhough_transform(p2, 15 * units::cm);

            bool flag_change = false;

            if (fabs(temp_dir1.angle(main_dir) - 3.1415926 / 2.) > fabs(temp_dir2.angle(main_dir) - 3.1415926 / 2.)) {
                if (fabs(temp_dir2.angle(main_dir) - 3.1415926 / 2.) > 80 / 180. * 3.1415926 &&
                    fabs(temp_dir1.angle(main_dir) - 3.1415926 / 2.) <
                        fabs(temp_dir2.angle(main_dir) - 3.1415926 / 2.) + 2.5 / 180. * 3.1415926) {
                }
                else {
                    flag_change = true;
                    start_wcpoint = independent_points.at(max_index);
                    main_dir = main_dir * -1;
                    max_index = min_index;
                }
            }

            if ((!flag_change) &&
                fabs(temp_dir1.angle(drift_dir_abs) - 3.1415926 / 2.) > fabs(temp_dir2.angle(drift_dir_abs) - 3.1415926 / 2.) &&
                fabs(temp_dir2.angle(drift_dir_abs) - 3.1415926 / 2.) / 3.1415926 * 180. < 10 &&
                fabs(temp_dir2.angle(main_dir) - 3.1415926 / 2.) / 3.1415926 * 180. < 80) {
                start_wcpoint = independent_points.at(max_index);
                main_dir = main_dir * -1;
                max_index = min_index;
            }

            if ((!flag_change) && fabs(temp_dir2.angle(drift_dir_abs) - 3.1415926 / 2.) < 1. / 180. * 3.1415926 &&
                fabs(temp_dir1.angle(drift_dir_abs) - 3.1415926 / 2.) > 3. / 180. * 3.1415926 &&
                fabs(temp_dir1.angle(main_dir) - 3.1415926 / 2.) / 3.1415926 * 180. > 70) {
                start_wcpoint = independent_points.at(max_index);
                main_dir = main_dir * -1;
                max_index = min_index;
            }
        }

        geo_point_t start_point(start_wcpoint.x(), start_wcpoint.y(), start_wcpoint.z());
        {
            // geo_point_t drift_dir_abs(1, 0, 0);
            dir = cluster->vhough_transform(start_point, 100 * units::cm);
            geo_point_t dir1 = cluster->vhough_transform(start_point, 30 * units::cm);
            if (dir.angle(dir1) > 20 * 3.1415926 / 180.) {
                if (fabs(dir.angle(drift_dir_abs) - 3.1415926 / 2.) < 5 * 3.1415926 / 180. ||
                    fabs(dir1.angle(drift_dir_abs) - 3.1415926 / 2.) < 5 * 3.1415926 / 180.) {
                    dir = cluster->vhough_transform(start_point, 200 * units::cm);
                }
                else {
                    dir = dir1;
                }
            }
        }
        dir = dir.norm();

        geo_point_t inv_dir = dir * (-1);
        start_wcpoint = cluster->get_furthest_wcpoint(start_wcpoint, inv_dir, 1 * units::cm, 0);
        end_wcpoint = cluster->get_furthest_wcpoint(start_wcpoint, dir);

   
        geo_point_t test_dir(end_wcpoint.x() - start_wcpoint.x(), end_wcpoint.y() - start_wcpoint.y(),
                             end_wcpoint.z() - start_wcpoint.z());
        start_wcpoint_idx = cluster->get_closest_point_index(start_wcpoint);
        end_wcpoint_idx = cluster->get_closest_point_index(end_wcpoint);
        if (fabs(test_dir.angle(drift_dir_abs) - 3.1415926 / 2.) < 2.5 * 3.1415926 / 180.) {
            cluster->adjust_wcpoints_parallel(start_wcpoint_idx, end_wcpoint_idx);
            start_wcpoint = cluster->point3d(start_wcpoint_idx);
            end_wcpoint = cluster->point3d(end_wcpoint_idx);
  
        }
    }
    if (sqrt(pow(start_wcpoint.x() - end_wcpoint.x(), 2) + pow(start_wcpoint.y() - end_wcpoint.y(), 2) +
             pow(start_wcpoint.z() - end_wcpoint.z(), 2)) < length / 3.) {
        // reverse the case ...
        start_wcpoint = independent_points.at(max_index);
        geo_point_t start_point(start_wcpoint.x(), start_wcpoint.y(), start_wcpoint.z());
        {
            // geo_point_t drift_dir_abs(1, 0, 0);
            dir = cluster->vhough_transform(start_point, 100 * units::cm);
            geo_point_t dir1 = cluster->vhough_transform(start_point, 30 * units::cm);
            if (dir.angle(dir1) > 20 * 3.1415926 / 180.) {
                if (fabs(dir.angle(drift_dir_abs) - 3.1415926 / 2.) < 5 * 3.1415926 / 180. ||
                    fabs(dir1.angle(drift_dir_abs) - 3.1415926 / 2.) < 5 * 3.1415926 / 180.) {
                    dir = cluster->vhough_transform(start_point, 200 * units::cm);
                }
                else {
                    dir = dir1;
                }
            }
        }
        dir = dir.norm();
        geo_point_t inv_dir = dir * (-1);
        start_wcpoint = cluster->get_furthest_wcpoint(start_wcpoint, inv_dir, 1 * units::cm, 0);
        end_wcpoint = cluster->get_furthest_wcpoint(start_wcpoint, dir);

        

        if (sqrt(pow(start_wcpoint.x() - end_wcpoint.x(), 2) + pow(start_wcpoint.y() - end_wcpoint.y(), 2) +
                 pow(start_wcpoint.z() - end_wcpoint.z(), 2)) < length / 3.) {  // reverse again ...
            start_wcpoint = end_wcpoint;
            geo_point_t start_point(start_wcpoint.x(), start_wcpoint.y(), start_wcpoint.z());
            {
                dir = cluster->vhough_transform(start_point, 100 * units::cm);
                geo_point_t dir1 = cluster->vhough_transform(start_point, 30 * units::cm);
                if (dir.angle(dir1) > 20 * 3.1415926 / 180.) {
                    if (fabs(dir.angle(drift_dir_abs) - 3.1415926 / 2.) < 5 * 3.1415926 / 180. ||
                        fabs(dir1.angle(drift_dir_abs) - 3.1415926 / 2.) < 5 * 3.1415926 / 180.) {
                        dir = cluster->vhough_transform(start_point, 200 * units::cm);
                    }
                    else {
                        dir = dir1;
                    }
                }
            }
            dir = dir.norm();
            geo_point_t inv_dir = dir * (-1);
            start_wcpoint = cluster->get_furthest_wcpoint(start_wcpoint, inv_dir, 1 * units::cm, 0);
            end_wcpoint = cluster->get_furthest_wcpoint(start_wcpoint, dir);
        }

    
      
        geo_point_t test_dir(end_wcpoint.x() - start_wcpoint.x(), end_wcpoint.y() - start_wcpoint.y(),
                             end_wcpoint.z() - start_wcpoint.z());
        start_wcpoint_idx = cluster->get_closest_point_index(start_wcpoint);
        end_wcpoint_idx = cluster->get_closest_point_index(end_wcpoint);
        if (fabs(test_dir.angle(drift_dir_abs) - 3.1415926 / 2.) < 2.5 * 3.1415926 / 180.) {
            cluster->adjust_wcpoints_parallel(start_wcpoint_idx, end_wcpoint_idx);
            start_wcpoint = cluster->point3d(start_wcpoint_idx);
            end_wcpoint = cluster->point3d(end_wcpoint_idx);
      
        }
    }
  
    // cluster->dijkstra_shortest_paths(dv, pcts, start_wcpoint_idx, use_ctpc);
    // cluster->cal_shortest_path(end_wcpoint_idx);
    // const auto& path_wcps = cluster->get_path_wcps();
    const auto& path_wcps = cluster->shortest_paths_graph(dv, pcts, use_ctpc).path(start_wcpoint_idx, end_wcpoint_idx);


    std::vector<bool> flag_u_pts, flag_v_pts, flag_w_pts;
    std::vector<bool> flag1_u_pts, flag1_v_pts, flag1_w_pts;
    std::vector<bool> flag2_u_pts, flag2_v_pts, flag2_w_pts;
    flag_u_pts.resize(cluster->npoints(), false);
    flag_v_pts.resize(cluster->npoints(), false);
    flag_w_pts.resize(cluster->npoints(), false);

    flag1_u_pts.resize(cluster->npoints(), false);
    flag1_v_pts.resize(cluster->npoints(), false);
    flag1_w_pts.resize(cluster->npoints(), false);

    flag2_u_pts.resize(cluster->npoints(), false);
    flag2_v_pts.resize(cluster->npoints(), false);
    flag2_w_pts.resize(cluster->npoints(), false);

    std::vector<geo_point_t> pts;

    // double acc_dis = 0;

    size_t prev_wcp_idx = path_wcps.front();
    for (auto it = path_wcps.begin(); it != path_wcps.end(); it++) {
        geo_point_t current_wcp = cluster->point3d(*it);
        geo_point_t prev_wcp = cluster->point3d(prev_wcp_idx);
        double dis = sqrt(pow(current_wcp.x() - prev_wcp.x(), 2) + pow(current_wcp.y() - prev_wcp.y(), 2) +
                          pow(current_wcp.z() - prev_wcp.z(), 2));
        // acc_dis += dis;
        // if (cluster->nchildren()==3449) std::cout << current_wcp << " " << dis/units::cm << " " << acc_dis/units::cm << std::endl;

        if (dis <= 1.0 * units::cm) {
            geo_point_t current_pt(current_wcp.x(), current_wcp.y(), current_wcp.z());
            pts.push_back(current_pt);
        }
        else {
            int num_points = int(dis / (1.0 * units::cm)) + 1;
            // double dis_seg = dis / num_points;
            for (int k = 0; k != num_points; k++) {
                geo_point_t current_pt(prev_wcp.x() + (k + 1.) / num_points * (current_wcp.x() - prev_wcp.x()),
                                       prev_wcp.y() + (k + 1.) / num_points * (current_wcp.y() - prev_wcp.y()),
                                       prev_wcp.z() + (k + 1.) / num_points * (current_wcp.z() - prev_wcp.z()));
                pts.push_back(current_pt);
            }
        }
        prev_wcp_idx = (*it);
    }
    for (const auto &pt : pts) {
        auto test_wpid = cluster->wpid(pt);
        if (test_wpid.apa()!=-1){
            af_temp_cloud.at(test_wpid.apa()).at(test_wpid.face())->add(pt);
        }
    }
    // if (flag_debug_porting) {
    //     std::cout << "temp_cloud->get_num_points() " << temp_cloud->get_num_points() << std::endl;
    // }

    const auto& winds = cluster->wire_indices();

    for (size_t j = 0; j != flag_u_pts.size(); j++) {
        geo_point_t test_p = cluster->point3d(j);
        auto test_wpid = cluster->wire_plane_id(j) ;
        std::pair<int, double> temp_results = af_temp_cloud.at(test_wpid.apa()).at(test_wpid.face())->get_closest_2d_dis(test_p, 0);
        double dis = temp_results.second;
        // if (flag_debug_porting && cluster->blob_with_point(j)->slice_index_min() == 8060) {
        //     std::cout << "get_closest_2d_dis(test_p, 0) " << test_p << " " << dis / units::cm << " cm " << (dis <= 2.4 * units::cm) << std::endl;
        //     auto closest_point = temp_cloud->point(0, temp_results.first);
        //     std::cout << "closest_point " << closest_point[0] << " " << closest_point[1] << std::endl;
        // }
        if (dis <= 1.5 * units::cm) {
            // flag_u_pts.at(j) = true;
            flag_u_pts.at(j) = true;
        }
        if (dis <= 2.4 * units::cm) {
            flag1_u_pts.at(j) = true;
        }
        else {
            auto& dead_u_index = af_dead_u_index.at(test_wpid.apa()).at(test_wpid.face());
            if (dead_u_index.find(winds[0][j]) != dead_u_index.end()) {
                // dead channels are corresponding to raw points
                if (cluster->point3d_raw(j).x() >= dead_u_index[winds[0][j]].first &&
                    cluster->point3d_raw(j).x() <= dead_u_index[winds[0][j]].second) {
                    if (dis < 10 * units::cm) flag1_u_pts.at(j) = true;
                    flag2_u_pts.at(j) = true;
                }
            }
        }
        temp_results = af_temp_cloud.at(test_wpid.apa()).at(test_wpid.face())->get_closest_2d_dis(test_p, 1);
        dis = temp_results.second;
        if (dis <= 1.5 * units::cm) {
            flag_v_pts.at(j) = true;
        }
        if (dis <= 2.4 * units::cm) {
            flag1_v_pts.at(j) = true;
        }
        else {
            auto& dead_v_index = af_dead_v_index.at(test_wpid.apa()).at(test_wpid.face());
            // dead channels are corresponding to raw points
            if (dead_v_index.find(winds[1][j]) != dead_v_index.end()) {
                if (cluster->point3d_raw(j).x() >= dead_v_index[winds[1][j]].first &&
                    cluster->point3d_raw(j).x() <= dead_v_index[winds[1][j]].second) {
                    if (dis < 10.0 * units::cm) flag1_v_pts.at(j) = true;
                    flag2_v_pts.at(j) = true;
                }
            }
        }
        temp_results = af_temp_cloud.at(test_wpid.apa()).at(test_wpid.face())->get_closest_2d_dis(test_p, 2);
        dis = temp_results.second;
        if (dis <= 1.5 * units::cm) {
            flag_w_pts.at(j) = true;
        }
        if (dis <= 2.4 * units::cm) {
            flag1_w_pts.at(j) = true;
        }
        else {
            auto& dead_w_index = af_dead_w_index.at(test_wpid.apa()).at(test_wpid.face());
            // dead channels are corresponding to raw points
            if (dead_w_index.find(winds[2][j]) != dead_w_index.end()) {
                if (cluster->point3d_raw(j).x() >= dead_w_index[winds[2][j]].first &&
                    cluster->point3d_raw(j).x() <= dead_w_index[winds[2][j]].second) {
                    if (dis < 10 * units::cm) flag1_w_pts.at(j) = true;
                    flag2_w_pts.at(j) = true;
                }
            }
        }
    }

    // special treatment of first and last point
    {
        auto wpid_front = cluster->wpid(pts.front());
        auto idx_front =  cluster->get_closest_point_index(pts.front());


        auto wpid_back = cluster->wpid(pts.back());
        auto idx_back =  cluster->get_closest_point_index(pts.back());

        std::vector<size_t> indices = cluster->get_closest_2d_index(cluster->point3d_raw(idx_front), 2.1 * units::cm, wpid_front.apa(), wpid_front.face(), 0);
        for (size_t k = 0; k != indices.size(); k++) {
            flag_u_pts.at(indices.at(k)) = true;
        }
        indices = cluster->get_closest_2d_index(cluster->point3d_raw(idx_front), 2.1 * units::cm, wpid_front.apa(), wpid_front.face(), 1);
        for (size_t k = 0; k != indices.size(); k++) {
            flag_v_pts.at(indices.at(k)) = true;
        }
        indices = cluster->get_closest_2d_index(cluster->point3d_raw(idx_front), 2.1 * units::cm, wpid_front.apa(), wpid_front.face(), 2);
        for (size_t k = 0; k != indices.size(); k++) {
            flag_w_pts.at(indices.at(k)) = true;
        }
        indices = cluster->get_closest_2d_index(cluster->point3d_raw(idx_back), 2.1 * units::cm, wpid_back.apa(), wpid_back.face(), 0);
        for (size_t k = 0; k != indices.size(); k++) {
            flag_u_pts.at(indices.at(k)) = true;
        }
        indices = cluster->get_closest_2d_index(cluster->point3d_raw(idx_back), 2.1 * units::cm, wpid_back.apa(), wpid_back.face(), 1);
        for (size_t k = 0; k != indices.size(); k++) {
            flag_v_pts.at(indices.at(k)) = true;
        }
        indices = cluster->get_closest_2d_index(cluster->point3d_raw(idx_back), 2.1 * units::cm, wpid_back.apa(), wpid_back.face(), 2);
        for (size_t k = 0; k != indices.size(); k++) {
            flag_w_pts.at(indices.at(k)) = true;
        }
    }

    // std::vector<Blob*>
    const auto &mcells = cluster->children();
    std::map<const Blob *, int> mcell_np_map, mcell_np_map1;
    for (auto it = mcells.begin(); it != mcells.end(); it++) {
        mcell_np_map[*it] = 0;
        mcell_np_map1[*it] = 0;
    }
    for (size_t j = 0; j != flag_u_pts.size(); j++) {
        const Blob* mcell = cluster->blob_with_point(j);
        
        if (flag_u_pts.at(j) && flag_v_pts.at(j) && flag1_w_pts.at(j) ||
            flag_u_pts.at(j) && flag_w_pts.at(j) && flag1_v_pts.at(j) ||
            flag_w_pts.at(j) && flag_v_pts.at(j) && flag1_u_pts.at(j)) {
            mcell_np_map[mcell]++;
        }

        if (flag_u_pts.at(j) && flag_v_pts.at(j) && (flag2_w_pts.at(j) || flag1_w_pts.at(j)) ||
            flag_u_pts.at(j) && flag_w_pts.at(j) && (flag2_v_pts.at(j) || flag1_v_pts.at(j)) ||
            flag_w_pts.at(j) && flag_v_pts.at(j) && (flag2_u_pts.at(j) || flag1_u_pts.at(j))) {
            mcell_np_map1[mcell]++;
        }
    }


    std::vector<Cluster *> final_clusters;

    std::vector<int> b2groupid(cluster->nchildren(), 0);
    std::set<int> groupids;

    for (size_t idx=0; idx < mcells.size(); idx++) {  
        Blob *mcell = mcells.at(idx);

        const size_t total_wires = mcell->u_wire_index_max() - mcell->u_wire_index_min() +
                             mcell->v_wire_index_max() - mcell->v_wire_index_min() +
                             mcell->w_wire_index_max() - mcell->w_wire_index_min();

        if (mcell_np_map[mcell] > 0.5 * mcell->nbpoints() ||
            (mcell_np_map[mcell] > 0.25 * mcell->nbpoints() && total_wires < 25)) {
            // cluster1->AddCell(mcell, mcell->GetTimeSlice());
            b2groupid[idx] = 0;
            groupids.insert(0);
        }
        else if (mcell_np_map1[mcell] >= 0.95 * mcell->nbpoints()) {
            // delete mcell;  // ghost cell ...
            b2groupid[idx] = -1; // to be deleted
            groupids.insert(-1);
        }
        else {
            // cluster2->AddCell(mcell, mcell->GetTimeSlice());
            b2groupid[idx] = 1;
            groupids.insert(1);
        }
    }
    auto scope_transform = cluster->get_scope_transform(scope);
    auto clusters_step0 = grouping->separate(cluster, b2groupid, true);
    // Apply the scope filter settings to all new clusters
    for (auto& [id, new_cluster] : clusters_step0) {
        new_cluster->set_scope_filter(scope, true);
        new_cluster->set_default_scope(scope);
        new_cluster->set_scope_transform(scope, scope_transform);
    }
    assert(cluster == nullptr);


    std::vector<Cluster*> other_clusters;

    if (clusters_step0.find(1) != clusters_step0.end()) {
        // other_clusters = Separate_2(clusters_step0[1], 5 * units::cm);
        const auto b2id = Separate_2(clusters_step0[1], scope, 5 * units::cm);
        auto other_clusters1 = grouping->separate(clusters_step0[1],b2id, true); // the cluster is now nullptr
        // Apply the scope filter settings to all new clusters
        for (auto& [id, new_cluster] : other_clusters1) {
            new_cluster->set_default_scope(scope);
            new_cluster->set_scope_filter(scope, true);
            new_cluster->set_scope_transform(scope, scope_transform);
        }
        assert(clusters_step0[1] == nullptr);

        for (auto it = other_clusters1.begin(); it != other_clusters1.end(); it++) {
            other_clusters.push_back(it->second);
        }
    }


//    if (clusters_step0.find(0) != clusters_step0.end()) {
    if (clusters_step0[0]->nchildren() > 0) {
        // merge some clusters from other_clusters to clusters_step0[0]
        {
            // cluster1->Create_point_cloud();
            // ToyPointCloud *cluster1_cloud = cluster1->get_point_cloud();
            std::vector<Cluster *> temp_merge_clusters;
            // check against other clusters
            for (size_t i = 0; i != other_clusters.size(); i++) {
                // other_clusters.at(i)->Create_point_cloud();
                std::tuple<int, int, double> temp_dis = other_clusters.at(i)->get_closest_points(*clusters_step0[0]);
                if (std::get<2>(temp_dis) < 0.5 * units::cm) {
                    double length_1 = other_clusters.at(i)->get_length();
                    geo_point_t p1(end_wcpoint.x(), end_wcpoint.y(), end_wcpoint.z());
                    double close_dis = other_clusters.at(i)->get_closest_dis(p1);

                    if (close_dis < 10 * units::cm && length_1 < 50 * units::cm) {
                        geo_point_t temp_dir1 = clusters_step0[0]->vhough_transform(p1, 15 * units::cm);
                        geo_point_t temp_dir2 = other_clusters.at(i)->vhough_transform(p1, 15 * units::cm);
                        if (temp_dir1.angle(temp_dir2) / 3.1415926 * 180. > 145 && length_1 < 30 * units::cm &&
                                close_dis < 3 * units::cm ||
                            fabs(temp_dir1.angle(drift_dir_abs) - 3.1415926 / 2.) / 3.1415926 * 180. < 3 &&
                                fabs(temp_dir2.angle(drift_dir_abs) - 3.1415926 / 2.) / 3.1415926 * 180. < 3) {
                            temp_merge_clusters.push_back(other_clusters.at(i));
                        }
                    }
                }
            }

            auto scope = clusters_step0[0]->get_default_scope();
            auto scope_transform = clusters_step0[0]->get_scope_transform(scope);

            for (auto temp_cluster : temp_merge_clusters) {
                other_clusters.erase(find(other_clusters.begin(),other_clusters.end(),temp_cluster));
                clusters_step0[0]->take_children(*temp_cluster, true);
                grouping->destroy_child(temp_cluster);
                assert(temp_cluster == nullptr);
            }
            clusters_step0[0]->set_default_scope(scope);
            clusters_step0[0]->set_scope_filter(scope, true);
            clusters_step0[0]->set_scope_transform(scope, scope_transform);

            final_clusters.push_back(clusters_step0[0]);
        }

        // ToyPointCloud *cluster1_cloud = cluster1->get_point_cloud();
        std::vector<Cluster *> saved_clusters;
        std::vector<Cluster *> to_be_merged_clusters;
        for (size_t i = 0; i != other_clusters.size(); i++) {
            // How to write???
            bool flag_save = false;
            double length_1 = other_clusters.at(i)->get_length();
      
            std::tuple<int, int, double> temp_dis = other_clusters.at(i)->get_closest_points(*clusters_step0[0]);
            if (length_1 < 30 * units::cm && std::get<2>(temp_dis) < 5 * units::cm) {
                int temp_total_points = other_clusters.at(i)->npoints();
                int temp_close_points = 0;
                for (int j = 0; j != other_clusters.at(i)->npoints(); j++) {
                    geo_point_t test_point = other_clusters.at(i)->point3d(j);
                    double close_dis = clusters_step0[0]->get_closest_dis(test_point);
                    if (close_dis < 10 * units::cm) {
                        temp_close_points++;
                    }
                }
                if (temp_close_points > 0.7 * temp_total_points) {
                    saved_clusters.push_back(other_clusters.at(i));
                    flag_save = true;
                }
            }
            else if (std::get<2>(temp_dis) < 2.5 * units::cm && length_1 >= 30 * units::cm) {
                int temp_total_points = other_clusters.at(i)->npoints();
                int temp_close_points = 0;
                for (int j = 0; j != other_clusters.at(i)->npoints(); j++) {
                    geo_point_t test_point = other_clusters.at(i)->point3d(j);
                    if (clusters_step0[0]->get_closest_dis(test_point) < 10 * units::cm) {
                        temp_close_points++;
                    }
                }
                if (temp_close_points > 0.85 * temp_total_points) {
                    saved_clusters.push_back(other_clusters.at(i));
                    flag_save = true;
                }
            }

            if (!flag_save) to_be_merged_clusters.push_back(other_clusters.at(i));
        }

        // add a protection
        std::vector<Cluster *> temp_save_clusters;

        for (size_t i = 0; i != saved_clusters.size(); i++) {
            Cluster *cluster1 = saved_clusters.at(i);
            if (cluster1->get_length() < 5 * units::cm) continue;
            // ToyPointCloud *cloud1 = cluster1->get_point_cloud();
            geo_point_t dir1(cluster1->get_pca().axis.at(0).x(), cluster1->get_pca().axis.at(0).y(),
                             cluster1->get_pca().axis.at(0).z());
            for (size_t j = 0; j != to_be_merged_clusters.size(); j++) {
                Cluster *cluster2 = to_be_merged_clusters.at(j);
                if (cluster2->get_length() < 10 * units::cm) continue;
                // ToyPointCloud *cloud2 = cluster2->get_point_cloud();
                geo_point_t dir2(cluster2->get_pca().axis.at(0).x(), cluster2->get_pca().axis.at(0).y(),
                                 cluster2->get_pca().axis.at(0).z());
                std::tuple<int, int, double> temp_dis = cluster1->get_closest_points(*cluster2);
                if (std::get<2>(temp_dis) < 15 * units::cm &&
                    fabs(dir1.angle(dir2) - 3.1415926 / 2.) / 3.1415926 * 180 > 75) {
                    //	std::cout << std::get<2>(temp_dis)/units::cm << " " <<  << std::endl;
                    temp_save_clusters.push_back(cluster1);
                    break;
                }
            }
        }
        for (size_t i = 0; i != temp_save_clusters.size(); i++) {
            Cluster *cluster1 = temp_save_clusters.at(i);
            to_be_merged_clusters.push_back(cluster1);
            saved_clusters.erase(find(saved_clusters.begin(), saved_clusters.end(), cluster1));
        }

        Cluster& cluster2 = grouping->make_child();
        cluster2.set_default_scope(scope);
        cluster2.set_scope_filter(scope, true);
        if (to_be_merged_clusters.size() > 0) cluster2.set_scope_transform(scope, to_be_merged_clusters[0]->get_scope_transform(scope));
        for (size_t i = 0; i != to_be_merged_clusters.size(); i++) {
            cluster2.take_children(*to_be_merged_clusters[i], true);
            grouping->destroy_child(to_be_merged_clusters[i]);
            assert(to_be_merged_clusters[i] == nullptr);
        }
        
        
        to_be_merged_clusters.clear();

        final_clusters.push_back(&cluster2);
        for (size_t i = 0; i != saved_clusters.size(); i++) {
            final_clusters.push_back(saved_clusters.at(i));
        }
        //saved_clusters.clear();
    }
    else {
        for (size_t i = 0; i != other_clusters.size(); i++) {
            final_clusters.push_back(other_clusters.at(i));
        }
    }
  
    return final_clusters;
}

#endif //_INDEV_

/// blob -> cluster_id
std::vector<int> WireCell::Clus::Facade::Separate_2(Cluster *cluster, 
                                                          const Tree::Scope& scope,
                                                          const double dis_cut)
{
    if (cluster->nchildren() == 0) {
        return std::vector<int>();
    }

    // std::cout << "Test: cluster has " << cluster->nchildren() << " blobs" << std::endl;
    auto& time_cells_set_map = cluster->time_blob_map();
    // Safe access to nested maps
    

    // std::cout << "Separate_2 nchildren: " << cluster->nchildren() << std::endl;
    const auto& mcells = cluster->children();

    // create graph for points between connected mcells, need to separate apa, face, and then ...
    std::map<int, std::map<int, std::vector<int> > > af_time_slices; // apa,face --> time slices 
    for (auto it = cluster->time_blob_map().begin(); it != cluster->time_blob_map().end(); it++) {
        int apa = it->first;
        for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++) {
            int face = it1->first;
            std::vector<int> time_slices_vec;
            for (auto it2 = it1->second.begin(); it2 != it1->second.end(); it2++) {
                time_slices_vec.push_back(it2->first);
            }
            af_time_slices[apa][face] = time_slices_vec;
        }
    }

    std::vector<std::pair<const Blob *, const Blob *>> connected_mcells;

    for (auto it = af_time_slices.begin(); it != af_time_slices.end(); it++) {
        int apa = it->first;
        for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++) {
            int face = it1->first;
            std::vector<int>& time_slices = it1->second;
 
            for (size_t i = 0; i != time_slices.size(); i++) {
                const BlobSet &mcells_set = time_cells_set_map.at(apa).at(face).at(time_slices.at(i));
                // std::cout << "time_slices.at(i)" << time_slices.at(i) << " mcells_set.size() " << mcells_set.size() << std::endl;

                // create graph for points in mcell inside the same time slice
                if (mcells_set.size() >= 2) {
                    for (auto it2 = mcells_set.begin(); it2 != mcells_set.end(); it2++) {
                        const Blob *mcell1 = *it2;
                        auto it2p = it2;
                        if (it2p != mcells_set.end()) {
                            it2p++;
                            for (auto it3 = it2p; it3 != mcells_set.end(); it3++) {
                                const Blob *mcell2 = *(it3);
                                if (mcell1->overlap_fast(*mcell2, 5)) connected_mcells.push_back(std::make_pair(mcell1, mcell2));
                            }
                        }
                    }
                }
                // create graph for points between connected mcells in adjacent time slices + 1, if not, + 2
                std::vector<BlobSet> vec_mcells_set;
                if (i + 1 < time_slices.size()) {
                    if (time_slices.at(i + 1) - time_slices.at(i) == (int)(1*cluster->grouping()->get_nticks_per_slice().at(apa).at(face))) {
                        vec_mcells_set.push_back(time_cells_set_map.at(apa).at(face).at(time_slices.at(i + 1)));
                        if (i + 2 < time_slices.size())
                            if (time_slices.at(i + 2) - time_slices.at(i) == (int)(2*cluster->grouping()->get_nticks_per_slice().at(apa).at(face)))
                                vec_mcells_set.push_back(time_cells_set_map.at(apa).at(face).at(time_slices.at(i + 2)));
                    }
                    else if (time_slices.at(i + 1) - time_slices.at(i) == (int)(2*cluster->grouping()->get_nticks_per_slice().at(apa).at(face))) {
                        vec_mcells_set.push_back(time_cells_set_map.at(apa).at(face).at(time_slices.at(i + 1)));
                    }
                }
                // std::cout << "time_slices.at(i)" << time_slices.at(i) << " vec_mcells_set.size() " << vec_mcells_set.size() << std::endl;
                bool flag = false;
                for (size_t j = 0; j != vec_mcells_set.size(); j++) {
                    if (flag) break;
                    BlobSet &next_mcells_set = vec_mcells_set.at(j);
                    for (auto it1 = mcells_set.begin(); it1 != mcells_set.end(); it1++) {
                        const Blob *mcell1 = (*it1);
                        for (auto it2 = next_mcells_set.begin(); it2 != next_mcells_set.end(); it2++) {
                            const Blob *mcell2 = (*it2);
                            if (mcell1->overlap_fast(*mcell2, 2)) {
                                flag = true;
                                connected_mcells.push_back(std::make_pair(mcell1, mcell2));
                            }
                        }
                    }
                }
            }
        }
    }

    // form ...

    const int N = mcells.size();
    Weighted::Graph graph(N);

    std::map<const Blob *, int> mcell_index_map;
    for (size_t i = 0; i != mcells.size(); i++) {
        Blob *curr_mcell = mcells.at(i);
        mcell_index_map[curr_mcell] = i;

        // auto v = vertex(i, graph);  // retrieve vertex descriptor
        // (graph)[v].ident = i;
    }

    for (auto it = connected_mcells.begin(); it != connected_mcells.end(); it++) {
        int index1 = mcell_index_map[it->first];
        int index2 = mcell_index_map[it->second];
        // auto edge = add_edge(index1, index2, graph);
        // if (edge.second) {
        //     (graph)[edge.first].dist = 1;
        // }
        /*auto edge =*/ add_edge(index1, index2, 1.0,graph);
    }

    {
        // std::string hack_pc_name = "3d";
        // std::vector<std::string> hack_coords = {"x", "y", "z"};
        // std::cout << "Separate_2: num_edges: " << num_edges(graph) << std::endl;
        std::vector<int> component(num_vertices(graph));
        const int num = connected_components(graph, &component[0]);

        if (num > 1) {
            std::vector<std::shared_ptr<Simple3DPointCloud>> pt_clouds;
            std::vector<std::vector<int>> vec_vec(num);
            for (int j = 0; j != num; j++) {
                pt_clouds.push_back(std::make_shared<Simple3DPointCloud>());
            }
            std::vector<int>::size_type i;
            for (i = 0; i != component.size(); ++i) {
                vec_vec.at(component[i]).push_back(i);
                Blob *mcell = mcells.at(i);
                for (const auto & pt : mcell->points(scope.pcname, scope.coords)) {
                    const std::vector<double> newpt = {pt.x(), pt.y(), pt.z()};
                    pt_clouds.at(component[i])->add(newpt);
                }
            }

            for (int j = 0; j != num; j++) {
                for (int k = j + 1; k != num; k++) {
                    std::tuple<int, int, double> temp_results = pt_clouds.at(j)->get_closest_points(*(pt_clouds.at(k)));
                    if (std::get<2>(temp_results) < dis_cut) {
                        int index1 = vec_vec[j].front();
                        int index2 = vec_vec[k].front();
                        // auto edge = add_edge(index1, index2, graph);
                        // if (edge.second) {
                        //     (graph)[edge.first].dist = 1;
                        // }
                        /*auto edge =*/ add_edge(index1, index2, 1.0,graph);
                    }
                }
            }
        }

        // std::cout << num << std::endl;
    }

  
    std::vector<int> component(num_vertices(graph));
    /*const int num =*/ connected_components(graph, &component[0]);
    return component;
}


