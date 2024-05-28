#include "WireCellImg/MultiAlgBlobClustering.h"
#include "WireCellImg/PointCloudFacade.h"
#include <WireCellImg/ClusteringFuncs.h>
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/ExecMon.h"
#include "WireCellAux/TensorDMpointtree.h"
#include "WireCellAux/TensorDMdataset.h"
#include "WireCellAux/TensorDMcommon.h"
#include "WireCellAux/SimpleTensorSet.h"

#include "WireCellUtil/Graph.h"


#include <fstream>

WIRECELL_FACTORY(MultiAlgBlobClustering, WireCell::Img::MultiAlgBlobClustering, WireCell::INamed,
                 WireCell::ITensorSetFilter, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Img;
using namespace WireCell::Aux;
using namespace WireCell::Aux::TensorDM;
using namespace WireCell::PointCloud::Facade;
using namespace WireCell::PointCloud::Tree;

MultiAlgBlobClustering::MultiAlgBlobClustering()
  : Aux::Logger("MultiAlgBlobClustering", "img")
{
}

void MultiAlgBlobClustering::configure(const WireCell::Configuration& cfg)
{
    m_inpath = get(cfg, "inpath", m_inpath);
    m_outpath = get(cfg, "outpath", m_outpath);
    m_bee_dir = get(cfg, "bee_dir", m_bee_dir);

    m_dead_live_overlap_offset = get(cfg, "dead_live_overlap_offset", m_dead_live_overlap_offset);
}

WireCell::Configuration MultiAlgBlobClustering::default_configuration() const
{
    Configuration cfg;
    cfg["inpath"] = m_inpath;
    cfg["outpath"] = m_outpath;
    cfg["bee_dir"] = m_bee_dir;

    cfg["dead_live_overlap_offset"] = m_dead_live_overlap_offset;
    return cfg;
}

namespace {
    void dump_bee(const Points::node_t& root, const std::string& fn)
    {
        using spdlog::debug;
        using WireCell::PointCloud::Facade::float_t;
        using WireCell::PointCloud::Facade::int_t;

        Configuration bee;
        bee["runNo"] = 0;
        bee["subRunNo"] = 0;
        bee["eventNo"] = 0;
        bee["geom"] = "uboone";
        bee["type"] = "cluster";

        std::vector<float_t> x;
        std::vector<float_t> y;
        std::vector<float_t> z;
        std::vector<float_t> q;
        std::vector<int_t> cluster_id;
        int_t cid = 0;
        for (const auto cnode : root.children()) {  // this is a loop through all clusters ...
            Scope scope = {"3d", {"x", "y", "z"}};
            const auto& sv = cnode->value.scoped_view(scope);

            const auto& spcs = sv.pcs();  // spcs 'contains' all blobs in this cluster ...
            // int npoints = 0;

            for (const auto& spc : spcs) {  // each little 3D pc --> (blobs)   spc represents x,y,z in a blob
                if (spc.get().get("x") == nullptr) {
                    debug("No x in point cloud, skip");
                    continue;
                }

                // assume others exist
                const auto& x_ = spc.get().get("x")->elements<float_t>();
                const auto& y_ = spc.get().get("y")->elements<float_t>();
                const auto& z_ = spc.get().get("z")->elements<float_t>();
                const size_t n = x_.size();
                x.insert(x.end(), x_.begin(), x_.end());  // Append x_ to x
                y.insert(y.end(), y_.begin(), y_.end());
                z.insert(z.end(), z_.begin(), z_.end());
                q.insert(q.end(), n, 1.0);
                cluster_id.insert(cluster_id.end(), n, cid);
                // npoints += n;
            }

            // spc.kd() // kdtree ...
            // const auto& skd = sv.kd();
            // std::cout << "xin6: " << sv.npoints() << " " << npoints << " " << spcs.size() << " " <<
            // skd.points().size() << std::endl;

            ++cid;
        }

        Json::Value json_x(Json::arrayValue);
        for (const auto& val : x) {
            json_x.append(val / units::cm);
        }
        bee["x"] = json_x;

        Json::Value json_y(Json::arrayValue);
        for (const auto& val : y) {
            json_y.append(val / units::cm);
        }
        bee["y"] = json_y;

        Json::Value json_z(Json::arrayValue);
        for (const auto& val : z) {
            json_z.append(val / units::cm);
        }
        bee["z"] = json_z;

        Json::Value json_q(Json::arrayValue);
        for (const auto& val : q) {
            json_q.append(val);
        }
        bee["q"] = json_q;

        Json::Value json_cluster_id(Json::arrayValue);
        for (const auto& val : cluster_id) {
            json_cluster_id.append(val);
        }
        bee["cluster_id"] = json_cluster_id;

        // Write cfg to file
        std::ofstream file(fn);
        if (file.is_open()) {
            Json::StreamWriterBuilder writer;
            writer["indentation"] = "    ";
            std::unique_ptr<Json::StreamWriter> jsonWriter(writer.newStreamWriter());
            jsonWriter->write(bee, &file);
            file.close();
        }
        else {
            raise<ValueError>("Failed to open file: " + fn);
        }
    }

    struct Point2D {
        Point2D(float x, float y)
          : x(x)
          , y(y)
        {
        }
        float x;
        float y;
    };

    bool anglular_less(const Point2D& a, const Point2D& b, const Point2D& center)
    {
        if (a.x - center.x >= 0 && b.x - center.x < 0) return true;
        if (a.x - center.x < 0 && b.x - center.x >= 0) return false;
        if (a.x - center.x == 0 && b.x - center.x == 0) {
            if (a.y - center.y >= 0 || b.y - center.y >= 0) return a.y > b.y;
            return b.y > a.y;
        }

        // compute the cross product of vectors (center -> a) x (center -> b)
        int det = (a.x - center.x) * (b.y - center.y) - (b.x - center.x) * (a.y - center.y);
        if (det < 0) return true;
        if (det > 0) return false;

        // points a and b are on the same line from the center
        // check which point is closer to the center
        int d1 = (a.x - center.x) * (a.x - center.x) + (a.y - center.y) * (a.y - center.y);
        int d2 = (b.x - center.x) * (b.x - center.x) + (b.y - center.y) * (b.y - center.y);
        return d1 > d2;
    }

    std::vector<Point2D> sort_angular(const std::vector<Point2D>& points)
    {
        if (points.size() < 3) {
            return points;
        }
        // Calculate the center of the points
        float sumX = 0.0;
        float sumY = 0.0;
        for (const auto& point : points) {
            sumX += point.x;
            sumY += point.y;
        }
        Point2D center(sumX / points.size(), sumY / points.size());

        std::vector<Point2D> sorted = points;
        std::sort(sorted.begin(), sorted.end(),
                  [&](const Point2D& a, const Point2D& b) { return anglular_less(a, b, center); });
        return sorted;
    }

    void dumpe_deadarea(const Points::node_t& root, const std::string& fn)
    {
        using WireCell::PointCloud::Facade::float_t;
        using WireCell::PointCloud::Facade::int_t;

        // Convert stringstream to JSON array
        Json::Value jdead;
        for (const auto cnode : root.children()) {
            for (const auto bnode : cnode->children()) {
                const auto& lpcs = bnode->value.local_pcs();
                const auto& pc_scalar = lpcs.at("corner");
                const auto& y = pc_scalar.get("y")->elements<float_t>();
                const auto& z = pc_scalar.get("z")->elements<float_t>();
                std::vector<Point2D> points;
                for (size_t i = 0; i < y.size(); ++i) {
                    points.push_back({(float)y[i], (float)z[i]});
                }
                auto sorted = sort_angular(points);
                Json::Value jarea(Json::arrayValue);
                for (const auto& point : sorted) {
                    Json::Value jpoint(Json::arrayValue);
                    jpoint.append(point.x / units::cm);
                    jpoint.append(point.y / units::cm);
                    jarea.append(jpoint);
                }
                jdead.append(jarea);
            }
        }

        // Output jsonArray to file
        std::ofstream file(fn);
        if (file.is_open()) {
            file << jdead.toStyledString();
            file.close();
        }
        else {
            std::cout << "Failed to open file: " << fn << std::endl;
        }
    }
}  // namespace

bool MultiAlgBlobClustering::operator()(const input_pointer& ints, output_pointer& outts)
{
    outts = nullptr;
    if (!ints) {
        log->debug("EOS at call {}", m_count++);
        return true;
    }

    bool flag_print = false;
    ExecMon em("starting MultiAlgBlobClustering");

    const int ident = ints->ident();
    std::string inpath = m_inpath;
    if (inpath.find("%") != std::string::npos) {
        inpath = String::format(inpath, ident);
    }

    const auto& intens = *ints->tensors();
    log->debug("Input {} tensors", intens.size());
    auto root_live = std::move(as_pctree(intens, inpath + "/live"));
    if (!root_live) {
        log->error("Failed to get point cloud tree from \"{}\"", inpath);
        return false;
    }
    log->debug("Got live pctree with {} children", root_live->nchildren());
    log->debug(em("got live pctree"));

    auto root_dead = as_pctree(intens, inpath + "/dead");
    if (!root_dead) {
        log->error("Failed to get point cloud tree from \"{}\"", inpath + "/dead");
        return false;
    }
    log->debug("Got dead pctree with {} children", root_dead->nchildren());
    log->debug(em("got dead pctree"));

    // BEE debug direct imaging output and dead blobs
    if (!m_bee_dir.empty()) {
        std::string sub_dir = String::format("%s/%d", m_bee_dir, ident);
        Persist::assuredir(sub_dir);
        dump_bee(*root_live.get(), String::format("%s/%d-img.json", sub_dir, ident));
        dumpe_deadarea(*root_dead.get(), String::format("%s/%d-channel-deadarea.json", sub_dir, ident));
    }
    log->debug(em("dump live+dead to bee"));

    if (flag_print) std::cout << em("Finish PC/Facade Conversion ") << std::endl;
    
    /// TODO: how to pass the parameters? for now, using default params
    WireCell::PointCloud::Facade::TPCParams tp;

    // Calculate the length of all the clusters and save them into a map
    std::map<const Cluster::const_pointer, double> cluster_length_map;
    std::set<Cluster::const_pointer > cluster_connected_dead;
    
    // initialize clusters ...
    live_clusters_t live_clusters;
    for (auto cnode : root_live->children()) {
        live_clusters.push_back(std::make_shared<Cluster>(cnode));
    }
    log->debug(em("make live clusters"));

    {  // ATTENTION, this block is just for debugging.
        // This is here just to trigger k-d tree building so that we can separate
        // that out its time/memory from those of the first algorithm
        size_t npts=0;
        for (const auto& cl : live_clusters) {
            cl->get_closest_point_mcell(geo_point_t(0,0,0));
            npts += cl->get_num_points();
        }
        log->debug("made {} live clusters with {} points", live_clusters.size(), npts);
        log->debug(em("make live k-d tree"));
    }

    // loop over all the clusters, and calculate length ...
    for (size_t ilive = 0; ilive < live_clusters.size(); ++ilive) {
        const auto& live = live_clusters[ilive];
        cluster_length_map[live] = live->get_length(tp);
        // std::cout << ilive << " xin " << live->get_length(tp)/units::cm << std::endl;
    }
    log->debug(em("get live lengths"));
    
    Cluster::const_vector dead_clusters;
    for (auto cnode : root_dead->children()) {
        dead_clusters.push_back(std::make_shared<Cluster>(cnode));
    }
    log->debug(em("make dead clusters"));

    {                           // trigger k-d tree building early for perf testing
        size_t npts=0;
        for (const auto& cl : dead_clusters) {
            cl->get_closest_point_mcell(geo_point_t(0,0,0));
            npts += cl->get_num_points();
        }
        log->debug("made {} dead clusters with {} points", dead_clusters.size(), npts);
        log->debug(em("make dead k-d tree"));
    }

    // dead_live
    clustering_live_dead(root_live, live_clusters, dead_clusters, cluster_length_map, cluster_connected_dead, tp,
                         m_dead_live_overlap_offset);
    log->debug(em("clustering_live_dead"));

    if (flag_print) std::cout << em("live_dead") << std::endl;
    // second function ...
    clustering_extend(root_live, live_clusters, cluster_length_map, cluster_connected_dead, tp, 4,60*units::cm,0,15*units::cm,1 );
    log->debug(em("clustering_extend"));
    if (flag_print) std::cout << em("first extend") << std::endl;
    
    // first round clustering
    clustering_regular(root_live, live_clusters, cluster_length_map,cluster_connected_dead,tp, 60*units::cm, false);
    log->debug(em("clustering_regular 1st"));
    if (flag_print) std::cout << em("1st regular") << std::endl;
    clustering_regular(root_live, live_clusters, cluster_length_map,cluster_connected_dead,tp, 30*units::cm, true); // do extension
    log->debug(em("clustering_regular 2nd"));
    if (flag_print) std::cout << em("2nd regular") << std::endl;

    
    //dedicated one dealing with parallel and prolonged track
    clustering_parallel_prolong(root_live, live_clusters, cluster_length_map,cluster_connected_dead,tp,35*units::cm);
    log->debug(em("clustering_parallel_prolong"));
    if (flag_print) std::cout << em("parallel prolong") << std::endl;
    
    //clustering close distance ones ... 
    clustering_close(root_live, live_clusters, cluster_length_map,cluster_connected_dead,tp, 1.2*units::cm);
    log->debug(em("clustering_close"));
    if (flag_print) std::cout << em("close") << std::endl;
    

    int num_try =3;
    // for very busy events do less ... 
    if (live_clusters.size() > 1100 ) num_try = 1;
    for (int i=0;i!= num_try ;i++){
      //extend the track ...
      // deal with prolong case
      clustering_extend(root_live,live_clusters, cluster_length_map,cluster_connected_dead,tp,1,150*units::cm,0);
      log->debug(em("clustering_extend prolong"));
      if (flag_print) std::cout << em("extend prolong") << std::endl;
      // deal with parallel case 
      clustering_extend(root_live,live_clusters, cluster_length_map,cluster_connected_dead,tp,2,30*units::cm,0);
      log->debug(em("clustering_extend parallel"));
      if (flag_print) std::cout << em("extend parallel") << std::endl;
      
      
      // extension regular case
      clustering_extend(root_live,live_clusters, cluster_length_map,cluster_connected_dead,tp,3,15*units::cm,0);
      log->debug(em("clustering_extend regular"));
      
      if (flag_print) std::cout << i << std::endl;
      
      if (flag_print) std::cout << em("extend regular") << std::endl;
      // extension ones connected to dead region ...
      if (i==0){
	clustering_extend(root_live,live_clusters, cluster_length_map,cluster_connected_dead,tp,4,60*units::cm,i);
      }else{
	clustering_extend(root_live,live_clusters, cluster_length_map,cluster_connected_dead,tp,4,35*units::cm,i);
      }
      log->debug(em("clustering_extend dead"));
      if (flag_print) std::cout << em("extend dead") << std::endl;
    }
        
    // BEE debug dead-live
    if (!m_bee_dir.empty()) {
        std::string sub_dir = String::format("%s/%d", m_bee_dir, ident);
        dump_bee(*root_live.get(), String::format("%s/%d-dead-live.json", sub_dir, ident));
    }
    log->debug(em("dump bee dead-live"));

    std::string outpath = m_outpath;
    if (outpath.find("%") != std::string::npos) {
        outpath = String::format(outpath, ident);
    }
    auto outtens = as_tensors(*root_live.get(), outpath + "/live");
    log->debug(em("as tensors live"));

    auto outtens_dead = as_tensors(*root_dead.get(), outpath + "/dead");
    log->debug(em("as tensors dead"));

    if (flag_print) std::cout << em("dump bee") << std::endl;
    
    // Merge
    outtens.insert(outtens.end(), outtens_dead.begin(), outtens_dead.end());
    log->debug("Total outtens {} tensors", outtens.size());
    outts = as_tensorset(outtens, ident);
    log->debug(em("as tensors set output"));

    live_clusters.clear();
    dead_clusters.clear();
    root_live = nullptr;
    root_dead = nullptr;

    log->debug(em("clear memory"));

    return true;
}
