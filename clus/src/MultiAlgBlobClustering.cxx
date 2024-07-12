#include "WireCellClus/MultiAlgBlobClustering.h"
#include "WireCellClus/Facade.h"
#include <WireCellClus/ClusteringFuncs.h>
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
  , m_bee_img("uboone", "img")
  , m_bee_ld("uboone", "live-dead")
  , m_bee_dead("channel-deadarea", 0.0, 3)
{
}

void MultiAlgBlobClustering::configure(const WireCell::Configuration& cfg)
{
    m_inpath = get(cfg, "inpath", m_inpath);
    m_outpath = get(cfg, "outpath", m_outpath);

    if (cfg.isMember("bee_dir")) {
        log->warn("the 'bee_dir' option is no longer supported, instead use 'bee_zip' to name a .zip file");
    }

    m_sink.reset(get<std::string>(cfg, "bee_zip", "mabc.zip"));
    m_save_deadarea = get(cfg, "save_deadarea", m_save_deadarea);

    m_dead_live_overlap_offset = get(cfg, "dead_live_overlap_offset", m_dead_live_overlap_offset);
    m_perf = get(cfg, "perf", m_perf);

    m_anode = Factory::find_tn<IAnodePlane>(cfg["anode"].asString());
}

WireCell::Configuration MultiAlgBlobClustering::default_configuration() const
{
    Configuration cfg;
    cfg["inpath"] = m_inpath;
    cfg["outpath"] = m_outpath;
    // cfg["bee_dir"] = m_bee_dir;
    cfg["bee_zip"] = "mabc.zip";
    cfg["save_deadarea"] = m_save_deadarea;

    cfg["dead_live_overlap_offset"] = m_dead_live_overlap_offset;
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
    flush(m_bee_img, ident);
    flush(m_bee_ld,  ident);
    if (m_save_deadarea && m_bee_dead.size()) {
        m_bee_dead.flush();
        m_sink.write(m_bee_dead);
        m_bee_dead.clear();
    }
    m_last_ident = ident;
}

// There are equivalent functions in Aux::Bee:: but the pc tree is subject to
// many schema and there is no "standard".  So we keep this dumper here, since
// it is here we know the pc tree schema.
static
void fill_bee_points(WireCell::Bee::Points& bpts, const Points::node_t& root)
{
    int clid = bpts.back_cluster_id();
    const double charge = 0;
    for (const auto cnode : root.children()) {  // this is a loop through all clusters ...
        ++clid;

        Scope scope = {"3d", {"x", "y", "z"}};
        const auto& sv = cnode->value.scoped_view(scope);

        const auto& spcs = sv.pcs();  // spcs 'contains' all blobs in this cluster ...

        for (const auto& spc : spcs) {  // each little 3D pc --> (blobs)   spc represents x,y,z in a blob
            auto x = spc.get().get("x")->elements<double>();
            auto y = spc.get().get("y")->elements<double>();
            auto z = spc.get().get("z")->elements<double>();
            const size_t size = x.size();
            // fixme: add to Bee::Points a method to append vector-like things...
            for (size_t ind = 0 ; ind<size; ++ind) {
                bpts.append(Point(x[ind], y[ind], z[ind]), charge, clid);
            }
        }
    }
}

static
void fill_bee_patches(WireCell::Bee::Patches& bee, const Points::node_t& root)
{
    int first_slice = -1;
    for (const auto cnode : root.children()) {
        for (const auto bnode : cnode->children()) {
            const auto& lpcs = bnode->value.local_pcs();

            const auto& pc_scalar = lpcs.at("scalar");
            int slice_index_min = pc_scalar.get("slice_index_min")->elements<int>()[0];
            if (first_slice < 0) {
                first_slice = slice_index_min;
            }
            if (slice_index_min != first_slice) continue;

            const auto& pc_corner = lpcs.at("corner");
            const auto& y = pc_corner.get("y")->elements<double>();
            const auto& z = pc_corner.get("z")->elements<double>();
            bee.append(y.begin(), y.end(), z.begin(), z.end());
        }
    }
}

namespace {

    struct Point2D {
        Point2D(float x, float y)
          : x(x)
          , y(y)
        {
        }
        float x;
        float y;
    };

    bool angular_less(const Point2D& a, const Point2D& b, const Point2D& center)
    {
        double angleA = std::atan2(a.y - center.y, a.x - center.x);
        double angleB = std::atan2(b.y - center.y, b.x - center.x);
        return angleA < angleB;
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
                  [&](const Point2D& a, const Point2D& b) { return angular_less(a, b, center); });
        return sorted;
    }
    std::vector<Point2D> unique(const std::vector<Point2D>& points, const float tolerance = 0.1)
    {
        auto less_with_tolerance = [&](const Point2D& a, const Point2D& b) {
            if (std::abs(a.x - b.x) > tolerance) return a.x < b.x;
            if (std::abs(a.y - b.y) > tolerance) return a.y < b.y;
            return false;
        };
        std::set<Point2D, decltype(less_with_tolerance)> unique_points(less_with_tolerance);
        for (const auto& point : points) {
            unique_points.insert(point);
        }
        std::vector<Point2D> unique_points_vec;
        for (const auto& point : unique_points) {
            unique_points_vec.push_back(point);
        }
        return unique_points_vec;
    }
#if 0
    bool valid(const std::vector<Point2D>& points, const float threshold = 1.)
    {
        if (points.size() < 3) return false;
        float min_x = points[0].x;
        float max_x = points[0].x;
        float min_y = points[0].y;
        float max_y = points[0].y;
        for (const auto& point : points) {
            if (point.x < min_x) min_x = point.x;
            if (point.x > max_x) max_x = point.x;
            if (point.y < min_y) min_y = point.y;
            if (point.y > max_y) max_y = point.y;
        }
        if (max_x - min_x < threshold && max_y - min_y < threshold) return false;
        return true;
    }
#endif

    void dumpe_deadarea(const Points::node_t& root, const std::string& fn)
    {
        using WireCell::PointCloud::Facade::float_t;
        using WireCell::PointCloud::Facade::int_t;

        // Convert stringstream to JSON array
        Json::Value jdead;
        for (const auto cnode : root.children()) {
            for (const auto bnode : cnode->children()) {
                const auto& lpcs = bnode->value.local_pcs();
                const auto& pc_scalar = lpcs.at("scalar");
                const auto& slice_index_min = pc_scalar.get("slice_index_min")->elements<int_t>()[0];
                if (slice_index_min != 0) continue;
                const auto& pc_corner = lpcs.at("corner");
                const auto& y = pc_corner.get("y")->elements<float_t>();
                const auto& z = pc_corner.get("z")->elements<float_t>();
                std::vector<Point2D> points;
                for (size_t i = 0; i < y.size(); ++i) {
                    points.push_back({(float) y[i] / (float) units::cm, (float) z[i] / (float) units::cm});
                }
                // Remove duplicate points with xxx cm tolerance
                auto unique_points = unique(points, 0.1);
                // if (!valid(unique_points, 1.0)) continue;
                if (unique_points.size() < 3) continue;
                auto sorted = sort_angular(unique_points);
                Json::Value jarea(Json::arrayValue);
                for (const auto& point : sorted) {
                    Json::Value jpoint(Json::arrayValue);
                    jpoint.append(point.x);
                    jpoint.append(point.y);
                    jarea.append(jpoint);
                }
                jdead.append(jarea);
            }
        }

        // Output jsonArray to file
        std::ofstream file(fn);
        if (file.is_open()) {
            Json::StreamWriterBuilder builder;
            builder.settings_["indentation"] = "    ";
            builder.settings_["precision"] = 6;  // significant digits
            // option available in jsoncpp 1.9.0
            // https://github.com/open-source-parsers/jsoncpp/blob/1.9.0/src/lib_json/json_writer.cpp#L128
            // builder.settings_["precisionType"] = "decimal";
            std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
            writer->write(jdead, &file);
            file.close();
        }
        else {
            std::cout << "Failed to open file: " << fn << std::endl;
        }
    }
}  // namespace

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
    if (m_last_ident < 0) {     // first time.
        reset_bee(ident, m_bee_img);
        reset_bee(ident, m_bee_ld);
        m_last_ident = ident;
    }
    else if (m_last_ident != ident) {
        flush(ident);
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
    grouping->set_anode(m_anode);
    perf("loaded live clusters");

    // log->debug("Got live pctree with {} children", root_live->nchildren());
    // log->debug(em("got live pctree"));
    auto root_dead = as_pctree(intens, inpath + "/dead");
    if (!root_dead) {
        log->error("Failed to get dead point cloud tree from \"{}\"", inpath + "/dead");
        raise<ValueError>("Failed to get dead point cloud tree from \"%s\"", inpath);
    }
    perf("loaded dead clusters");

    // FIXME: we do not yet have a replacement for dumpe_deadarea() in WireCell::Bee
    // BEE debug direct imaging output and dead blobs
    // if (!m_bee_dir.empty()) {
    //     std::string sub_dir = String::format("%s/%d", m_bee_dir, ident);
    //     Persist::assuredir(sub_dir);
    //     dump_bee(*root_live.get(), String::format("%s/%d-img.json", sub_dir, ident));
    //     if (m_save_deadarea) {
    //         dumpe_deadarea(*root_dead.get(), String::format("%s/%d-channel-deadarea.json", sub_dir, ident));
    //     }
    //     perf("loaded dump live clusters to bee");
    // }
    fill_bee_points(m_bee_img, *root_live.get());
    perf("loaded dump live clusters to bee");
    if (m_save_deadarea) {
        fill_bee_patches(m_bee_dead, *root_dead.get());
        perf("loaded dump dead regions to bee");
    }
    log->debug("will {} {} dead patches", m_save_deadarea ? "save" : "not save", m_bee_dead.size());

    cluster_set_t cluster_connected_dead;

    // initialize clusters ...
    Grouping& live_grouping = *root_live->value.facade<Grouping>();
    Grouping& dead_grouping = *root_dead->value.facade<Grouping>();

    //perf.dump("original live clusters", live_grouping, false, false);
    //perf.dump("original dead clusters", dead_grouping, false, false);

    perf.dump("pre clustering", live_grouping);

    // dead_live
    clustering_live_dead(live_grouping, dead_grouping, cluster_connected_dead, m_dead_live_overlap_offset);
    perf.dump("clustering live-dead", live_grouping);

    // second function ...
    clustering_extend(live_grouping, cluster_connected_dead, 4, 60 * units::cm, 0, 15 * units::cm, 1);
    perf.dump("clustering extend", live_grouping);

    // first round clustering
    clustering_regular(live_grouping, cluster_connected_dead, 60 * units::cm, false);
    perf.dump("clustering regular no extension", live_grouping);

    clustering_regular(live_grouping, cluster_connected_dead, 30 * units::cm, true);  // do extension
    perf.dump("clustering regular with extension", live_grouping);

    // dedicated one dealing with parallel and prolonged track
    clustering_parallel_prolong(live_grouping, cluster_connected_dead, 35 * units::cm);
    perf.dump("clustering parallel prolong", live_grouping);

    // clustering close distance ones ...
    clustering_close(live_grouping, cluster_connected_dead, 1.2 * units::cm);
    perf.dump("clustering close", live_grouping);

    int num_try = 3;
    // for very busy events do less ...
    if (live_grouping.nchildren() > 1100) num_try = 1;
    for (int i = 0; i != num_try; i++) {
        // extend the track ...

        // deal with prolong case
        clustering_extend(live_grouping, cluster_connected_dead, 1, 150 * units::cm, 0);
        perf.dump("clustering extend 1", live_grouping);

        // deal with parallel case
        clustering_extend(live_grouping, cluster_connected_dead, 2, 30 * units::cm, 0);
        perf.dump("clustering extend 2", live_grouping);

        // extension regular case
        clustering_extend(live_grouping, cluster_connected_dead, 3, 15 * units::cm, 0);
        perf.dump("clustering extend 3", live_grouping);

        // extension ones connected to dead region ...
        if (i == 0) {
            clustering_extend(live_grouping, cluster_connected_dead, 4, 60 * units::cm, i);
        }
        else {
            clustering_extend(live_grouping, cluster_connected_dead, 4, 35 * units::cm, i);
        }
        perf.dump("clustering extend 4", live_grouping);
    }

    /// PLACEHOLDER: just to test the function
    std::map<int, std::pair<double, double>> dead_u_index;
    std::map<int, std::pair<double, double>> dead_v_index;
    std::map<int, std::pair<double, double>> dead_w_index;
    clustering_separate(live_grouping, dead_u_index, dead_v_index, dead_w_index);
    perf.dump("clustering_separate", live_grouping);

    // BEE debug dead-live
    // if (!m_bee_dir.empty()) {
    //     std::string sub_dir = String::format("%s/%d", m_bee_dir, ident);
    //     dump_bee(*root_live.get(), String::format("%s/%d-dead-live.json", sub_dir, ident));
    //     perf("dump live clusters to bee");
    // }
    fill_bee_points(m_bee_ld, *root_live.get());
    perf("dump live clusters to bee");

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
