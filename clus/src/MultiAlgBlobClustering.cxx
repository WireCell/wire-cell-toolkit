#include "WireCellClus/MultiAlgBlobClustering.h"
#include "WireCellClus/Facade.h"
#include <WireCellClus/ClusteringFuncs.h>
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/ExecMon.h"
#include "WireCellUtil/String.h"
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
  , m_bee_ld("uboone", "clustering")
  , m_bee_dead("channel-deadarea", 1*units::mm, 3) // tolerance, minpts
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
    // m_x_boundary_low_limit = get(cfg, "x_boundary_low_limit", m_x_boundary_low_limit);
    // m_x_boundary_high_limit = get(cfg, "x_boundary_high_limit", m_x_boundary_high_limit);

    m_func_cfgs = cfg["func_cfgs"];


    m_perf = get(cfg, "perf", m_perf);

    m_anode = Factory::find_tn<IAnodePlane>(cfg["anode"].asString());

    m_face = get<int>(cfg, "face", 0);

    m_bee_img.detector(get<std::string>(cfg, "bee_detector", "uboone"));
    m_bee_img.algorithm(String::format("%s-%d-%d", m_bee_img.algorithm().c_str(), m_anode->ident(), m_face));
    m_bee_ld.detector(get<std::string>(cfg, "bee_detector", "uboone"));
    m_bee_ld.algorithm(String::format("%s-%d-%d", m_bee_ld.algorithm().c_str(), m_anode->ident(), m_face));

    m_geomhelper = Factory::find_tn<IClusGeomHelper>(cfg["geom_helper"].asString());
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
        if (m_use_config_rse) {
            // Set RSE in the sink
            m_sink.set_rse(m_runNo, m_subRunNo, m_eventNo);
        }
        // Use default behavior
        reset_bee(ident, m_bee_img);
        reset_bee(ident, m_bee_ld);
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
    grouping->set_anode(m_anode);
    grouping->set_params(m_geomhelper->get_params(m_anode->ident(), m_face));
    perf("loaded live clusters");

    // log->debug("Got live pctree with {} children", root_live->nchildren());
    // log->debug(em("got live pctree"));
    log->debug("as_pctree from \"{}\"", inpath + "/dead");
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

    std::map<int, std::pair<double, double>>& dead_u_index = live_grouping.get_dead_winds(0, 0);
    std::map<int, std::pair<double, double>>& dead_v_index = live_grouping.get_dead_winds(0, 1);
    std::map<int, std::pair<double, double>>& dead_w_index = live_grouping.get_dead_winds(0, 2);
    log->debug("dead_u_index size {}", dead_u_index.size());
    log->debug("dead_v_index size {}", dead_v_index.size());
    log->debug("dead_w_index size {}", dead_w_index.size());

#define __HIDE__

    for (const auto& func_cfg : m_func_cfgs) {
        std::cout << "func_cfg: " << func_cfg << std::endl;
        auto func = getClusteringFunction(func_cfg);
        func(live_grouping, dead_grouping, cluster_connected_dead);
        perf.dump(func_cfg["name"].asString(), live_grouping);
    }
#ifndef __HIDE__
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

    // log->debug("clustering_separate nclusters {}", live_grouping.nchildren());
    clustering_separate(live_grouping, true);
    // log->debug("clustering_separate nclusters {}", live_grouping.nchildren());
    perf.dump("clustering_separate", live_grouping);

    const auto &tp = live_grouping.get_params();
    auto global_point_cloud = std::make_shared<DynamicPointCloud>(tp.angle_u, tp.angle_v, tp.angle_w);
    for (const Cluster *cluster : live_grouping.children()) {
        global_point_cloud->add_points(cluster, 0);
    }
    std::vector<double> cluster_lengths;
    for (const Cluster *cluster : live_grouping.children()) {
        cluster_lengths.push_back(cluster->get_length());
    }
    std::cout << "large clusters 10mm " << std::count_if(cluster_lengths.begin(), cluster_lengths.end(), [](double x) { return x > 10*units::mm; }) << std::endl;
    std::cout << "large clusters 5cm " << std::count_if(cluster_lengths.begin(), cluster_lengths.end(), [](double x) { return x > 5*units::cm; }) << std::endl;
    clustering_connect1(live_grouping, global_point_cloud, dead_u_index, dead_v_index, dead_w_index);
    cluster_lengths.clear();
    for (const Cluster *cluster : live_grouping.children()) {
        cluster_lengths.push_back(cluster->get_length());
    }
    std::cout << "large clusters 10mm " << std::count_if(cluster_lengths.begin(), cluster_lengths.end(), [](double x) { return x > 10*units::mm; }) << std::endl;
    std::cout << "large clusters 5cm " << std::count_if(cluster_lengths.begin(), cluster_lengths.end(), [](double x) { return x > 5*units::cm; }) << std::endl;
    perf.dump("clustering_connect1", live_grouping);

    clustering_deghost(live_grouping, dead_u_index, dead_v_index, dead_w_index, true);
    perf.dump("clustering_deghost", live_grouping);

    clustering_examine_x_boundary(live_grouping);
    perf.dump("clustering_examine_x_boundary", live_grouping);

    clustering_protect_overclustering(live_grouping);
    perf.dump("clustering_protect_overclustering", live_grouping);
#endif

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
