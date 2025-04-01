#include "WireCellClus/MultiAlgBlobClustering.h"
#include "WireCellClus/Facade.h"
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

    m_func_cfgs = cfg["func_cfgs"];


    m_perf = get(cfg, "perf", m_perf);

    for (const auto& aname : cfg["anodes"]) {
        auto anode = Factory::find_tn<IAnodePlane>(aname.asString());
        m_anodes.push_back(anode);
    }

    m_dv = Factory::find_tn<IDetectorVolumes>(cfg["detector_volumes"].asString());

    m_face = get<int>(cfg, "face", 0);

    m_bee_img.detector(get<std::string>(cfg, "bee_detector", "uboone"));
    m_bee_img.algorithm(String::format("%s-%d-%d", m_bee_img.algorithm().c_str(), m_anodes.front()->ident(), m_face));
    log->debug("m_bee_img.algorithm: {}", m_bee_img.algorithm());
    m_bee_ld.detector(get<std::string>(cfg, "bee_detector", "uboone"));
    m_bee_ld.algorithm(String::format("%s-%d-%d", m_bee_ld.algorithm().c_str(), m_anodes.front()->ident(), m_face));
    log->debug("m_bee_ld.algorithm: {}", m_bee_ld.algorithm());

    // m_geomhelper = Factory::find_tn<IClusGeomHelper>(cfg["geom_helper"].asString());

    m_dump_json = get<bool>(cfg, "dump_json", false);
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
        for (const auto bnode : cnode->children()) {  // blobs ...
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
    log->debug("loading tensor set ident={} (last={})", ident, m_last_ident);
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
    dead_grouping.set_anodes(m_anodes);
    dead_grouping.set_detector_volumes(m_dv);
    // dead_grouping.set_params(m_geomhelper->get_params(m_anodes.front()->ident(), m_face));
    

    //perf.dump("original live clusters", live_grouping, false, false);
    //perf.dump("original dead clusters", dead_grouping, false, false);

    perf.dump("pre clustering", live_grouping);

  

    for (const auto& func_cfg : m_func_cfgs) {
        // std::cout << "func_cfg: " << func_cfg << std::endl;
        auto func = getClusteringFunction(func_cfg);

        func(live_grouping, dead_grouping, cluster_connected_dead);

        perf.dump(func_cfg["name"].asString(), live_grouping);
    }

    fill_bee_points(m_bee_ld, *root_live.get());
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
