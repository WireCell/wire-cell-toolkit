#include "WireCellClus/MultiAlgBlobClustering.h"
#include "WireCellClus/Facade_Summary.h"
#include "WireCellClus/PRSegment.h"
#include "WireCellClus/PRVertex.h"
#include "WireCellClus/PRShower.h"
#include "WireCellClus/PRShowerFunctions.h"  // doc pr/84 r2: shower_get_closest_point
#include "WireCellClus/PRSegmentFunctions.h" // doc pr/93 r4: segment_orphan_confident_track
#include "WireCellClus/TrackFitting.h"


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
#include "WireCellUtil/GraphTools.h"

#include <chrono>
#include <algorithm>
#include <map>
#include <set>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cstdio>

WIRECELL_FACTORY(MultiAlgBlobClustering, WireCell::Clus::MultiAlgBlobClustering, WireCell::INamed,
                 WireCell::ITensorSetFilter, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Aux;
using namespace WireCell::Aux::TensorDM;
using namespace WireCell::Clus::Facade;
using namespace WireCell::PointCloud::Tree;
using WireCell::GraphTools::mir;

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

static std::string format_flag_names(const std::set<std::string>& flag_names)
{
    std::ostringstream ss;
    ss << "[";
    bool first = true;
    for (const auto& flag_name : flag_names) {
        if (!first) {
            ss << ",";
        }
        ss << flag_name;
        first = false;
    }
    ss << "]";
    return ss.str();
}

// NOT static: the larwirecell QLMatching plugin links against this symbol,
// forward-declaring it as WireCell::Clus::Facade::normalize_cluster_flags
// (no installed header carries it), so define it in that namespace.
namespace WireCell::Clus::Facade {
void normalize_cluster_flags(Grouping& grouping, Log::logptr_t log, const std::string& grouping_name, int ident)
{
    std::set<std::string> flag_names;
    for (const auto* cluster : grouping.children()) {
        for (const auto& flag_name : cluster->flag_names()) {
            flag_names.insert(flag_name);
        }
    }

    SPDLOG_LOGGER_DEBUG(log, "normalize_cluster_flags: ident={} grouping={} nclusters={} all_flags={}",
                        ident, grouping_name, grouping.children().size(), format_flag_names(flag_names));

    if (flag_names.empty()) {
        return;
    }

    size_t nmissing = 0;
    for (auto* cluster : grouping.children()) {
        const auto cluster_flags = cluster->flag_names();
        const std::set<std::string> cluster_flag_set(cluster_flags.begin(), cluster_flags.end());

        for (const auto& flag_name : flag_names) {
            if (cluster_flag_set.count(flag_name)) {
                continue;
            }
            cluster->set_flag(flag_name, 0);
            ++nmissing;
        }
    }
    SPDLOG_LOGGER_DEBUG(log, "normalize_cluster_flags: ident={} grouping={} added={} missing flag values",
                        ident, grouping_name, nmissing);
}
}  // namespace WireCell::Clus::Facade


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
        SPDLOG_LOGGER_DEBUG(log, "the 'bee_dir' option is no longer supported, instead use 'bee_zip' to name a .zip file");
    }
    std::string bee_zip = get<std::string>(cfg, "bee_zip", "mabc.zip");
    // Add new configuration option for initial index
    m_initial_index = get<int>(cfg, "initial_index", m_initial_index);

    // Optional shared Bee sink: when "bee_sink" names an IBeeSink component, all
    // Bee writes go to that single shared zip (at an explicit per-event index)
    // and m_sink is left unused, so a multi-node chain emits one .zip.
    std::string bee_sink_tn = get<std::string>(cfg, "bee_sink", "");
    if (!bee_sink_tn.empty()) {
        m_shared_sink = Factory::find_tn<IBeeSink>(bee_sink_tn);
        m_shared_sink->acquire();
        m_use_shared_sink = true;
        m_bee_event_index = m_initial_index;
        log->debug("using shared Bee sink {} (own bee_zip disabled)", bee_sink_tn);
    }
    else {
        //std::cout << "Xin: " << m_initial_index << " " << bee_zip << std::endl;
        m_bee_zip = bee_zip;
        // A '%' means one zip per event, opened lazily by ensure_own_sink().
        m_bee_zip_templated = m_bee_zip.find('%') != std::string::npos;
        if (!m_bee_zip_templated) {
            m_sink.reset(bee_zip, m_initial_index);  // Use the new reset with initial index
        }
        else {
            log->debug("own Bee zip, one per event, name template {}", m_bee_zip);
        }
    }

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
    // Take the event number from the per-event tensor ident (raw; run/subrun = 0).
    m_rse_from_ident = get(cfg, "rse_from_ident", m_rse_from_ident);

    // Same, but keep the configured run/subrun -- with an optional per-ident
    // override table, because a group of events can span several runs.
    m_event_from_ident = get(cfg, "event_from_ident", m_event_from_ident);
    if (cfg.isMember("rse_map")) {
        const auto& rm = cfg["rse_map"];
        for (const auto& key : rm.getMemberNames()) {
            const auto& v = rm[key];
            if (!v.isArray() || v.size() < 2) continue;
            // Keys come from a jsonnet-built object, so a non-numeric one is
            // a config bug, not user input -- say which key rather than let
            // stoi throw an unlabelled invalid_argument out of configure().
            try {
                m_rse_map[std::stoi(key)] = {v[0].asInt(), v[1].asInt()};
            }
            catch (const std::exception&) {
                THROW(ValueError() << errmsg{"MultiAlgBlobClustering: rse_map key is not an integer: " + key});
            }
        }
        log->debug("rse_map with {} entries", m_rse_map.size());
    }

    m_grouping2file_prefix = get(cfg, "grouping2file_prefix", m_grouping2file_prefix);

    m_save_deadarea = get(cfg, "save_deadarea", m_save_deadarea);
    m_save_real_cluster_id = get(cfg, "save_real_cluster_id", m_save_real_cluster_id);
    m_save_assoc_cluster_id = get(cfg, "save_assoc_cluster_id", m_save_assoc_cluster_id);
    m_real_cluster_id_global = get(cfg, "real_cluster_id_global", m_real_cluster_id_global);
    m_dead_area_version = get(cfg, "dead_area_version", m_dead_area_version);

    m_save_opflash = get(cfg, "save_opflash", m_save_opflash);
    if (m_save_opflash) {
        m_bee_flash = Bee::Flashes(get<std::string>(cfg, "bee_detector", "uboone"), "op");
    }
    m_bee_flash_per_flash = get(cfg, "bee_flash_per_flash", m_bee_flash_per_flash);
    m_bee_flash_pred_min = get(cfg, "bee_flash_pred_min", m_bee_flash_pred_min);
    m_flash_group_window = get(cfg, "flash_group_window", m_flash_group_window);
    m_flash_group_greedy = get(cfg, "flash_group_greedy", m_flash_group_greedy);

    m_dead_live_overlap_offset = get(cfg, "dead_live_overlap_offset", m_dead_live_overlap_offset);

    for (auto jtn : cfg["pipeline"]) {
        std::string tn = jtn.asString();
        SPDLOG_LOGGER_DEBUG(log, "configuring clustering method: {}", tn);
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
            bpc.dQdx_scale = get<double>(bps, "dQdx_scale", 1.0);
            bpc.dQdx_offset = get<double>(bps, "dQdx_offset", 0.0);
            bpc.use_associate_points = get<bool>(bps, "use_associate_points", false);
            bpc.use_graph_vertices = get<bool>(bps, "use_graph_vertices", false);
            bpc.require_pr_graph   = get<bool>(bps, "require_pr_graph", false);
            // Prototype-parity options; absent => false => byte-identical legacy output.
            bpc.particle_ids = get<bool>(bps, "particle_ids", false);
            bpc.pseudo_shower_track_paint = get<bool>(bps, "pseudo_shower_track_paint", false);
            bpc.include_vertex_points = get<bool>(bps, "include_vertex_points", false);
            // Dump this set at the pre-clustering point (like the special "img"
            // set) even if its name isn't "img".  Absent => false => legacy
            // name-based routing => byte-identical.
            bpc.prepipeline = get<bool>(bps, "prepipeline", false);

            // Optional drift-side / APA grouping (additive; absent -> unchanged)
            if (bps.isMember("apa_groups")) {
                for (const auto& grp : bps["apa_groups"]) {
                    ApaGroup ag;
                    ag.name = get<std::string>(grp, "name", "");
                    if (grp.isMember("apas")) {
                        for (const auto& a : grp["apas"]) {
                            ag.apas.insert(a.asInt());
                        }
                    }
                    bpc.apa_groups.push_back(ag);
                }
            }

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
            }else if (!bpc.apa_groups.empty()) { // one bee instance per APA group
                for (const auto& grp : bpc.apa_groups) {
                    m_bee_points[bpc.name].by_group[grp.name] =
                        Bee::Points(bpc.detector,
                                    String::format("%s-%s", bpc.algorithm.c_str(), grp.name.c_str()));
                }
            }else{
                m_bee_points[bpc.name].global.detector(bpc.detector);
                m_bee_points[bpc.name].global.algorithm(String::format("%s-global", bpc.name));
                // std::cout << "Test: Global: " << m_bee_points[bpc.name].global.algorithm() << std::endl;
            }
            
            SPDLOG_LOGGER_DEBUG(log, "Configured bee points set: {}, algorithm: {}, individual: {}", 
                        bpc.name, bpc.algorithm, bpc.individual ? "true" : "false");
        }
    } 

    // Configure particle-flow bee output (mc tree)
    if (cfg.isMember("bee_pf")) {
        for (const auto& pf : cfg["bee_pf"]) {
            BeePFConfig pfc;
            pfc.name     = get<std::string>(pf, "name",     "mc");
            pfc.visitor  = get<std::string>(pf, "visitor",  "");
            pfc.grouping = get<std::string>(pf, "grouping", "live");
            // Prototype-parity options; absent => legacy output, byte-identical.
            pfc.prototype_names = get<bool>(pf, "prototype_names", false);
            pfc.em_ke_min = get<double>(pf, "em_ke_min", 0.0);
            pfc.np_ke_min = get<double>(pf, "np_ke_min", 0.0);
            // doc pr/34 §10 port-fidelity knobs; absent => legacy, byte-identical.
            pfc.pf_track_main_cluster_only = get<bool>(pf, "pf_track_main_cluster_only", false);
            pfc.pf_track_bridged_clusters = get<bool>(pf, "pf_track_bridged_clusters", false);  // doc pr/40 round 9 B2
            pfc.pf_shower_vertex_barrier = get<bool>(pf, "pf_shower_vertex_barrier", false);
            pfc.pf_shower_parent_precedence = get<bool>(pf, "pf_shower_parent_precedence", false);
            pfc.pf_pi0_node_per_id = get<bool>(pf, "pf_pi0_node_per_id", false);
            pfc.pf_pdg_name_prototype_fallback = get<bool>(pf, "pf_pdg_name_prototype_fallback", false);
            // doc pr/38 Round 4; absent => legacy flat orphan roots, byte-identical.
            pfc.pf_orphan_track_parentage = get<bool>(pf, "pf_orphan_track_parentage", false);
            // doc pr/65 round 3; absent => legacy fabricated orphan roots, byte-identical.
            pfc.pf_orphan_audit_only = get<bool>(pf, "pf_orphan_audit_only", false);
            // doc pr/84 round 2; absent => legacy pseudo-parent rendering, byte-identical.
            pfc.pf_direct_when_touching = get<bool>(pf, "pf_direct_when_touching", false);
            pfc.pf_touch_max = get<double>(pf, "pf_touch_max", pfc.pf_touch_max);
            // doc 77 round 1 (2026-08-24): pf_touch_cross_main/_max removed --
            // zero movers on all 7 census candidates, F1.0 probe failure
            // (pr/84 sec 607/622).  See sbnd_xin/docs/77_knob-ledger.tsv.
            pfc.pf_pseudo_gap_from_main = get<bool>(pf, "pf_pseudo_gap_from_main", false);
            pfc.pf_unique_node_ids = get<bool>(pf, "pf_unique_node_ids", false);
            // doc pr/92; absent => legacy (dropped-satellite set unread), byte-identical.
            pfc.pf_drop_stray_satellites = get<bool>(pf, "pf_drop_stray_satellites", false);
            // doc pr/93 round 4; absent => legacy audit-only orphans / legacy
            // shower-view vertex precedence, byte-identical.
            pfc.pf_orphan_confident_track = get<bool>(pf, "pf_orphan_confident_track", false);
            pfc.pf_orphan_track_min = get<double>(pf, "pf_orphan_track_min", pfc.pf_orphan_track_min);
            // doc pr/123 round 2; absent => legacy (guard-freed tracks not re-rooted), byte-identical.
            pfc.pf_orphan_guard_freed = get<bool>(pf, "pf_orphan_guard_freed", false);
            pfc.pf_track_owns_loose_vertex = get<bool>(pf, "pf_track_owns_loose_vertex", false);
            m_bee_pf_configs.push_back(pfc);
            m_bee_pf_trees[pfc.name] = Bee::ParticleTree(pfc.name);
            SPDLOG_LOGGER_DEBUG(log, "Configured bee_pf: name={} visitor={}", pfc.name, pfc.visitor);
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
                // dead_area_version=2 emits the wire-cell-bee3 v2 wrapper with
                // tpc=apa so the slab lands on the correct anode face.  Use -1
                // for v1 (legacy bare-array JSON, default).
                int tpc = (m_dead_area_version >= 2) ? apa : -1;
                m_bee_dead_patches[apa].insert({face,Bee::Patches(name, 1*units::mm, 3, tpc)}); // Same parameters as the global one
            }
        }
    }

    // Optional dead-area drift-side grouping (mirrors the clustering grouping).
    if (cfg.isMember("dead_apa_groups")) {
        for (const auto& grp : cfg["dead_apa_groups"]) {
            ApaGroup ag;
            ag.name = get<std::string>(grp, "name", "");
            if (grp.isMember("apas")) {
                for (const auto& a : grp["apas"]) {
                    ag.apas.insert(a.asInt());
                }
            }
            m_dead_apa_groups.push_back(ag);
        }
    }
    if (m_save_deadarea && !m_dead_apa_groups.empty()) {
        for (const auto& ag : m_dead_apa_groups) {
            std::string name = String::format("channel-deadarea-%s", ag.name.c_str());
            // All APAs in a group share the same anode-x / drift direction, so
            // any member's tpc index places the v2 slab correctly; use the first.
            int tpc = -1;
            if (m_dead_area_version >= 2 && !ag.apas.empty()) tpc = *ag.apas.begin();
            m_bee_dead_groups.insert({ag.name, Bee::Patches(name, 1*units::mm, 3, tpc)});
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
    cfg["save_real_cluster_id"] = m_save_real_cluster_id;
    cfg["save_assoc_cluster_id"] = m_save_assoc_cluster_id;
    cfg["real_cluster_id_global"] = m_real_cluster_id_global;
    cfg["bee_flash_pred_min"] = m_bee_flash_pred_min;

    // Add the new parameter to default configuration
    cfg["initial_index"] = m_initial_index;

    cfg["dead_live_overlap_offset"] = m_dead_live_overlap_offset;

    cfg["use_config_rse"] = false;  // By default, don't use configured RSE
    cfg["runNo"] = m_runNo;
    cfg["subRunNo"] = m_subRunNo;
    cfg["eventNo"] = m_eventNo;
    cfg["rse_from_ident"] = m_rse_from_ident;
    cfg["event_from_ident"] = m_event_from_ident;

    return cfg;
}

void MultiAlgBlobClustering::finalize()
{
    flush();
    if (m_use_shared_sink) {
        // Shared zip is reference-counted: it is closed on the last release(),
        // independent of which node finalizes last.
        m_shared_sink->release();
    }
    else if (!m_bee_zip_templated || m_bee_zip_open_evt >= 0) {
        m_sink.close();
    }
}

void MultiAlgBlobClustering::apply_event_from_ident(int ident)
{
    m_eventNo = ident;
    auto it = m_rse_map.find(ident);
    if (it != m_rse_map.end()) {
        m_runNo = it->second.first;
        m_subRunNo = it->second.second;
    }
    // else: keep the configured run/subrun -- correct whenever the whole group
    // is one (run, subrun), and the honest fallback when the caller supplied no
    // table.
    if (!m_use_shared_sink) {
        m_sink.set_rse(m_runNo, m_subRunNo, m_eventNo);
    }
}

void MultiAlgBlobClustering::ensure_own_sink()
{
    if (!m_bee_zip_templated) return;
    if (m_eventNo == m_bee_zip_open_evt) return;
    if (m_bee_zip_open_evt >= 0) {
        m_sink.close();
    }
    const std::string name = String::format(m_bee_zip, m_eventNo);
    m_sink.reset(name, m_initial_index);
    m_sink.set_rse(m_runNo, m_subRunNo, m_eventNo);
    m_bee_zip_open_evt = m_eventNo;
    log->debug("own Bee zip {} (event {})", name, m_eventNo);
}

size_t MultiAlgBlobClustering::write_obj(const WireCell::Bee::Object& obj)
{
    if (m_use_shared_sink) {
        return m_shared_sink->write(obj, m_bee_event_index, m_runNo, m_subRunNo, m_eventNo);
    }
    ensure_own_sink();
    return m_sink.write(obj);
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

    write_obj(bpts);
    reset_bee(ident, bpts);
}

void MultiAlgBlobClustering::flush(int ident)
{
    // Nothing has been loaded since the last flush -- this is the second,
    // redundant flush that arrives when EOS is followed by finalize().  Return
    // rather than run the whole empty sweep: every write below is already
    // no-op, but the unconditional ++m_bee_event_index at the bottom is not,
    // and in a multi-event run the Bee index must count event boundaries, not
    // flush calls.
    if (m_last_ident < 0 && ident < 0) {
        return;
    }
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
        bool grouped = (it != m_bee_points_configs.end()) ? !it->apa_groups.empty() : false;

        if (individual) {
            // Write individual bee points
            for (auto& [anode_id, face_map] : apa_bpts.by_apa_face) {
                for (auto& [face, bpts] : face_map) {
                    if (!bpts.empty()) {
                        write_obj(bpts);
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
        } else if (grouped) {
            // Write one bee instance per APA group
            for (auto& [gname, bpts] : apa_bpts.by_group) {
                if (!bpts.empty()) {
                    write_obj(bpts);
                    int run = 0, evt = 0;
                    if (ident > 0) {
                        run = (ident >> 16) & 0x7fff;
                        evt = (ident) & 0xffff;
                    }
                    bpts.reset(evt, 0, run);
                }
            }
        } else {
            // Write global bee points
            if (!apa_bpts.global.empty()) {
                write_obj(apa_bpts.global);
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

        if (!m_dead_apa_groups.empty()) {
            // Flush one dead-area patch set per drift-side group
            for (auto& [gname, patches] : m_bee_dead_groups) {
                if (patches.size()) {
                    patches.flush();
                    write_obj(patches);
                    patches.clear();
                }
            }
        } else {
            // Flush individual patches
            for (auto& [apa, face_map] : m_bee_dead_patches) {
                for (auto& [face, patches] : face_map) {
                    if (patches.size()) {
                        patches.flush();
                        write_obj(patches);
                        patches.clear();
                    }
                }
            }
        }
    }

    // Flush the optical flash / charge-light "op" display.
    if (m_save_opflash && !m_bee_flash.empty()) {
        write_obj(m_bee_flash);
        int run = 0, evt = 0;
        if (ident > 0) {
            run = (ident >> 16) & 0x7fff;
            evt = (ident) & 0xffff;
        }
        m_bee_flash.reset(evt, 0, run);
    }

    // Flush particle-flow mc trees
    for (auto& [name, tree] : m_bee_pf_trees) {
        if (!tree.empty()) {
            write_obj(tree);
            tree.reset();
        }
    }

    // One flush(int) call == one event boundary (or EOS).  In shared-sink mode
    // advance the explicit per-event index so the next event's objects land at
    // the next Bee index.  Incrementing unconditionally (not only when this
    // node wrote something) keeps all nodes' indices aligned by event ordinal
    // even when a given node has an empty event.
    if (m_use_shared_sink) {
        ++m_bee_event_index;
    }

    m_last_ident = ident;
}



// Helper function remains the same as in the previous response

// real_cluster_id_global (doc 53): collapse the TWO ident numbering epochs the
// "real_cluster_id" array mixes into one event-wide epoch.
//
// merge_clusters() records each member's ident() as of the moment it runs
// (ClusteringFuncs.cxx), the save-time fill-in writes the CURRENT ident, and
// Grouping::enumerate_idents() re-runs after every visitor -- so the two are
// different dense 1..N numberings and nothing in the array says which one a row
// belongs to.  Measured on the SBND d52ron 30-event set: 31% of distinct values
// named two different clusters.
//
// Representative rows (real_cluster_main != 0) take the cluster's own current
// ident; each other pre-merge group takes a fresh id above the largest ident in
// the grouping.  So:
//   - group MEMBERSHIP is untouched => every consumer that compares rows within
//     one cluster (ClusteringUnmergeBundle's split, TaggerCheckTGM
//     main_component_mode="real") is unchanged;
//   - a cluster nothing merged is rewritten to exactly the values it had;
//   - the value becomes a valid event-wide key, and a merged cluster's
//     representative rows carry that cluster's own id.
//
// WHERE THIS RUNS MATTERS.  It is called once per event right after the
// clustering pipeline and BEFORE the Bee fills, so the Bee per-blob label
// (fill_bee_points_from_cluster's real_clid) and the saved pctree carry the SAME
// ids.  It used to sit in the tensor-save block, which is after every
// fill_bee_points() -- that fixed the tarball and left the Bee zips on the old
// two-epoch labels.  Per-visitor Bee dumps (trace_bee) are snapshots taken
// mid-pipeline and necessarily keep whatever ids existed at that step.
//
// Deliberately NOT gated on save_real_cluster_id: examine_bundles writes the
// array in memory whether or not the tarball is saved, so the Bee labels are
// wrong either way.  Detectors that never write it (PDHD, PDVD -- their
// examine_bundles is disabled) see a structural no-op.
//
// Deterministic: clusters in children() (tree) order, ids handed out in
// ascending old-value order.  No pointer-keyed iteration.  Fails open on a
// row-count mismatch, as every other provenance carve does.
void MultiAlgBlobClustering::restamp_real_cluster_id(Grouping& grouping) const
{
    int next = 0;
    for (const Cluster* cluster : grouping.children()) {
        next = std::max(next, cluster->ident());
    }
    ++next;
    size_t nstamped = 0;
    for (Cluster* cluster : grouping.children()) {
        if (!cluster->has_pcarray<int>("real_cluster_id", "perblob")) continue;
        if (!cluster->has_pcarray<int>("real_cluster_main", "perblob")) continue;
        const size_t nb = cluster->nchildren();
        auto rid = cluster->get_pcarray<int>("real_cluster_id", "perblob");
        auto rmain = cluster->get_pcarray<int>("real_cluster_main", "perblob");
        if (rid.size() != nb || rmain.size() != nb) {
            log->warn("real_cluster_id_global: cluster {} has {}/{} rows for {} "
                      "blobs, not re-stamped", cluster->ident(),
                      rid.size(), rmain.size(), nb);
            continue;
        }
        std::map<int, int> remap;          // old id -> fresh id, ascending
        for (size_t i = 0; i < nb; ++i) {
            if (rmain[i] == 0) remap.emplace(rid[i], 0);
        }
        for (auto& [old_id, fresh] : remap) {
            (void) old_id;
            fresh = next++;
        }
        std::vector<int> out(nb);
        for (size_t i = 0; i < nb; ++i) {
            out[i] = rmain[i] != 0 ? cluster->ident() : remap.at(rid[i]);
        }
        cluster->put_pcarray(out, "real_cluster_id", "perblob");
        ++nstamped;
    }
    if (nstamped) {
        log->debug("real_cluster_id_global: re-stamped {} cluster(s) into one epoch, "
                   "ids 1..{}", nstamped, next - 1);
    }
}

void MultiAlgBlobClustering::fill_bee_points(const std::string& name, const Grouping& grouping)
{
    // std::cout << "Test: " << name << " " << grouping.wpids().size() << std::endl;
    
    if (m_bee_points.find(name) == m_bee_points.end()) {
        SPDLOG_LOGGER_WARN(log, "Bee points set '{}' not found, skipping", name);
        return;
    }
    
    auto& apa_bpts = m_bee_points[name];
    
    // Find the configuration for this name
    auto it = std::find_if(m_bee_points_configs.begin(), m_bee_points_configs.end(),
                          [&name](const BeePointsConfig& cfg) { return cfg.name == name; });
    
    if (it == m_bee_points_configs.end()) {
        SPDLOG_LOGGER_WARN(log, "Configuration for bee points set '{}' not found, skipping", name);
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
        for (auto& [gname, bpts] : apa_bpts.by_group) {
            bpts.rse(m_runNo, m_subRunNo, m_eventNo);
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
        for (auto& [gname, bpts] : apa_bpts.by_group) {
            bpts.reset(evt, 0, run);
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
                        fill_bee_points_from_cluster(it2->second, *cluster, config.pcname, config.coords, config.filter, config.dQdx_scale, config.dQdx_offset);
                    }
                }
            }
        }
    }else if (!config.apa_groups.empty()){ // route each cluster to its APA group
        for (const auto* cluster : grouping.children()) {
            // Collect the APA(s) this cluster's blobs live on.
            std::set<int> capas;
            for (const auto& w : cluster->wpids_blob_set()) {
                capas.insert(w.apa());
            }
            // Route the whole cluster to the first group that owns one of its
            // APAs (a rare cross-group cluster lands wholly in one group).
            for (const auto& grp : config.apa_groups) {
                bool match = false;
                for (int a : capas) { if (grp.apas.count(a)) { match = true; break; } }
                if (!match) continue;
                auto it = apa_bpts.by_group.find(grp.name);
                if (it != apa_bpts.by_group.end()) {
                    fill_bee_points_from_cluster(it->second, *cluster, config.pcname, config.coords, config.filter, config.dQdx_scale, config.dQdx_offset);
                }
                break;
            }
        }
    }else{ // fill in the global
        // std::cout << "Test: " << name << " " << grouping.wpids().size() << " " << grouping.nchildren() << std::endl;

        for (const auto* cluster : grouping.children()) {
            fill_bee_points_from_cluster(apa_bpts.global, *cluster, config.pcname, config.coords, config.filter, config.dQdx_scale, config.dQdx_offset);
        }
    }
}

// Fill bee points from PRGraph track trajectories
void MultiAlgBlobClustering::fill_bee_points_from_pr_graph(const std::string& name, const Grouping& grouping,
                                                           std::shared_ptr<WireCell::Clus::TrackFitting> tf_in,
                                                           bool do_reset)
{
    if (m_bee_points.find(name) == m_bee_points.end()) {
        SPDLOG_LOGGER_WARN(log, "Bee points set '{}' not found for PR graph, skipping", name);
        return;
    }

    auto& apa_bpts = m_bee_points[name];

    // Find the configuration for this name
    auto it = std::find_if(m_bee_points_configs.begin(), m_bee_points_configs.end(),
                          [&name](const BeePointsConfig& cfg) { return cfg.name == name; });

    if (it == m_bee_points_configs.end()) {
        SPDLOG_LOGGER_WARN(log, "Configuration for bee points set '{}' not found, skipping", name);
        return;
    }

    const auto& config = *it;

    // Reset RSE values for all points objects.  Skipped for the 2nd..Nth
    // bundle of a per-bundle sequence -- reset() clears the accumulated points.
    if (do_reset) {
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
    }

    // doc pr/94 Phase 4b: render the caller's fitter when given one, and take
    // the graph from THAT fitter.  grouping.get_pr_graph() is by definition
    // m_track_fitting->get_graph() (Facade_Grouping.cxx:76-79) = the unnamed
    // slot = bundle 0, so reading it here emitted bundle 0's trajectories once
    // per bundle instead of each bundle's own -- i.e. every candidate after
    // the first contributed no track_fit/shower_track points at all (SBND
    // 18255/18625, owner Bee scan 2026-08-19).  Null tf_in reproduces the
    // legacy single-candidate resolution exactly.
    auto tf_sel = tf_in ? tf_in : grouping.get_track_fitting();
    if (!tf_sel) {
        SPDLOG_LOGGER_WARN(log, "No TrackFitting in grouping for bee points set '{}'", name);
        return;
    }
    auto pr_graph = tf_sel->get_graph();
    if (!pr_graph) {
        SPDLOG_LOGGER_WARN(log, "No PR graph found in grouping for bee points set '{}'", name);
        return;
    }

    SPDLOG_LOGGER_TRACE(log, "Filling bee points '{}' from PR graph with {} vertices and {} edges",
               name, boost::num_vertices(*pr_graph), boost::num_edges(*pr_graph));

    // Build segment → shower map for shower_track mode so each point gets
    // the shower's ID as cluster_id (all points from the same shower share
    // the same color in Bee).
    std::map<PR::SegmentPtr, PR::ShowerPtr, PR::SegmentIndexCmp> seg_to_shower;
    if (config.use_associate_points) {
        // Same fitter the graph came from: a shower list from bundle 0 would
        // classify bundle i's segments as tracks (charge 0) at random.
        auto tf = tf_sel;
        if (tf) {
            for (const auto& shower : tf->get_showers()) {
                PR::IndexedVertexSet sv; PR::IndexedSegmentSet ss;
                shower->fill_sets(sv, ss, /*flag_exclude_start_segment=*/false);
                for (const auto& seg : ss) {
                    seg_to_shower[seg] = shower;
                }
            }
        }
    }

    // Iterate through all segments (edges) in the graph.  Stable edge-index
    // order: boost::edges on setS iterates in pointer order, which reorders
    // the emitted Bee point arrays run-to-run.
    int segment_count = 0;
    for (auto edge_desc : PR::ordered_edges(*pr_graph)) {
        const auto& edge_bundle = (*pr_graph)[edge_desc];
        auto segment = edge_bundle.segment;

        if (!segment) continue;

        // Encode ID as cluster_id * 1000 + segment graph index for global uniqueness
        const int cluster_id = segment->cluster() ? segment->cluster()->get_cluster_id() : 0;
        const int encoded_id = cluster_id * 1000 + static_cast<int>(segment->get_graph_index());

        if (config.use_associate_points) {
            // --- shower_track mode: use associated points, charge from shower membership ---
            // Shower membership is the authoritative shower-vs-track answer: segments
            // absorbed from other clusters may not have kShowerTrajectory/kShowerTopology
            // flags or pdg=11 updated, but they are correctly recorded in seg_to_shower
            // by the clustering step. Fall back to per-segment flags only for segments
            // that are not part of any shower (standalone shower-like segments).
            auto shower_it = seg_to_shower.find(segment);
            bool is_shower = (shower_it != seg_to_shower.end()) ||
                segment->flags_any(PR::SegmentFlags::kShowerTrajectory) ||
                segment->flags_any(PR::SegmentFlags::kShowerTopology) ||
                (segment->has_particle_info() && std::abs(segment->particle_info()->pdg()) == 11);
            // A long-muon pseudo-shower (Shower::get_particle_type() == +-13)
            // is displayed as a muon by the PF tree (make_shower_leaf reads
            // the same cached type); paint its points as track for
            // consistency.  Overrides all disjuncts: the PF verdict is the
            // shower's cached type.  Doc sbnd_xin/docs/pr/45.
            if (config.pseudo_shower_track_paint && shower_it != seg_to_shower.end() &&
                std::abs(shower_it->second->get_particle_type()) == 13) {
                is_shower = false;
            }
            const double charge = is_shower ? 15000.0 : 0.0;

            auto dpc = segment->dpcloud("associate_points");
            if (!dpc) {
                // doc pr/55 (2026-08-09): the Family-C "phantom segment" mechanism --
                // a segment that DOES get a non-empty fits() (so it draws a
                // trajectory in the track_fit layer) but has no associate_points
                // cloud at all, so it contributes 0 points to shower_track.
                // Sentinel log only, no behavior change.
                SPDLOG_LOGGER_DEBUG(log,
                    "pr55 shower_track layer: segment {} (cluster {}) has no "
                    "associate_points dpcloud -- contributes 0 points to '{}'",
                    encoded_id, cluster_id, name);
                segment_count++;
                continue;
            }
            if (dpc->npoints() == 0) {
                SPDLOG_LOGGER_DEBUG(log,
                    "pr55 shower_track layer: segment {} (cluster {}) has an empty "
                    "associate_points dpcloud (0 points) -- contributes 0 points to '{}'",
                    encoded_id, cluster_id, name);
            }
            // Use the shower's start-segment encoded ID as cluster_id when the
            // segment belongs to a shower (mirrors seg_display_id in fill_bee_pf_tree:
            // cluster_id * 1000 + seg_id), so all points from the same shower share
            // the same ID in Bee.
            const int shower_cluster_id = [&]() -> int {
                if (shower_it == seg_to_shower.end()) {
                    // particle_ids: prototype per-particle convention
                    // (NeutrinoID::fill_point_info else-branch): a segment not
                    // absorbed into a shower is its own particle,
                    // cluster_id*1000 + seg id.  Legacy (off): collapse to the
                    // plain cluster_id.
                    if (config.particle_ids) {
                        int sid = segment->id();
                        if (sid < 0) sid = static_cast<int>(segment->get_graph_index());
                        return cluster_id * 1000 + sid;
                    }
                    return cluster_id;
                }
                auto start_seg = shower_it->second->start_segment();
                if (!start_seg) return cluster_id;
                int sid = start_seg->id();
                if (sid < 0) sid = static_cast<int>(start_seg->get_graph_index());
                const auto* cl = start_seg->cluster();
                return cl ? cl->get_cluster_id() * 1000 + sid : sid;
            }();
            const size_t ndp = dpc->npoints();
            for (size_t ip = 0; ip < ndp; ++ip) {
                apa_bpts.global.append(dpc->point3d(ip), charge, cluster_id, shower_cluster_id);
            }
        } else {
            // --- default mode: use fitted points with dQdx scale/offset ---
            const auto& fits = segment->fits();

            if (fits.empty()) {
                // doc pr/55 (2026-08-09): previously silent -- this is the exact
                // mechanism behind a "phantom segment" (fit-layer contributes
                // nothing for a segment that DOES have associate points in the
                // shower_track layer, since that branch runs independently
                // above). Sentinel log only, no behavior change.
                SPDLOG_LOGGER_DEBUG(log,
                    "pr55 track_fit layer: segment {} (cluster {}) reached the Bee dump "
                    "with an empty fits() -- contributes 0 points to '{}'",
                    encoded_id, cluster_id, name);
                segment_count++;
                continue;
            }

            SPDLOG_LOGGER_TRACE(log, "Segment {} has {} fitted points", encoded_id, fits.size());

            for (const auto& fit : fits) {
                if (!fit.valid()) continue;

                const auto& point = fit.point;
                double charge = fit.dQ;
                charge = charge * config.dQdx_scale + config.dQdx_offset;
                if (charge < 0) charge = 0;

                if (config.individual) {
                    if (fit.paf.first >= 0 && fit.paf.second >= 0) {
                        int apa = fit.paf.first;
                        int face = fit.paf.second;
                        auto it_apa = apa_bpts.by_apa_face.find(apa);
                        if (it_apa != apa_bpts.by_apa_face.end()) {
                            auto it_face = it_apa->second.find(face);
                            if (it_face != it_apa->second.end()) {
                                it_face->second.append(point, charge, cluster_id, encoded_id);
                            }
                        }
                    }
                } else {
                    apa_bpts.global.append(point, charge, cluster_id, encoded_id);
                }
            }
        }

        segment_count++;
    }

    // include_vertex_points: append the PR-graph vertex fit points to the
    // fitted-point (track_fit) dump, real_cluster_id = -1, matching the
    // prototype's fill_skeleton_info_magnify vertex rows.
    if (config.include_vertex_points && !config.use_associate_points) {
        int vtx_count = 0;
        for (auto node_desc : PR::ordered_nodes(*pr_graph)) {
            const auto& node_bundle = (*pr_graph)[node_desc];
            auto vertex = node_bundle.vertex;
            if (!vertex) continue;
            const auto& fit = vertex->fit();
            if (!fit.valid()) continue;

            const int cluster_id = vertex->cluster() ? vertex->cluster()->get_cluster_id() : 0;
            double charge = fit.dQ * config.dQdx_scale + config.dQdx_offset;
            if (charge < 0) charge = 0;

            if (config.individual) {
                if (fit.paf.first >= 0 && fit.paf.second >= 0) {
                    auto it_apa = apa_bpts.by_apa_face.find(fit.paf.first);
                    if (it_apa != apa_bpts.by_apa_face.end()) {
                        auto it_face = it_apa->second.find(fit.paf.second);
                        if (it_face != it_apa->second.end()) {
                            it_face->second.append(fit.point, charge, cluster_id, -1);
                        }
                    }
                }
            } else {
                apa_bpts.global.append(fit.point, charge, cluster_id, -1);
            }
            ++vtx_count;
        }
        SPDLOG_LOGGER_TRACE(log, "Appended {} vertex fit points to bee set '{}'", vtx_count, name);
    }

    SPDLOG_LOGGER_TRACE(log, "Filled bee points '{}' from {} segments", name, segment_count);
}


void MultiAlgBlobClustering::fill_bee_vertices_from_pr_graph(const std::string& name, const Facade::Grouping& grouping,
                                                             std::shared_ptr<WireCell::Clus::TrackFitting> tf_in,
                                                             bool do_reset)
{
    if (m_bee_points.find(name) == m_bee_points.end()) {
        SPDLOG_LOGGER_WARN(log, "Bee points set '{}' not found for graph vertices, skipping", name);
        return;
    }

    auto& apa_bpts = m_bee_points[name];

    // Reset RSE.  See fill_bee_points_from_pr_graph for why do_reset exists.
    if (do_reset) {
        if (m_use_config_rse) {
            apa_bpts.global.rse(m_runNo, m_subRunNo, m_eventNo);
        } else {
            int run = 0, evt = 0;
            if (m_last_ident > 0) {
                run = (m_last_ident >> 16) & 0x7fff;
                evt = (m_last_ident) & 0xffff;
            }
            apa_bpts.global.reset(evt, 0, run);
        }
    }

    // doc pr/94 Phase 4b: the caller's fitter and ITS graph, not the unnamed
    // slot -- otherwise every bundle re-emits bundle 0's vertices.
    auto tf_sel = tf_in ? tf_in : grouping.get_track_fitting();
    if (!tf_sel) {
        SPDLOG_LOGGER_WARN(log, "No TrackFitting in grouping for vertices bee set '{}'", name);
        return;
    }
    auto pr_graph = tf_sel->get_graph();
    if (!pr_graph) {
        SPDLOG_LOGGER_WARN(log, "No PR graph found in grouping for vertices bee set '{}'", name);
        return;
    }

    int vertex_count = 0;
    for (auto node_desc : PR::ordered_nodes(*pr_graph)) {
        const auto& node_bundle = (*pr_graph)[node_desc];
        auto vertex = node_bundle.vertex;
        if (!vertex) { ++vertex_count; continue; }

        // Encode ID as cluster_id * 1000 + vertex graph index for global uniqueness
        const int cluster_id = vertex->cluster() ? vertex->cluster()->get_cluster_id() : 0;
        const int encoded_id = cluster_id * 1000 + static_cast<int>(vertex->get_graph_index());

        const WireCell::Point& point = vertex->fit().point;
        const double charge = vertex->flags_any(PR::VertexFlags::kNeutrinoVertex) ? 15000.0 : 0.0;
        apa_bpts.global.append(point, charge, cluster_id, encoded_id);
        ++vertex_count;
    }

    SPDLOG_LOGGER_TRACE(log, "Filled bee vertices '{}' from {} vertices", name, vertex_count);
}


// Helper: map PDG code to a short display name for the Bee particle tree.
// prototype=true follows the prototype's TDatabasePDG naming
// (bee/WCReader.cc PDGName: "proton"/"neutron", not "p"/"n").
// proto_fallback=true (doc pr/34 F5) fills the table's gaps the way the
// prototype does: pi0 by name, nuclei decoded from the 10LZZZAAAI code
// (bee/WCReader.cc:529-547), everything else as the PDG number -- instead of
// collapsing to "particle".
static std::string pf_pdg_to_name(int pdg, bool prototype = false, bool proto_fallback = false)
{
    switch (pdg) {
        case  11: return "e-";
        case -11: return "e+";
        case  13: return "mu-";
        case -13: return "mu+";
        case  22: return "gamma";
        case  211: return "pi+";
        case -211: return "pi-";
        case 2212: return prototype ? "proton" : "p";
        case 2112: return prototype ? "neutron" : "n";
        case  321: return "K+";
        case -321: return "K-";
        default:   break;
    }
    if (!proto_fallback) return "particle";
    if (pdg == 111) return "pi0";   // TDatabasePDG knows it; the table above did not
    if (pdg > 1000000000) {          // nuclei: 10LZZZAAAI (WCReader.cc:533-535)
        const int z = (pdg - 1000000000) / 10000;
        const int a = (pdg - 1000000000 - z * 10000) / 10;
        const char* elem = z == 18 ? "Ar" : z == 17 ? "Cl" : z == 19 ? "Ca"
                         : z == 16 ? "S"  : z == 15 ? "P"  : z == 14 ? "Si"
                         : z == 1  ? "H"  : z == 2  ? "He" : nullptr;
        if (elem) return std::string(elem) + "-" + std::to_string(a);
        return std::to_string(pdg);  // WCReader returns the bare code for unknown Z
    }
    return std::to_string(pdg);
}


// Hierarchical particle-flow tree matching the prototype "mc" Bee JSON format.
//
// Output is a bare JSON array of jsTree nodes, each node:
//   { "id":N, "text":"name  KE MeV",
//     "data":{"start":[x,y,z],"end":[x,y,z]},
//     "children":[...] }
// Leaf nodes (no children) additionally carry "icon":"jstree-file".
//
// Algorithm (mirrors prototype NeutrinoID::fill_particle_tree):
//   1. BFS from main_vertex through non-shower track segments, establishing
//      parent-child relationships among segments and recording which track
//      segment "arrived at" each vertex (vtx_incoming_seg map).
//   2. Disconnected track segments (not reachable from main_vertex) are
//      collected as additional root-level nodes so nothing is lost.
//   3. Each shower is attached under its parent track segment according to
//      start_connection_type:
//        type 1 (direct)   – nested directly as a leaf child.
//        type 2/3 (indirect/gap) – an intermediate pseudo-gamma leaf is
//                                  inserted first, then the shower under it.
//   4. Node IDs follow the prototype convention: cluster_id*1000 + seg_id.
void MultiAlgBlobClustering::fill_bee_pf_tree(const BeePFConfig& cfg,
                                               const Facade::Grouping& grouping,
                                               bool flag_print,
                                               std::shared_ptr<WireCell::Clus::TrackFitting> tf_in,
                                               std::set<int>* shared_used_ids,
                                               Configuration* out_particles)
{
    // Debug dump of the PF-tree assembly.  Previously forced on (unconditional
    // stdout spam); now opt-in via env var.  Log/stdout only -- no effect on
    // the emitted mc.json bytes.
    flag_print = flag_print || (std::getenv("WCT_BEE_PF_PRINT") != nullptr);

    auto map_it = m_bee_pf_trees.find(cfg.name);
    if (map_it == m_bee_pf_trees.end()) {
        SPDLOG_LOGGER_WARN(log, "bee_pf tree storage '{}' not found", cfg.name);
        return;
    }
    auto& tree = map_it->second;

    // doc pr/94 Phase 4: render the caller's fitter when given one.  Resolving
    // the unnamed slot implicitly is correct only while there is exactly one
    // candidate; with per-bundle fitters the unnamed slot is always bundle 0,
    // so every later bundle's flow would silently render as a repeat of it.
    auto tf = tf_in ? tf_in : grouping.get_track_fitting();
    if (!tf) return;

    // ...and take the GRAPH from that same fitter.  Grouping::get_pr_graph()
    // is defined as m_track_fitting->get_graph() (Facade_Grouping.cxx:76-79),
    // i.e. the unnamed slot again -- so reading it here would have walked
    // bundle 0's graph from bundle i's vertex.  Caught by the doc pr/94 §10.1
    // sync check on NCpi0 evt 18625, whose second bundle reconstructed a real
    // 1498 MeV candidate and emitted no Bee node at all.
    auto pr_graph = tf->get_graph();
    if (!pr_graph) return;

    auto main_vertex = tf->get_main_vertex();
    if (!main_vertex) {
        SPDLOG_LOGGER_DEBUG(log, "fill_bee_pf_tree '{}': no main vertex, skipping", cfg.name);
        return;
    }
    const auto* main_cluster = main_vertex->cluster();
    // F1 (doc pr/34 §10.2): the prototype's track loop keeps only segments in
    // the main vertex's cluster, compared by cluster ID (NeutrinoID.cxx:1488).
    // Ident, not pointer: split products inherit the parent's ident
    // (Facade_Grouping.cxx:182), so ident-equality is the weaker, correct test.
    auto same_cluster = [&](const PR::SegmentPtr& s) {
        const auto* c = s->cluster();
        return main_cluster && c && c->get_cluster_id() == main_cluster->get_cluster_id();
    };
    // doc sbnd_xin/docs/pr/40 round 9 B2: clusters graph-connected to the
    // main cluster by an nv_bridge_track bridge segment.  Widens ONLY the
    // two BFS gates below (the orphan pools stay main-cluster-only); knob
    // off, or no bridge fired => the set is empty => byte-identical.
    const auto& bridged_ids = tf->get_bridged_cluster_ids();
    auto bridged_cluster = [&](const PR::SegmentPtr& s) {
        if (!cfg.pf_track_bridged_clusters || bridged_ids.empty()) return false;
        const auto* c = s->cluster();
        return c && bridged_ids.count(c->get_cluster_id()) > 0;
    };

    const auto& showers            = tf->get_showers();
    const auto& pi0_showers        = tf->get_pi0_showers();
    const auto& map_shower_pio_id  = tf->get_map_shower_pio_id();
    const auto& map_pio_id_mass    = tf->get_map_pio_id_mass();
    // doc pr/92: satellites dropped from the kine tree; mirror the drop
    // here so PF and Enu describe the same particle set.  pi0 conjunct is
    // pure defense -- the kine side never drops pi0-paired showers.
    const auto& dropped_sat_ids = tf->get_dropped_satellite_shower_ids();
    auto sat_dropped = [&](const PR::ShowerPtr& sh) {
        return cfg.pf_drop_stray_satellites && !dropped_sat_ids.empty() &&
               !pi0_showers.count(sh) &&
               dropped_sat_ids.count(sh->get_shower_id()) > 0;
    };
    PR::IndexedSegmentSet conn4_skip_segs;

    // --- Vertex → node-descriptor map ---
    std::map<PR::VertexPtr, PR::node_descriptor, PR::VertexIndexCmp> vtx_to_nd;
    for (auto nd : PR::ordered_nodes(*pr_graph)) {
        if (auto vtx = (*pr_graph)[nd].vertex) vtx_to_nd[vtx] = nd;
    }
    if (!vtx_to_nd.count(main_vertex)) return;

    // --- Segment → shower map; collect all shower segments ---
    std::map<PR::SegmentPtr, PR::ShowerPtr, PR::SegmentIndexCmp> seg_to_shower;
    PR::IndexedSegmentSet shower_segs;
    PR::IndexedVertexSet shower_vtxs;   // F2: every vertex in any shower's view
    for (const auto& shower : showers) {
        PR::IndexedVertexSet sv; PR::IndexedSegmentSet ss;
        // doc pr/38: with the F2 barrier on, sv excludes the shower's START
        // vertex -- the prototype's barrier source map_vtx_segs never holds
        // it (`if (vtx == start_vertex) continue;`, WCShower.cxx:547,
        // :708-716, :733-745), so a main-track attachment junction stays
        // traversable and only shower-INTERIOR vertices block.  Blocking the
        // junction silently dropped every track segment behind it (SBND
        // 18255-219295 proton seg 15001 behind a conn-2 attachment;
        // 18255-489330 proton seg 4044 behind a conn-1 attachment).  Barrier
        // off: sv is unused, legacy path byte-identical.
        shower->fill_sets(sv, ss, /*flag_exclude_start_segment=*/false,
                          /*exclude_start_vertex=*/cfg.pf_shower_vertex_barrier);
        for (const auto& seg : ss) { seg_to_shower[seg] = shower; shower_segs.insert(seg); }
        shower_vtxs.insert(sv.begin(), sv.end());

        auto [_, conn_type] = shower->get_start_vertex_and_type();
        if (conn_type == 4) {
            for (const auto& seg : ss) conn4_skip_segs.insert(seg);
            if (auto start_seg = shower->start_segment()) conn4_skip_segs.insert(start_seg);
        }
    }

    // --- BFS from main_vertex through track-only (non-shower) segments ---
    //   seg_parent[S] = nullptr  →  S is a root (direct daughter of neutrino vtx)
    //   seg_parent[S] = P        →  S is a child of P
    //   seg_endpoints[S] = {near_vtx, far_vtx}  (near = toward main_vertex)
    //   vtx_incoming_seg[V] = S  →  S is the segment that first reached V from main_vertex
    std::map<PR::SegmentPtr, PR::SegmentPtr,  PR::SegmentIndexCmp> seg_parent;
    std::map<PR::SegmentPtr, std::vector<PR::SegmentPtr>, PR::SegmentIndexCmp> seg_children;
    std::map<PR::SegmentPtr, std::pair<PR::VertexPtr,PR::VertexPtr>, PR::SegmentIndexCmp> seg_endpoints;
    std::map<PR::VertexPtr,  PR::SegmentPtr,  PR::VertexIndexCmp>  vtx_incoming_seg;

    PR::IndexedVertexSet visited_vtxs;
    PR::IndexedSegmentSet used_segs = shower_segs;   // pre-mark showers as visited
    // F2 (doc pr/34 §10.3): the prototype also pre-seeds used_vertices from
    // every shower (NeutrinoID.cxx:1597-1602), so the track BFS never expands
    // THROUGH a shower vertex.  Pop-time barrier only: the segment that
    // reaches such a vertex is still claimed, matching the prototype's
    // curr_vtx check.  No extra filter at the seed/expansion guards.
    if (cfg.pf_shower_vertex_barrier) visited_vtxs.insert(shower_vtxs.begin(), shower_vtxs.end());

    visited_vtxs.insert(main_vertex);
    std::vector<std::pair<PR::VertexPtr, PR::SegmentPtr>> bfs_cur;

    // Seed BFS: all non-shower edges adjacent to main_vertex.  Stable
    // edge-index order so the PF-tree child order is run-reproducible.
    for (auto edesc : PR::sorted_out_edges(vtx_to_nd.at(main_vertex), *pr_graph)) {
        auto seg = (*pr_graph)[edesc].segment;
        if (!seg || used_segs.count(seg) || conn4_skip_segs.count(seg) ||
            (cfg.pf_track_main_cluster_only && !same_cluster(seg) && !bridged_cluster(seg))) continue;  // doc pr/40 round 9 B2
        auto far = PR::find_other_vertex(*pr_graph, seg, main_vertex);
        if (!far) continue;
        used_segs.insert(seg);
        seg_parent[seg]    = nullptr;
        seg_endpoints[seg] = {main_vertex, far};
        vtx_incoming_seg[far] = seg;
        bfs_cur.push_back({far, seg});
    }

    while (!bfs_cur.empty()) {
        std::vector<std::pair<PR::VertexPtr, PR::SegmentPtr>> bfs_next;
        for (auto& [cur_vtx, inc_seg] : bfs_cur) {
            if (visited_vtxs.count(cur_vtx)) continue;
            visited_vtxs.insert(cur_vtx);
            auto nd_it = vtx_to_nd.find(cur_vtx);
            if (nd_it == vtx_to_nd.end()) continue;
            for (auto edesc : PR::sorted_out_edges(nd_it->second, *pr_graph)) {
                auto seg = (*pr_graph)[edesc].segment;
                if (!seg || used_segs.count(seg) || conn4_skip_segs.count(seg) ||
                    (cfg.pf_track_main_cluster_only && !same_cluster(seg) && !bridged_cluster(seg))) continue;  // doc pr/40 round 9 B2
                auto far = PR::find_other_vertex(*pr_graph, seg, cur_vtx);
                if (!far) continue;
                used_segs.insert(seg);
                seg_parent[seg]    = inc_seg;
                seg_endpoints[seg] = {cur_vtx, far};
                seg_children[inc_seg].push_back(seg);
                if (!vtx_incoming_seg.count(far)) vtx_incoming_seg[far] = seg;
                bfs_next.push_back({far, seg});
            }
        }
        bfs_cur = std::move(bfs_next);
    }

    // doc pr/93 round 4 (pf_track_owns_loose_vertex): condition (a) of the
    // guard below must mean "a REAL track segment was walked here from the
    // main vertex".  vtx_incoming_seg is extended DURING the fixed point
    // below (non-root branch) with shower-derived parents, so snapshot the
    // BFS-only key set here.  Knob off => empty, no allocation.
    PR::IndexedVertexSet track_bfs_vtxs;
    if (cfg.pf_track_owns_loose_vertex) {
        for (const auto& [v, s] : vtx_incoming_seg) track_bfs_vtxs.insert(v);
    }

    // // Log disconnected non-shower track segments (not added to particle flow).
    // for (auto edge_desc : mir(boost::edges(*pr_graph))) {
    //     auto seg = (*pr_graph)[edge_desc].segment;
    //     if (!seg || seg_to_shower.count(seg) || seg_parent.count(seg)) continue;
    //     auto [va, vb] = PR::find_vertices(*pr_graph, seg);
    //     const int cluster_id = seg->cluster() ? seg->cluster()->get_cluster_id() : -1;
    //     const int graph_idx  = static_cast<int>(seg->get_graph_index());
    //     const auto& fits = seg->fits();
    //     const double length  = fits.empty() ? -1.0
    //         : PR::walk_length(fits.begin(), fits.end(),
    //                           [](const PR::Fit& f) -> WireCell::Point { return f.point; }) / units::cm;
    //     std::string pi_name = "?";
    //     if (seg->has_particle_info()) {
    //         pi_name = seg->particle_info()->name();
    //     }
    //     std::cout << "[fill_bee_pf_tree] DISCONNECTED track seg"
    //               << "  cluster=" << cluster_id
    //               << "  graph_idx=" << graph_idx
    //               << "  encoded_id=" << (cluster_id >= 0 ? cluster_id * 1000 + graph_idx : graph_idx)
    //               << "  length_cm=" << std::fixed << std::setprecision(2) << length
    //               << "  particle=" << pi_name
    //               << "  has_va=" << (va ? 1 : 0) << " " << (va ? va->fit().point : WireCell::Point(0,0,0)) << " " <<  va->wcpt().point << " " << seg->fits().size()  << " " << seg->fits()[0].point << " " << seg->fits()[1].point
    //               << "  has_vb=" << (vb ? 1 : 0) << " " << (vb ? vb->fit().point : WireCell::Point(0,0,0)) << " " << vb->wcpt().point 
    //               << "\n";
    // }

    // --- Extend vtx_incoming_seg through shower vertex sets (mirrors prototype) ---
    // The prototype guarantees a shower's start_vtx is always picked from
    // (main-cluster vertices ∪ existing shower vertices) during clustering, so it is
    // always reachable at fill time.  We replicate that guarantee here by propagating
    // vtx_incoming_seg into every vertex belonging to an already-resolved shower, then
    // repeating to a fixed point so that showers nested inside other showers resolve too.
    //
    // Two flavors of "resolved":
    //   a) start_vtx == main_vertex (or in root_reachable_vtxs)
    //      → shower hangs from root; add its vertices to root_reachable_vtxs
    //   b) start_vtx in vtx_incoming_seg
    //      → shower hangs from a track segment; extend vtx_incoming_seg with its vertices
    PR::IndexedVertexSet root_reachable_vtxs;
    std::map<PR::VertexPtr, PR::ShowerPtr, PR::VertexIndexCmp> vtx_to_parent_shower;
    {
        bool any_added = true;
        while (any_added) {
            any_added = false;
            for (const auto& shower : showers) {
                auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();
                if (conn_type == 4) continue;
                if (!start_vtx) continue;

                const bool at_main  = (start_vtx == main_vertex);
                const bool at_root  = root_reachable_vtxs.count(start_vtx) > 0;
                const bool at_track = vtx_incoming_seg.count(start_vtx) > 0;
                if (!at_main && !at_root && !at_track) continue;

                PR::SegmentPtr parent_seg = at_track ? vtx_incoming_seg.at(start_vtx) : nullptr;
                PR::ShowerPtr parent_shower = (at_main || at_root) ? shower : nullptr;

                PR::IndexedVertexSet sv; PR::IndexedSegmentSet ss;
                shower->fill_sets(sv, ss, /*flag_exclude_start_segment=*/false);
                // doc pr/93 round 4 (pf_track_owns_loose_vertex): vertices
                // that are STRUCTURE (an endpoint of a member segment) as
                // opposed to loose point-cloud association.  A shower's view
                // can hold a vertex none of whose incident segments it owns
                // -- the F12 absorb guard (shower_absorb_track_guard)
                // add_vertex()es the frontier BEFORE refusing the segment
                // beyond it, and add_shower()/add_segment(seg,true) can do
                // the same.  Deliberately EXCLUDES start_vtx: the root
                // branch's own guards already prevent self-stamping, and
                // adding it would make the attachment loop's defensive
                // self-test load-bearing.
                PR::IndexedVertexSet struct_vtxs;
                if (cfg.pf_track_owns_loose_vertex) {
                    for (const auto& mseg : ss) {
                        auto [ma, mb] = PR::find_vertices(*pr_graph, mseg);
                        if (ma && ma->descriptor_valid()) struct_vtxs.insert(ma);
                        if (mb && mb->descriptor_valid()) struct_vtxs.insert(mb);
                    }
                }
                for (const auto& vtx : sv) {
                    if (vtx == main_vertex) continue;
                    if (at_main || at_root) {
                        // hangs from root → add to root_reachable_vtxs
                        if (!root_reachable_vtxs.count(vtx)) {
                            // A root-anchored shower's vertex SET is a loose
                            // association; the track BFS above reached its
                            // vertices by walking real track segments out of
                            // the neutrino vertex.  Where the two disagree
                            // this branch silently wins, and F3's
                            // parent-shower precedence (:1419) then hangs
                            // everything anchored there under the shower
                            // instead of under the track.
                            //
                            // doc sbnd_xin/docs/pr/74 round 4: that is exactly
                            // wrong for a Michel.  18255-506746 -- once K6
                            // demotes seg 21048 to a stopping muon, the track
                            // BFS reaches its far vertex 21037, but the
                            // neighbouring 102 MeV shower's set also contains
                            // 21037 and claimed it, so the 64 MeV Michel
                            // rendered as a daughter of that shower rather
                            // than of the muon that produced it (measured:
                            // "PROPAGATE-OVER-TRACK vtx_gidx=37
                            // claimed_by_shower_ke=102.128
                            // over_incoming_seg_gidx=48").
                            //
                            // Keyed on kMuonStemGuard, which ONLY the
                            // default-OFF K6 pass ever sets -- so with K6 off
                            // no vertex is ever protected and this is
                            // byte-identical.  Narrow on purpose: the general
                            // "track BFS beats shower set" rule would
                            // restructure the tree at every vertex where the
                            // two disagree, which is a separate change with
                            // its own census.
                            auto vis = vtx_incoming_seg.find(vtx);
                            const bool track_owns_via_michel_stem =
                                vis != vtx_incoming_seg.end() && vis->second &&
                                vis->second->flags_any(PR::SegmentFlags::kMuonStemGuard);
                            // doc pr/93 round 4 (pf_track_owns_loose_vertex):
                            // the general "track BFS beats shower set" rule
                            // the comment above deferred, restricted to the
                            // LOOSE-association case: skip the claim when
                            // BOTH (a) the real track BFS walked a segment to
                            // this vertex (snapshot set -- the live map gains
                            // shower-derived entries during this fixed
                            // point), and (b) the vertex is not an endpoint
                            // of ANY member segment of this shower.  SBND
                            // 18264-69314: the 151.9cm muon's far endpoint
                            // (deg 2) claimed by the 595 MeV root shower
                            // whose nearest member is 35cm away, stealing the
                            // muon's own 67 MeV conn-1 daughter.  NOT a
                            // superset of the michel term (that one also
                            // protects structural vertices); both kept.
                            // C++ default false => byte-identical.
                            const bool track_owns_loose =
                                cfg.pf_track_owns_loose_vertex &&
                                track_bfs_vtxs.count(vtx) > 0 &&
                                struct_vtxs.count(vtx) == 0;
                            if (flag_print && vis != vtx_incoming_seg.end()) {
                                std::cout << "[fill_bee_pf_tree] PROPAGATE-OVER-TRACK"
                                          << "  vtx_gidx=" << vtx->get_graph_index()
                                          << "  claimed_by_shower_ke=" << shower->get_kine_best()/units::MeV
                                          << "  over_incoming_seg_gidx=" << vis->second->get_graph_index()
                                          << "  michel_stem_protected=" << (track_owns_via_michel_stem ? 1 : 0)
                                          << "  loose_protected=" << (track_owns_loose ? 1 : 0)
                                          << "\n";
                            }
                            if (track_owns_via_michel_stem || track_owns_loose) {
                                if (track_owns_loose && !track_owns_via_michel_stem) {
                                    SPDLOG_LOGGER_DEBUG(log,
                                        "pf_track_owns_loose_vertex: vtx_gidx={} kept by track "
                                        "seg_gidx={} (shower_id={} ke_mev={:.1f} nseg={} claim=loose)",
                                        vtx->get_graph_index(),
                                        (vis != vtx_incoming_seg.end() && vis->second)
                                            ? static_cast<long>(vis->second->get_graph_index()) : -1L,
                                        shower->get_shower_id(),
                                        shower->get_kine_best() / units::MeV,
                                        shower->get_num_segments());
                                }
                                continue;
                            }
                            root_reachable_vtxs.insert(vtx);
                            vtx_to_parent_shower[vtx] = parent_shower;
                            any_added = true;
                        }
                    } else {
                        // hangs from a track segment → extend vtx_incoming_seg
                        if (!vtx_incoming_seg.count(vtx) && !root_reachable_vtxs.count(vtx)) {
                            vtx_incoming_seg[vtx] = parent_seg;
                            // F3 (doc pr/34 §10.4): also record the parent
                            // SHOWER, exactly as the root branch does -- the
                            // half-populated map is what hangs a nested shower
                            // under the track segment instead of its shower.
                            // Inside the guard: the shower's own start vertex
                            // is already in vtx_incoming_seg, so no
                            // self-parenting.
                            if (cfg.pf_shower_parent_precedence) vtx_to_parent_shower[vtx] = shower;
                            any_added = true;
                        }
                    }
                }
            }
        }
    }

    // --- Attach each shower to its parent (track segment, other shower, or root) ---
    // type 1 (direct):    nested directly under parent as shower leaf
    // type 2/3 (indirect): a pseudo-gamma node is inserted between parent and shower
    // shower_parent_vtx[shower] = the connection vertex used by that shower
    using ShowerSegMap = std::map<PR::SegmentPtr,
                                  std::vector<std::pair<PR::ShowerPtr,PR::VertexPtr>>,
                                  PR::SegmentIndexCmp>;
    using ShowerShowerMap = std::map<PR::ShowerPtr,
                                     std::vector<std::pair<PR::ShowerPtr,PR::VertexPtr>>,
                                     PR::ShowerIndexCmp>;
    ShowerSegMap seg_direct_showers;    // seg → [(shower, conn_vtx)]
    ShowerSegMap seg_indirect_showers;  // seg → [(shower, conn_vtx)]
    ShowerShowerMap shower_direct_showers;    // shower → [(child_shower, conn_vtx)]
    ShowerShowerMap shower_indirect_showers;  // shower → [(child_shower, conn_vtx)]
    std::vector<std::pair<PR::ShowerPtr,PR::VertexPtr>> root_direct_showers;
    std::vector<std::pair<PR::ShowerPtr,PR::VertexPtr>> root_indirect_showers;

    // doc pr/38 Round 4 (pf_orphan_track_parentage): orphan TRACK segments
    // anchored at a vertex inside a shower's view attach as children of that
    // shower's displayed leaf.  Filled by the anchoring pass below (after the
    // shower-attachment loop), read by make_shower_leaf.  Empty when the knob
    // is off => legacy output.
    std::map<PR::ShowerPtr, std::vector<PR::SegmentPtr>, PR::ShowerIndexCmp> shower_child_segs;

    for (const auto& shower : showers) {
        auto [start_vtx, conn_type] = shower->get_start_vertex_and_type();

        // type 4 = "not clearly connected"; skip entirely (prototype behaviour)
        if (conn_type == 4) {
            if (flag_print) {
                auto start_seg = shower->start_segment();
                const auto* cl = start_seg ? start_seg->cluster() : nullptr;
                std::cout << "[fill_bee_pf_tree] SKIP shower (conn_type=4)"
                          << "  pdg=" << shower->get_particle_type()
                          << "  ke=" << shower->get_kine_best() / units::MeV << " MeV"
                          << "  cluster=" << (cl ? std::to_string(cl->get_cluster_id()) : "?")
                          << " nsegments=" << shower->get_num_segments() 
                          << "\n";
            }
            continue;
        }

        // doc pr/92: satellite dropped from the kine tree -- skip here too so
        // the PF tree and kine_reco_Enu describe the same particle set.  The
        // shower never enters any root/seg/shower pool, so neither
        // append_showers nor the pseudo-carrier path can resurrect it.
        if (sat_dropped(shower)) {
            if (flag_print) {
                auto start_seg = shower->start_segment();
                const auto* cl = start_seg ? start_seg->cluster() : nullptr;
                std::cout << "[fill_bee_pf_tree] SKIP shower (stray satellite, pr/92)"
                          << "  conn_type=" << conn_type
                          << "  pdg=" << shower->get_particle_type()
                          << "  ke=" << shower->get_kine_best() / units::MeV << " MeV"
                          << "  cluster=" << (cl ? std::to_string(cl->get_cluster_id()) : "?")
                          << " nsegments=" << shower->get_num_segments()
                          << "\n";
            }
            continue;
        }

        bool direct = (conn_type == 1);

        if (!start_vtx || start_vtx == main_vertex) {
            if (flag_print) {
                auto start_seg = shower->start_segment();
                const auto* cl = start_seg ? start_seg->cluster() : nullptr;
                WireCell::Point sp = shower->get_start_point();
                WireCell::Point ep = shower->get_end_point();
                std::cout << "[fill_bee_pf_tree] ROOT shower"
                          << "  conn_type=" << conn_type
                          << "  reason=" << (!start_vtx ? "null_start_vtx" : "start_vtx==main_vertex")
                          << "  pdg=" << shower->get_particle_type()
                          << "  ke=" << shower->get_kine_best() / units::MeV << " MeV"
                          << "  cluster=" << (cl ? std::to_string(cl->get_cluster_id()) : "?")
                          << "  start=(" << sp.x()/units::cm << "," << sp.y()/units::cm << "," << sp.z()/units::cm << ") cm"
                          << "  end=(" << ep.x()/units::cm << "," << ep.y()/units::cm << "," << ep.z()/units::cm << ") cm"
                        << " nsegments=" << shower->get_num_segments() 
                          << "\n";
            }
            auto& vec = direct ? root_direct_showers : root_indirect_showers;
            vec.push_back({shower, main_vertex});
        } else {
            // F3 (doc pr/34 §10.4b): the prototype resolves a shower's parent
            // via map_vertex_in_shower FIRST (NeutrinoID.cxx:1655/1680/1720);
            // with the map fully populated by F3a, test the parent shower
            // before the incoming track segment.  The self test is defensive
            // only -- the F3a/root-branch write guards already exclude a
            // shower's own start vertex.
            auto ps_it = cfg.pf_shower_parent_precedence
                       ? vtx_to_parent_shower.find(start_vtx) : vtx_to_parent_shower.end();
            if (ps_it != vtx_to_parent_shower.end() && ps_it->second && ps_it->second != shower &&
                !sat_dropped(ps_it->second)) {   // pr/92: a dropped parent never renders; fall through
                PR::ShowerPtr parent_shower = ps_it->second;
                if (flag_print) {
                    std::cout << "[fill_bee_pf_tree] SHOWER-attached shower (parent-shower precedence)"
                              << "  conn_type=" << conn_type
                              << "  parent_shower_pdg=" << parent_shower->get_particle_type()
                              << "  pdg=" << shower->get_particle_type()
                              << "  ke=" << shower->get_kine_best() / units::MeV << " MeV"
                              << "\n";
                }
                auto& mp = direct ? shower_direct_showers : shower_indirect_showers;
                mp[parent_shower].push_back({shower, start_vtx});
                continue;   // resolved; skip the legacy incoming-segment path
            }
            auto it = vtx_incoming_seg.find(start_vtx);
            if (it != vtx_incoming_seg.end()) {
                if (flag_print) {
                    auto start_seg = shower->start_segment();
                    const auto* cl = start_seg ? start_seg->cluster() : nullptr;
                    WireCell::Point sp = shower->get_start_point();
                    WireCell::Point ep = shower->get_end_point();
                    const int parent_seg_id = (it->second->id() >= 0)
                                                ? it->second->id()
                                                : static_cast<int>(it->second->get_graph_index());
                    std::cout << "[fill_bee_pf_tree] SEGMENT-attached shower"
                              << "  conn_type=" << conn_type
                              << "  parent_seg=" << parent_seg_id
                              << "  pdg=" << shower->get_particle_type()
                              << "  ke=" << shower->get_kine_best() / units::MeV << " MeV"
                              << "  cluster=" << (cl ? std::to_string(cl->get_cluster_id()) : "?")
                              << "  start=(" << sp.x()/units::cm << "," << sp.y()/units::cm << "," << sp.z()/units::cm << ") cm"
                              << " end=(" << ep.x()/units::cm << "," << ep.y()/units::cm << "," << ep.z()/units::cm << ") cm"
                              << " nsegments=" << shower->get_num_segments() 
                              << "\n";
                }
                auto& mp = direct ? seg_direct_showers : seg_indirect_showers;
                mp[it->second].push_back({shower, start_vtx});
            } else {
                // start_vtx not reachable from main_vertex → check for parent shower
                if (root_reachable_vtxs.count(start_vtx)) {
                    // start_vtx is inside a root-level shower → attach to that parent shower
                    auto parent_shower_it = vtx_to_parent_shower.find(start_vtx);
                    if (parent_shower_it != vtx_to_parent_shower.end() &&
                        !sat_dropped(parent_shower_it->second)) {   // pr/92: dropped parent -> root fallback
                        PR::ShowerPtr parent_shower = parent_shower_it->second;
                        if (flag_print) {
                            auto start_seg = shower->start_segment();
                            const auto* cl = start_seg ? start_seg->cluster() : nullptr;
                            WireCell::Point sp = shower->get_start_point();
                            WireCell::Point ep = shower->get_end_point();
                            std::cout << "[fill_bee_pf_tree] SHOWER-attached shower (via parent shower vtx)"
                                      << "  conn_type=" << conn_type
                                      << "  parent_shower_pdg=" << parent_shower->get_particle_type()
                                      << "  pdg=" << shower->get_particle_type()
                                      << "  ke=" << shower->get_kine_best() / units::MeV << " MeV"
                                      << "  cluster=" << (cl ? std::to_string(cl->get_cluster_id()) : "?")
                                      << "  start=(" << sp.x()/units::cm << "," << sp.y()/units::cm << "," << sp.z()/units::cm << ") cm"
                                      << "  end=(" << ep.x()/units::cm << "," << ep.y()/units::cm << "," << ep.z()/units::cm << ") cm"
                                        << " nsegments=" << shower->get_num_segments()
                                      << "\n";
                        }
                        auto& mp = direct ? shower_direct_showers : shower_indirect_showers;
                        mp[parent_shower].push_back({shower, start_vtx});
                    } else {
                        if (flag_print) {
                            auto start_seg = shower->start_segment();
                            const auto* cl = start_seg ? start_seg->cluster() : nullptr;
                            WireCell::Point sp = shower->get_start_point();
                            WireCell::Point ep = shower->get_end_point();
                            std::cout << "[fill_bee_pf_tree] ROOT shower (via root-reachable shower vtx, no parent found)"
                                      << "  conn_type=" << conn_type
                                      << "  pdg=" << shower->get_particle_type()
                                      << "  ke=" << shower->get_kine_best() / units::MeV << " MeV"
                                      << "  cluster=" << (cl ? std::to_string(cl->get_cluster_id()) : "?")
                                      << "  start=(" << sp.x()/units::cm << "," << sp.y()/units::cm << "," << sp.z()/units::cm << ") cm"
                                      << "  end=(" << ep.x()/units::cm << "," << ep.y()/units::cm << "," << ep.z()/units::cm << ") cm"
                                        << " nsegments=" << shower->get_num_segments()
                                      << "\n";
                        }
                        auto& vec = direct ? root_direct_showers : root_indirect_showers;
                        vec.push_back({shower, main_vertex});
                    }
                } else {
                    // start_vtx truly isolated from main_vertex → fallback to root
                    if (flag_print) {
                        auto start_seg = shower->start_segment();
                        const auto* cl = start_seg ? start_seg->cluster() : nullptr;
                        WireCell::Point sp = shower->get_start_point();
                        WireCell::Point ep = shower->get_end_point();

                        std::cout << "[fill_bee_pf_tree] ROOT shower (fallback: start_vtx not in BFS tree)"
                                  << "  conn_type=" << conn_type
                                  << "  pdg=" << shower->get_particle_type()
                                  << "  ke=" << shower->get_kine_best() / units::MeV << " MeV"
                                  << "  cluster=" << (cl ? std::to_string(cl->get_cluster_id()) : "?")
                                  << "  start=(" << sp.x()/units::cm << "," << sp.y()/units::cm << "," << sp.z()/units::cm << ") cm"
                                    << "  end=(" << ep.x()/units::cm << "," << ep.y()/units::cm << "," << ep.z()/units::cm << ") cm"
                                        << " nsegments=" << shower->get_num_segments()
                                  << "\n";
                    }
                    auto& vec = direct ? root_direct_showers : root_indirect_showers;
                    // doc pr/84 r2 F2 (= pr/84 P3): with the knob on, anchor the
                    // pseudo carrier at the MAIN vertex so a remote association
                    // draws its real gap.  Legacy passes the shower's own start
                    // vertex, which collapses the carrier to zero length
                    // (gstart==gend in append_pseudo_shower) no matter how far
                    // away the shower is.
                    vec.push_back({shower, cfg.pf_pseudo_gap_from_main ? main_vertex : start_vtx});
                }
            }
        }
    }

    // --- Helpers ---
    auto get_vtx_pt = [](PR::VertexPtr v) -> WireCell::Point {
        return v->fit().valid() ? v->fit().point : v->wcpt().point;
    };

    // ID following prototype convention: cluster_id * 1000 + seg_id
    auto seg_display_id = [](PR::SegmentPtr seg) -> int {
        int sid = seg->id();
        if (sid < 0) sid = static_cast<int>(seg->get_graph_index());
        const auto* cl = seg->cluster();
        return cl ? cl->get_cluster_id() * 1000 + sid : sid;
    };

    // doc pr/38 Round 4 (pf_orphan_track_parentage): graph-faithful parentage
    // for barrier-orphaned track segments.  The flat safety net below emits
    // every BFS-unreached segment as a parentless, childless root -- even when
    // the graph chains it off a claimed track (pi+ -> proton at a shared
    // vertex) or off a shower's interior vertex (a guard-excluded muon
    // continuing an EM arm; SBND 18255 evt 142421 segs 7011/7012/7018).  The
    // stage-1 fixed point above already computed the correct parent for those
    // vertices (vtx_incoming_seg / vtx_to_parent_shower) and threw it away.
    // When on: TRACK anchor first (either endpoint carries a claimed incoming
    // segment -> insert into seg_parent/seg_children so build_seg_node's
    // recursion renders it), else SHOWER anchor (endpoint inside a shower's
    // view -> child of that shower's displayed leaf via shower_child_segs;
    // deliberately NOT put in seg_parent, whose null entries the root loop
    // emits).  An anchored orphan extends vtx_incoming_seg at its far vertex
    // so orphan-of-orphan chains anchor in a later round of the fixed point.
    // Orphans with no anchor at all fall to the flat net exactly as before.
    // This state is UNREACHABLE in the prototype (no shower_absorb_track
    // guard there -- such tracks are absorbed into the shower), so this is a
    // designed divergence; see porting_dictionary.  Off => pass skipped,
    // byte-identical legacy output.
    if (cfg.pf_shower_vertex_barrier && cfg.pf_orphan_track_parentage) {
        // P1: orphan pool -- IDENTICAL selection to the flat safety net below.
        std::vector<PR::SegmentPtr> orphan_pool;
        for (auto edesc : mir(boost::edges(*pr_graph))) {
            auto seg = (*pr_graph)[edesc].segment;
            if (!seg || used_segs.count(seg) || conn4_skip_segs.count(seg)) continue;
            if (!same_cluster(seg)) continue;      // prototype NeutrinoID.cxx:1488
            if (seg->dirsign() == 0) continue;     // prototype NeutrinoID.cxx:1215
            if (seg->fits().empty()) continue;
            orphan_pool.push_back(seg);
        }
        // Deterministic order: display id with graph-index tie-break (distinct
        // segments can compare EQUAL under SegmentIndexCmp, PRGraph.h:315).
        std::sort(orphan_pool.begin(), orphan_pool.end(),
                  [&](const PR::SegmentPtr& a, const PR::SegmentPtr& b) {
                      const int ida = seg_display_id(a), idb = seg_display_id(b);
                      if (ida != idb) return ida < idb;
                      return a->get_graph_index() < b->get_graph_index();
                  });
        // P2: fixed-point anchoring over the sorted pool.
        PR::IndexedSegmentSet claimed;
        auto anchor_common = [&](const PR::SegmentPtr& seg, PR::VertexPtr near_vtx,
                                 const std::string& parent_text) {
            auto far_vtx = PR::find_other_vertex(*pr_graph, seg, near_vtx);
            seg_endpoints[seg] = {near_vtx, far_vtx};
            used_segs.insert(seg);   // bars the flat net below from re-emitting
            claimed.insert(seg);
            if (far_vtx && !vtx_incoming_seg.count(far_vtx)) vtx_incoming_seg[far_vtx] = seg;
            if (flag_print) {
                std::cout << "[fill_bee_pf_tree] ANCHOR orphan seg=" << seg_display_id(seg)
                          << " -> parent=" << parent_text << "\n";
            }
        };
        bool changed = true;
        while (changed) {
            changed = false;
            for (const auto& seg : orphan_pool) {
                if (claimed.count(seg)) continue;
                auto [va, vb] = PR::find_vertices(*pr_graph, seg);
                PR::VertexPtr track_anchor = nullptr;
                if (va && vtx_incoming_seg.count(va)) track_anchor = va;
                else if (vb && vtx_incoming_seg.count(vb)) track_anchor = vb;
                if (track_anchor) {
                    auto parent = vtx_incoming_seg.at(track_anchor);
                    seg_parent[seg] = parent;
                    seg_children[parent].push_back(seg);
                    anchor_common(seg, track_anchor, std::to_string(seg_display_id(parent)));
                    changed = true;
                    continue;
                }
                PR::VertexPtr shower_anchor = nullptr;
                if (va && vtx_to_parent_shower.count(va)) shower_anchor = va;
                else if (vb && vtx_to_parent_shower.count(vb)) shower_anchor = vb;
                if (shower_anchor) {
                    auto psh = vtx_to_parent_shower.at(shower_anchor);
                    shower_child_segs[psh].push_back(seg);
                    auto pss = psh->start_segment();
                    anchor_common(seg, shower_anchor,
                                  std::string("shower:") + (pss ? std::to_string(seg_display_id(pss)) : "?"));
                    changed = true;
                    continue;
                }
            }
        }
    }

    int next_id = 1;  // fallback counter for nodes without a natural ID

    // doc pr/84 round 3 G1 (pf_unique_node_ids): jsTree keys its model by node
    // id, so a repeated id is invalid input -- see the knob's docstring.  The
    // id addresses nothing outside the tree (bee3 mc.js draws from
    // data.start/data.end), so a collision is resolved by re-issuing from a
    // range above the natural `cluster_id*1000 + seg_id` space.  Every firing
    // is logged: with shower_dedup_start_seg on there should be none.
    // doc pr/94 Phase 4: when the caller supplies the set, ids stay unique
    // ACROSS bundles -- a per-call set would restart the reissue range at
    // 1000000 for every bundle and collide by construction.
    std::set<int> local_used_node_ids;
    std::set<int>& used_node_ids = shared_used_ids ? *shared_used_ids : local_used_node_ids;

    auto make_node = [&](int id,
                         const std::string& text,
                         const WireCell::Point& start,
                         const WireCell::Point& end) -> Configuration {
        if (cfg.pf_unique_node_ids && !used_node_ids.insert(id).second) {
            int fresh = 1000000;
            while (!used_node_ids.insert(fresh).second) ++fresh;
            SPDLOG_LOGGER_DEBUG(log, "pr84 pf_id_collision: id={} reissued={} text='{}'",
                                id, fresh, text);
            id = fresh;
        }
        Configuration node;
        node["id"]   = id;
        node["text"] = text;
        Configuration dj;
        dj["start"][0] = start.x() / units::cm;
        dj["start"][1] = start.y() / units::cm;
        dj["start"][2] = start.z() / units::cm;
        dj["end"][0]   = end.x() / units::cm;
        dj["end"][1]   = end.y() / units::cm;
        dj["end"][2]   = end.z() / units::cm;
        node["data"] = dj;
        node["children"] = Json::arrayValue;
        return node;
    };

    auto format_mev = [&cfg](double energy) -> std::string {
        // prototype_names: integer MeV like the prototype's
        // WCReader::MCJSON ("int e = KE(...)*1000").  Legacy: "%.2f".
        if (cfg.prototype_names) {
            return std::to_string(static_cast<int>(energy / units::MeV));
        }
        char buf[32];
        std::snprintf(buf, sizeof(buf), "%.2f", energy / units::MeV);
        return std::string(buf);
    };

    // KeepMC-style pruning (prototype bee/WCReader.cc KeepMC): drop a LEAF
    // node whose kinetic energy is below the configured floor.  Nodes with
    // surviving children are always kept so the flow hierarchy stays intact.
    // Thresholds default 0 => keep everything (byte-identical legacy output).
    auto keep_node = [&cfg](int pdg, double ke, const Configuration& node) -> bool {
        if (cfg.em_ke_min <= 0.0 && cfg.np_ke_min <= 0.0) return true;
        if (!node["children"].empty()) return true;
        const int apdg = std::abs(pdg);
        if (apdg == 22 || apdg == 11) return ke >= cfg.em_ke_min;
        if (apdg == 2112 || apdg == 2212 || apdg > 1000000000) return ke >= cfg.np_ke_min;
        return true;
    };

    // Forward declare as std::function to handle mutual recursion with make_shower_leaf
    std::function<void(Configuration&, const std::vector<std::pair<PR::ShowerPtr,PR::VertexPtr>>&,
                       const std::vector<std::pair<PR::ShowerPtr,PR::VertexPtr>>&, PR::VertexPtr)> append_showers;
    // Forward declare (assigned at its original definition site below) so
    // make_shower_leaf can render shower_child_segs orphan children through
    // it (pf_orphan_track_parentage).  Declaration mechanics only: the first
    // call happens at assembly time, well after the assignment.
    std::function<Configuration(PR::SegmentPtr)> build_seg_node;

    auto make_shower_leaf = [&](PR::ShowerPtr shower) -> Configuration {
        const int pdg = shower->get_particle_type();
        const std::string ke = format_mev(shower->get_kine_best());
        auto [svtx, sconn] = shower->get_start_vertex_and_type();
        auto start_seg = shower->start_segment();
        int id = start_seg ? seg_display_id(start_seg) : (next_id++);
        auto node = make_node(id,
                              pf_pdg_to_name(pdg, cfg.prototype_names, cfg.pf_pdg_name_prototype_fallback) + "  " + ke + " MeV",
                              shower->get_start_point(), shower->get_end_point());
        if (flag_print) {
            const auto* cl = start_seg ? start_seg->cluster() : nullptr;
            std::cout << "[fill_bee_pf_tree] ADD shower-leaf"
                      << "  id=" << id
                      << "  shower_id=" << shower->get_shower_id()
                      << "  pdg=" << pdg
                      << "  conn_type=" << sconn
                      << "  ke=" << ke << " MeV"
                      << "  cluster=" << (cl ? std::to_string(cl->get_cluster_id()) : "?")
                      << "  has_start_vtx=" << (svtx ? 1 : 0)
                      << "\n";
        }
        
        // Append any showers attached to this shower
        static const std::vector<std::pair<PR::ShowerPtr,PR::VertexPtr>> empty_showers;
        auto di_it  = shower_direct_showers.find(shower);
        auto ind_it = shower_indirect_showers.find(shower);
        append_showers(node["children"],
                       di_it  != shower_direct_showers.end()   ? di_it->second  : empty_showers,
                       ind_it != shower_indirect_showers.end() ? ind_it->second : empty_showers,
                       svtx);

        // doc pr/38 Round 4 (pf_orphan_track_parentage): orphan TRACK
        // segments anchored inside this shower's view render as its children
        // (single funnel point -- covers direct leaves, pseudo-gamma wrappers
        // and pi0-grouped leaves alike).  Same KeepMC convention as
        // build_seg_node's child recursion; each child renders its own
        // seg_children chain.  Map empty when the knob is off.
        if (cfg.pf_orphan_track_parentage) {
            auto cs_it = shower_child_segs.find(shower);
            if (cs_it != shower_child_segs.end()) {
                for (const auto& child : cs_it->second) {
                    auto child_node = build_seg_node(child);
                    const int    cpdg = child->has_particle_info() ? child->particle_info()->pdg() : 0;
                    const double cke  = child->has_particle_info() ? child->particle_info()->kinetic_energy() : 0.0;
                    if (!keep_node(cpdg, cke, child_node)) continue;
                    node["children"].append(child_node);
                }
            }
        }

        if (node["children"].empty()) node["icon"] = "jstree-file";
        return node;
    };

    // Helper: insert one pseudo-particle node + shower-leaf into a parent children array.
    // Mirrors prototype fill_psuedo_reco_tree: pseudo PDG is gamma (22) for EM showers,
    // neutron (2112) for all others (e.g. isolated proton activities that were not
    // absorbed into an EMShower — the neutral carrier is assumed to be an unseen neutron).
    auto append_pseudo_shower = [&](Configuration& parent_children, PR::ShowerPtr sh, PR::VertexPtr conn_vtx) {
        const int pdg = (std::abs(sh->get_particle_type()) == 11 ||
                         std::abs(sh->get_particle_type()) == 22) ? 22 : 2112;
        const std::string pname = pf_pdg_to_name(pdg, cfg.prototype_names, cfg.pf_pdg_name_prototype_fallback);
        PR::VertexPtr cv = conn_vtx;
        WireCell::Point gstart = cv ? get_vtx_pt(cv) : sh->get_start_point();
        WireCell::Point gend   = sh->get_start_point();
        const std::string pseudo_ke = format_mev(sh->get_kine_best());
        const int pseudo_id = next_id++;
        auto pseudo = make_node(pseudo_id, pname + "  " + pseudo_ke + " MeV", gstart, gend);
        if (flag_print) {
            std::cout << "[fill_bee_pf_tree] ADD pseudo-" << pname
                      << "  id=" << pseudo_id
                      << "  ke=" << pseudo_ke << " MeV"
                      << "  child_pdg=" << sh->get_particle_type()
                      << "\n";
        }
        auto leaf = make_shower_leaf(sh);
        // KeepMC pruning: a pseudo carrier with its only shower dropped is
        // itself dropped (nothing left to carry).
        if (!keep_node(sh->get_particle_type(), sh->get_kine_best(), leaf)) return;
        pseudo["children"].append(leaf);
        if (pseudo["children"].empty()) pseudo["icon"] = "jstree-file";
        parent_children.append(pseudo);
    };

    // doc pr/84 r2 F1 (pf_direct_when_touching): a conn-2/3 shower whose
    // fitted charge comes within pf_touch_max of the main vertex is a graph
    // artifact ("the BFS could not walk there"), not a neutral daughter --
    // render it as a direct leaf.  Distance deliberately excludes the pr/84
    // sec 4 remote-association population (min 4.91 cm > 3 cm default),
    // which must KEEP its carrier (see F2).  pi0 daughters never reach this
    // test -- their carrier is correct.  doc 77 round 1 (2026-08-24): rung 2
    // (pf_touch_cross_main/_max) removed -- zero movers, F1.0 probe failure
    // (pr/84 sec 607/622).  See sbnd_xin/docs/77_knob-ledger.tsv.
    auto effectively_touching = [&](PR::ShowerPtr sh) -> bool {
        if (!cfg.pf_direct_when_touching || !main_vertex) return false;
        const int conn = sh->get_start_vertex_and_type().second;
        if (conn != 2 && conn != 3) return false;
        const double d_fit = shower_get_closest_point(*sh, get_vtx_pt(main_vertex), "fit").first;
        if (d_fit < 0) return false;  // no fit cloud: fail safe, keep the carrier
        return d_fit <= cfg.pf_touch_max;
    };

    // Append all showers (direct + indirect via pseudo-gamma) into a children array,
    // given the connection vertex for the indirect case.
    // F4 (doc pr/34 §10.5): one pi0 node per pi0 id.  The prototype memoizes
    // the pi0 node on a member (map_pio_id_saved_pair, NeutrinoID.cxx:1326/
    // :1361), so a pi0 whose daughters hang under two parents still renders
    // once; the toolkit's invocation-local grouping renders it once PER
    // PARENT.  A jsTree node has exactly one parent, so the merged node needs
    // a single home: the HIGHEST-ENERGY daughter's parent (owner decision
    // 2026-08-04 -- deliberately NOT the prototype's first-writer-wins).
    std::map<int, std::vector<std::pair<PR::ShowerPtr, PR::VertexPtr>>> pi0_all_groups;
    std::map<int, PR::ShowerPtr> pi0_home_daughter;
    if (cfg.pf_pi0_node_per_id) {
        auto collect = [&](const std::vector<std::pair<PR::ShowerPtr, PR::VertexPtr>>& vec) {
            for (const auto& pr : vec)
                if (pi0_showers.count(pr.first))
                    pi0_all_groups[map_shower_pio_id.at(pr.first)].push_back(pr);
        };
        collect(root_direct_showers);
        collect(root_indirect_showers);
        for (const auto& [seg, vec] : seg_direct_showers) collect(vec);
        for (const auto& [seg, vec] : seg_indirect_showers) collect(vec);
        for (const auto& [sh, vec] : shower_direct_showers) collect(vec);
        for (const auto& [sh, vec] : shower_indirect_showers) collect(vec);
        for (const auto& [pi0_id, group] : pi0_all_groups) {
            PR::ShowerPtr best = nullptr;
            for (const auto& [sh, cv] : group)
                if (!best || sh->get_kine_best() > best->get_kine_best()) best = sh;
            pi0_home_daughter[pi0_id] = best;   // ties: first in collection order
        }
    }

    // Pi0 showers are grouped by pi0_id and rendered as: pi0 node → gamma → shower_leaf.
    append_showers = [&](Configuration& children,
                               const std::vector<std::pair<PR::ShowerPtr,PR::VertexPtr>>& direct,
                               const std::vector<std::pair<PR::ShowerPtr,PR::VertexPtr>>& indirect,
                               PR::VertexPtr fallback_conn_vtx) {
        // --- Non-pi0 direct showers ---
        for (auto& [sh, _] : direct) {
            if (pi0_showers.count(sh)) continue;
            auto leaf = make_shower_leaf(sh);
            if (!keep_node(sh->get_particle_type(), sh->get_kine_best(), leaf)) continue;
            children.append(leaf);
        }
        // --- Non-pi0 indirect showers (pseudo-gamma) ---
        for (auto& [sh, conn_vtx] : indirect) {
            if (pi0_showers.count(sh)) continue;
            // doc pr/84 r2 F1: vertex-touching shower renders directly, in the
            // same children array the pseudo carrier would have landed in --
            // parent linkage and sibling order otherwise unchanged.
            if (effectively_touching(sh)) {
                auto leaf = make_shower_leaf(sh);
                if (!keep_node(sh->get_particle_type(), sh->get_kine_best(), leaf)) continue;
                if (flag_print) {
                    std::cout << "[fill_bee_pf_tree] SUPPRESS pseudo (pr84 touching)"
                              << "  pdg=" << sh->get_particle_type()
                              << "  ke=" << sh->get_kine_best() / units::MeV << " MeV"
                              << "  conn_type=" << sh->get_start_vertex_and_type().second
                              << "\n";
                }
                children.append(leaf);
                continue;
            }
            PR::VertexPtr cv = conn_vtx ? conn_vtx : fallback_conn_vtx;
            append_pseudo_shower(children, sh, cv);
        }

        // --- Pi0 showers: group by pi0_id, emit one pi0 node per group ---
        // Collect from both direct and indirect lists.
        std::map<int, std::vector<std::pair<PR::ShowerPtr, PR::VertexPtr>>> pi0_groups;
        for (auto& [sh, cv] : direct)   if (pi0_showers.count(sh)) pi0_groups[map_shower_pio_id.at(sh)].push_back({sh, cv});
        for (auto& [sh, cv] : indirect) if (pi0_showers.count(sh)) pi0_groups[map_shower_pio_id.at(sh)].push_back({sh, cv});

        for (auto& [pi0_id, local_group] : pi0_groups) {
            // F4: only the invocation at the home parent (the highest-energy
            // daughter's) emits the pi0 node, and it carries ALL daughters;
            // other parents skip theirs.
            const bool merged = cfg.pf_pi0_node_per_id && pi0_all_groups.count(pi0_id);
            if (merged) {
                const auto& home = pi0_home_daughter.at(pi0_id);
                bool has_home = false;
                for (auto& [sh, cv] : local_group) if (sh == home) { has_home = true; break; }
                if (!has_home) continue;
            }
            const auto& group = merged ? pi0_all_groups.at(pi0_id) : local_group;
            auto mass_it = map_pio_id_mass.find(pi0_id);
            const double pi0_ke = (mass_it != map_pio_id_mass.end()) ? mass_it->second.first : 0.0;
            // Pi0 sits at the connection vertex (point particle: start == end)
            PR::VertexPtr conn_vtx = group[0].second ? group[0].second : fallback_conn_vtx;
            if (merged) {
                // seat the merged node at its home daughter's connection vertex
                for (const auto& [sh, cv] : group)
                    if (sh == pi0_home_daughter.at(pi0_id)) { conn_vtx = cv ? cv : fallback_conn_vtx; break; }
            }
            WireCell::Point pi0_pt = conn_vtx ? get_vtx_pt(conn_vtx) : WireCell::Point(0,0,0);
            const int pi0_node_id = next_id++;
            auto pi0_node = make_node(pi0_node_id, "pi0  " + format_mev(pi0_ke) + " MeV", pi0_pt, pi0_pt);
            if (flag_print) {
                std::cout << "[fill_bee_pf_tree] ADD pi0-node"
                          << "  id=" << pi0_node_id
                          << "  pi0_id=" << pi0_id
                          << "  ke=" << format_mev(pi0_ke) << " MeV"
                          << "  n_showers=" << group.size()
                          << "\n";
            }
            for (auto& [sh, cv] : group) {
                PR::VertexPtr gcv = cv ? cv : fallback_conn_vtx;
                append_pseudo_shower(pi0_node["children"], sh, gcv);
            }
            // KeepMC pruning: skip a pi0 whose daughter showers were all dropped.
            if (pi0_node["children"].empty() &&
                (cfg.em_ke_min > 0.0 || cfg.np_ke_min > 0.0)) continue;
            if (pi0_node["children"].empty()) pi0_node["icon"] = "jstree-file";
            children.append(pi0_node);
        }
    };

    // Build the full JSON subtree for a track segment (recursive).
    // (declared as std::function above, next to append_showers)
    build_seg_node =
        [&](PR::SegmentPtr seg) -> Configuration {

        // F5 (doc pr/34 §10.6): the prototype's no-PID node text is PDGName(0)
        // -- TDatabasePDG does not know 0, so the literal "0" -- with the KE of
        // an all-zero 4-vector: "0  0 MeV", not "particle  0 MeV".
        std::string pname = cfg.pf_pdg_name_prototype_fallback ? "0" : "particle";
        std::string ke_str = cfg.prototype_names ? "0" : "0.00";
        if (seg->has_particle_info()) {
            auto pi = seg->particle_info();
            // prototype_names: TDatabasePDG-style short names ("mu-",
            // "proton") like the prototype mc.json; legacy: the
            // ParticleDataSet long name ("muon", "proton").
            pname  = cfg.prototype_names ? pf_pdg_to_name(pi->pdg(), true, cfg.pf_pdg_name_prototype_fallback) : pi->name();
            ke_str = format_mev(pi->kinetic_energy());
        }

        auto ep_it = seg_endpoints.find(seg);
        WireCell::Point start_pt, end_pt;
        PR::VertexPtr far_vtx = nullptr;
        if (ep_it != seg_endpoints.end()) {
            start_pt = ep_it->second.first  ? get_vtx_pt(ep_it->second.first)  : WireCell::Point(0,0,0);
            end_pt   = ep_it->second.second ? get_vtx_pt(ep_it->second.second) : start_pt;
            far_vtx  = ep_it->second.second;
        }

        if (flag_print) {
            const auto* seg_cluster = seg->cluster();
            // Ident comparison, matching F1's guard (§10.2): the old pointer
            // test was stricter and over-reported out-of-main-cluster segments.
            const bool in_main_cluster = (main_cluster && seg_cluster &&
                seg_cluster->get_cluster_id() == main_cluster->get_cluster_id());
            const bool is_shower_seg = (seg_to_shower.count(seg) > 0);
            auto pit = seg_parent.find(seg);
            auto parent_seg_dbg = (pit != seg_parent.end()) ? pit->second : nullptr;
            const std::string parent_text = parent_seg_dbg ? std::to_string(seg_display_id(parent_seg_dbg)) : "ROOT";
            std::cout << "[fill_bee_pf_tree] ADD track-node"
                      << "  seg=" << seg_display_id(seg)
                      << "  parent=" << parent_text
                      << "  name=" << pname
                      << "  ke=" << ke_str << " MeV"
                      << "  cluster=" << (seg_cluster ? std::to_string(seg_cluster->get_cluster_id()) : "?")
                      << "  in_main_cluster=" << (in_main_cluster ? 1 : 0)
                      << "  is_shower_seg=" << (is_shower_seg ? 1 : 0)
                      << "  start=(" << start_pt.x()/units::cm << "," << start_pt.y()/units::cm << "," << start_pt.z()/units::cm << ") cm"
                      << "  end=(" << end_pt.x()/units::cm << "," << end_pt.y()/units::cm << "," << end_pt.z()/units::cm << ") cm"
                      << "\n";
        }

        auto node = make_node(seg_display_id(seg),
                              pname + "  " + ke_str + " MeV",
                              start_pt, end_pt);

        // Attach showers connected at the far-end vertex of this segment
        static const std::vector<std::pair<PR::ShowerPtr,PR::VertexPtr>> empty_showers;
        auto di_it  = seg_direct_showers.find(seg);
        auto ind_it = seg_indirect_showers.find(seg);
        append_showers(node["children"],
                       di_it  != seg_direct_showers.end()   ? di_it->second  : empty_showers,
                       ind_it != seg_indirect_showers.end() ? ind_it->second : empty_showers,
                       far_vtx);

        // Recurse into child track segments
        auto ch_it = seg_children.find(seg);
        if (ch_it != seg_children.end()) {
            for (const auto& child : ch_it->second) {
                if (conn4_skip_segs.count(child)) continue;
                auto child_node = build_seg_node(child);
                const int cpdg = child->has_particle_info() ? child->particle_info()->pdg() : 0;
                const double cke = child->has_particle_info() ? child->particle_info()->kinetic_energy() : 0.0;
                if (!keep_node(cpdg, cke, child_node)) continue;
                node["children"].append(child_node);
            }
        }

        if (node["children"].empty()) node["icon"] = "jstree-file";
        return node;
    };

    // --- Assemble top-level particles array ---
    Configuration particles = Json::arrayValue;

    // Showers attached directly to the neutrino vertex
    append_showers(particles, root_direct_showers, root_indirect_showers, main_vertex);

    // Root track segments (direct daughters of the neutrino vertex).
    // Disconnected segments (orphaned fragments unreachable from main_vertex)
    // are now skipped entirely to avoid adding zero-energy orphaned particles.
    for (auto& [seg, parent] : seg_parent) {
        if (parent != nullptr) continue;   // skip non-roots
        if (conn4_skip_segs.count(seg)) continue;
        auto root_node = build_seg_node(seg);
        const int rpdg = seg->has_particle_info() ? seg->particle_info()->pdg() : 0;
        const double rke = seg->has_particle_info() ? seg->particle_info()->kinetic_energy() : 0.0;
        if (!keep_node(rpdg, rke, root_node)) continue;
        particles.append(root_node);
    }

    // doc pr/38: orphan safety net, gated under pf_shower_vertex_barrier as
    // the no-silent-drop complement of the barrier.  The prototype's flat
    // loop gives EVERY non-shower main-cluster segment a node with
    // mc_mother=0 before the mother-assignment BFS runs
    // (NeutrinoID.cxx:1485-1489), so a segment the BFS never reaches stays in
    // its tree as a root-level node; the BFS-built tree above silently
    // dropped it.  Emitted as root-level leaves with the prototype's own node
    // conventions: endpoints from the segment's fit points oriented by
    // dirsign (fill_reco_tree, NeutrinoID.cxx:1217-1239), dirsign==0 segments
    // not plotted (:1215), KeepMC floors as for every other node.  Emission
    // order is sorted by encoded id -- the prototype's is pointer-map order,
    // which is not reproducible; a stable order is chosen deliberately.
    // doc pr/65 round 3 (rung 3): audit-only variant of the safety net below.
    // Keep the VISIBILITY the net was added for, drop the FABRICATION: name
    // every still-unclaimed segment in the log -- WITHOUT the display filters,
    // so dirsign==0 / empty-fit / KeepMC-floored segments (which today vanish
    // from the tree silently) are named too -- and append no node.  The
    // owner's requirement: an unclaimed orphan must never become a root-level
    // PF particle.  Load-bearing only with rung 1 (shower_absorb_unreachable_main)
    // claiming the fragments first.
    if (cfg.pf_shower_vertex_barrier && cfg.pf_orphan_audit_only) {
        std::vector<PR::SegmentPtr> unclaimed;
        for (auto edesc : mir(boost::edges(*pr_graph))) {
            auto seg = (*pr_graph)[edesc].segment;
            if (!seg || used_segs.count(seg) || conn4_skip_segs.count(seg)) continue;
            if (!same_cluster(seg)) continue;      // prototype NeutrinoID.cxx:1488
            unclaimed.push_back(seg);
        }
        std::sort(unclaimed.begin(), unclaimed.end(),
                  [&](const PR::SegmentPtr& a, const PR::SegmentPtr& b) {
                      return seg_display_id(a) < seg_display_id(b);
                  });
        int n_emitted = 0;
        for (const auto& seg : unclaimed) {
            const auto* cl = seg->cluster();
            const bool haspi = seg->has_particle_info() && seg->particle_info();
            SPDLOG_LOGGER_INFO(log,
                "pr65 pf-orphan-audit: unclaimed seg={} cluster={} pdg={} ke_mev={:.2f} nfits={} dirsign={}",
                seg_display_id(seg),
                cl ? std::to_string(cl->get_cluster_id()) : "?",
                haspi ? seg->particle_info()->pdg() : 0,
                haspi ? seg->particle_info()->kinetic_energy() / units::MeV : 0.0,
                seg->fits().size(), seg->dirsign());
            // doc pr/93 round 4 (pf_orphan_confident_track, pr/65 rung 4):
            // emit a root node for the narrow confident-long-straight-track
            // class only -- e.g. SBND 18255-315167's 150.7cm proton, freed
            // from shower membership by shower_cone_absorb_guard but graph-
            // disconnected from the main vertex.  Node construction mirrors
            // the (production-disabled) flat net below, including its
            // dirsign/fit display filters and the KeepMC floors.  Audit
            // lines above stay untouched for every segment.  Knob off =>
            // this block never runs => byte-identical.
            if (cfg.pf_orphan_confident_track &&
                seg->dirsign() != 0 && !seg->fits().empty() &&
                PR::segment_orphan_confident_track(seg, cfg.pf_orphan_track_min)) {
                auto pi = seg->particle_info();
                const std::string pname = cfg.prototype_names
                    ? pf_pdg_to_name(pi->pdg(), true, cfg.pf_pdg_name_prototype_fallback)
                    : pi->name();
                const std::string ke_str = format_mev(pi->kinetic_energy());
                const auto& fits = seg->fits();
                const WireCell::Point& p_front = fits.front().point;
                const WireCell::Point& p_back  = fits.back().point;
                const bool fwd = (seg->dirsign() == 1);
                auto node = make_node(seg_display_id(seg),
                                      pname + "  " + ke_str + " MeV",
                                      fwd ? p_front : p_back,
                                      fwd ? p_back : p_front);
                node["icon"] = "jstree-file";
                if (keep_node(pi->pdg(), pi->kinetic_energy(), node)) {
                    particles.append(node);
                    ++n_emitted;
                    SPDLOG_LOGGER_INFO(log,
                        "pr93 pf-orphan-confident-track: EMIT root seg={} cluster={} pdg={} "
                        "ke_mev={:.2f} len_cm={:.1f} dirsign={}",
                        seg_display_id(seg),
                        cl ? std::to_string(cl->get_cluster_id()) : "?",
                        pi->pdg(), pi->kinetic_energy() / units::MeV,
                        PR::segment_track_length(seg) / units::cm, seg->dirsign());
                }
            }
        }
        SPDLOG_LOGGER_INFO(log,
            "pr65 pf-orphan-audit: {} unclaimed segment(s), no PF node fabricated (pf_orphan_audit_only)",
            unclaimed.size());
        if (cfg.pf_orphan_confident_track && n_emitted) {
            SPDLOG_LOGGER_INFO(log,
                "pr93 pf-orphan-confident-track: {} of them emitted as root track node(s)", n_emitted);
        }
    }
    else if (cfg.pf_shower_vertex_barrier) {
        std::vector<PR::SegmentPtr> orphans;
        for (auto edesc : mir(boost::edges(*pr_graph))) {
            auto seg = (*pr_graph)[edesc].segment;
            if (!seg || used_segs.count(seg) || conn4_skip_segs.count(seg)) continue;
            if (!same_cluster(seg)) continue;      // prototype NeutrinoID.cxx:1488
            if (seg->dirsign() == 0) continue;     // prototype NeutrinoID.cxx:1215
            if (seg->fits().empty()) continue;
            orphans.push_back(seg);
        }
        std::sort(orphans.begin(), orphans.end(),
                  [&](const PR::SegmentPtr& a, const PR::SegmentPtr& b) {
                      return seg_display_id(a) < seg_display_id(b);
                  });
        for (const auto& seg : orphans) {
            std::string pname = cfg.pf_pdg_name_prototype_fallback ? "0" : "particle";
            std::string ke_str = cfg.prototype_names ? "0" : "0.00";
            int opdg = 0;
            double oke = 0.0;
            if (seg->has_particle_info()) {
                auto pi = seg->particle_info();
                pname  = cfg.prototype_names
                       ? pf_pdg_to_name(pi->pdg(), true, cfg.pf_pdg_name_prototype_fallback)
                       : pi->name();
                ke_str = format_mev(pi->kinetic_energy());
                opdg = pi->pdg();
                oke  = pi->kinetic_energy();
            }
            const auto& fits = seg->fits();
            const WireCell::Point& p_front = fits.front().point;
            const WireCell::Point& p_back  = fits.back().point;
            const bool fwd = (seg->dirsign() == 1);
            auto node = make_node(seg_display_id(seg),
                                  pname + "  " + ke_str + " MeV",
                                  fwd ? p_front : p_back,
                                  fwd ? p_back : p_front);
            node["icon"] = "jstree-file";
            if (!keep_node(opdg, oke, node)) continue;
            if (flag_print) {
                const auto* cl = seg->cluster();
                std::cout << "[fill_bee_pf_tree] ADD orphan-track-root"
                          << "  seg=" << seg_display_id(seg)
                          << "  name=" << pname
                          << "  ke=" << ke_str << " MeV"
                          << "  cluster=" << (cl ? std::to_string(cl->get_cluster_id()) : "?")
                          << "\n";
            }
            particles.append(node);
        }
    }

    // doc pr/123 round 2 (pf_orphan_guard_freed): a track the pass4
    // long-track guard declined (SegmentFlags::kPass4GuardFreed) that no
    // shower and no BFS claimed gets a root PF node -- SBND 18255-171572's
    // 125cm muon: correctly expelled from the EM shower, but cross-cluster
    // and score-100-sentinel PID, so both the pr/93 confident-track class
    // and the main-cluster audit scope miss it and it vanished from the PF
    // tree.  Node construction mirrors the pr/93 emission (dirsign/fit
    // display filters, KeepMC floors).  The flag is the predicate: nothing
    // outside the guard's own decline set is touched.  Knob off => this
    // block never runs => byte-identical.
    if (cfg.pf_orphan_guard_freed) {
        std::vector<PR::SegmentPtr> freed;
        for (auto edesc : mir(boost::edges(*pr_graph))) {
            auto seg = (*pr_graph)[edesc].segment;
            if (!seg || used_segs.count(seg) || conn4_skip_segs.count(seg)) continue;
            if (!seg->flags_any(PR::SegmentFlags::kPass4GuardFreed)) continue;
            if (seg->dirsign() == 0) continue;
            if (seg->fits().empty()) continue;
            if (!seg->has_particle_info() || !seg->particle_info()) continue;
            freed.push_back(seg);
        }
        std::sort(freed.begin(), freed.end(),
                  [&](const PR::SegmentPtr& a, const PR::SegmentPtr& b) {
                      return seg_display_id(a) < seg_display_id(b);
                  });
        for (const auto& seg : freed) {
            auto pi = seg->particle_info();
            const std::string pname = cfg.prototype_names
                ? pf_pdg_to_name(pi->pdg(), true, cfg.pf_pdg_name_prototype_fallback)
                : pi->name();
            const std::string ke_str = format_mev(pi->kinetic_energy());
            const auto& fits = seg->fits();
            const WireCell::Point& p_front = fits.front().point;
            const WireCell::Point& p_back  = fits.back().point;
            const bool fwd = (seg->dirsign() == 1);
            auto node = make_node(seg_display_id(seg),
                                  pname + "  " + ke_str + " MeV",
                                  fwd ? p_front : p_back,
                                  fwd ? p_back : p_front);
            node["icon"] = "jstree-file";
            if (!keep_node(pi->pdg(), pi->kinetic_energy(), node)) continue;
            particles.append(node);
            const auto* cl = seg->cluster();
            SPDLOG_LOGGER_INFO(log,
                "pr123 pf-orphan-guard-freed: EMIT root seg={} cluster={} pdg={} "
                "ke_mev={:.2f} len_cm={:.1f}",
                seg_display_id(seg),
                cl ? std::to_string(cl->get_cluster_id()) : "?",
                pi->pdg(), pi->kinetic_energy() / units::MeV,
                PR::segment_track_length(seg) / units::cm);
        }
    }

    if (out_particles) {
        // doc pr/94 Phase 4: accumulate.  set_particles() REPLACES the array
        // (Bee.cxx:549-551), so a second call would erase the first bundle.
        // Wrap this bundle's roots under one synthetic node -- the Bee "mc"
        // layer is already a bare JSON forest, so an extra root needs no
        // format change -- and let the caller set the concatenation once.
        const auto& mvp = main_vertex->fit().valid() ? main_vertex->fit().point
                                                     : main_vertex->wcpt().point;
        std::string label = "nu";
        {
            const auto& ti = tf->get_tagger_info();
            if (ti.nu_index >= 0) {
                label = "nu " + std::to_string(ti.nu_index)
                      + " (gid " + std::to_string(ti.matched_flash_gid)
                      + ", cluster " + std::to_string(ti.cluster_id) + ")";
            }
        }
        auto root = make_node(1000000 + static_cast<int>(out_particles->size()), label, mvp, mvp);
        root["children"] = particles;
        out_particles->append(root);
        SPDLOG_LOGGER_TRACE(log, "fill_bee_pf_tree '{}': bundle root '{}' with {} top-level particles",
                            cfg.name, label, particles.size());
    }
    else {
        tree.set_particles(particles);
        SPDLOG_LOGGER_TRACE(log, "fill_bee_pf_tree '{}': {} top-level particles",
                            cfg.name, particles.size());
    }
}


// Helper function to fill bee points from a single cluster
void MultiAlgBlobClustering::fill_bee_points_from_cluster(
    Bee::Points& bpts, const Cluster& cluster,
    const std::string& pcname, const std::vector<std::string>& coords, int filter,
    double dQdx_scale, double dQdx_offset)
{
    int clid = cluster.get_cluster_id(); //bpts.back_cluster_id() + 1;

    // std::cout << "Test: " << bpts.size() << " " << bpts.back_cluster_id() << " " <<  clid << std::endl;

    if (pcname == "stm_fit"){
        // STM fit-trajectory layer (TaggerCheckSTM save_stm_fit knob): the
        // per-point fitted dQ is encoded into q with the same
        // dQdx_scale/dQdx_offset convention as the PRGraph track_fit layer.
        // Reachable only when a bee_points_sets entry names this PC, so
        // existing configs are byte-identical.
        auto& fit_pc = cluster.get_pc(pcname);
        if (fit_pc.empty()) {
            return;
        }
        const auto& fx = fit_pc.get(coords.at(0))->elements<double>();
        const auto& fy = fit_pc.get(coords.at(1))->elements<double>();
        const auto& fz = fit_pc.get(coords.at(2))->elements<double>();
        const auto& fdQ = fit_pc.get("dQ")->elements<double>();
        for (size_t i = 0; i < fx.size(); ++i) {
            bpts.append(Point(fx[i], fy[i], fz[i]), fdQ[i]*dQdx_scale + dQdx_offset, clid, 0);
        }
        return;
    }

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

        // std::cout << "Steiner Test: " << x_coords.size() << " " << y_coords.size() << " " << z_coords.size() << std::endl;

         for (size_t i = 0; i < x_coords.size(); ++i) {
            // Create point from steiner point cloud
            Point vtx(x_coords[i], y_coords[i], z_coords[i]);

            // Get the point index from the default scope
            auto point_index = cluster.get_closest_point_index(vtx);
            
            auto charge_result = cluster.calc_charge_wcp(point_index, 4000, true);
            double point_charge = charge_result.second; // Extract the charge value from the pair

            if (flag_steiner_terminal[i]) {
                bpts.append(Point(x_coords[i], y_coords[i], z_coords[i]), point_charge, clid, 1);  // terminals  ... 
            }else{
                bpts.append(Point(x_coords[i], y_coords[i], z_coords[i]), point_charge, clid, 0); // non-terminals ...
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

            // Original (pre-flash-merge) per-blob cluster ids, if present (written
            // by examine_bundles' flash-time merge into "real_cluster_id"/"perblob").
            // One id per blob, in the same blob order as spcs.  When absent (any
            // non-flash-merged cluster / other layers / detectors) the Bee
            // real_cluster_id stays == clid, so output is backward-compatible.
            std::vector<int> orig_ids;
            if (cluster.has_pcarray<int>("real_cluster_id", "perblob")) {
                auto sp = cluster.get_pcarray<int>("real_cluster_id", "perblob");
                orig_ids.assign(sp.begin(), sp.end());
            }
            const bool use_orig = (orig_ids.size() == spcs.size());

            // Dead-channel threshold: uncertainty > 1e10 flags a dead wire
            // (matches PointTreeBuilding m_dead_threshold and Facade_Cluster::is_wire_dead).
            const double dead_threshold = 1e10;

            // For each scoped point cloud (one per blob), compute per-point wire charge.
            // Each point gets q = mean(Q_U, Q_V, Q_W) over the non-dead planes for the
            // specific wires that define this 3D point — matching the prototype formula.
            // The global cluster point index is obtained via a KD-tree nearest-neighbour
            // lookup (exact match: every scoped-view point is also in the cluster PC).
            for (size_t spc_idx = 0; spc_idx < spcs.size(); ++spc_idx) {
                const auto& spc = spcs[spc_idx];
                auto x = spc.get().get(coords[0])->elements<double>();
                auto y = spc.get().get(coords[1])->elements<double>();
                auto z = spc.get().get(coords[2])->elements<double>();

                // Per-blob original cluster id (falls back to the merged clid).
                const int real_clid = use_orig ? orig_ids[spc_idx] : clid;

                const size_t size = x.size();
                for (size_t ind = 0; ind < size; ++ind) {
                    // Resolve global point index via spatial lookup.
                    // charge_value() caches all per-plane charge vectors on first call,
                    // so subsequent lookups are O(1) vector reads.
                    const WireCell::Point pt(x[ind], y[ind], z[ind]);
                    const size_t pt_idx = cluster.get_closest_point_index(pt);

                    // Per-plane charge mean (prototype formula), excluding dead planes.
                    double sum = 0.0;
                    int nplanes = 0;
                    for (int plane = 0; plane < 3; ++plane) {
                        if (!cluster.is_wire_dead(pt_idx, plane, dead_threshold)) {
                            sum += cluster.charge_value(pt_idx, plane);
                            ++nplanes;
                        }
                    }
                    const double point_charge = (nplanes > 0) ? sum / nplanes : 0.0;

                    bpts.append(pt, point_charge, clid, real_clid);
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
        SPDLOG_LOGGER_WARN(log, "Cannot access node for cluster");
        return;
    }
    
    // Iterate through child nodes (blobs)
    for (const auto* bnode : cluster_node->children()) {
        auto wpid = bnode->value.facade<Blob>()->wpid();
        int apa = wpid.apa();
        int face = wpid.face();


        // Select the destination patches: a drift-side group bucket when dead
        // grouping is enabled, otherwise the per-(apa,face) bucket.
        Bee::Patches* patches = nullptr;
        if (!m_dead_apa_groups.empty()) {
            for (const auto& ag : m_dead_apa_groups) {
                if (ag.apas.count(apa)) {
                    auto it = m_bee_dead_groups.find(ag.name);
                    if (it != m_bee_dead_groups.end()) patches = &it->second;
                    break;
                }
            }
        } else {
            auto it_apa = m_bee_dead_patches.find(apa);
            if (it_apa != m_bee_dead_patches.end()) {
                auto it_face = it_apa->second.find(face);
                if (it_face != it_apa->second.end()) patches = &it_face->second;
            }
        }
        if (!patches) continue;

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
        patches->append(y.begin(), y.end(), z.begin(), z.end());
    }
}


// Group the root "opflash" flashes (both TPC sides) by their time with the
// given window and stash a per-row "group" array (parallel to gid/time/ch/pe)
// on the root opflash PC.  No-op when the window is non-positive or the opflash
// arrays are missing, so it is off by default and adds nothing to the output.
// Run pre-pipeline: the array then survives the whole pipeline and is read by
// fill_bee_flashes (and is available to later steps for reuse).
// Group the two per-side flashes of a cathode-crosser for the Bee display.
// A crosser scintillates independently in each drift volume and is matched as
// two cluster halves, one per side, each carrying its own matched flash.  We
// pair ONLY such matched flashes ("good bundles"): within each window-coincident
// neighborhood keep the single closest cross-side pair (one flash per side,
// minimal |dt|) and leave everything else ungrouped.  This replaces the old
// blind single-linkage chaining, which at a wide (1 us) window merged unrelated
// busy-event flashes.  window<=0 => no grouping (bit-identical ungrouped dump).
// greedy=false: one closest cross-side pair per coincidence neighborhood (the
// rest solo).  greedy=true: disjoint pairs, repeatedly taking the next-closest
// available cross-side pair, so two distinct close crossers both group.
static void store_flash_groups(WireCell::Clus::Facade::Grouping& grouping, double window, bool greedy)
{
    if (window <= 0) return;

    auto& lpcs = grouping.local_pcs();
    auto it = lpcs.find("opflash");
    if (it == lpcs.end()) return;
    auto& ds = it->second;
    auto a_gid = ds.get("gid");
    auto a_time = ds.get("time");
    auto a_apa = ds.get("apa");           // physical drift side (0/1)
    if (!a_gid || !a_time) return;
    const auto gid = a_gid->elements<int>();
    const auto time = a_time->elements<double>();
    const size_t nrow = gid.size();
    if (nrow == 0) return;
    const bool have_apa = (a_apa != nullptr);
    std::vector<int> apa_col;
    if (have_apa) { auto ac = a_apa->elements<int>(); apa_col.assign(ac.begin(), ac.end()); }

    // Eligible flashes = those matched to a charge cluster with predicted light
    // >= 100 (the same criterion that yields a "matching: N" row in the Bee op
    // display).  Only these participate in grouping.
    std::set<int> matched_gids;
    for (const auto* cluster : grouping.children()) {
        const int mgid = cluster->get_scalar<int>("matched_flash_gid", -1);
        if (mgid < 0) continue;
        if (!cluster->has_pcarray<double>("pe", "flashpred")) continue;
        auto pred = cluster->get_pcarray<double>("pe", "flashpred");
        double pred_tot = 0; for (double v : pred) pred_tot += v;
        if (pred_tot < 100) continue;
        matched_gids.insert(mgid);
    }

    // Unique flashes (first row per gid): time + drift side.
    struct F { int gid; double t; int side; };
    std::map<int, F> uniq;
    for (size_t i = 0; i < nrow; ++i) {
        if (uniq.find(gid[i]) == uniq.end())
            uniq[gid[i]] = F{gid[i], time[i], have_apa ? apa_col[i] : -1};
    }

    // Eligible (matched, known side) flashes, time-sorted (tie-break gid).
    std::vector<F> elig;
    for (const auto& kv : uniq)
        if (matched_gids.count(kv.first) && kv.second.side >= 0) elig.push_back(kv.second);
    std::sort(elig.begin(), elig.end(), [](const F& a, const F& b) {
        if (a.t != b.t) return a.t < b.t;
        return a.gid < b.gid;
    });

    std::map<int, int> group_of;
    int next_group = 0;
    std::set<int> paired;
    auto make_pair = [&](size_t a, size_t b) {
        const int g = ++next_group;
        group_of[elig[a].gid] = g; group_of[elig[b].gid] = g;
        paired.insert(elig[a].gid); paired.insert(elig[b].gid);
    };
    if (greedy) {
        // Greedy disjoint pairs: take the next-closest available cross-side pair
        // (|dt| <= window) until none remain; each flash used at most once.
        struct Cand { double dt; size_t a, b; };
        std::vector<Cand> cands;
        for (size_t a = 0; a < elig.size(); ++a)
            for (size_t b = a + 1; b < elig.size(); ++b) {
                const double dt = elig[b].t - elig[a].t;   // elig is time-sorted
                if (dt > window) break;                    // no later b can fit
                if (elig[a].side != elig[b].side) cands.push_back({dt, a, b});
            }
        std::sort(cands.begin(), cands.end(), [&](const Cand& x, const Cand& y) {
            if (x.dt != y.dt) return x.dt < y.dt;
            if (elig[x.a].gid != elig[y.a].gid) return elig[x.a].gid < elig[y.a].gid;
            return elig[x.b].gid < elig[y.b].gid;
        });
        std::set<size_t> used;
        for (const auto& c : cands) {
            if (used.count(c.a) || used.count(c.b)) continue;
            make_pair(c.a, c.b);
            used.insert(c.a); used.insert(c.b);
        }
    } else {
        // Single-linkage neighborhoods over eligible flashes; in each, keep only
        // the closest cross-side pair (one per side, |dt| <= window).
        size_t i = 0;
        while (i < elig.size()) {
            size_t j = i + 1;
            while (j < elig.size() && (elig[j].t - elig[j - 1].t) <= window) ++j;
            double best = -1; int bi = -1, bj = -1;
            for (size_t a = i; a < j; ++a)
                for (size_t b = a + 1; b < j; ++b) {
                    if (elig[a].side == elig[b].side) continue;
                    double dt = elig[b].t - elig[a].t; if (dt < 0) dt = -dt;
                    if (dt <= window && (bi < 0 || dt < best)) { best = dt; bi = (int) a; bj = (int) b; }
                }
            if (bi >= 0) make_pair((size_t) bi, (size_t) bj);
            i = j;
        }
    }

    // Everyone else gets a unique singleton id so the viewer never merges them.
    for (const auto& kv : uniq)
        if (!paired.count(kv.first)) group_of[kv.first] = ++next_group;

    // Per-row group array, aligned to the opflash rows.
    std::vector<int> group(nrow);
    for (size_t i = 0; i < nrow; ++i) group[i] = group_of[gid[i]];
    grouping.put_pcarray<int>(group, "group", "opflash");
}

void MultiAlgBlobClustering::fill_bee_flashes(const WireCell::Clus::Facade::Grouping& grouping)
{
    // Run/sub/event numbers (same convention as the points/patches dumps).
    if (m_use_config_rse) {
        m_bee_flash.rse(m_runNo, m_subRunNo, m_eventNo);
    } else {
        int run = 0, evt = 0;
        if (m_last_ident > 0) { run = (m_last_ident >> 16) & 0x7fff; evt = (m_last_ident) & 0xffff; }
        m_bee_flash.reset(evt, 0, run);
    }

    // The self-contained per-flash optical display PC on the (merged) root,
    // written by QLMatching and carried across the per-APA -> all-APA merge.
    // One row per (flash, channel): gid, time (raw ns), ch, pe.
    const auto& lpcs = grouping.local_pcs();
    auto it = lpcs.find("opflash");
    if (it == lpcs.end()) return;   // no flashes attached
    const auto& ds = it->second;
    auto a_gid = ds.get("gid");
    auto a_time = ds.get("time");
    auto a_ch = ds.get("ch");
    auto a_pe = ds.get("pe");
    if (!a_gid || !a_time || !a_ch || !a_pe) return;
    const auto gid = a_gid->elements<int>();
    const auto time = a_time->elements<double>();
    const auto ch = a_ch->elements<int>();
    const auto pe = a_pe->elements<double>();
    const size_t nrow = gid.size();

    // Per-flash physical drift side (TPC 0 = low-x, TPC 1 = high-x), written by
    // QLMatching from the flash's measured light pattern.  Use it when present;
    // fall back to the gid-derived side for older dumps that lack the column.
    // The gid alone encodes the processing node's anode, not the lit volume, so
    // on the merged root it tags every flash with one side (e.g. all "TPC0").
    auto a_apa = ds.get("apa");
    const bool have_apa = (a_apa != nullptr);
    std::vector<int> apa_col;
    if (have_apa) {
        auto ac = a_apa->elements<int>();
        apa_col.assign(ac.begin(), ac.end());
    }

    // Optional matched cross-side flash-pair grouping, stored on the root
    // by store_flash_groups (flash_group_window>0).  Absent => no grouping
    // is emitted and the op JSON stays bit-identical to the ungrouped case.
    auto a_group = ds.get("group");
    const bool have_group = (a_group != nullptr);
    std::vector<int> group_col;
    if (have_group) {
        auto gc = a_group->elements<int>();
        group_col.assign(gc.begin(), gc.end());
    }

    // Optional per-side-clock flash time (QLMatching "time1": the flash time on
    // input-1's/top charge clock, ns; "time" is input-0's/bottom).  Present only
    // when per-input trigger_offsets are configured (PDVD BDE/TDE); absent =>
    // op_t1 is not emitted and the op JSON stays bit-identical (PDHD/SBND).
    auto a_time1 = ds.get("time1");
    const bool have_time1 = (a_time1 != nullptr);
    std::vector<double> time1_col;
    if (have_time1) {
        auto tc = a_time1->elements<double>();
        time1_col.assign(tc.begin(), tc.end());
    }

    // Group rows by global flash id (first-seen order) into per-flash time +
    // dense per-channel measured PE.
    std::vector<int> flash_order;
    std::map<int, double> flash_time;
    std::map<int, double> flash_time1;               // gid -> input-1-clock time
    std::map<int, int> flash_group;                  // gid -> flash-group id
    std::map<int, int> flash_apa;                    // gid -> physical drift side
    std::map<int, std::map<int, double>> flash_pe;   // gid -> (ch -> pe)
    for (size_t i = 0; i < nrow; ++i) {
        const int g = gid[i];
        if (flash_pe.find(g) == flash_pe.end()) {
            flash_order.push_back(g);
            flash_time[g] = time[i];
            if (have_time1) flash_time1[g] = time1_col[i];
            if (have_group) flash_group[g] = group_col[i];
            if (have_apa) flash_apa[g] = apa_col[i];
        }
        flash_pe[g][ch[i]] = pe[i];
    }

    // Order flashes by ascending flash time so the Bee viewer steps through them
    // low->high (it walks the op arrays by index and does no sorting of its own).
    // stable_sort keeps the original first-seen order among equal-time flashes
    // (e.g. the two sides of a cross-side TPC0/TPC1 pair).
    std::stable_sort(flash_order.begin(), flash_order.end(),
                     [&](int a, int b) { return flash_time[a] < flash_time[b]; });

    // Matched clusters: predicted per-channel PE keyed by global flash id, with
    // the same total-predicted-light >= m_bee_flash_pred_min filter as the
    // legacy dump_light (default 100 PE; doc pr/94 sec 9.9 -- a genuine match
    // under the cut is drawn as "no flash match", which is a display artifact
    // and not a statement about the matching).
    // cluster_id is the cluster's own id, identical to the "img" charge dump
    // enumeration (this runs at the same pre-pipeline point), so the Bee viewer
    // associates each flash to the same physical charge cluster.
    std::map<int, std::vector<std::pair<int, std::vector<double>>>> matched;
    for (const auto* cluster : grouping.children()) {
        const int mgid = cluster->get_scalar<int>("matched_flash_gid", -1);
        if (mgid < 0) continue;
        if (!cluster->has_pcarray<double>("pe", "flashpred")) continue;
        auto pred_span = cluster->get_pcarray<double>("pe", "flashpred");
        std::vector<double> pred(pred_span.begin(), pred_span.end());
        double pred_tot = 0;
        for (double v : pred) pred_tot += v;
        // doc pr/94 round 3 (owner: "why is this piece shown as non-matched in
        // Bee?").  Dump every genuinely matched cluster together with the
        // predicted light that decides whether the display keeps it, at THIS
        // stage so the ids printed are exactly the ones the "img"/"op" JSON
        // carries (enumerate_idents re-issues ids after every visitor, so the
        // QLMatching log's per-run idents are a different epoch and cannot be
        // bridged after the fact).  Diagnostic only; gated so it costs nothing
        // and changes nothing when unset.
        if (std::getenv("WCT_OPDUMP_DEBUG")) {
            log->info("op-dump debug: cluster {} matched_flash_gid={} pred_tot={:.3f} PE "
                      "L={:.2f} cm nblobs={} kept_by_display={}",
                      cluster->get_cluster_id(), mgid, pred_tot,
                      cluster->get_length() / units::cm, cluster->nchildren(),
                      pred_tot >= 100);
        }
        if (pred_tot < m_bee_flash_pred_min) continue;
        matched[mgid].push_back({cluster->get_cluster_id(), std::move(pred)});
    }

    // Emit ALL flashes (matched + unmatched). A matched flash emits one row per
    // matched cluster; an unmatched flash emits one row with empty cluster_id.
    // The flash's drift side (apa) comes from the "apa" column when present
    // (QLMatching derives it from the measured light pattern); otherwise fall
    // back to the legacy gid encoding gid = anode_ident*kFlashGidStride + idx,
    // so apa = gid / kFlashGidStride (correct only for single-face anodes).
    constexpr int kFlashGidStride = 1000000;
    std::vector<int> appended_groups;   // one per appended row, same order
    std::vector<double> appended_t1;    // ditto, input-1-clock time (us)
    for (const int g : flash_order) {
        const int apa = have_apa ? flash_apa[g] : (g / kFlashGidStride);
        const int grp = have_group ? flash_group[g] : -1;
        const double t1_us = have_time1 ? flash_time1[g] * 1e-3 : 0.0;   // ns -> us
        int maxch = -1;
        for (const auto& cv : flash_pe[g]) if (cv.first > maxch) maxch = cv.first;
        std::vector<double> pes(maxch + 1, 0.0);
        double peTotal = 0;
        for (const auto& cv : flash_pe[g]) { pes[cv.first] = cv.second; peTotal += cv.second; }
        const double t_us = flash_time[g] * 1e-3;   // ns -> us (matches dump_light)

        auto mit = matched.find(g);
        if (mit != matched.end() && !mit->second.empty()) {
            if (m_bee_flash_per_flash) {
                // One row per flash: all matched cluster ids together, and the
                // predicted light is the element-wise sum over the matched
                // clusters (each cluster contributes its own predicted PE).
                std::vector<int> cids;
                std::vector<double> pred_sum;
                for (const auto& cp : mit->second) {
                    cids.push_back(cp.first);
                    if (pred_sum.size() < cp.second.size()) pred_sum.resize(cp.second.size(), 0.0);
                    for (size_t k = 0; k < cp.second.size(); ++k) pred_sum[k] += cp.second[k];
                }
                m_bee_flash.append(t_us, pes, peTotal, cids, pred_sum, apa);
                appended_groups.push_back(grp);
                appended_t1.push_back(t1_us);
            } else {
                for (const auto& cp : mit->second) {
                    m_bee_flash.append(t_us, pes, peTotal, std::vector<int>{cp.first}, cp.second, apa);
                    appended_groups.push_back(grp);
                    appended_t1.push_back(t1_us);
                }
            }
        } else {
            m_bee_flash.append(t_us, pes, peTotal, std::vector<int>{}, std::vector<double>{}, apa);
            appended_groups.push_back(grp);
            appended_t1.push_back(t1_us);
        }
    }

    // Attach the per-row flash-group array only when grouping was computed, so
    // the ungrouped output is unchanged.
    if (have_group) m_bee_flash.set_groups(appended_groups);
    // Attach the per-row input-1-clock time only when the opflash PC carries it
    // (per-input trigger_offsets, PDVD), so other detectors' op JSON is unchanged.
    if (have_time1) m_bee_flash.set_t1(appended_t1);
}

struct Perf {
    using Clock = std::chrono::steady_clock;
    using MS    = std::chrono::duration<double, std::milli>;

    bool enable;
    Log::logptr_t log;
    ExecMon em;
    Clock::time_point t_start;
    Clock::time_point t_last;

    Perf(bool e, Log::logptr_t l, const std::string& t = "starting MultiAlgBlobClustering")
      : enable(e)
      , log(l)
      , em(t)
    {
        t_start = t_last = Clock::now();
    }

    ~Perf()
    {
        if (!enable) return;
        SPDLOG_LOGGER_DEBUG(log, "MultiAlgBlobClustering performance summary:\n{}", em.summary());
    }

    void operator()(const std::string& ctx)
    {
        if (!enable) return;
        auto now = Clock::now();
        SPDLOG_LOGGER_DEBUG(log, "MABC timing: {} took {} ms (cumulative {} ms)",
                            ctx, MS(now - t_last).count(), MS(now - t_start).count());
        t_last = now;
        em(ctx);
    }

    void dump(const std::string& ctx, const Ensemble& ensemble, bool shallow = true, bool mon = true)
    {
        if (!enable) return;
        if (mon) (*this)(ctx);

        SPDLOG_LOGGER_TRACE(log, "{} ensemble with {} groupings:", ctx, ensemble.nchildren());

        for (const auto* grouping : ensemble.children()) {

            {
                auto name = grouping->get_name();
                size_t npoints_total = 0;
                size_t nzero = 0;
                size_t count = 0;
                for (const auto* cluster : grouping->children()) {
                    int n = cluster->npoints();
                    if (n == 0) {
                        ++nzero;
                    }
                    npoints_total += n;
                    // SPDLOG_LOGGER_DEBUG(log, "loaded cluster {} with {} points out of {}", count, n, npoints_total);
                    ++count;
                    // std::cout << "Xin: " << name << " loaded cluster " << count << " with " << n << "points and " << cluster->nchildren() << "blobs" << std::endl;
                }

                

                SPDLOG_LOGGER_TRACE(log, "\tgrouping \"{}\": {}, {} points and {} clusters with no points",
                           name, *grouping, npoints_total, nzero);
                (void)count;
            }

            if (shallow) continue;

            auto children = grouping->children();  // copy
            sort_clusters(children);
            size_t count = 0;
            for (const auto* cluster : children) {
                bool sane = cluster->sanity(log);
                SPDLOG_LOGGER_TRACE(log, "\t\tcluster {} {} sane:{}", count++, *cluster, sane);
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
        SPDLOG_LOGGER_WARN(log, "No pc-tree at tensor datapath {}, making empty", path);
        ensemble.make_grouping(name);
    }
        
    Grouping* grouping = ensemble.with_name(name).at(0);
    if (!grouping) {
        raise<KeyError>("failed to make grouping node %s at %s", name, path);
    }

    grouping->enumerate_idents();
    grouping->set_anodes(m_anodes);
    grouping->set_detector_volumes(m_dv);
    check_perblob_provenance(*grouping->node(), "load:" + path);
    return *grouping;
}

bool MultiAlgBlobClustering::operator()(const input_pointer& ints, output_pointer& outts)
{
    outts = nullptr;
    if (!ints) {
        flush();
        SPDLOG_LOGGER_DEBUG(log, "EOS at call {}", m_count++);
        return true;
    }

    Perf perf{m_perf, log};

    const int ident = ints->ident();
    SPDLOG_LOGGER_DEBUG(log, "loading tensor set ident={} (last={})", ident, m_last_ident);
    if (m_last_ident < 0) {     // first time.
        if (m_rse_from_ident) {
            // The tensor ident already carries the real event id (raw, unmasked);
            // run/subrun are not available in this chain.
            m_runNo = 0; m_subRunNo = 0; m_eventNo = ident;
            if (!m_use_shared_sink) m_sink.set_rse(m_runNo, m_subRunNo, m_eventNo);
        }
        else if (m_event_from_ident) {
            apply_event_from_ident(ident);
        }
        else if (m_use_config_rse && !m_use_shared_sink) {
            // Set RSE in the sink (shared-sink mode passes RSE per write_obj).
            m_sink.set_rse(m_runNo, m_subRunNo, m_eventNo);
        }
        // Use default behavior
        // reset_bee(ident, m_bee_img);
        // reset_bee(ident, m_bee_ld);
        m_last_ident = ident;
    }
    else if (m_last_ident != ident) {
        flush(ident);   // writes the previous event with its still-current RSE
        if (m_rse_from_ident) {
            m_runNo = 0; m_subRunNo = 0; m_eventNo = ident;
            if (!m_use_shared_sink) m_sink.set_rse(m_runNo, m_subRunNo, m_eventNo);
        }
        else if (m_event_from_ident) {
            apply_event_from_ident(ident);
        }
        else if (m_use_config_rse) {
            // Update event number for next event
            m_eventNo++;
            // Update RSE in sink (shared-sink mode passes RSE per write_obj).
            if (!m_use_shared_sink) {
                m_sink.set_rse(m_runNo, m_subRunNo, m_eventNo);
            }
        }
    }
    // else do nothing when ident is unchanged.


    Points::node_t root;
    Ensemble& ensemble = *root.value.facade<Ensemble>();
    // Tell the visitors which event this is.  A visitor that writes one file
    // per event (SbndPrMagnifyTrackingVisitor, PrDisplayDump) has no other way
    // to know: its own runNo/eventNo come from configure() and are constant for
    // the whole process.
    ensemble.set_ident(ident);
    // Publish the RSE this node resolved (config, auto-increment, ident or
    // rse_map) so the per-event writers downstream stamp THIS event rather than
    // their own configure-time constant.  Only when a multi-event mode is on:
    // in the one-event-per-process production path the visitors' own numbers
    // are already right and must keep being used, byte for byte.
    if (m_rse_from_ident || m_event_from_ident) {
        ensemble.set_rse(m_runNo, m_subRunNo, m_eventNo);
    }

    for (const auto& gname : m_groupings) {
        const auto datapath = inpath(gname, ident);
        load_grouping(ensemble, gname, datapath, ints);
        perf.dump("loaded " + gname, ensemble);
    }    

    if (m_save_deadarea) {
        // Fill patches from the dead grouping (not "live" — that was a
        // regression from the ensemble-facade refactor; the result was
        // empty Bee::Patches and no channel-deadarea-*.json in the
        // mabc-*.zip output).
        auto gs = ensemble.with_name("dead");
        if (gs.size()) {
            fill_bee_patches_from_grouping(*gs[0]);
            perf("dump dead regions to bee");
        }
    }

    perf.dump("pre clustering", ensemble);

    for (const auto& config : m_bee_points_configs) {
        if (config.name != "img" && !config.prepipeline) {
            continue;
        }
        auto gs = ensemble.with_name("live");
        if (gs.empty()) {
            continue;
        }
        fill_bee_points(config.name, *gs[0]);
    }

    // Dump the optical flash / charge-light "op" display at the SAME point as
    // the "img" charge dump, BEFORE the clustering pipeline runs: at this point
    // the live clusters are exactly the per-APA matched clusters (their
    // matched-flash association + predicted PE intact, cluster ids == the "img"
    // enumeration). The pipeline below re-clusters/merges them, after which the
    // 1:1 cluster<->flash mapping no longer holds.
    if (m_save_opflash) {
        auto gs = ensemble.with_name("live");
        if (gs.size()) {
            // Stash the ±window TPC0/TPC1 flash grouping on the root opflash PC
            // first, so the op dump (and every later pipeline step) can read it.
            store_flash_groups(*gs[0], m_flash_group_window, m_flash_group_greedy);
            fill_bee_flashes(*gs[0]);
            perf("dump op flashes to bee");
        }
    }

    perf.dump("start clustering", ensemble);

    // THE MAIN LOOP
    for (const auto& cmeth : m_pipeline) {
        cmeth.meth->visit(ensemble);
        perf.dump(cmeth.name, ensemble);

        for (auto* grouping : ensemble.children()) {
            grouping->enumerate_idents(m_clusters_id_order);
        }
        {
            auto gs = ensemble.with_name("live");
            if (gs.size()) check_perblob_provenance(*gs[0]->node(), "post:" + cmeth.name);
        }

        // Dump bee points right after specific visitor runs
        for (const auto& config : m_bee_points_configs) {
            if (config.name == "img" || config.prepipeline) continue;
            if (config.visitor.empty() || config.visitor != cmeth.name) continue;

            auto gs = ensemble.with_name(config.grouping);
            if (gs.empty()) {
                continue;
            }

            // Check if this visitor produced a PRGraph that we should save
            auto pr_graph = gs[0]->get_pr_graph();

            // std::cout << "Test: Visitor: " << cmeth.name << " Grouping: " << config.grouping << " " << pr_graph << std::endl;

            if (pr_graph) {
                // doc pr/94 Phase 4b: the point layers (track_fit,
                // shower_track, vertices) render every per-bundle candidate,
                // the way the "mc" particle-flow layer already does.  The
                // "nu<i>" slots exist only in per-bundle mode; with none
                // present this is exactly the single legacy call, and the
                // pr_graph test above still gates it correctly because
                // candidate 0 always also publishes to the unnamed slot.
                std::vector<std::shared_ptr<WireCell::Clus::TrackFitting>> nu_tfs;
                for (int i = 0;; ++i) {
                    auto tfi = gs[0]->get_track_fitting("nu" + std::to_string(i));
                    if (!tfi) break;
                    nu_tfs.push_back(tfi);
                }
                if (nu_tfs.empty()) nu_tfs.push_back(nullptr);   // legacy: unnamed slot
                for (size_t i = 0; i < nu_tfs.size(); ++i) {
                    // reset ONLY on the first pass, else bundle i wipes i-1.
                    const bool first = (i == 0);
                    if (config.use_graph_vertices) {
                        fill_bee_vertices_from_pr_graph(config.name, *gs[0], nu_tfs[i], first);
                    } else {
                        // Fill bee points from PRGraph (for track trajectories)
                        fill_bee_points_from_pr_graph(config.name, *gs[0], nu_tfs[i], first);
                    }
                }
                // std::cout << "Filled bee points from PR graph for visitor: " << cmeth.name << " grouping: " << config.grouping << std::endl;
            } else if (config.require_pr_graph) {
                // PR-output set with no PR graph: leave it EMPTY.  The generic
                // dump below would fill a track_fit/shower_track/vertices layer
                // with the whole clustering point set -- and in RAW coordinates,
                // since those sets declare coords x/y/z (PR fit points are
                // already T0-corrected).  SBND evt 18255/52195 with the bundle
                // veto on: 30127 points per layer, x out to -234 cm, three
                // byte-identical copies (sbnd_xin/docs/pr/3 sec. 9).
                SPDLOG_LOGGER_DEBUG(log, "bee points set '{}': visitor {} produced no PR graph; leaving the set empty (require_pr_graph)",
                                    config.name, cmeth.name);
            } else {
                // Fill bee points from clusters normally
                fill_bee_points(config.name, *gs[0]);
                // std::cout << "Filled bee points from clusters for visitor: " << cmeth.name << " grouping: " << config.grouping << std::endl;
            }
        }

        // Particle-flow dump triggered by the same visitor
        for (const auto& pf_cfg : m_bee_pf_configs) {
            if (pf_cfg.visitor.empty() || pf_cfg.visitor != cmeth.name) continue;
            auto pf_gs = ensemble.with_name(pf_cfg.grouping);
            if (pf_gs.empty()) continue;
            const auto& pf_grouping = *pf_gs[0];
            auto tf = pf_grouping.get_track_fitting();
            if (!tf) continue;
            // doc pr/94 Phase 4: render every per-bundle candidate.  The
            // "nu<i>" named slots exist only in per-bundle mode, so with none
            // present this is exactly the single legacy call, byte-identical.
            std::vector<std::shared_ptr<WireCell::Clus::TrackFitting>> nu_tfs;
            for (int i = 0;; ++i) {
                auto tfi = pf_grouping.get_track_fitting("nu" + std::to_string(i));
                if (!tfi) break;
                nu_tfs.push_back(tfi);
            }
            if (nu_tfs.empty()) {
                fill_bee_pf_tree(pf_cfg, pf_grouping);
            }
            else {
                Configuration all = Json::arrayValue;
                std::set<int> used_ids;
                for (const auto& tfi : nu_tfs) {
                    fill_bee_pf_tree(pf_cfg, pf_grouping, false, tfi, &used_ids, &all);
                }
                auto pf_it = m_bee_pf_trees.find(pf_cfg.name);
                if (pf_it != m_bee_pf_trees.end()) pf_it->second.set_particles(all);
                SPDLOG_LOGGER_DEBUG(log, "fill_bee_pf_tree '{}': nu_per_bundle wrote {} bundle root(s)",
                                    pf_cfg.name, all.size());
            }
        }
    }

    //
    // At this point, the ensemble may have more or fewer groupings just "live"
    // and "dead" including no groupings at all.  But for now, we assume the
    // original "live" and "dead" still exist and with their original facades.
    // Famous last words....
    //
    

    // Collapse real_cluster_id into one ident epoch BEFORE any Bee fill below and
    // before the tensor save, so the Bee per-blob labels and the saved pctree
    // agree (doc 53).  No-op where the array was never written.
    if (m_real_cluster_id_global) {
        for (const auto& gname : m_groupings) {
            auto gs = ensemble.with_name(gname);
            if (gs.empty()) continue;
            restamp_real_cluster_id(*gs[0]);
        }
    }

    // Fill all configured bee points sets (except those with visitor-specific handling)
    for (const auto& config : m_bee_points_configs) {
        if(config.name == "img" || config.prepipeline) continue;

        // Skip configs with visitor specified - they were already handled in the visitor loop
        if (!config.visitor.empty()) continue;

        auto gs = ensemble.with_name(config.grouping);
        if (gs.empty()) {
            continue;
        }
        fill_bee_points(config.name, *gs[0]);

    }
    perf("dump live clusters to bee");

    if (m_grouping2file_prefix.size()) {
        std::string fname = String::format("%s-%d.npz", m_grouping2file_prefix, m_count);
        auto live = ensemble.with_name("live");
        grouping2file(*live[0], fname);
    }
    auto grouping_names = ensemble.names();

    if (m_dump_json) {
        for (const auto& name : grouping_names) {
            auto gs = ensemble.with_name(name);
            Persist::dump(String::format("%s-summary-%d.json", name, ident),
                          json_summary(*gs[0]), true);
        }
    }

    SPDLOG_LOGGER_DEBUG(log, "Produce pctrees with {} groupings", grouping_names.size());
    
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
        normalize_cluster_flags(grouping, log, name, ident);
        // Homogenize the "perblob" key set so the flash-merge provenance
        // survives serialization (see m_save_real_cluster_id in the header).
        // Done after normalize_cluster_flags so the fill-in "main_cluster"
        // flag values are the final, normalized ones.
        // The isolated grouping's pair (doc 52) must be homogenized FIRST, and it
        // needs a wider gate than the flash pair's.  The invariant that
        // Dataset::append enforces is "every cluster that has a 'perblob' PC has
        // the SAME keys in it" -- it raises when an incoming per-blob dataset is
        // missing a key the target already has
        // (util/src/PointCloudDataset.cxx:261); a cluster with no "perblob" PC at
        // all is simply absent and fine.  Gating on "isolated" the way the flash
        // loop does is not sufficient here, because the assoc pair already exists
        // when the all-APA switch_scope runs, so its carve puts it on the
        // out-of-volume shards too -- and examine_bundles then skips those (not in
        // scope), so they never receive "isolated".  Measured symptom: the whole
        // per-APA -> all-APA handoff threw "missing keys in append: 3 missing:
        // isolated real_cluster_id real_cluster_main".
        // So: any cluster that has a "perblob" PC gets the full key set.  -1 in
        // "isolated" is that array's documented "this blob is in the main
        // sub-component" value, which is exactly what a cluster nothing merged is.
        // Running before the flash loop lets that loop's "isolated" gate then
        // cover these clusters too.
        if (m_save_assoc_cluster_id) {
            for (Cluster* cluster : grouping.children()) {
                const auto& lpcs = cluster->value().local_pcs();
                if (lpcs.find("perblob") == lpcs.end()) continue;
                const size_t nb = cluster->nchildren();
                if (!cluster->has_pcarray<int>("isolated", "perblob")) {
                    cluster->put_pcarray(std::vector<int>(nb, -1), "isolated", "perblob");
                }
                if (!cluster->has_pcarray<int>("assoc_cluster_id", "perblob")) {
                    cluster->put_pcarray(std::vector<int>(nb, cluster->ident()),
                                         "assoc_cluster_id", "perblob");
                }
                // A cluster the isolated grouping never merged is a main, not an
                // associated fragment: all 1.  This is the "absent provenance =>
                // main" sentinel the un-merge needs to keep crossers whole.
                if (!cluster->has_pcarray<int>("assoc_cluster_main", "perblob")) {
                    cluster->put_pcarray(std::vector<int>(nb, 1),
                                         "assoc_cluster_main", "perblob");
                }
            }
        }
        if (m_save_real_cluster_id) {
            for (Cluster* cluster : grouping.children()) {
                if (!cluster->has_pcarray<int>("isolated", "perblob")) continue;
                const size_t nb = cluster->nchildren();
                if (!cluster->has_pcarray<int>("real_cluster_id", "perblob")) {
                    cluster->put_pcarray(std::vector<int>(nb, cluster->ident()),
                                         "real_cluster_id", "perblob");
                }
                // An unmerged cluster is its own representative: all 1.
                if (!cluster->has_pcarray<int>("real_cluster_main", "perblob")) {
                    cluster->put_pcarray(std::vector<int>(nb, 1),
                                         "real_cluster_main", "perblob");
                }
            }
            // real_cluster_id_global's re-stamp is NOT here: it must run before
            // the Bee fills so the pctree and the Bee zip carry the SAME ids
            // (doc 53).  See restamp_real_cluster_id(), called right after the
            // clustering pipeline.
        }
        // "real_cluster_was_main" (doc pr/20 Part I P1) needs the same
        // key-homogeneity fill-in as the pair above, but is deliberately
        // PRESENCE-triggered instead of getting a knob of its own: the writer is
        // ClusteringExamineBundles' save_bundle_main_provenance, and two
        // independent booleans would admit the one configuration that throws
        // (writer on, fill-in off => a "perblob" PC whose key set differs
        // between clusters => Dataset::append raises on the next merge, see
        // aux/src/TensorDMpointtree.cxx).  Fill iff somebody actually wrote it.
        {
            bool any_wasmain = false;
            for (Cluster* cluster : grouping.children()) {
                if (cluster->has_pcarray<int>("real_cluster_was_main", "perblob")) {
                    any_wasmain = true;
                    break;
                }
            }
            if (any_wasmain) {
                // Gate on "has a perblob PC at all", which is the exact scope of
                // the invariant append enforces -- wider than the flash pair's
                // "isolated" gate, for the reason spelled out above it.
                for (Cluster* cluster : grouping.children()) {
                    const auto& lpcs = cluster->value().local_pcs();
                    if (lpcs.find("perblob") == lpcs.end()) continue;
                    if (cluster->has_pcarray<int>("real_cluster_was_main", "perblob")) continue;
                    // A cluster the flash merge never touched still speaks for
                    // itself: all 1, the same "absent provenance => own main"
                    // sentinel real_cluster_main uses.
                    cluster->put_pcarray(std::vector<int>(cluster->nchildren(), 1),
                                         "real_cluster_was_main", "perblob");
                }
            }
        }
        // "nu_band_veto_role" (ClusteringNeutrino's record_band_veto, doc
        // pr/66) needs the same key-homogeneity fill-in as the block above,
        // and for the same reason: presence-triggered rather than a knob of
        // its own, so writer-on/fill-in-off can never happen and leave a
        // "perblob" PC whose key set differs between clusters (Dataset::append
        // raises at the next merge).  Fill iff somebody actually wrote it.
        {
            bool any_bandveto = false;
            for (Cluster* cluster : grouping.children()) {
                if (cluster->has_pcarray<int>("nu_band_veto_role", "perblob")) {
                    any_bandveto = true;
                    break;
                }
            }
            if (any_bandveto) {
                // Gate on "has a perblob PC at all" -- same wide gate as
                // real_cluster_was_main above, for the same reason.
                for (Cluster* cluster : grouping.children()) {
                    const auto& lpcs = cluster->value().local_pcs();
                    if (lpcs.find("perblob") == lpcs.end()) continue;
                    if (cluster->has_pcarray<int>("nu_band_veto_role", "perblob")) continue;
                    // 0 = unmarked, the array's own documented value -- and
                    // exactly what PointTreeMerging's normalize_pctree_local_pcs
                    // would zero-fill anyway, so the two agree.
                    cluster->put_pcarray(std::vector<int>(cluster->nchildren(), 0),
                                         "nu_band_veto_role", "perblob");
                }
            }
        }
        auto node = ensemble.remove_child(grouping);
        check_perblob_provenance(*node, "save:" + outpath(name, ident));
        auto tens = as_tensors(*node, outpath(name, ident));
        // Serialize -> deserialize self-test: if "save" is clean but "rtrip"
        // is not, the corruption is inside the TensorDM round trip itself.
        if (std::getenv("WCT_PROV_CHECK")) {
            auto rt = Aux::TensorDM::as_pctree(tens, outpath(name, ident));
            if (rt) check_perblob_provenance(*rt, "rtrip:" + outpath(name, ident));
        }
        outtens.insert(outtens.end(), tens.begin(), tens.end());
        SPDLOG_LOGGER_DEBUG(log, "Produce {} tensors for grouping {}", tens.size(), name);
    }
    outts = as_tensorset(outtens, ident);

    perf("done");

    return true;
}
