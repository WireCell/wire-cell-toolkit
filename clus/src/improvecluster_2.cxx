// ImproveCluster_2 - Second level cluster improvement 
//
// This class inherits from ImproveCluster_1 and provides additional
// cluster improvement functionality, building upon the Steiner tree
// enhancements from the first level.

#include "improvecluster_1.h"  // Include the ImproveCluster_1 header
#include "SteinerGrapher.h"
#include "RetileDump.h"
#include "WireCellClus/ResampleLive.h"
#include "WireCellClus/BeamParticleFunctions.h"
#include "WireCellClus/GapFill.h"
#include "WireCellAux/CascadeRun.h"
#include "WireCellIface/ITensorForward.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Units.h"
#include <atomic>
#include <chrono>
#include <climits>
#include <cmath>
#include <map>
#include <set>

#include <vector>

namespace WireCell::Clus {

    class ImproveCluster_2 : public ImproveCluster_1 {

    public:

        ImproveCluster_2();
        virtual ~ImproveCluster_2();

        // IConfigurable API - extend the base configuration
        void configure(const WireCell::Configuration& config) override;
        virtual Configuration default_configuration() const override;

        // IPCTreeMutate API - override to add second level improvements
        virtual std::unique_ptr<node_t> mutate(node_t& node) const override;

    private:

        // doc pdvd/31 round 6, the owner's Q5: this component ran the Steiner
        // terminal finder at the C++ default 4000 e regardless of what
        // CreateSteinerGraph was configured with, and the value was unreachable
        // from jsonnet -- so PDVD carried two terminal thresholds in one stage
        // (500 here, 4000 there).  The knob makes them syncable from one place.
        // Default 4000 = the historical value = byte-identical when the key is
        // absent (every detector today).
        double m_steiner_terminal_charge{4000.0};

        // doc pdvd/113: how much of the retile the Steiner copy gets.  The
        // retile (inherited from MicroBooNE, which has far more dead channels)
        // fabricates activity before re-tiling; each mode below removes one
        // more piece of that fabrication.  The original cluster's basic_pid
        // path is computed in every mode, so the graph it caches on the source
        // cluster (TrackFitting, PRSegmentFunctions read it) is unchanged.
        //   "full"      -- the historical retile (C++ default => a config
        //                  without the key is byte-identical);
        //   "no_paint"  -- no hack_activity_improved (neither the orig- nor the
        //                  temp-path disc painting); dead channels and good
        //                  charge within 20 cm are still added and re-tiled;
        //   "footprint" -- also no dead / nearby-charge extension: the cluster's
        //                  own blob footprints are re-tiled;
        //   "none"      -- no tiling at all: every original blob is re-sampled
        //                  in place with the configured sampler (charge_stepped
        //                  in PDHD/PDVD production), activity rebuilt from the
        //                  grouping's CTPC the way ClusteringResampleLive does.
        std::string m_retile_mode{"full"};

        // doc pdvd/129 phase 1.  dump_dir (C++ default "" = off: nothing is
        // written and the retile is byte-identical).  Set => for every cluster
        // retiled in "full" mode the final pass's tiled blobs are written
        // BEFORE remove_bad_blobs, as an imaging cluster archive plus a JSON
        // side file (RetileDump.h).  Study only; the returned cluster is
        // unchanged.  The reduced retile modes do not dump.
        std::string m_dump_dir{""};
        mutable std::atomic<int> m_dump_serial{0};

        // The non-"full" modes.  orig_cluster is the reinitialized source.
        std::unique_ptr<node_t> mutate_reduced(Cluster& orig_cluster) const;

        // doc pdvd/130 phase 3: deghosting inside the retile.  Config object
        // "retile_deghost" (absent = off, nothing below runs, output byte-
        // identical).  Full mode, final pass, per face: tile, sample, filter 1
        // (the own-blob comparison, setting `filter`), then the cascade model
        // (Aux::Cascade::run_cascade with the slice sums of the retile's
        // slices) on the surviving blobs, the Steiner repair with the length
        // pricing, and the kept cells sampled in place of the survivors.  The
        // source cluster is tagged (scalar retile_deghost_writeback = 1) so
        // CreateSteinerGraph replaces its blobs by the copy's.
        //   filter  "run" (the members: today's rule), "vote" (group vote
        //           only), "painted_run" (run bound only on blobs with a plane
        //           without measured charge), "none" (no filter 1)
        //   scope   "bundle" (the beam flash bundle, beam_particle_pick_bundle
        //           over the grouping's mains in [beam_window_low, high)),
        //           "all" (every cluster retiled: the simulation gate)
        //   levels, cut_*, guard, max_level_nodes, repair, repair_*, policy,
        //   charge_scale, uncer_cut, nthreads, final_guard, iso_fallback,
        //   ident_base: as CascadeDeghosting
        //   bridge_len_cost, bridge_max: the length pricing of the repair
        //   dump_dir: the model's input / output cells (RetileDump, "rdin" /
        //           "rdout"), for the replay-equality gate
        struct RetileDeghost {
            bool on{false};
            std::string filter{"run"};
            std::string scope{"bundle"};
            double beam_window_low{0.0}, beam_window_high{0.0};
            std::vector<Aux::Cascade::LevelSpec> levels;
            Aux::Cascade::RunParams rp;
            double charge_scale{0.25};
            double uncer_cut{1e11};
            double bridge_len_cost{0.0};
            double bridge_max{0.0};
            std::string dump_dir;
            // doc pdvd/131: mode "retile" (doc 130: tile the 20 cm extension and the painting, filter 1,
            // the model) or "gapfill": the cluster's own blobs are the fixed base (with -beam-deghost the
            // imaging-level cascade cells); no extension, no painting, no filter 1; the gaps of the base
            // along the cluster's path are bridged by one-wire corridors tiled at cell resolution; wires
            // owned by other clusters' original blobs (owner map of the retile's region) are never
            // admitted; the model sees the residual charge (measured minus what others' blobs explain) and
            // decides the corridor cells, the Steiner repair picks the chain.
            std::string mode{"retile"};
            bool residual{true};             // residual charge (owner ruling: doc 131 step D)
            int corridor_wires{1};           // corridor half-width in wires per plane (ruling 3)
            int gap_min_slices{2};           // a gap spans at least this many slices
            double gap_max{20 * units::cm};  // longer gaps are left to the Steiner graph
            double reach{20 * units::cm};    // the owner map's region = the base's bounds widened by this
        };
        RetileDeghost m_rd;
        // doc pdvd/131 gapfill, one face: the base cells and the corridor candidates as IBlobs
        struct GapFace {
            std::vector<IBlob::pointer> base, cand;
            size_t npath{0}, nuncovered{0}, ngaps{0}, nskipped{0}, nowner{0}, ncand_tiled{0};
        };
        GapFace rd_gapfill_face(const Cluster& orig, const std::vector<size_t>& path, int apa, int face) const;
        mutable std::atomic<int> m_rd_call{0};
        mutable const Grouping* m_rd_pick_grouping{nullptr};   // compared, never dereferenced
        mutable int m_rd_pick_gid{-1};
        bool rd_in_scope(const Cluster& cluster) const;
        void configure_retile_deghost(const WireCell::Configuration& jrd);

    };

} // namespace WireCell::Clus

WIRECELL_FACTORY(ImproveCluster_2, WireCell::Clus::ImproveCluster_2,
                 WireCell::INamed, WireCell::IConfigurable, WireCell::IPCTreeMutate)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;
using namespace WireCell::PointCloud::Tree;

// Segregate this weird choice for namespace.
namespace WCF = WireCell::Clus::Facade;

// Nick name for less typing.
namespace WRG = WireCell::RayGrid;
namespace WireCell::Clus {

    ImproveCluster_2::ImproveCluster_2()
    {
    }

    ImproveCluster_2::~ImproveCluster_2() 
    {
    }

    void ImproveCluster_2::configure(const WireCell::Configuration& cfg)
    {
        // doc pdvd/113.  Absent key => "full" => byte-identical.  Validated
        // before the base configure (which needs live detector volumes) so a
        // typo fails on its own message.
        m_retile_mode = get<std::string>(cfg, "retile_mode", m_retile_mode);
        m_dump_dir = get<std::string>(cfg, "dump_dir", m_dump_dir);
        if (m_retile_mode != "full" && m_retile_mode != "no_paint" &&
            m_retile_mode != "footprint" && m_retile_mode != "none") {
            raise<ValueError>("ImproveCluster_2: unknown retile_mode '%s' "
                              "(expected full, no_paint, footprint or none)", m_retile_mode.c_str());
        }

        // doc pdvd/130 phase 3: absent key => off
        if (cfg.isMember("retile_deghost") && cfg["retile_deghost"].isObject()) {
            configure_retile_deghost(cfg["retile_deghost"]);
        }

        // Configure base class first
        ImproveCluster_1::configure(cfg);

        // doc pdvd/31 round 6 (owner Q5).  Same key name as CreateSteinerGraph's
        // so one jsonnet value can feed both and they cannot drift apart.
        m_steiner_terminal_charge =
            get(cfg, "terminal_charge_threshold", m_steiner_terminal_charge);
    }

    void ImproveCluster_2::configure_retile_deghost(const WireCell::Configuration& j)
    {
        auto& rd = m_rd;
        rd.on = true;
        rd.filter = get<std::string>(j, "filter", rd.filter);
        if (rd.filter != "run" && rd.filter != "vote" && rd.filter != "painted_run" && rd.filter != "none") {
            raise<ValueError>("ImproveCluster_2 retile_deghost: unknown filter '%s' (run, vote, painted_run, none)", rd.filter.c_str());
        }
        rd.scope = get<std::string>(j, "scope", rd.scope);
        if (rd.scope != "bundle" && rd.scope != "all") {
            raise<ValueError>("ImproveCluster_2 retile_deghost: unknown scope '%s' (bundle, all)", rd.scope.c_str());
        }
        rd.beam_window_low = get(j, "beam_window_low", rd.beam_window_low);
        rd.beam_window_high = get(j, "beam_window_high", rd.beam_window_high);
        if (rd.scope == "bundle" && !(rd.beam_window_low < rd.beam_window_high)) {
            raise<ValueError>("ImproveCluster_2 retile_deghost: scope bundle needs beam_window_low < beam_window_high");
        }
        auto& rp = rd.rp;
        rp.cut.max_depth = get(j, "cut_max_depth", rp.cut.max_depth);
        rp.cut.min_length = get(j, "cut_min_length", rp.cut.min_length);
        rp.cut.nudge = get(j, "cut_nudge", rp.cut.nudge);
        rp.guard = get(j, "guard", rp.guard);
        rp.max_level_nodes = get(j, "max_level_nodes", rp.max_level_nodes);
        rp.repair = get(j, "repair", rp.repair);
        rp.steiner.p_term = get(j, "repair_p_term", rp.steiner.p_term);
        rp.steiner.q_floor = get(j, "repair_q_floor", rp.steiner.q_floor);
        rp.steiner.budget = get(j, "repair_budget", 0.0);
        rp.policy = get(j, "policy", rp.policy);
        rp.nthreads = std::max(1, get(j, "nthreads", rp.nthreads));
        rp.final_guard = get(j, "final_guard", rp.final_guard);
        rp.ident_base = get(j, "ident_base", rp.ident_base);
        rp.dump_dir = get<std::string>(j, "cascade_dump_dir", "");
        const auto& jiso = j["iso_fallback"];
        rp.iso = jiso.isObject();
        if (rp.iso) {
            rp.iso_params.nmin = get(jiso, "nmin", rp.iso_params.nmin);
            rp.iso_params.mmin = get(jiso, "mmin", rp.iso_params.mmin);
            rp.iso_params.amin = get(jiso, "amin", rp.iso_params.amin);
            rp.iso_params.t_keep = get(jiso, "t_keep", rp.iso_params.t_keep);
            rp.iso_params.amb_lo = get(jiso, "amb_lo", rp.iso_params.amb_lo);
            rp.iso_params.amb_hi = get(jiso, "amb_hi", rp.iso_params.amb_hi);
        }
        rd.charge_scale = get(j, "charge_scale", rd.charge_scale);
        rd.uncer_cut = get(j, "uncer_cut", rd.uncer_cut);
        rd.bridge_len_cost = get(j, "bridge_len_cost", rd.bridge_len_cost);
        rd.bridge_max = get(j, "bridge_max", rd.bridge_max);
        rp.steiner.len_cost = rd.bridge_len_cost;
        rp.steiner.max_bridge_len = rd.bridge_max;
        rd.dump_dir = get<std::string>(j, "dump_dir", "");
        // doc pdvd/131
        rd.mode = get<std::string>(j, "mode", rd.mode);
        if (rd.mode != "retile" && rd.mode != "gapfill") {
            raise<ValueError>("ImproveCluster_2 retile_deghost: unknown mode '%s' (retile, gapfill)", rd.mode.c_str());
        }
        rd.residual = get(j, "residual", rd.residual);
        rd.corridor_wires = std::max(0, get(j, "corridor_wires", rd.corridor_wires));
        rd.gap_min_slices = std::max(1, get(j, "gap_min_slices", rd.gap_min_slices));
        rd.gap_max = get(j, "gap_max", rd.gap_max);
        rd.reach = get(j, "reach", rd.reach);
        rd.levels.clear();
        for (const auto& jl : j["levels"]) {
            Aux::Cascade::LevelSpec lc;
            lc.width = get<int>(jl, "width", 0);
            lc.name = get<std::string>(jl, "forward", "");
            lc.superwire = get<int>(jl, "superwire", 1);
            lc.threshold = get<double>(jl, "threshold", 0.0);
            if (lc.name.empty()) raise<ValueError>("ImproveCluster_2 retile_deghost: every level needs a forward");
            if (rd.levels.empty() != (lc.width <= 0)) {
                raise<ValueError>("ImproveCluster_2 retile_deghost: level 0 is uncut (width 0), every later level has a width");
            }
            lc.forward = Factory::find_tn<ITensorForward>(lc.name);
            rd.levels.push_back(lc);
        }
        if (rd.levels.empty()) raise<ValueError>("ImproveCluster_2 retile_deghost: no levels");
        if (rp.steiner.budget <= 0) rp.steiner.budget = Aux::Cascade::default_repair_budget(rd.levels.back().threshold);
        std::string desc;
        for (const auto& lc : rd.levels) {
            desc += fmt::format(" [width={} k={} thr={:.4f} {}]", lc.width, lc.superwire, lc.threshold, lc.name);
        }
        SPDLOG_LOGGER_INFO(log, "retile_deghost on: filter={} scope={} window=[{:.3f}, {:.3f}) us levels:{} repair={} budget={:.4f} "
                           "bridge_len_cost={} bridge_max={:.1f} cm charge_scale={} nthreads={} dump_dir='{}'",
                           rd.filter, rd.scope, rd.beam_window_low / units::us, rd.beam_window_high / units::us, desc,
                           rp.repair, rp.steiner.budget, rd.bridge_len_cost, rd.bridge_max / units::cm, rd.charge_scale,
                           rp.nthreads, rd.dump_dir);
        if (rd.mode == "gapfill") {
            SPDLOG_LOGGER_INFO(log, "retile_deghost mode=gapfill (doc pdvd/131): residual={} corridor_wires={} gap_min_slices={} "
                               "gap_max={:.1f} cm reach={:.1f} cm", rd.residual, rd.corridor_wires, rd.gap_min_slices,
                               rd.gap_max / units::cm, rd.reach / units::cm);
        }
    }

    bool ImproveCluster_2::rd_in_scope(const Cluster& cluster) const
    {
        if (m_rd.scope == "all") return true;
        const Grouping* g = cluster.grouping();
        if (!g) return false;
        if (g != m_rd_pick_grouping) {
            // the beam bundle, as CheckBeamParticle / ClusteringBeamDeghost pick it; once per grouping
            std::vector<PR::BeamBundleRow> rows;
            for (const auto* c : g->children()) {
                if (!c->get_flag(Flags::main_cluster)) continue;
                rows.push_back({c->get_cluster_id(), c->get_scalar<int>("matched_flash_gid", -1),
                                c->get_cluster_t0(), c->get_length()});
            }
            std::set<int> gids;
            for (const auto& r : rows) {
                if (r.gid >= 0 && r.t0 >= m_rd.beam_window_low && r.t0 < m_rd.beam_window_high) gids.insert(r.gid);
            }
            std::vector<PR::BeamFlashInfo> flashes;
            for (int gid : gids) {
                auto fl = g->flash_by_gid(gid);
                flashes.push_back({gid, fl ? fl.value() : -1.0, static_cast<bool>(fl)});
            }
            const auto pick = PR::beam_particle_pick_bundle(rows, flashes, m_rd.beam_window_low, m_rd.beam_window_high);
            m_rd_pick_grouping = g;
            m_rd_pick_gid = pick.gid;
            SPDLOG_LOGGER_INFO(log, "retile_deghost: {} main(s), {} in window over {} bundle(s); beam bundle gid {} ({})",
                               rows.size(), pick.n_in_window, pick.n_gids, pick.gid, pick.why);
        }
        return m_rd_pick_gid >= 0 && cluster.get_scalar<int>("matched_flash_gid", -1) == m_rd_pick_gid;
    }

    Configuration ImproveCluster_2::default_configuration() const
    {
        Configuration cfg = ImproveCluster_1::default_configuration();

        cfg["terminal_charge_threshold"] = m_steiner_terminal_charge;
        cfg["retile_mode"] = m_retile_mode;
        cfg["dump_dir"] = m_dump_dir;

        return cfg;
    }

    std::unique_ptr<ImproveCluster_2::node_t> ImproveCluster_2::mutate(node_t& node) const
    {
        using Clock = std::chrono::steady_clock;
        using MS = std::chrono::duration<double, std::milli>;
        auto t_mutate_start = Clock::now();
        auto t0 = Clock::now();

        // get the original cluster
        auto* orig_cluster = reinitialize(node);
        SPDLOG_LOGGER_TRACE(log, "timing: reinitialize took {} ms", MS(Clock::now()-t0).count());

        // First: get the shortest path from the original cluster
        // Create a SteinerGrapher instance with the cluster
        // You'll need to provide appropriate configuration
        Steiner::Grapher::Config grapher_config;
        // Configure as needed - you may need to access member variables
        // that provide detector volumes and point cloud transform sets
        grapher_config.dv = m_dv;     // From NeedDV mixin
        grapher_config.pcts = m_pcts; // From NeedPCTS mixin
        grapher_config.perf = m_verbose; // toggle timing printouts with verbose=true
        // doc pdvd/31 round 6 (owner Q5): the one place this stage's terminal
        // threshold is set.  Both graphers below share this config object, so
        // the two internal find_steiner_terminals passes stay in step.
        grapher_config.terminal_charge_threshold = m_steiner_terminal_charge;

        // Create the Steiner::Grapher instance
        t0 = Clock::now();
        Steiner::Grapher orig_steiner_grapher(*orig_cluster, grapher_config, log);
        auto& orig_graph = orig_steiner_grapher.get_graph("basic_pid"); // this is good for the original cluster
        SPDLOG_LOGGER_TRACE(log, "timing: get_graph(basic_pid) took {} ms", MS(Clock::now()-t0).count());
        SPDLOG_LOGGER_TRACE(log, "Orig Graph vertices: {}, edges: {}", boost::num_vertices(orig_graph), boost::num_edges(orig_graph));

        // Establish same blob steiner edges
        t0 = Clock::now();
        orig_steiner_grapher.establish_same_blob_steiner_edges("basic_pid", true);
        std::vector<size_t> orig_path_point_indices;
        {  
            auto pair_points = orig_cluster->get_two_boundary_wcps();
            auto first_index  =   orig_cluster->get_closest_point_index(pair_points.first);
            auto second_index =   orig_cluster->get_closest_point_index(pair_points.second);
            orig_path_point_indices = orig_cluster->graph_algorithms("basic_pid").shortest_path(first_index, second_index);
        }
        SPDLOG_LOGGER_TRACE(log, "timing: establish+shortest_path(basic_pid) took {} ms", MS(Clock::now()-t0).count());
        SPDLOG_LOGGER_TRACE(log, "Orig shortest path indices: {} ; Graph vertices: {}, edges: {}", orig_path_point_indices.size(), boost::num_vertices(orig_graph), boost::num_edges(orig_graph));
        
        t0 = Clock::now();
        orig_steiner_grapher.remove_same_blob_steiner_edges("basic_pid");
        SPDLOG_LOGGER_TRACE(log, "timing: remove_same_blob_steiner_edges(basic_pid) took {} ms", MS(Clock::now()-t0).count());
        SPDLOG_LOGGER_TRACE(log, "Orig Graph vertices: {}, edges: {}", boost::num_vertices(orig_graph), boost::num_edges(orig_graph));

        // doc pdvd/113: everything below this point is the retile proper.
        if (m_retile_mode != "full") {
            return mutate_reduced(*orig_cluster);
        }

        // Second, make a temp_cluster based on the original cluster via ImproveCluster_1
        SPDLOG_LOGGER_TRACE(log, "Grouping {} {}", m_grouping->get_name(), m_grouping->children().size());

        t0 = Clock::now();
        auto temp_node = ImproveCluster_1::mutate(node);
        auto temp_cluster_1 = temp_node->value.facade<Cluster>();
        auto& temp_cluster = m_grouping->make_child();
        temp_cluster.take_children(*temp_cluster_1);  // Move all blobs from improved cluster
        temp_cluster.from(*orig_cluster);
        // doc pdvd/25 M3: a retile can yield NO points (PDVD run 039252 evt 4,
        // a protect_bundle fragment); everything below indexes point 0 and
        // threw std::out_of_range.  Hand the (empty) ImproveCluster_1 node
        // back unchanged; the caller skips such a cluster.
        if (temp_cluster.npoints() == 0) {
            SPDLOG_LOGGER_WARN(log, "ImproveCluster_2: retiled cluster {} has no points; skipping the improvement", orig_cluster->ident());
            auto* temp_cluster_ptr0 = &temp_cluster;
            m_grouping->destroy_child(temp_cluster_ptr0, true);
            return temp_node;
        }
        SPDLOG_LOGGER_TRACE(log, "timing: ImproveCluster_1::mutate took {} ms", MS(Clock::now()-t0).count());
        SPDLOG_LOGGER_TRACE(log, "Grouping {} {}", m_grouping->get_name(), m_grouping->children().size());

        t0 = Clock::now();
        Steiner::Grapher temp_steiner_grapher(temp_cluster, grapher_config, log);
        
        // this requires CTPC and ref_point cloud of original cluster
        auto& temp_graph = temp_cluster.find_graph("ctpc_ref_pid", *orig_cluster, m_dv, m_pcts);
        SPDLOG_LOGGER_TRACE(log, "timing: Grapher+find_graph(ctpc_ref_pid) took {} ms", MS(Clock::now()-t0).count());
        // (temp_steiner_grapher used below for establish/remove_same_blob_steiner_edges)


        SPDLOG_LOGGER_TRACE(log, "Temp Graph vertices: {}, edges: {}", boost::num_vertices(temp_graph), boost::num_edges(temp_graph));
        t0 = Clock::now();
        temp_steiner_grapher.establish_same_blob_steiner_edges("ctpc_ref_pid", false);
        std::vector<size_t> temp_path_point_indices;
        {  
            auto pair_points = temp_cluster.get_two_boundary_wcps();
            auto first_index  =   temp_cluster.get_closest_point_index(pair_points.first);
            auto second_index =   temp_cluster.get_closest_point_index(pair_points.second);
            temp_path_point_indices = temp_cluster.graph_algorithms("ctpc_ref_pid").shortest_path(first_index, second_index);
        }
        SPDLOG_LOGGER_TRACE(log, "timing: establish+shortest_path(ctpc_ref_pid) took {} ms", MS(Clock::now()-t0).count());
        SPDLOG_LOGGER_TRACE(log, "Temp shortest path indices: {} ; Graph vertices: {}, edges: {}", temp_path_point_indices.size(), boost::num_vertices(temp_graph), boost::num_edges(temp_graph));
        t0 = Clock::now();
        temp_steiner_grapher.remove_same_blob_steiner_edges("ctpc_ref_pid");
        SPDLOG_LOGGER_TRACE(log, "timing: remove_same_blob_steiner_edges(ctpc_ref_pid) took {} ms", MS(Clock::now()-t0).count());
        SPDLOG_LOGGER_TRACE(log, "Temp Graph vertices: {}, edges: {}", boost::num_vertices(temp_graph), boost::num_edges(temp_graph));


        // star to construct a new cluster
        const auto wpid_set = orig_cluster->wpids_blob_set();

        // make a new node from the existing grouping
        auto& new_cluster = m_grouping->make_child(); // make a new cluster inside 

        // doc pdvd/129: study dump, off unless dump_dir is set
        const bool dumping = !m_dump_dir.empty();
        std::map<int, std::vector<RetileDump::Face>> dump_faces;   // apa -> faces

        // doc pdvd/130 phase 3: the deghosting step, for this cluster or not
        const bool rd_cluster = m_rd.on && rd_in_scope(*orig_cluster);
        const bool rd_dump = rd_cluster && !m_rd.dump_dir.empty();
        const bool rd_gap = rd_cluster && m_rd.mode == "gapfill";   // doc pdvd/131
        struct RdFace {
            int face{0};
            double tick{0};
            std::vector<double> angles;
            size_t ntiled{0}, nsampled{0};
            size_t nbase{0}, ngaps{0}, ncand{0}, nskipped{0}, nuncovered{0}, npath{0}, nowner{0};   // gapfill
            RetileDump::Face dump;   // iblobs = this face's model input
        };
        struct RdApa {
            std::vector<const Blob*> blobs;        // the survivors (to be replaced by the kept cells)
            std::vector<IBlob::pointer> in;        // their tiled blobs, face by face, in tiling order
            std::vector<IBlob::pointer> base;      // gapfill: the fixed base cells (a prefix-free subset of `in`)
            std::vector<RdFace> faces;
        };
        std::map<int, RdApa> rd_apa;
        const int rd_call = rd_cluster ? m_rd_call++ : -1;

        for (auto it = wpid_set.begin(); it != wpid_set.end(); ++it) {
            int apa = it->apa();
            int face = it->face();
            const auto& angles = m_wpid_angles.at(*it);

            if (rd_gap) {
                // doc pdvd/131 gapfill: this face gets no extension, no painting and no filter 1; the base
                // cells and the corridor candidates go to the per-anode model below.
                t0 = Clock::now();
                auto gf = rd_gapfill_face(*orig_cluster, orig_path_point_indices, apa, face);
                auto& ra = rd_apa[apa];
                ra.faces.emplace_back();
                auto& rf = ra.faces.back();
                rf.face = face;
                rf.tick = m_grouping->get_tick().at(apa).at(face);
                rf.angles = angles;
                rf.dump.face = face;
                rf.dump.tick_span = m_grouping->get_nticks_per_slice().at(apa).at(face);
                rf.nbase = gf.base.size(); rf.ngaps = gf.ngaps; rf.ncand = gf.cand.size(); rf.nskipped = gf.nskipped;
                rf.nuncovered = gf.nuncovered; rf.npath = gf.npath; rf.nowner = gf.nowner;
                rf.ntiled = rf.nsampled = gf.base.size() + gf.cand.size();
                for (const auto& b : gf.base) { ra.in.push_back(b); ra.base.push_back(b); rf.dump.iblobs.push_back(b); }
                for (const auto& b : gf.cand) { ra.in.push_back(b); rf.dump.iblobs.push_back(b); }
                if (rd_dump) {
                    for (const Blob* b : orig_cluster->children()) {
                        if (b->wpid().apa() != apa || b->wpid().face() != face) continue;
                        rf.dump.orig.push_back({b->slice_index_min(), b->u_wire_index_min(), b->u_wire_index_max(),
                                                b->v_wire_index_min(), b->v_wire_index_max(),
                                                b->w_wire_index_min(), b->w_wire_index_max()});
                    }
                }
                SPDLOG_LOGGER_TRACE(log, "timing: rd_gapfill_face (apa={},face={}) took {} ms", apa, face, MS(Clock::now()-t0).count());
                continue;
            }

            std::map<std::pair<int, int>, std::vector<WRG::measure_t> > map_slices_measures;
            
            // get original activities ...
            t0 = Clock::now();
            get_activity_improved(*orig_cluster, map_slices_measures, apa, face);
            std::set<RetileDump::cell_t> dump_before;
            if (dumping || rd_dump) dump_before = RetileDump::sentinel_cells(map_slices_measures);
            SPDLOG_LOGGER_TRACE(log, "timing: get_activity_improved (apa={},face={}) took {} ms", apa, face, MS(Clock::now()-t0).count());

            // doc pdhd/08 stage 1.  One cell-provenance set shared by the two
            // calls below: they paint into the SAME map_slices_measures, so
            // without this the temp pass reports a bridge retracing the orig
            // pass's ghost as "nothing new" and the ghost reads as harmless.
            // Allocated only when the census is on (it can hold millions of
            // cells); nullptr in production.
            bridge_cells_t bridge_cells;
            bridge_cells_t* bcp = m_bad_blob_report ? &bridge_cells : nullptr;

            // hack activity according to original cluster
            t0 = Clock::now();
            hack_activity_improved(*orig_cluster, map_slices_measures, orig_path_point_indices, apa, face, "orig", bcp); // may need more args
            SPDLOG_LOGGER_TRACE(log, "timing: hack_activity(orig) (apa={},face={}) took {} ms", apa, face, MS(Clock::now()-t0).count());

            // hack activities according to the new cluster
            t0 = Clock::now();
            hack_activity_improved(temp_cluster, map_slices_measures, temp_path_point_indices, apa, face, "temp", bcp); // may need more args
            SPDLOG_LOGGER_TRACE(log, "timing: hack_activity(temp) (apa={},face={}) took {} ms", apa, face, MS(Clock::now()-t0).count());

            // Step 3.
            t0 = Clock::now();
            auto iblobs = make_iblobs_improved(map_slices_measures, apa, face);
            SPDLOG_LOGGER_TRACE(log, "timing: make_iblobs_improved (apa={},face={}) took {} ms", apa, face, MS(Clock::now()-t0).count());
            SPDLOG_LOGGER_TRACE(log, "new cluster {} iblobs for apa {} face {}", iblobs.size(), apa, face);

            auto niblobs = iblobs.size();
            RetileDump::Face* dface = nullptr;
            if (dumping) {
                auto& dfs = dump_faces[apa];
                dfs.emplace_back();
                dface = &dfs.back();
                dface->face = face;
                dface->iblobs = iblobs;
                dface->sampled.assign(niblobs, 0);
                dface->cells = RetileDump::classify(dump_before, map_slices_measures);
                if (!map_slices_measures.empty()) {
                    dface->tick_span = map_slices_measures.begin()->first.second - map_slices_measures.begin()->first.first;
                }
                for (const Blob* b : orig_cluster->children()) {
                    if (b->wpid().apa() != apa || b->wpid().face() != face) continue;
                    dface->orig.push_back({b->slice_index_min(), b->u_wire_index_min(), b->u_wire_index_max(),
                                           b->v_wire_index_min(), b->v_wire_index_max(),
                                           b->w_wire_index_min(), b->w_wire_index_max()});
                }
            }
            // start to sampling points 
            int npoints = 0;
            t0 = Clock::now();
            std::map<const Blob*, int> rd_bind;   // doc pdvd/130: sampled blob -> iblob index (lookup only)
            for (size_t bind=0; bind<niblobs; ++bind) {
          
                const IBlob::pointer iblob = iblobs[bind];
                auto sampler = m_samplers.at(apa).at(face);
                const double tick = m_grouping->get_tick().at(apa).at(face);

                auto pcs = Aux::sample_live(sampler, iblob, angles, tick, bind);
                // DO NOT EXTEND FURTHER! see #426, #430

                if (pcs["3d"].size()==0) continue; // no points ...
                if (dface) dface->sampled[bind] = 1;
                // Access 3D coordinates
                auto pc3d = pcs["3d"];  // Get the 3D point cloud dataset
                auto x_coords = pc3d.get("x")->elements<double>();  // Get X coordinates
                // auto y_coords = pc3d.get("y")->elements<double>();  // Get Y coordinates  
                // auto z_coords = pc3d.get("z")->elements<double>();  // Get Z coordinates
                // auto ucharge_val = pc3d.get("ucharge_val")->elements<double>();  // Get U charge
                // auto vcharge_val = pc3d.get("vcharge_val")->elements<double>();  // Get V charge
                // auto wcharge_val = pc3d.get("wcharge_val")->elements<double>();  // Get W charge
                // auto ucharge_err = pc3d.get("ucharge_unc")->elements<double>();  // Get U charge error
                // auto vcharge_err = pc3d.get("vcharge_unc")->elements<double>();  // Get V charge error
                // auto wcharge_err = pc3d.get("wcharge_unc")->elements<double>();  // Get W charge error

                // std::cout << "ImproveCluster_1 PCS: " << pcs.size() << " " 
                //           << pcs["3d"].size() << " " 
                //           << x_coords.size() << std::endl;     

                npoints +=x_coords.size();
                if (pcs.empty()) {
                    SPDLOG_LOGGER_TRACE(log, "skipping blob {} with no points", iblob->ident());
                    continue;
                }
                auto* bnode = new_cluster.node()->insert(Tree::Points(std::move(pcs)));
                if (rd_cluster) rd_bind[bnode->value.facade<Blob>()] = (int) bind;
            }
            SPDLOG_LOGGER_TRACE(log, "timing: sample_live loop (apa={},face={}) took {} ms", apa, face, MS(Clock::now()-t0).count());
            SPDLOG_LOGGER_TRACE(log, "{} points sampled for apa {} face {} Blobs {}", npoints, apa, face, niblobs);

            // remove bad blobs
            t0 = Clock::now();
            if (map_slices_measures.empty()) continue; // no tiled blobs for this face
            int tick_span = map_slices_measures.begin()->first.second -  map_slices_measures.begin()->first.first;
            std::vector<const Blob*> blobs_to_remove;
            // doc pdvd/130: a blob with a plane whose every wire in its range carries the sentinel (dead or painted)
            const double rd_tick = m_grouping->get_tick().at(apa).at(face);
            auto rd_eligible = [&](const IBlob::pointer& ib) {
                const auto& sl = ib->slice();
                const int ta = (int) std::lround(sl->start() / rd_tick);
                const int tb = (int) std::lround((sl->start() + sl->span()) / rd_tick);
                auto it = map_slices_measures.find({ta, tb});
                if (it == map_slices_measures.end()) return false;
                for (const auto& strip : ib->shape().strips()) {
                    if (strip.layer < 2 || strip.layer >= (int) it->second.size()) continue;
                    const auto& m = it->second[strip.layer];
                    bool any = false, all = true;
                    for (int w = strip.bounds.first; w < strip.bounds.second && w < (int) m.size(); ++w) {
                        any = true;
                        if (!RetileDump::is_sentinel(m[w])) { all = false; break; }
                    }
                    if (any && all) return true;
                }
                return false;
            };
            if (!rd_cluster) {
                blobs_to_remove = remove_bad_blobs(*orig_cluster, new_cluster, tick_span, apa, face);
            }
            else if (m_rd.filter != "none") {
                BadBlobOpts opts;
                opts.max_run = (m_rd.filter == "vote") ? 0.0 : m_bad_blob_max_run;
                opts.run_merge = (m_rd.filter == "run") ? m_bad_blob_run_merge : 0.0;
                std::map<const Blob*, bool> elig;
                if (m_rd.filter == "painted_run") {
                    for (const Blob* b : new_cluster.children()) {
                        auto it = rd_bind.find(b);
                        if (it != rd_bind.end()) elig[b] = rd_eligible(iblobs[it->second]);
                    }
                    opts.run_eligible = &elig;
                }
                blobs_to_remove = remove_bad_blobs(*orig_cluster, new_cluster, tick_span, apa, face, &opts);
            }
            // doc pdvd/40 round 3 census: point count before/after the removal
            // (bad_blob_report only) -- the proof that a removed blob's points
            // leave the retiled cluster the Steiner build sees.
            // Counted over the blob children, NOT via Cluster::npoints(): that
            // memo lives in the ClusterCache, which child insert/remove does
            // not invalidate (see remove_bad_blobs), so it would read stale.
            auto count_pts = [&]() { size_t n = 0; for (const Blob* b : new_cluster.children()) n += b->nbpoints(); return n; };
            const size_t npts_before = m_bad_blob_report ? count_pts() : 0;
            for (const Blob* blob : blobs_to_remove) {
                if (dface) {
                    dface->removed.push_back({blob->slice_index_min(), blob->u_wire_index_min(), blob->u_wire_index_max(),
                                              blob->v_wire_index_min(), blob->v_wire_index_max(),
                                              blob->w_wire_index_min(), blob->w_wire_index_max()});
                }
                Blob& b = const_cast<Blob&>(*blob);
                new_cluster.remove_child(b);
            }
            if (m_bad_blob_report) {
                SPDLOG_LOGGER_DEBUG(log, "BADBLOBRM ident={} apa={} face={} removed={} npts {} -> {} nblobs {}",
                                    orig_cluster->ident(), apa, face, blobs_to_remove.size(), npts_before,
                                    count_pts(), new_cluster.children().size());
            }
            SPDLOG_LOGGER_TRACE(log, "timing: remove_bad_blobs (apa={},face={}) took {} ms", apa, face, MS(Clock::now()-t0).count());
            SPDLOG_LOGGER_TRACE(log, "{} blobs removed for apa {} face {} remaining {}", blobs_to_remove.size(), apa, face, new_cluster.children().size());

            // doc pdvd/130 phase 3: collect this face's survivors for the model, which runs once per anode with both
            // faces in one graph (as the imaging node and the training do)
            if (rd_cluster) {
                auto& ra = rd_apa[apa];
                std::vector<int> binds;
                for (const Blob* b : new_cluster.children()) {
                    if (b->wpid().apa() != apa || b->wpid().face() != face) continue;
                    auto it = rd_bind.find(b);
                    if (it == rd_bind.end()) continue;
                    ra.blobs.push_back(b);
                    binds.push_back(it->second);
                }
                std::sort(binds.begin(), binds.end());   // tiling order: slice time, then blob
                ra.faces.emplace_back();
                auto& rf = ra.faces.back();
                rf.face = face;
                rf.tick = rd_tick;
                rf.angles = angles;
                rf.dump.face = face;
                rf.dump.tick_span = tick_span;
                for (int k : binds) {
                    ra.in.push_back(iblobs[k]);
                    rf.dump.iblobs.push_back(iblobs[k]);
                }
                rf.nsampled = rd_bind.size();
                rf.ntiled = niblobs;
                if (rd_dump) {
                    rf.dump.cells = RetileDump::classify(dump_before, map_slices_measures);
                    for (const Blob* b : orig_cluster->children()) {
                        if (b->wpid().apa() != apa || b->wpid().face() != face) continue;
                        rf.dump.orig.push_back({b->slice_index_min(), b->u_wire_index_min(), b->u_wire_index_max(),
                                                b->v_wire_index_min(), b->v_wire_index_max(),
                                                b->w_wire_index_min(), b->w_wire_index_max()});
                    }
                }
            }
        }

        if (dumping) {
            const int serial = m_dump_serial++;
            for (const auto& [dapa, dfs] : dump_faces) {
                RetileDump::write(m_dump_dir, orig_cluster->ident(), serial, dapa, dfs);
            }
        }
        // doc pdvd/130 phase 3: filter 2 per anode, the cascade model on the survivors of both faces, then the kept
        // cells sampled in their place
        if (rd_cluster) {
            const auto drift_speeds = m_grouping->get_drift_speed();
            for (auto& [apa, ra] : rd_apa) {   // int-keyed: ascending
                if (ra.in.empty()) continue;
                auto t_rd = Clock::now();
                std::vector<ISlice::pointer> slices;
                std::set<const ISlice*> seen;   // dedup only
                for (const auto& b : ra.in) if (seen.insert(b->slice().get()).second) slices.push_back(b->slice());
                const auto sc = Aux::Cascade::make_slice_charge(slices, m_rd.charge_scale, m_rd.uncer_cut);
                auto rp = m_rd.rp;
                const auto& ds = drift_speeds.at(apa);
                auto rd_center = [&](const IBlob::pointer& b) {
                    const double v = ds.at(b->face()->which());
                    Vector c(b->slice()->start() * v, 0, 0);
                    const auto& corners = b->shape().corners();
                    if (corners.empty()) return c;
                    const auto& coords = b->face()->raygrid();
                    Vector yz(0, 0, 0);
                    for (const auto& cr : corners) yz += coords.ray_crossing(cr.first, cr.second);
                    yz = yz * (1.0 / corners.size());
                    return Vector(c.x(), yz.y(), yz.z());
                };
                if (m_rd.bridge_len_cost != 0.0 || m_rd.bridge_max > 0.0) {
                    rp.edge_length = [&](const IBlob::pointer& a, const IBlob::pointer& b) {
                        return (rd_center(a) - rd_center(b)).magnitude();
                    };
                }
                const auto res = Aux::Cascade::run_cascade(ra.in, sc, m_rd.levels, rp, orig_cluster->ident(), log, rd_call);
                for (const Blob* blob : ra.blobs) {
                    Blob& b = const_cast<Blob&>(*blob);
                    new_cluster.remove_child(b);
                }
                // doc pdvd/131 gapfill: the base is fixed.  Every base cell is sampled whatever the model said of it;
                // of the model's kept cells only those that are not base cells (the corridor cells) are added.
                using rd_key_t = std::array<int, 8>;
                auto rd_key = [&](const IBlob::pointer& b) {
                    rd_key_t k{};
                    k[0] = b->face()->which();
                    const double tk = m_grouping->get_tick().at(apa).at(b->face()->which());
                    k[1] = (int) std::lround(b->slice()->start() / tk);
                    for (const auto& strip : b->shape().strips()) {
                        if (strip.layer < 2 || strip.layer > 4) continue;
                        k[2 + 2 * (strip.layer - 2)] = strip.bounds.first;
                        k[3 + 2 * (strip.layer - 2)] = strip.bounds.second;
                    }
                    return k;
                };
                std::set<rd_key_t> base_keys;
                if (rd_gap) for (const auto& b : ra.base) base_keys.insert(rd_key(b));
                std::vector<IBlob::pointer> final_cells;   // what is sampled, in order
                if (rd_gap) final_cells = ra.base;
                size_t ncand_kept = 0;
                for (const auto& cell : res.kept) {
                    if (rd_gap && base_keys.count(rd_key(cell))) continue;
                    final_cells.push_back(cell);
                    ++ncand_kept;
                }
                size_t nsampled = 0;
                std::map<int, size_t> kept_by_face;
                for (size_t k = 0; k < final_cells.size(); ++k) {
                    const auto& cell = final_cells[k];
                    const int face = cell->face()->which();
                    ++kept_by_face[face];
                    const RdFace* rf = nullptr;
                    for (const auto& f : ra.faces) if (f.face == face) rf = &f;
                    if (!rf) continue;
                    auto sampler = m_samplers.at(apa).at(face);
                    auto pcs = Aux::sample_live(sampler, cell, rf->angles, rf->tick, (int) k);
                    if (pcs["3d"].size() == 0) continue;
                    ++nsampled;
                    new_cluster.node()->insert(Tree::Points(std::move(pcs)));
                }
                std::string per_face;
                for (const auto& f : ra.faces) {
                    per_face += fmt::format(" [face {}: tiled {} sampled {} after filter 1 {} kept {}]", f.face, f.ntiled, f.nsampled,
                                            f.dump.iblobs.size(), kept_by_face[f.face]);
                }
                if (rd_gap) {
                    std::string gf_face;
                    for (const auto& f : ra.faces) {
                        gf_face += fmt::format(" [face {}: base {} path {} uncovered {} gaps {} candidates {} skipped {} owner {}]",
                                               f.face, f.nbase, f.npath, f.nuncovered, f.ngaps, f.ncand, f.nskipped, f.nowner);
                    }
                    SPDLOG_LOGGER_INFO(log, "retile_gapfill ident={} apa={}: base {} + candidates {} -> candidates kept {} (model kept {}, "
                                       "threshold {}, weak dropped {}, bridges {}, added {}) sampled {}{} ({:.1f} ms)",
                                       orig_cluster->ident(), apa, ra.base.size(), ra.in.size() - ra.base.size(), ncand_kept,
                                       res.kept.size(), res.nkeep_thr, res.srep.nweak_cells, res.srep.nbridges, res.srep.nadded,
                                       nsampled, gf_face, MS(Clock::now() - t_rd).count());
                }
                else {
                SPDLOG_LOGGER_INFO(log, "retile_deghost ident={} apa={} filter={}: model in {} -> kept {} (threshold {}, weak dropped {}, "
                                   "bridges {}, added {}) sampled {}{} ({:.1f} ms)",
                                   orig_cluster->ident(), apa, m_rd.filter, ra.in.size(), res.kept.size(), res.nkeep_thr,
                                   res.srep.nweak_cells, res.srep.nbridges, res.srep.nadded, nsampled, per_face,
                                   MS(Clock::now() - t_rd).count());
                }
                if (rd_dump) {
                    std::vector<RetileDump::Face> fin, fout;
                    for (const auto& f : ra.faces) {
                        fin.push_back(f.dump);
                        fin.back().sampled.assign(f.dump.iblobs.size(), 1);
                        // gapfill: the base cells come first (flag 1), then the corridor candidates (flag 2)
                        if (rd_gap) for (size_t i = f.nbase; i < fin.back().sampled.size(); ++i) fin.back().sampled[i] = 2;
                        fout.push_back(f.dump);
                        fout.back().iblobs.clear();
                        for (const auto& cell : final_cells) if (cell->face()->which() == f.face) fout.back().iblobs.push_back(cell);
                        fout.back().sampled.assign(fout.back().iblobs.size(), 1);
                    }
                    RetileDump::write(m_rd.dump_dir, orig_cluster->ident(), rd_call, apa, fin, "rdin");
                    RetileDump::write(m_rd.dump_dir, orig_cluster->ident(), rd_call, apa, fout, "rdout");
                }
            }
            orig_cluster->set_scalar<int>("retile_deghost_writeback", 1);   // CreateSteinerGraph replaces the blobs
        }

        // Remove this cluster from the grouping
        t0 = Clock::now();
        auto* temp_cluster_ptr = &temp_cluster;
        m_grouping->destroy_child(temp_cluster_ptr, true);
        SPDLOG_LOGGER_TRACE(log, "timing: destroy_child(temp_cluster) took {} ms", MS(Clock::now()-t0).count());
        SPDLOG_LOGGER_TRACE(log, "Grouping {} {}", m_grouping->get_name(), m_grouping->children().size());

        auto& default_scope = orig_cluster->get_default_scope();
        auto& raw_scope = orig_cluster->get_raw_scope();

        SPDLOG_LOGGER_TRACE(log, "Scope: {} {}", default_scope.hash(), raw_scope.hash());
        if (default_scope.hash()!=raw_scope.hash()){
            t0 = Clock::now();
            auto correction_name = orig_cluster->get_scope_transform(default_scope);
            // add_corrected_points builds x_t0cor from get_cluster_t0(), but a
            // freshly-retiled cluster's T0 defaults to 0 and was only copied by
            // from() AFTERWARD -- so the correction applied a ZERO drift shift,
            // leaving x_t0cor equal to the raw drift-x.  Harmless for in-time
            // clusters, but for an out-of-time cosmic it is off by v_drift*|T0|
            // (~164 cm at T0 = -1051 us), which pushed steiner boundary points
            // past the anode and gave TaggerCheckFC/STM/Neutrino spurious
            // out-of-fiducial "exit" verdicts.  Set the real T0 first.
            new_cluster.set_cluster_t0(orig_cluster->get_cluster_t0());
            new_cluster.add_corrected_points(m_pcts, correction_name);
            new_cluster.from(*orig_cluster); // copy remaining state from original cluster
            SPDLOG_LOGGER_TRACE(log, "timing: add_corrected_points took {} ms", MS(Clock::now()-t0).count());
        }

        // auto retiled_node = new_cluster.node();

        SPDLOG_LOGGER_TRACE(log, "timing: mutate() TOTAL took {} ms", MS(Clock::now()-t_mutate_start).count());
        return m_grouping->remove_child(new_cluster);

    }

    // doc pdvd/131.  One face of the "gapfill" mode: the base (the cluster's own blobs as cells), the owner map
    // of the region, the path's coverage and gaps, the corridor candidates tiled at cell resolution.
    ImproveCluster_2::GapFace ImproveCluster_2::rd_gapfill_face(const Cluster& orig, const std::vector<size_t>& path,
                                                                int apa, int face) const
    {
        namespace RL = WireCell::Clus::ResampleLive;
        namespace GF = WireCell::Clus::GapFill;
        GapFace out;
        const auto& iface = m_face.at(apa).at(face);
        const auto& ianode = m_anode.at(apa);
        const double tick = m_grouping->get_tick().at(apa).at(face);
        const int tick_span = m_grouping->get_nticks_per_slice().at(apa).at(face);
        const double vdrift = m_grouping->get_drift_speed().at(apa).at(face);
        const auto& coords = iface->raygrid();
        const int which = iface->which();
        const auto& pinfo = m_plane_infos.at(apa).at(face);
        std::array<int, 3> nw{};
        for (int p = 0; p < 3; ++p) nw[p] = pinfo[p].total_wires;

        // The base: the cluster's own blobs on this face, as index boxes, in children order.
        std::vector<GF::Box> boxes;
        std::vector<const Blob*> bblobs;
        int smin_all = INT_MAX, smax_all = INT_MIN;
        std::array<int, 3> wlo{INT_MAX, INT_MAX, INT_MAX}, whi{INT_MIN, INT_MIN, INT_MIN};
        for (const Blob* b : orig.children()) {
            if (b->wpid().apa() != apa || b->wpid().face() != face) continue;
            GF::Box bx;
            bx.smin = b->slice_index_min();
            bx.smax = b->slice_index_max();
            bx.w = {std::make_pair(b->u_wire_index_min(), b->u_wire_index_max()),
                    std::make_pair(b->v_wire_index_min(), b->v_wire_index_max()),
                    std::make_pair(b->w_wire_index_min(), b->w_wire_index_max())};
            boxes.push_back(bx);
            bblobs.push_back(b);
            smin_all = std::min(smin_all, bx.smin);
            smax_all = std::max(smax_all, bx.smax);
            for (int p = 0; p < 3; ++p) { wlo[p] = std::min(wlo[p], bx.w[p].first); whi[p] = std::max(whi[p], bx.w[p].second); }
        }
        if (boxes.empty()) return out;
        std::map<int, std::vector<size_t>> boxes_by_slice;   // smin -> box indices
        for (size_t i = 0; i < boxes.size(); ++i) boxes_by_slice[boxes[i].smin].push_back(i);

        // The region: the base's bounds widened by the reach (ticks along the drift, wires per plane).
        const int reach_ticks = (vdrift > 0 && tick > 0) ? (int) std::ceil(m_rd.reach / (vdrift * tick)) : 0;
        const int rs0 = smin_all - reach_ticks, rs1 = smax_all + reach_ticks;
        std::array<int, 3> rlo{}, rhi{};
        for (int p = 0; p < 3; ++p) {
            const double y0 = m_grouping->convert_time_wire_2Dpoint(smin_all, 0, apa, face, p).second;
            const double y1 = m_grouping->convert_time_wire_2Dpoint(smin_all, 1, apa, face, p).second;
            const double pitch = std::abs(y1 - y0);
            const int rw = pitch > 0 ? (int) std::ceil(m_rd.reach / pitch) : 0;
            rlo[p] = std::max(0, wlo[p] - rw);
            rhi[p] = std::min(nw[p], whi[p] + rw);
        }

        // The owner map (doc 131 step B): other clusters' original blobs with a slice inside the region.  The
        // bundle-mates (same matched flash) are not "other".  Blobs are visited in a sorted order (the per-slice
        // sets are pointer-keyed).
        GF::OwnerMap owner;
        const int gid = (m_rd.scope == "bundle") ? m_rd_pick_gid : -1;
        for (const Cluster* c : m_grouping->children()) {
            if (c == &orig) continue;
            if (gid >= 0 && c->get_scalar<int>("matched_flash_gid", -1) == gid) continue;
            const auto& tbm = c->time_blob_map();
            auto ia = tbm.find(apa);
            if (ia == tbm.end()) continue;
            auto ifc = ia->second.find(face);
            if (ifc == ia->second.end()) continue;
            std::vector<const Blob*> obs;
            for (auto it = ifc->second.lower_bound(rs0); it != ifc->second.end() && it->first < rs1; ++it) {
                for (const Blob* b : it->second) obs.push_back(b);
            }
            std::sort(obs.begin(), obs.end(), [](const Blob* a, const Blob* b) {
                return std::make_tuple(a->slice_index_min(), a->u_wire_index_min(), a->v_wire_index_min(), a->w_wire_index_min(),
                                       a->u_wire_index_max(), a->v_wire_index_max(), a->w_wire_index_max())
                     < std::make_tuple(b->slice_index_min(), b->u_wire_index_min(), b->v_wire_index_min(), b->w_wire_index_min(),
                                       b->u_wire_index_max(), b->v_wire_index_max(), b->w_wire_index_max());
            });
            for (const Blob* b : obs) {
                const std::array<std::pair<int, int>, 3> r = {std::make_pair(b->u_wire_index_min(), b->u_wire_index_max()),
                                                              std::make_pair(b->v_wire_index_min(), b->v_wire_index_max()),
                                                              std::make_pair(b->w_wire_index_min(), b->w_wire_index_max())};
                bool inside = false;
                for (int p = 0; p < 3; ++p) if (r[p].second > rlo[p] && r[p].first < rhi[p]) inside = true;
                if (!inside) continue;
                const double q = b->charge();
                for (int s = b->slice_index_min(); s < b->slice_index_max(); s += tick_span) {
                    for (int p = 0; p < 3; ++p) owner.add(p, s, r[p].first, r[p].second, std::isfinite(q) ? q : 0.0);
                }
            }
        }
        out.nowner = owner.size();

        // The base cells as IBlobs: activity as retile_mode "none" (CTPC row, else dead registry, else absent;
        // bounds +- 2 wires), the live charge reduced to the residual (step D); one SimpleSlice per slice start.
        const int wire_margin = 2;
        std::map<int, ISlice::map_t> act;
        auto live_q = [&](int p, int s, int w, double& q, double& err) {
            const auto* row = m_grouping->wire_charge_row(apa, face, p, s);
            if (!row) return false;
            auto it = row->find(w);
            if (it == row->end()) return false;
            q = it->second.first;
            err = it->second.second;
            if (m_rd.residual) q = GF::residual(q, owner.other_share(p, s, w));
            return true;
        };
        for (const GF::Box& bx : boxes) {
            auto& a = act[bx.smin];
            for (int p = 0; p < 3; ++p) {
                const auto& wires = iface->planes()[p]->wires();
                const int lo = std::max(0, bx.w[p].first - wire_margin);
                const int hi = std::min((int) wires.size() - 1, bx.w[p].second + wire_margin);
                auto live = [&](int w, double& q, double& err) { return live_q(p, bx.smin, w, q, err); };
                auto dead = [&](int w) { return m_grouping->is_wire_dead(apa, face, p, w, bx.smin); };
                for (const auto& wv : RL::compose_activity(lo, hi, live, dead)) {
                    auto ich = ianode->channel(wires[wv.wire]->channel());
                    if (!ich) continue;
                    a[ich] = ISlice::value_t(wv.charge, wv.error);
                }
            }
        }
        std::map<int, ISlice::pointer> slice_of;
        for (const auto& [s, a] : act) {
            slice_of[s] = std::make_shared<Aux::SimpleSlice>(nullptr, s, s * tick, tick_span * tick, a);
        }
        int ident = 0;
        for (const GF::Box& bx : boxes) {
            const RL::plane_bounds_t bounds = {bx.w[0], bx.w[1], bx.w[2]};
            RayGrid::Blob shape = RL::shape_from_bounds(coords, bounds);
            out.base.push_back(std::make_shared<Aux::SimpleBlob>(ident++, 0.0f, 0.0f, shape, slice_of.at(bx.smin), iface));
        }

        // The path on this face, resampled at 0.3 cm as the painting does (hack_activity_improved), each sample
        // converted to the aligned slice start and the three wire indices.
        const double step = 0.3 * units::cm;
        std::vector<GF::PathPt> pts;
        std::vector<geo_point_t> ppts;
        std::vector<bool> pon;
        {
            bool have = false;
            geo_point_t prev;
            bool prev_on = false;
            for (size_t idx : path) {
                const geo_point_t pt = orig.point3d_raw(idx);
                const auto wp = orig.wire_plane_id(idx);
                const bool on = (wp.apa() == apa && wp.face() == face);
                if (!have) {
                    ppts.push_back(pt); pon.push_back(on);
                    have = true; prev = pt; prev_on = on;
                    continue;
                }
                const double dis = (pt - prev).magnitude();
                if (dis < step) {
                    ppts.push_back(pt); pon.push_back(on);
                }
                else {
                    const int ncount = (int) (dis / step) + 1;
                    for (int i = 0; i < ncount; ++i) {
                        ppts.push_back(prev + (pt - prev) * ((i + 1.0) / ncount));
                        pon.push_back(on && prev_on);
                    }
                    pon.back() = on;
                }
                prev = pt; prev_on = on;
            }
        }
        double s_acc = 0;
        for (size_t i = 0; i < ppts.size(); ++i) {
            if (i) s_acc += (ppts[i] - ppts[i - 1]).magnitude();
            GF::PathPt pp;
            pp.s = s_acc;
            pp.on_face = pon[i];
            if (pp.on_face) {
                int t0 = 0;
                for (int p = 0; p < 3; ++p) {
                    auto [tind, wind] = m_grouping->convert_3Dpoint_time_ch(ppts[i], apa, which, p);
                    if (p == 0) t0 = tind;
                    pp.wire[p] = wind;
                }
                pp.slice = (int) std::lround(t0 * 1.0 / tick_span) * tick_span;
            }
            pts.push_back(pp);
        }
        out.npath = pts.size();
        std::vector<bool> covered(pts.size(), false);
        for (size_t i = 0; i < pts.size(); ++i) {
            if (!pts[i].on_face) continue;
            // boxes whose slice start lies within one slice of the sample
            for (auto it = boxes_by_slice.lower_bound(pts[i].slice - 2 * tick_span); it != boxes_by_slice.end() && it->first <= pts[i].slice + tick_span; ++it) {
                for (size_t bi : it->second) {
                    if (GF::covers(boxes[bi], pts[i].slice, pts[i].wire, 1, tick_span)) { covered[i] = true; break; }
                }
                if (covered[i]) break;
            }
            if (!covered[i]) ++out.nuncovered;
        }
        const auto gaps = GF::find_gaps(pts, covered, m_rd.gap_min_slices, m_rd.gap_max);
        out.ngaps = gaps.size();
        if (gaps.empty()) return out;

        // The corridors (step C): around every sample of a gap, +- corridor_wires per plane in the sample's slice;
        // wires owned by another cluster are never admitted; live wires carry the residual charge; dead wires and,
        // in at most one plane, live wires without charge are bridge wires (the tiling sentinel).
        std::map<std::pair<int, int>, std::vector<WRG::measure_t>> msm;
        for (const auto& g : gaps) {
            for (size_t i = g.first; i <= g.last; ++i) {
                const auto& pt = pts[i];
                if (!pt.on_face) continue;
                std::array<int, 3> nlive{};
                std::array<std::vector<std::pair<int, double>>, 3> adm;   // (wire, value) per plane; value 1e-3 = sentinel
                for (int p = 0; p < 3; ++p) {
                    for (int w = pt.wire[p] - m_rd.corridor_wires; w <= pt.wire[p] + m_rd.corridor_wires; ++w) {
                        if (w < 0 || w >= nw[p]) continue;
                        if (owner.is_other(p, pt.slice, w)) continue;
                        double q = 0, err = 0;
                        if (live_q(p, pt.slice, w, q, err) && q > 0) {
                            adm[p].emplace_back(w, q);
                            ++nlive[p];
                        }
                        else {
                            adm[p].emplace_back(w, 1e-3);   // dead, or live without charge: a bridge wire
                        }
                    }
                }
                if (GF::missing_planes(nlive) >= 2) { ++out.nskipped; continue; }
                auto& measures = msm[{pt.slice, pt.slice + tick_span}];
                if (measures.empty()) {
                    measures.resize(5);
                    measures[0].push_back(1);
                    measures[1].push_back(1);
                    for (int p = 0; p < 3; ++p) measures[2 + p].resize(nw[p], 0);
                }
                for (int p = 0; p < 3; ++p) {
                    for (const auto& [w, v] : adm[p]) {
                        auto& m = measures[2 + p][w];
                        if (v == 1e-3) { if (m == 0) m = 1e-3; }   // a live value already there wins
                        else m = v;
                    }
                }
            }
        }
        auto cand = make_iblobs_improved(msm, apa, face);
        out.ncand_tiled = cand.size();
        // candidates inside a base cell (same slice, every range within) are not new
        for (const auto& cb : cand) {
            const int s = (int) std::lround(cb->slice()->start() / tick);
            std::array<std::pair<int, int>, 3> r{};
            for (const auto& strip : cb->shape().strips()) {
                if (strip.layer < 2 || strip.layer > 4) continue;
                r[strip.layer - 2] = strip.bounds;
            }
            bool inside = false;
            auto it = boxes_by_slice.find(s);
            if (it != boxes_by_slice.end()) {
                for (size_t bi : it->second) {
                    const auto& bx = boxes[bi];
                    bool in = true;
                    for (int p = 0; p < 3; ++p) if (r[p].first < bx.w[p].first || r[p].second > bx.w[p].second) in = false;
                    if (in) { inside = true; break; }
                }
            }
            if (!inside) out.cand.push_back(cb);
        }
        return out;
    }

    // doc pdvd/113.  The reduced retiles.  Each mirrors the "full" path's
    // bookkeeping (new child of the grouping, T0-corrected points, removed from
    // the grouping on return) so CreateSteinerGraph consumes the result exactly
    // as it consumes a full retile.
    std::unique_ptr<ImproveCluster_2::node_t> ImproveCluster_2::mutate_reduced(Cluster& orig_cluster) const
    {
        namespace RL = WireCell::Clus::ResampleLive;

        auto& new_cluster = m_grouping->make_child();

        if (m_retile_mode == "none") {
            // No tiling.  The blob loop is duplicated from ClusteringResampleLive::visit
            // (clustering_resample_live.cxx, untouched): shape from the blob's wire
            // bounds, activity over bounds +- 2 wires from the grouping's CTPC row
            // (live), else the dead-wind registry (dead), else absent.  Children are
            // taken in the source cluster's order.
            const int wire_margin = 2;
            const auto ticks = m_grouping->get_tick();
            int blob_counter = 0;
            for (const Blob* fblob : orig_cluster.children()) {
                const WirePlaneId wpid = fblob->wpid();
                const int apa = wpid.apa();
                const int face = wpid.face();
                const auto& iface = m_face.at(apa).at(face);
                const auto& ianode = m_anode.at(apa);
                const auto& sampler = m_samplers.at(apa).at(face);
                const double tick = ticks.at(apa).at(face);
                const int smin = fblob->slice_index_min();
                const int smax = fblob->slice_index_max();
                const RL::plane_bounds_t bounds = {
                    std::make_pair(fblob->u_wire_index_min(), fblob->u_wire_index_max()),
                    std::make_pair(fblob->v_wire_index_min(), fblob->v_wire_index_max()),
                    std::make_pair(fblob->w_wire_index_min(), fblob->w_wire_index_max()),
                };

                ISlice::map_t activity;
                for (int p = 0; p < 3; ++p) {
                    const auto& wires = iface->planes()[p]->wires();
                    const int nwires = (int) wires.size();
                    const int lo = std::max(0, bounds[p].first - wire_margin);
                    const int hi = std::min(nwires - 1, bounds[p].second + wire_margin);
                    const auto* row = m_grouping->wire_charge_row(apa, face, p, smin);
                    auto live = [&](int w, double& q, double& err) {
                        if (!row) return false;
                        auto it = row->find(w);
                        if (it == row->end()) return false;
                        q = it->second.first;
                        err = it->second.second;
                        return true;
                    };
                    auto dead = [&](int w) { return m_grouping->is_wire_dead(apa, face, p, w, smin); };
                    for (const auto& wv : RL::compose_activity(lo, hi, live, dead)) {
                        // Channel by IDENT through the anode: PDHD's wrapped induction
                        // wires are absent from IWirePlane::channels() (doc pdvd/31).
                        auto ich = ianode->channel(wires[wv.wire]->channel());
                        if (!ich) continue;
                        activity[ich] = ISlice::value_t(wv.charge, wv.error);
                    }
                }

                auto islice = std::make_shared<Aux::SimpleSlice>(nullptr, smin, smin * tick, (smax - smin) * tick, activity);
                WireCell::RayGrid::Blob shape = RL::shape_from_bounds(iface->raygrid(), bounds);
                auto iblob = std::make_shared<Aux::SimpleBlob>(blob_counter, 0.0f, 0.0f, shape, islice, iface);
                auto pcs = Aux::sample_live(sampler, iblob, m_wpid_angles.at(wpid), tick, blob_counter);
                ++blob_counter;
                if (pcs["3d"].size() == 0) continue;   // as the full path: a blob with no points is dropped
                new_cluster.node()->insert(Tree::Points(std::move(pcs)));
            }
        }
        else {
            // "no_paint" / "footprint": the full path's tiling, sampling and
            // remove_bad_blobs, without the two hack_activity_improved calls (and so
            // without the temp cluster, whose only use is the second call's path).
            const bool extend = (m_retile_mode == "no_paint");
            const auto wpid_set = orig_cluster.wpids_blob_set();
            for (auto it = wpid_set.begin(); it != wpid_set.end(); ++it) {
                int apa = it->apa();
                int face = it->face();
                const auto& angles = m_wpid_angles.at(*it);

                std::map<std::pair<int, int>, std::vector<WRG::measure_t> > map_slices_measures;
                get_activity_improved(orig_cluster, map_slices_measures, apa, face, extend);

                auto iblobs = make_iblobs_improved(map_slices_measures, apa, face);
                const size_t niblobs = iblobs.size();
                for (size_t bind = 0; bind < niblobs; ++bind) {
                    const IBlob::pointer iblob = iblobs[bind];
                    auto sampler = m_samplers.at(apa).at(face);
                    const double tick = m_grouping->get_tick().at(apa).at(face);
                    auto pcs = Aux::sample_live(sampler, iblob, angles, tick, bind);
                    if (pcs["3d"].size() == 0) continue;
                    new_cluster.node()->insert(Tree::Points(std::move(pcs)));
                }

                if (map_slices_measures.empty()) continue;
                int tick_span = map_slices_measures.begin()->first.second - map_slices_measures.begin()->first.first;
                auto blobs_to_remove = remove_bad_blobs(orig_cluster, new_cluster, tick_span, apa, face);
                for (const Blob* blob : blobs_to_remove) {
                    Blob& b = const_cast<Blob&>(*blob);
                    new_cluster.remove_child(b);
                }
            }
        }

        // Log-only census: one line per mutate.  Point counts and raw-coordinate
        // sums over the blob "3d" PCs (not Cluster::npoints(), whose cache child
        // insert/remove does not invalidate).  With a 'stepped' sampler configured
        // as the clustering job's, retile_mode "none" must reproduce the source
        // cloud: pts and sums equal.
        auto census = [](const Cluster& c, size_t& npts, double& sx, double& sy, double& sz) {
            npts = 0; sx = sy = sz = 0;
            for (const Blob* b : c.children()) {
                const auto& lpcs = b->node()->value.local_pcs();
                auto pit = lpcs.find("3d");
                if (pit == lpcs.end()) continue;
                const auto& ds = pit->second;
                auto ax = ds.get("x"); auto ay = ds.get("y"); auto az = ds.get("z");
                if (!ax || !ay || !az) continue;
                for (double v : ax->elements<double>()) sx += v;
                for (double v : ay->elements<double>()) sy += v;
                for (double v : az->elements<double>()) sz += v;
                npts += ds.size_major();
            }
        };
        size_t n_orig = 0, n_new = 0;
        double ox = 0, oy = 0, oz = 0, nx = 0, ny = 0, nz = 0;
        census(orig_cluster, n_orig, ox, oy, oz);
        census(new_cluster, n_new, nx, ny, nz);
        SPDLOG_LOGGER_DEBUG(log, "RETILEMODE mode={} ident={} blobs_orig={} blobs_new={} pts_orig={} pts_new={} "
                            "sum_orig=({:.17g},{:.17g},{:.17g}) sum_new=({:.17g},{:.17g},{:.17g})",
                            m_retile_mode, orig_cluster.ident(), orig_cluster.children().size(),
                            new_cluster.children().size(), n_orig, n_new, ox, oy, oz, nx, ny, nz);

        // Same T0 handling as the full path (see the comment there).
        auto& default_scope = orig_cluster.get_default_scope();
        auto& raw_scope = orig_cluster.get_raw_scope();
        if (default_scope.hash() != raw_scope.hash()) {
            auto correction_name = orig_cluster.get_scope_transform(default_scope);
            new_cluster.set_cluster_t0(orig_cluster.get_cluster_t0());
            new_cluster.add_corrected_points(m_pcts, correction_name);
            new_cluster.from(orig_cluster);
        }

        return m_grouping->remove_child(new_cluster);
    }

} // namespace WireCell::Clus
