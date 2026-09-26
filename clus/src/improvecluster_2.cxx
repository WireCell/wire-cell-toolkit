// ImproveCluster_2 - Second level cluster improvement 
//
// This class inherits from ImproveCluster_1 and provides additional
// cluster improvement functionality, building upon the Steiner tree
// enhancements from the first level.

#include "improvecluster_1.h"  // Include the ImproveCluster_1 header
#include "SteinerGrapher.h"
#include "WireCellClus/ResampleLive.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"
#include <chrono>

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

        // The non-"full" modes.  orig_cluster is the reinitialized source.
        std::unique_ptr<node_t> mutate_reduced(Cluster& orig_cluster) const;

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
        if (m_retile_mode != "full" && m_retile_mode != "no_paint" &&
            m_retile_mode != "footprint" && m_retile_mode != "none") {
            raise<ValueError>("ImproveCluster_2: unknown retile_mode '%s' "
                              "(expected full, no_paint, footprint or none)", m_retile_mode.c_str());
        }

        // Configure base class first
        ImproveCluster_1::configure(cfg);

        // doc pdvd/31 round 6 (owner Q5).  Same key name as CreateSteinerGraph's
        // so one jsonnet value can feed both and they cannot drift apart.
        m_steiner_terminal_charge =
            get(cfg, "terminal_charge_threshold", m_steiner_terminal_charge);
    }

    Configuration ImproveCluster_2::default_configuration() const
    {
        Configuration cfg = ImproveCluster_1::default_configuration();

        cfg["terminal_charge_threshold"] = m_steiner_terminal_charge;
        cfg["retile_mode"] = m_retile_mode;

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

        for (auto it = wpid_set.begin(); it != wpid_set.end(); ++it) {
            int apa = it->apa();
            int face = it->face();
            const auto& angles = m_wpid_angles.at(*it);

            std::map<std::pair<int, int>, std::vector<WRG::measure_t> > map_slices_measures;
            
            // get original activities ...
            t0 = Clock::now();
            get_activity_improved(*orig_cluster, map_slices_measures, apa, face);
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
            // start to sampling points 
            int npoints = 0;
            t0 = Clock::now();
            for (size_t bind=0; bind<niblobs; ++bind) {
          
                const IBlob::pointer iblob = iblobs[bind];
                auto sampler = m_samplers.at(apa).at(face);
                const double tick = m_grouping->get_tick().at(apa).at(face);

                auto pcs = Aux::sample_live(sampler, iblob, angles, tick, bind);
                // DO NOT EXTEND FURTHER! see #426, #430

                if (pcs["3d"].size()==0) continue; // no points ...
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
                new_cluster.node()->insert(Tree::Points(std::move(pcs)));
            }
            SPDLOG_LOGGER_TRACE(log, "timing: sample_live loop (apa={},face={}) took {} ms", apa, face, MS(Clock::now()-t0).count());
            SPDLOG_LOGGER_TRACE(log, "{} points sampled for apa {} face {} Blobs {}", npoints, apa, face, niblobs);

            // remove bad blobs
            t0 = Clock::now();
            if (map_slices_measures.empty()) continue; // no tiled blobs for this face
            int tick_span = map_slices_measures.begin()->first.second -  map_slices_measures.begin()->first.first;
            auto blobs_to_remove = remove_bad_blobs(*orig_cluster, new_cluster, tick_span, apa, face);
            // doc pdvd/40 round 3 census: point count before/after the removal
            // (bad_blob_report only) -- the proof that a removed blob's points
            // leave the retiled cluster the Steiner build sees.
            // Counted over the blob children, NOT via Cluster::npoints(): that
            // memo lives in the ClusterCache, which child insert/remove does
            // not invalidate (see remove_bad_blobs), so it would read stale.
            auto count_pts = [&]() { size_t n = 0; for (const Blob* b : new_cluster.children()) n += b->nbpoints(); return n; };
            const size_t npts_before = m_bad_blob_report ? count_pts() : 0;
            for (const Blob* blob : blobs_to_remove) {
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
