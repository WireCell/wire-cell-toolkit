#ifndef WIRECELLCLUS_TRACKFITTING_H
#define WIRECELLCLUS_TRACKFITTING_H

#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellUtil/Logging.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellClus/PRShower.h"
#include "WireCellClus/NeutrinoTaggerInfo.h"
#include "WireCellClus/PRVertexScoreboard.h"

#include <Eigen/IterativeLinearSolvers>
#include <unordered_map>
#include <unordered_set>


namespace WireCell::Clus {

    /**
     * Dedicated TrackFitting class that can be instantiated and used by 
     * other ensemble visitors without needing to be configured as a component.
     * 
     * This class encapsulates track fitting algorithms that can work on
     * individual clusters or collections of clusters.
     */
    class TrackFitting {
    public:
        
        enum class FittingType {
            Single,
            Multiple
        };

        /**
         * Structure to hold all track fitting parameters in one place
         */
        struct Parameters {
            // Diffusion coefficients (LArTPC standard values)
            double DL = 6.4* pow(units::cm,2)/units::second;                    // m²/s, longitudinal diffusion
            double DT = 9.8* pow(units::cm,2)/units::second;                    // m²/s, transverse diffusion

            // Software filter effects (wire dimension broadening)
            double col_sigma_w_T = 0.188060 * 3*units::mm * 0.2; // Collection plane, units: wire pitch
            double ind_sigma_u_T = 0.402993 * 3*units::mm * 0.3; // U induction plane
            double ind_sigma_v_T = 0.402993 * 3*units::mm * 0.5; // V induction plane

            // Uncertainty parameters
            double rel_uncer_ind = 0.075;          // Relative uncertainty for induction planes
            double rel_uncer_col = 0.05;           // Relative uncertainty for collection plane
            double add_uncer_ind = 0.0;            // Additional uncertainty for induction
            double add_uncer_col = 300.0;          // Additional uncertainty for collection
            
            // Longitudinal filter effects (time dimension)
            double add_sigma_L = 1.428249 *0.5505*units::mm /0.5;   
             
            // Additional useful parameters for charge err estimation ...
            double rel_charge_uncer = 0.1; // 10% 
            double add_charge_uncer = 600; // electrons

            double default_charge_th = 100;
            double default_charge_err = 1000; 

            double scaling_quality_th = 0.5;
            double scaling_ratio = 0.05;

            double area_ratio1 = 1.8*units::mm;
            double area_ratio2 = 1.7;

            double skip_default_ratio_1 = 0.25;
            double skip_ratio_cut = 0.97;
            double skip_ratio_1_cut = 0.75;

            double skip_angle_cut_1 = 160;
            double skip_angle_cut_2 = 90;
            double skip_angle_cut_3 = 45;
            double skip_dis_cut = 0.5*units::cm;

            // doc pr/28 S17: on an isochronous cluster (small blob-center drift-x
            // extent relative to its y/z footprint, the iso_band_like measure --
            // clustering_neutrino.cxx) the multi-track charge veto's revert
            // (skip_trajectory_point's `p = ps_point`) has no distance-based safety
            // margin -- see porting_dictionary.md.  C++ default -1 = off: the revert
            // is unconditional, matching the prototype (trajectory_fit.h:745-750)
            // and the toolkit's own pre-S17 behaviour (23bd6783).  When >= 0, WCT
            // internal length units (units::cm=10, so 20 cm = 200.0 -- the JSON
            // loader in TaggerCheckNeutrino/TaggerCheckSTM does no unit conversion,
            // same as every other skip_* knob here), abstain from the revert (keep
            // the fitted point) for any point whose segment's cluster has
            // blob-center drift-x extent below this cut.
            double skip_revert_iso_xext_cut = -1;

            // doc pr/49: cross-cluster projection-ghost filter on the 2D
            // charge association (18255-57441: a physically unrelated
            // cluster's charge aliases with the fitted track in the V view
            // only and detours the fit).  C++ default -1 = off: legacy path,
            // byte-identical.  When >= 0, a candidate 2D cell whose charge is
            // live (flag 1), which is NOT covered by the fitted cluster's OWN
            // blobs but IS covered by an OUT-OF-SCOPE cluster's (per-plane
            // exact wire interval x slice interval from
            // Cluster::time_blob_map()) is down-weighted in the fit.
            // "Out-of-scope" (round 3, owner decision 2026-08-08): a cluster
            // with NO segment in the current fit -- clusters fitted together
            // in the same PR graph are aware of each other and their shared
            // projections are legitimate, so they never count as foreign
            // (see m_cov_fit_scope).  Cells covered by nobody (own-track
            // charge just past the tiled envelope, dead-channel single-view
            // charge with no 3D image) are kept at full weight, so events
            // with no out-of-scope overlap are untouched.  The value is the
            // wire/slice tolerance in cells (0 = strict; the 57441
            // contamination sits ONE cell away, so tolerance >= 1 re-admits
            // it).  Dead-derived cells (flag 0), rescue-injected anchors and
            // end/vertex points (flag_end_point) are exempt.
            double fit_blob_coverage = -1;

            // doc pr/49: OPTIONAL 3D gate on the foreign test above -- when
            // > 0, an out-of-scope cluster's claim only triggers the
            // deweight if that cluster is FARTHER than this from the 3D
            // point being fit.  Round 3 default 0 = DISABLED (scope-only,
            // owner decision 2026-08-08: any out-of-scope claim deweights
            // regardless of distance; graph-scope membership replaced the 15
            // cm far-gate as the "not in the fitting range" criterion).
            // WCT internal length units (units::cm = 10; the JSON loader
            // does no unit conversion).  Inert while fit_blob_coverage < 0.
            double fit_blob_coverage_ghost_dis = 0;

            // doc pr/49: the least-squares weight multiplier applied in
            // fit_point to foreign-ghost cells (owner: "deweight the foreign
            // charge in the fitting, but not disable them completely" -- a
            // dead-channel region can leave good single-view charge with no
            // 3D image that the fit must still use).  1.0 would keep full
            // weight (filter becomes a no-op); 0 would be a hard drop.
            // Inert while fit_blob_coverage < 0.
            double fit_blob_coverage_weight = 0.1;

            // doc sbnd_xin/docs/pr/67 -- LOG-ONLY probe (0 = off = no lines =
            // byte-identical).  examine_end_ps_vec is the primary END trimmer:
            // it pops points off the front and back of a trajectory while
            // Grouping::is_good_point(p_raw, apa, face, 0.2 cm, 0, 0) is false.
            // That makes it the direct mechanism behind the owner's second
            // hypothesis for the pr/67 cases -- "it is also possible that we
            // have the situation covered, but then the track trajectory was
            // removed somehow".  Today it removes points silently, so a tip
            // that was fitted and then amputated looks identical to one that
            // was never reached.  Reported as a double for the existing
            // set_parameter(name, value) plumbing.
            double traj_cover_probe = 0;

            // doc sbnd_xin/docs/pr/107 -- dQ/dx-fit point retention (prototype
            // parity).  do_multi_tracking runs a THIRD form_map_graph pass
            // right before dQ_dx_multi_fit that the prototype does not have
            // (pr/28 T4: PR::Fit::reset() clears the fit indices, the
            // prototype's reset_fit_prop() is a resize that keeps them).
            // That pass re-applies form_map_graph's "store the point only if
            // its U+V+W quantity > 0" rule to the FINAL trajectory, so with
            // fit_exclusion on every interior point whose cells were all
            // stripped by update_association (cells equidistant from two
            // prongs -- the last ~1 cm before a multi-prong junction) is
            // DELETED from the trajectory before the dQ/dx fit, from the
            // segment's fit point cloud and from the DL vertex net's input.
            // The prototype fits dQ/dx on every trajectory point (its dQ/dx
            // fit, like the toolkit's, never reads the association maps).
            // > 0: the pre-dQ/dx pass keeps every point (zero-quantity
            // points get an index and an empty association, exactly what the
            // prototype's resize leaves).  Trajectory rounds 1-2 are
            // untouched.  0 = legacy drop = byte-identical.  Reported as a
            // double for the set_parameter(name, value) plumbing.
            double dqdx_fit_keep_all_points = 0;

            // doc pdvd/30 round 2 (PDVD 039252/2 evt 298595, cluster 86).
            // organize_segments_path_3rd takes its input path from
            // segment->fits() whenever that is merely *non-empty*, and only
            // falls back to the raw wcpts() when it is exactly empty.  A
            // segment whose fits() has collapsed to its two endpoint vertices
            // therefore re-enters the pass as a 2-point straight chord, gets
            // resampled into a straight interpolation, and never again
            // consults the (bent) wcpts() path still held in the same object.
            // Measured on evt 298595: a 128-wcpt, ~64 cm arm rendered as a
            // chord with max perpendicular deviation 0.001 cm sitting a median
            // 3.45 cm off the trajectory TaggerCheckSTM fits on the same
            // charge.  > 0: when fits() carries no shape (<= 2 points) while
            // wcpts() carries substantially more, prefer wcpts().  The
            // prototype needs no such test -- its ProtoSegment ctor seeds
            // fit_pt_vec (= get_point_vec(), the same field, see
            // clus/docs/porting/porting_dictionary.md:218-219) from the FULL
            // path_wcps, so its equivalent input is never degenerate by
            // construction.  0 = legacy = byte-identical.  Applies to
            // organize_segments_path_3rd only; organize_segments_path_2nd
            // (:1716) has the identical pattern but is deliberately left
            // untouched for want of measurement (doc pdvd/30 sec "Scope").
            double traj_degenerate_wcpts_fallback = 0;

            double default_dQ_dx = 5000;

            double end_point_factor=0.6;
            double mid_point_factor=0.9;
            int nlevel=3;
            double charge_cut=2000;
            
            double low_dis_limit = 1.2*units::cm;            // cm, lower distance limit for point organization
            double end_point_limit = 0.6*units::cm;          // cm, extension distance for end points
            double time_tick_cut = 20;   //            //  tick cut for point association

            // addition parameters
            double share_charge_err = 8000;
            double min_drift_time = 50*units::us;
            double search_range = 10; // wires, or time slices (not ticks)

            double dead_ind_weight = 0.3;
            double dead_col_weight = 0.9;
            double close_ind_weight = 0.15;
            double close_col_weight = 0.45;
            double overlap_th = 0.5;
            double dx_norm_length = 0.6*units::cm;
            double lambda= 0.0005;

            double div_sigma = 0.6*units::cm;            
        };
 
        /**
         * Constructor
         * @param fitting_type The type of fitting to perform (single or multiple tracks)
         */
        explicit TrackFitting(FittingType fitting_type = FittingType::Single);   
        virtual ~TrackFitting();

        /**
         * Set the fitting type
         * @param fitting_type The new fitting type to use
         */
        void set_fitting_type(FittingType fitting_type) { m_fitting_type = fitting_type; }

        /**
         * Enable/disable per-step timing output (printed to stdout)
         */
        void set_perf(bool perf) { m_perf = perf; }
        bool get_perf() const { return m_perf; }

        /**
         * Get the current fitting type
         * @return The current fitting type
         */
        FittingType get_fitting_type() const { return m_fitting_type; }

        // Parameter management methods
        
        /**
         * Get read-only access to current parameters
         */
        const Parameters& get_parameters() const { return m_params; }
        
        /**
         * Set new parameters (replaces all current parameters)
         */
        void set_parameters(const Parameters& params) { m_params = params; }
        
        /**
         * Set specific parameter by name
         */
        void set_parameter(const std::string& name, double value);
        
        /**
         * Get specific parameter by name
         */
        double get_parameter(const std::string& name) const;

        // single track fitting utilizes the segments ... 
        void add_segment(std::shared_ptr<PR::Segment> segment);
        /**
         * Get the set of segments currently stored in this TrackFitting instance.
         * @return Set of shared pointers to PR::Segment
         */
        std::set<std::shared_ptr<PR::Segment>> get_segments() const { return m_segments; }
        void clear_segments();
 
        // multi-track fitting utilized the Graph ... 
        void add_graph(std::shared_ptr<PR::Graph> graph);
        std::shared_ptr<PR::Graph> get_graph() const { return m_graph; }

        /// Store / retrieve the identified neutrino interaction vertex.
        void set_main_vertex(PR::VertexPtr v) { m_main_vertex = v; }
        PR::VertexPtr get_main_vertex() const { return m_main_vertex; }

        /// Store / retrieve the full set of reconstructed showers.
        void set_showers(PR::IndexedShowerSet showers) { m_showers = std::move(showers); }
        const PR::IndexedShowerSet& get_showers() const { return m_showers; }

        /// doc sbnd_xin/docs/pr/40 round 9 B2 -- clusters connected to the
        /// main cluster by an nv_bridge_track bridge segment.  Written by
        /// PatternAlgorithms::nv_bridge_track (cleared at each
        /// shower_clustering_with_nv entry); read by fill_bee_pf_tree's
        /// pf_track_bridged_clusters gate.  Empty when the bridge knob is
        /// off => downstream predicates are byte-identical to legacy.
        void clear_bridged_cluster_ids() { m_bridged_cluster_ids.clear(); }
        void add_bridged_cluster_id(int id) { m_bridged_cluster_ids.insert(id); }
        const std::set<int>& get_bridged_cluster_ids() const { return m_bridged_cluster_ids; }

        /// doc sbnd_xin/docs/pr/92 -- stray satellite showers dropped from the
        /// kinematics tree by fill_kine_tree's kine_drop_stray_satellites gate
        /// (shower ids, PR::Shower::get_shower_id()).  Stashed unconditionally
        /// per event by TaggerCheckNeutrino (replace semantics: empty when the
        /// knob is off or no vertex was found, so no cross-event carryover);
        /// read by fill_bee_pf_tree's pf_drop_stray_satellites gate.
        void set_dropped_satellite_shower_ids(std::set<int> ids) { m_dropped_satellite_shower_ids = std::move(ids); }
        const std::set<int>& get_dropped_satellite_shower_ids() const { return m_dropped_satellite_shower_ids; }

        /// Store / retrieve pi0 identification results from TaggerCheckNeutrino.
        void set_pi0_data(PR::IndexedShowerSet pi0_showers,
                          PR::ShowerIntMap map_shower_pio_id,
                          std::map<int, std::vector<PR::ShowerPtr>> map_pio_id_showers,
                          std::map<int, std::pair<double, int>> map_pio_id_mass)
        {
            m_pi0_showers        = std::move(pi0_showers);
            m_map_shower_pio_id  = std::move(map_shower_pio_id);
            m_map_pio_id_showers = std::move(map_pio_id_showers);
            m_map_pio_id_mass    = std::move(map_pio_id_mass);
        }
        const PR::IndexedShowerSet& get_pi0_showers() const { return m_pi0_showers; }
        const PR::ShowerIntMap& get_map_shower_pio_id() const { return m_map_shower_pio_id; }
        const std::map<int, std::vector<PR::ShowerPtr>>& get_map_pio_id_showers() const { return m_map_pio_id_showers; }
        const std::map<int, std::pair<double, int>>& get_map_pio_id_mass() const { return m_map_pio_id_mass; }

        /// Store / retrieve reconstructed neutrino kinematics (filled by TaggerCheckNeutrino).
        void set_kine_info(PR::KineInfo ki)  { m_kine_info = std::move(ki); }
        const PR::KineInfo& get_kine_info()  const { return m_kine_info; }

        /// Store / retrieve BDT input features (filled by TaggerCheckNeutrino).
        void set_tagger_info(PR::TaggerInfo ti) { m_tagger_info = std::move(ti); }
        const PR::TaggerInfo& get_tagger_info() const { return m_tagger_info; }
        PR::TaggerInfo& get_tagger_info_mutable() { return m_tagger_info; }

        /// Store / retrieve the per-event vertex scoreboard (doc
        /// sbnd_xin/docs/pr/75).  Stashed by TaggerCheckNeutrino beside
        /// set_tagger_info, i.e. AFTER snap_main_vertex_to_kink and the final
        /// improve_vertex, and read by PrDisplayDump.  Empty (filled==false)
        /// unless the `vertex_scoreboard` knob was on -- read that as "no
        /// scoreboard taken", never as "no candidates".
        void set_vertex_scoreboard(PR::VertexScoreboard vsb) { m_vertex_scoreboard = std::move(vsb); }
        const PR::VertexScoreboard& get_vertex_scoreboard() const { return m_vertex_scoreboard; }

        void clear_graph();

        /// Drop every piece of state that belongs to one event.
        ///
        /// A visitor that owns a TrackFitting as a member (TaggerCheckSTM,
        /// TaggerCheckNeutrino) keeps that fitter for the whole process, but
        /// almost everything the fitter caches -- m_grouping, the cluster and
        /// blob sets, global_rb_map, the charge maps -- points into the
        /// per-event Points tree that MultiAlgBlobClustering::operator() builds
        /// as a local and destroys when the event ends.  Streaming a second
        /// event through the same process without this call is a use-after-free
        /// (SBND mcp1k 48367 then 48895: SIGSEGV in form_point_association).
        ///
        /// Call it once per event, at the top of visit() -- never per
        /// neutrino candidate.  Within one event the candidate loop
        /// deliberately reuses the member fitter for candidate 0 and builds a
        /// fresh one for the rest; this call adds freshness ACROSS events
        /// without removing it WITHIN one.
        ///
        /// Inert in a one-event process: at the first visit() there is nothing
        /// to drop, so the legacy per-event job is byte-identical.
        ///
        /// doc 76 round 3: it must drop the state that gives a WRONG ANSWER as
        /// well as the state that CRASHES.  The ident-keyed memo
        /// m_cluster_xext_cache is the one that moved a physics number -- see
        /// the comment in the body.  Set WCT_TF_RESET_CENSUS=1 to have it log
        /// what the previous event actually left behind.
        void reset_for_new_event();

        void add_cluster(std::shared_ptr<Facade::Cluster> cluster);

        /// Pre-load all clusters at once and call prepare_data() a single time.
        /// Call this before starting pattern recognition so that subsequent
        /// do_multi_tracking calls can use flag_force_load_data=false.
        void preload_clusters(const std::vector<Facade::Cluster*>& clusters);

        // collect charge
        void prepare_data();

        // Fill the global readout map
        void fill_global_rb_map();

        /**
         * Organize original path from segment points with distance limits
         * @param segment Pointer to PR::Segment containing the path points
         * @param low_dis_limit Lower distance limit for point organization
         * @param end_point_limit Extension distance for end points
         * @return Vector of organized 3D points
         */
        std::vector<WireCell::Point> organize_orig_path(std::shared_ptr<PR::Segment> segment, double low_dis_limit=1.2*units::cm, double end_point_limit=0.6*units::cm);

        std::vector<WireCell::Point> examine_end_ps_vec(std::shared_ptr<PR::Segment> segment, const std::vector<WireCell::Point>& pts, bool flag_start, bool flag_end);

        void organize_ps_path(std::shared_ptr<PR::Segment> segment, std::vector<WireCell::Point>& pts, double low_dis_limit, double end_point_limit);

        // use the m_graph to organize ...
        void organize_segments_path(double low_dis_limit, double end_point_limit);
        // use m_graph, after first round of fitting
        void organize_segments_path_2nd(double low_dis_limit, double end_point_limit);
        // use m_graph, after second round of fitting 
        void organize_segments_path_3rd(double step_size);

    private:
        // Helper functions for organize_segments_path methods
        
        /**
         * Check and reset vertices that are too close together
         * @param edge_range The range of edges (segments) to check
         */
        void check_and_reset_close_vertices();
        
        /**
         * Get segment vertices in correct order (start, end)
         * @param segment The segment to process
         * @param ed The edge descriptor
         * @param start_v Output: pointer to start vertex
         * @param end_v Output: pointer to end vertex
         * @param vd1 Output: descriptor for vertex 1
         * @param vd2 Output: descriptor for vertex 2
         * @return true if successful, false otherwise
         */
        bool get_ordered_segment_vertices(
            std::shared_ptr<PR::Segment> segment,
            const PR::edge_descriptor& ed,
            std::shared_ptr<PR::Vertex>& start_v,
            std::shared_ptr<PR::Vertex>& end_v,
            PR::node_descriptor& vd1,
            PR::node_descriptor& vd2
        );
        
        /**
         * Generate 2D projections and create fit vector from 3D points
         * @param segment The segment for which to generate fits
         * @param pts The 3D points
         * @return Vector of Fit objects with 3D points and 2D projections
         */
        std::vector<PR::Fit> generate_fits_with_projections(
            std::shared_ptr<PR::Segment> segment,
            const std::vector<WireCell::Point>& pts
        );

    public:
        
        

                /// Internal coordinate (can be more complex)
        struct Coord2D {
            int apa, face, time, wire, channel;
            WirePlaneLayer_t plane;  // Additional internal information

            Coord2D(int a, int f, int t, int w, int c, WirePlaneLayer_t p)
                : apa(a), face(f), time(t), wire(w), channel(c), plane(p) {}

            bool operator<(const Coord2D& other) const {
                if (apa != other.apa) return apa < other.apa;
                if (face != other.face) return face < other.face;
                if (time != other.time) return time < other.time;
                if (wire != other.wire) return wire < other.wire;
                if (channel != other.channel) return channel < other.channel;
                return plane < other.plane;
            }
        };

        /// Per-plane data for 3D points (exactly matches prototype)
        struct PlaneData {
            std::set<Coord2D> associated_2d_points;
            // doc pr/49 (fit_blob_coverage knob): subset of
            // associated_2d_points classified as foreign-ghost cells (live,
            // outside the fitted cluster's own blob coverage, inside an
            // out-of-scope cluster's -- round 3: one with no segment in the
            // current fit).  They stay in the fit -- a
            // dead-channel region can leave good single-view charge with no
            // 3D image, which the fit must still use -- but fit_point scales
            // their least-squares weight by fit_blob_coverage_weight.
            // Empty on the legacy path (knob off).
            std::set<Coord2D> deweighted_2d_points;
            double quantity;

            PlaneData() : quantity(0.0) {}
        };

        /// 3D point with per-plane associations (corrected structure)
        struct Point3DInfo {
            std::map<WirePlaneLayer_t, PlaneData> plane_data;
            
            const PlaneData& get_plane_data(WirePlaneLayer_t plane) const {
                static PlaneData empty;
                auto it = plane_data.find(plane);
                return (it != plane_data.end()) ? it->second : empty;
            }
            
            void set_plane_data(WirePlaneLayer_t plane, const PlaneData& data) {
                plane_data[plane] = data;
            }
        };

        struct CoordReadout {
            int apa, time, channel;

            CoordReadout(int a, int t, int c)
            : apa(a), time(t), channel(c) {}

            bool operator<(const CoordReadout& other) const {
            if (apa != other.apa) return apa < other.apa;
            if (time != other.time) return time < other.time;
            return channel < other.channel;
            }

            bool operator==(const CoordReadout& other) const {
                return apa == other.apa && time == other.time && channel == other.channel;
            }
        };

        struct CoordReadoutHash {
            size_t operator()(const CoordReadout& k) const {
                size_t h = std::hash<int>{}(k.apa);
                h ^= std::hash<int>{}(k.time)    + 0x9e3779b9 + (h << 6) + (h >> 2);
                h ^= std::hash<int>{}(k.channel) + 0x9e3779b9 + (h << 6) + (h >> 2);
                return h;
            }
        };


        /// Simple charge measurement (in ternal interface)
        struct ChargeMeasurement {
            double charge, charge_err;
            int flag;

            ChargeMeasurement(double q = 0.0, double qe = 0.0, int f = 0)
                : charge(q), charge_err(qe), flag(f) {}
        };

        /// Fitted 2D charge result: measured + predicted + cluster association
        struct FittedCharge2D {
            double charge;        // original measurement charge
            double charge_err;    // original measurement uncertainty
            double pred_charge;   // predicted charge (un-whitened, same units as charge)
            int flag;             // 0=dead, 1=live, 2=bad
            /// Clusters owning this cell.  Ident-ordered, not pointer-ordered:
            /// the only current consumer (PrDisplayDump::dump_proj) sorts the
            /// ids itself, so this changes no output today -- it removes the
            /// trap for the next consumer that iterates it directly.
            std::set<Facade::Cluster*, PR::ClusterPtrCmp> clusters;
        };

        using WireTime = std::pair<int, int>;            // (wire_index, time_slice)
        using APAFacePlane = std::tuple<int, int, int>;   // (apa, face, plane);

        /// doc pr/109 §8: one captured snapshot of m_fitted_charge_2d, for one
        /// cluster's dQ/dx fit.  Kept per fit (not merged) so a writer can emit
        /// the prediction that THIS cluster's fit produced -- see
        /// m_cluster_fitted_charge_2d for why the merged map cannot.
        struct ClusterFitted2D {
            Facade::Cluster* cluster{nullptr};
            int ident{-1};
            std::map<APAFacePlane, std::map<WireTime, FittedCharge2D>> cells;
        };

        // Fill fitted 2D charge results after dQ/dx fitting
        void fill_fitted_charge_2d(
            const std::map<CoordReadout, std::pair<ChargeMeasurement, std::set<Coord2D>>>& map_U,
            const std::map<CoordReadout, std::pair<ChargeMeasurement, std::set<Coord2D>>>& map_V,
            const std::map<CoordReadout, std::pair<ChargeMeasurement, std::set<Coord2D>>>& map_W,
            const Eigen::VectorXd& pred_u, const Eigen::VectorXd& pred_v, const Eigen::VectorXd& pred_w,
            double rel_uncer_ind, double rel_uncer_col,
            double add_uncer_ind, double add_uncer_col);

        /// Merge every per-cluster snapshot captured inside fill_fitted_charge_2d
        /// into m_fitted_charge_2d so the flat map covers every cluster fit this
        /// event.  Call once per event after all do_multi_tracking() calls are done
        /// (e.g. at the end of TaggerCheckNeutrino::visit, before set_track_fitting).
        void assemble_fitted_charge_2d();

        /// Merge an externally-accumulated fitted-charge map into
        /// m_fitted_charge_2d (last-writer-wins on overlap, same semantics as
        /// assemble_fitted_charge_2d).  Used by TaggerCheckSTM's save_stm_fit
        /// dump to hand per-pass snapshots to a fresh holder fitter that
        /// downstream writers (SbndMagnifyTrackingVisitor) read from the
        /// grouping's named "stm" slot.
        void merge_fitted_charge_2d(const std::map<APAFacePlane, std::map<WireTime, FittedCharge2D>>& other);

        // point associations
        void form_point_association(std::shared_ptr<PR::Segment> segment, WireCell::Point &p, PlaneData& temp_2dut, PlaneData& temp_2dvt, PlaneData& temp_2dwt, double dis_cut, int nlevel, double time_tick_cut );

        void examine_point_association(std::shared_ptr<PR::Segment> segment, WireCell::Point &p, PlaneData& temp_2dut, PlaneData& temp_2dvt, PlaneData& temp_2dwt, bool flag_end_point = false, double charge_cut = 2000);

        /// doc pr/49 (fit_blob_coverage knob): true iff any of the cluster's
        /// OWN blobs covers (wire, time) in the given plane -- wire in the
        /// blob's half-open per-plane interval [min-tol, max+tol) and time
        /// (tick) in [slice_index_min-tol_ticks, slice_index_max+tol_ticks).
        /// Interval search over Cluster::time_blob_map() keys, so it accepts
        /// both blob-aligned times and the floor-quantized ticks of the
        /// Steiner/fallback candidate branches without an alignment
        /// assumption.  Order-invariant existential test (OR over blobs), so
        /// iterating the pointer-keyed BlobSet cannot affect the result.
        bool is_cell_covered_by_own_blobs(const Facade::Cluster* cluster, int apa, int face,
                                          int plane, int wire, int time,
                                          int tol_cells, int nticks_per_slice) const;

        /// doc pr/49 (round 3, scope-aware): true iff some OUT-OF-SCOPE
        /// cluster covers the cell per the same test.  Out-of-scope = not
        /// `cluster` itself AND not in m_cov_fit_scope (clusters owning a
        /// segment in the current fit -- "fitted together" clusters are
        /// aware of each other and their shared projections are legitimate,
        /// so they never count as foreign).  ghost_dis > 0 additionally
        /// requires `other->get_closest_dis(p) > ghost_dis` (optional extra
        /// gate; round-3 default 0 = scope-only).  A cell covered by nobody
        /// (own-track charge spilling just past the tiled envelope, or
        /// dead-channel single-view charge with no 3D image) is kept at
        /// full weight.  Order-invariant OR over clusters and blobs.
        /// `claimant` (optional, doc pr/50 diagnostics) receives the first
        /// covering out-of-scope cluster; the result is unchanged.
        bool is_cell_covered_by_foreign_blobs(const Facade::Grouping* grouping,
                                              const Facade::Cluster* cluster,
                                              const WireCell::Point& p, double ghost_dis,
                                              int apa, int face,
                                              int plane, int wire, int time,
                                              int tol_cells, int nticks_per_slice,
                                              const Facade::Cluster** claimant = nullptr) const;

        /// doc pr/49 round 3: rebuild m_cov_fit_scope = clusters owning a
        /// segment in the current fit.  Walks m_graph's edges when a graph
        /// is present (fresh walk, NOT the stale m_cluster_edges -- the
        /// single-tracking path never rebuilds that map) and adds `seg`'s
        /// cluster (the one segment being fitted on the form_map path;
        /// nullptr on the form_map_graph path, where the graph supplies
        /// everything).  Called once per form_map/form_map_graph invocation,
        /// only while the fit_blob_coverage knob is on.  doc pr/50: the same
        /// walk also refreshes m_cov_vtx_info (graph vertex positions +
        /// degrees) for the deweight sentinel diagnostics.
        void rebuild_cov_fit_scope(const std::shared_ptr<PR::Segment>& seg);
        /// doc pr/98 perf: the keep/strip decision depends only on (cell,
        /// segment) -- never on the fit point -- and every segment cloud is
        /// static while ONE segment's fit points are processed (its own
        /// set_fit_associate_vec runs after its loop; other segments' ran
        /// before or run after).  Consecutive fit points re-claim the same
        /// cells, so form_map_graph passes a per-segment decision cache.
        /// nullptr disables caching (identical decisions either way).
        /// Keyed by the packed (apa,face,plane,wire,time) identity of a
        /// Coord2D (see exclusion_cache_key in TrackFitting.cxx) -- an
        /// unordered_map on one integer beats std::map's Coord2D
        /// operator< chain in the pr/98 profile.
        using ExclusionDecisionCache = std::unordered_map<uint64_t, bool>;
        void update_association(std::shared_ptr<PR::Segment> segment,
                                const std::vector<std::shared_ptr<PR::Segment>>& all_segments,
                                PlaneData& temp_2dut, PlaneData& temp_2dvt, PlaneData& temp_2dwt,
                                ExclusionDecisionCache* decision_cache = nullptr);

        void form_map(std::vector<std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>>& ptss, double end_point_factor=0.6, double mid_point_factor=0.9, int nlevel=3, double time_tick_cut=20, double charge_cut=2000);
        void form_map_graph(bool flag_exclusion, double end_point_factor=0.6, double mid_point_factor=0.9, int nlevel=3, double time_tick_cut=20, double charge_cut=2000);

        // track trajectory fitting // should fit all APA ...
        void trajectory_fit(std::vector<std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>>& pss_vec, int charge_div_method = 1, double div_sigma = 0.6*units::cm);
        WireCell::Point fit_point(WireCell::Point& init_p, int i, std::shared_ptr<PR::Segment> segment,std::map<std::pair<int, int>, std::map<std::tuple<int, int, int>, double>>& map_Udiv_fac, std::map<std::pair<int, int>, std::map<std::tuple<int, int, int>, double>>& map_Vdiv_fac, std::map<std::pair<int, int>, std::map<std::tuple<int, int, int>, double>>& map_Wdiv_fac, double offset_t, double slope_x, double offset_u, double slope_yu, double slope_zu, double offset_v, double slope_yv, double slope_zv, double offset_w, double slope_yw, double slope_zw);
        void multi_trajectory_fit(int charge_div_method = 1, double div_sigma = 0.6*units::cm);

        // examine trajectory ...
        // doc pr/28 T3: init_indices carries the GLOBAL m_3d_to_2d key of each
        // point (prototype's init_indices, multi_track_fitting.h:429).  The
        // per-segment loop position is NOT that key in the multi-track path.
        std::vector<WireCell::Point> examine_segment_trajectory(std::shared_ptr<PR::Segment> segment, std::vector<WireCell::Point>& final_ps_vec, std::vector<WireCell::Point>& init_ps_vec, const std::vector<int>& init_indices);
        bool skip_trajectory_point(WireCell::Point& p, std::pair<int,int>& apa_face, int i, int index, std::vector<std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>>& pss_vec,  std::vector<std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>>& fine_tracking_path);

        // prepare for dQ/dx fitting
        double cal_gaus_integral(int tbin, int wbin, double t_center, double t_sigma, 
                                       double w_center, double w_sigma, int flag, double nsigma, int cur_ntime_ticks);

        double cal_gaus_integral_seg(int tbin, int wbin, std::vector<double>& t_centers, std::vector<double>& t_sigmas, std::vector<double>& w_centers, std::vector<double>& w_sigmas, std::vector<double>& weights, int flag, double nsigma, int cur_ntime_ticks);

        void update_dQ_dx_data();
        void recover_original_charge_data();

        /**
         * Calculate compact matrix analysis for wire plane sharing
         * 
         * This function analyzes the sharing patterns between 2D measurements and 3D positions
         * to compute overlap ratios and adjust weight matrix coefficients. It processes sparse
         * matrices representing the relationship between 2D wire measurements and 3D positions.
         * 
         * @param weight_matrix Reference to sparse weight matrix (MW, MV, or MU) to be modified
         * @param response_matrix_transpose Transposed response matrix (RWT, RVT, or RUT)  
         * @param n_2d_measurements Number of 2D measurements (wire/time points)
         * @param n_3d_positions Number of 3D positions
         * @param cut_position Threshold for wire sharing cut (default 2.0)
         * @return Vector of pairs containing overlap ratios for each 3D position
         *         Each pair contains (previous_neighbor_ratio, next_neighbor_ratio)
         */
        std::vector<std::pair<double, double>> calculate_compact_matrix(Eigen::SparseMatrix<double>& weight_matrix, const Eigen::SparseMatrix<double>& response_matrix_transpose, int n_2d_measurements, int n_3d_positions, double cut_position = 2.0);
        std::vector<std::vector<double>> calculate_compact_matrix_multi(std::vector<std::vector<int> >& connected_vec,Eigen::SparseMatrix<double>& weight_matrix, const Eigen::SparseMatrix<double>& response_matrix_transpose, int n_2d_measurements, int n_3d_positions, double cut_position = 2.0);

        void dQ_dx_fill(double dis_end_point_ext=0.45*units::cm);

        void dQ_dx_fit(double dis_end_point_ext=0.45*units::cm, bool flag_dQ_dx_fit_reg=true);
        void dQ_dx_multi_fit(double dis_end_point_ext=0.45*units::cm, bool flag_dQ_dx_fit_reg=true);

        void do_single_tracking(std::shared_ptr<PR::Segment> segment, bool flag_dQ_dx_fit_reg= true, bool flag_dQ_dx_fit= true, bool flag_force_load_data = false, bool flag_hack = false, Facade::Cluster* cluster_filter = nullptr);
        void do_multi_tracking(bool flag_dQ_dx_fit_reg= true, bool flag_dQ_dx_fit= true, bool flag_force_load_data = false, bool flag_exclusion =false, bool flag_hack = false, Facade::Cluster* cluster_filter = nullptr);


        

        /**  
         * Get anode for a specific APA identifier
         * @param apa_ident APA identifier (typically same as APA number)
         * @return Pointer to IAnodePlane, or nullptr if not found
         */
        IAnodePlane::pointer get_anode(int apa_ident = 0) const;

        /**
         * Get all available anodes from the grouping
         * @return Map of APA identifier to anode pointer
         */
        std::map<int, IAnodePlane::pointer> get_all_anodes() const;

        /**
         * Get the grouping associated with this TrackFitting instance
         * @return Pointer to Grouping, or nullptr if not set
         */
        Facade::Grouping* grouping() const { return m_grouping; }

        /**
         * Get channel number for a specific wire location
         * Uses hybrid caching for optimal performance
         * @param apa APA number
         * @param face Face number (0 or 1)
         * @param plane Plane index (0=U, 1=V, 2=W typically)  
         * @param wire Wire index within the plane
         * @return Channel number, or -1 if invalid
         */
        int get_channel_for_wire(int apa, int face, int plane, int wire) const;

        /**
         * Get all wires that belong to a specific channel
         * @param apa APA number
         * @param channel_number Channel identifier
         * @return Vector of wire information (face, plane, wire_index)
         */
        std::vector<std::tuple<int, int, int>> get_wires_for_channel(int apa, int channel_number) const;

        // map_apa_ch_plane_wires: (apa,channel) -> vector of (face, plane, wire)
        void collect_2D_charge(std::map<CoordReadout, ChargeMeasurement>& charge_2d_u, std::map<CoordReadout, ChargeMeasurement>& charge_2d_v, std::map<CoordReadout, ChargeMeasurement>& charge_2d_w, std::map<std::pair<int, int>, std::vector<std::tuple<int, int, int>>>& map_apa_ch_plane_wires);
        /**
         * Clear all caches (useful for memory management)
         */
        void clear_cache() const;

        /**
         * Get cache statistics for monitoring/debugging
         */
        struct CacheStats {
            size_t hot_planes_count;
            size_t cold_entries_count;
            size_t total_lookups;
            size_t hot_hits;
            size_t cold_hits;
            double hit_rate() const { 
                return total_lookups > 0 ? (double)(hot_hits + cold_hits) / total_lookups : 0.0; 
            }
        };
        CacheStats get_cache_stats() const;

        /**
         * Inherit pre-built geometry and cluster charge data from a parent fitter.
         *
         * Copies the wire-plane geometry, wire-channel cache, and the already-computed
         * charge data for @p cluster from @p src into this fitter.  After this call:
         *   - BuildGeometry() will NOT be called again (m_grouping is set)
         *   - prepare_data() will skip @p cluster (it is pre-populated in m_loaded_clusters)
         *
         * Intended for lightweight child fitters (e.g. in compare_main_vertices_all_showers)
         * that need to fit a single temporary segment without re-loading all cluster blobs.
         */
        void inherit_from(const TrackFitting& src, Facade::Cluster* cluster);

        /**
         * Set the detector volume for this TrackFitting instance
         * @param dv Pointer to IDetectorVolumes
         */
        void set_detector_volume(IDetectorVolumes::pointer dv) { m_dv = dv; }
        
        /**
         * Set the PCTransformSet for coordinate transformations
         * @param pcts Pointer to PCTransformSet interface
         */
        void set_pc_transforms(IPCTransformSet::pointer pcts) { m_pcts = pcts; }

        /**
         * Get the current detector volumes
         * @return Pointer to detector volumes interface
         */
        IDetectorVolumes::pointer get_detector_volume() const { return m_dv; }

        /**
         * Get the current PCTransformSet
         * @return Pointer to PCTransformSet interface
         */
        IPCTransformSet::pointer get_pc_transforms() const { return m_pcts; }

        std::vector<std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>> get_fine_tracking_path() const { return fine_tracking_path; }
        std::vector<double> get_dQ() const { return dQ; }
        std::vector<double> get_dx() const { return dx; }
        std::vector<double> get_pu() const { return pu; }
        std::vector<double> get_pv() const { return pv; }
        std::vector<double> get_pw() const { return pw; }
        std::vector<double> get_pt() const { return pt; }
        std::vector<std::pair<int,int>> get_paf() const {return paf;}

        /// Which dx rule applies to trajectory point i of n in dQ_dx_fit, given
        /// the per-point (apa, face) list: 0 = first point of an (apa,face) run
        /// (extrapolate backward from i+1), 1 = last point of a run (extrapolate
        /// forward from i-1), 2 = interior point, 3 = isolated.  A run boundary
        /// at the very LAST point is a "last" point, not a "first" one: the
        /// former predicate read i+1 there and threw std::out_of_range
        /// (doc pdvd/25 M3: PDVD trajectories cross 16 (anode,face) volumes).
        /// Identical to the former branch order whenever i+1 exists.
        static int dqdx_path_point_role(int i, int n, const std::vector<std::pair<int, int>>& paf);
        std::vector<double> get_reduced_chi2() const { return reduced_chi2; }

        // Measured 2D charge data access
        const std::unordered_map<CoordReadout, ChargeMeasurement, CoordReadoutHash>& get_charge_data() const { return m_charge_data; }

        // Fitted 2D charge data organized by (apa, face, plane) -> (wire, time)
        const std::map<APAFacePlane, std::map<WireTime, FittedCharge2D>>& get_fitted_charge_2d() const { return m_fitted_charge_2d; }

        /// Restrict the next fit (and the snapshot it captures) to one
        /// cluster.  do_multi_tracking()/do_single_tracking() set this from
        /// their cluster_filter argument; exposed so a test can drive
        /// fill_fitted_charge_2d() the way a real fit does.
        void set_cluster_filter(Facade::Cluster* cluster) { m_cluster_filter = cluster; }

        /// doc pr/109 §8: the per-cluster snapshots, in capture order.  Unlike
        /// get_fitted_charge_2d() these are NOT merged, so pred_charge is the
        /// value the named cluster's own fit produced.
        const std::vector<ClusterFitted2D>& get_cluster_fitted_charge_2d() const { return m_cluster_fitted_charge_2d; }

        /**
         * Get geometry information for wire plane offsets
         * @return Map of WirePlaneId to tuple (offset_t, offset_u, offset_v, offset_w)
         */
        const std::map<WirePlaneId, std::tuple<double, double, double, double>>& get_wpid_offsets() const { return wpid_offsets; }

        /**
         * Get geometry information for wire plane slopes
         * @return Map of WirePlaneId to tuple (slope_t, slope_yu_zu, slope_yv_zv, slope_yw_zw)
         */
        const std::map<WirePlaneId, std::tuple<double, std::pair<double, double>, std::pair<double, double>, std::pair<double, double>>>& get_wpid_slopes() const { return wpid_slopes; }

    private:
        // doc pr/98: test seam.  Lets the update_association doctest inject
        // synthetic wpid_offsets/wpid_slopes without a detector-volumes
        // fixture (BuildGeometry needs live IDetectorVolumes).  Declared
        // here only; defined in clus/test/doctest_update_association.cxx.
        friend struct TrackFittingTestHarness;

         // Core parameters - centralized storage
        Parameters m_params;

        // doc pr/28 S17: memoized cluster ident() -> blob-center drift-x extent
        // (cm), for skip_revert_iso_xext_cut.  Lookup/insert only, never
        // iterated -- no ordering hazard (CLAUDE.md only forbids ITERATING
        // pointer-keyed containers; this is int-keyed regardless).  Lives for
        // this TrackFitting instance (one multi-track fitting pass).
        std::unordered_map<int, double> m_cluster_xext_cache;

        // Helper method to get parameter value or default
        double get_param_or_default(double param_value, double default_value) const {
            return (param_value < 0) ? default_value : param_value;
        }

        bool m_perf{false};  // if true, print per-step timing to stdout
        FittingType m_fitting_type;
        IDetectorVolumes::pointer m_dv{nullptr};  
        IPCTransformSet::pointer m_pcts{nullptr};          // PC Transform Set
        
        // cluster and grouping, CTPC is from m_grouping ...
        // Ordered by cluster id (not pointer): prepare_data() iterates these and
        // FP-accumulates shared dead-channel charge, so iteration order must be
        // content-stable across runs.
        Facade::Grouping* m_grouping{nullptr};
        std::set<Facade::Cluster*, PR::ClusterPtrCmp> m_clusters;
        std::set<Facade::Cluster*, PR::ClusterPtrCmp> m_loaded_clusters;  ///< Clusters whose charge data has been loaded into m_charge_data
        bool m_charge_data_dirty{true};                ///< True when m_clusters has clusters not yet in m_charge_data
        Facade::Cluster* m_cluster_filter{nullptr};    ///< If non-null, restrict fitting to segments of this cluster
        // doc pr/107: set by do_multi_tracking around its pre-dQ/dx
        // form_map_graph call only (never on the trajectory rounds); while
        // true form_map_graph stores zero-quantity interior points instead
        // of dropping them.  Always false while the knob is 0.
        bool m_keep_zero_quantity_points{false};
        // doc pr/108 stage dump (debug only, env WCT_TRAJ_DUMP=<path>; unset => no
        // code path): per trajectory round, every interior point's association
        // counts before/after exclusion, live-plane quantities, kept flag, and
        // the fitted positions after each round.  Mirrors WCP_TRAJ_DUMP in the
        // prototype so the two can be diffed stage by stage.
        FILE* m_traj_dump{nullptr};
        int   m_traj_dump_call{0};
        int   m_traj_dump_stage{0};
        void  traj_dump_fits(const char* tag);

        // doc pr/109 sec 9 exclusion-decision dump (debug only, env
        // WCT_EXCL_DUMP=<path>; unset => no code path).  update_association's
        // keep rule is "strictly closest of all siblings, or within 0.3 cm",
        // but it discards the two distances the decision turns on and
        // early-breaks out of the competitor scan as soon as one competitor
        // ties or beats the segment.  With the dump on, the scan runs to
        // completion so the TRUE nearest-competitor distance is recorded --
        // the decision is unchanged, because
        //     drop  <=>  exists other with dis <= min_dis_track
        //           <=>  min-over-others <= min_dis_track,
        // which is what the full scan evaluates.  One line per fresh
        // (segment, cell) decision; cache hits are not re-dumped, so each
        // decision appears exactly once.  Lets the recovery curve of a
        // proposed tie margin or a larger keep floor be computed offline,
        // from one run per detector, without building either knob.
        FILE* m_excl_dump{nullptr};
        WireCell::Point m_excl_dump_pt;   // the fit point whose association is being trimmed
        bool  exclusion_keep_dump(const std::shared_ptr<PR::Segment>& segment,
                                  const std::vector<std::shared_ptr<PR::Segment>>& all_segments,
                                  const WireCell::Point& test_point, const Coord2D& coord,
                                  int apa, int face, int plane, double min_dis_track);

        // doc pr/49 round 3 (fit_blob_coverage knob): clusters owning a
        // segment in the current fit -- the "fitting scope" whose members
        // never count as foreign in is_cell_covered_by_foreign_blobs.
        // Rebuilt by rebuild_cov_fit_scope() at the top of each
        // form_map/form_map_graph call while the knob is on; membership-only
        // (.count()), never iterated, so pointer keying cannot affect
        // determinism.  NOT m_clusters: that set is polluted by
        // preload_clusters() with charge-cache-only clusters that own no
        // segment.
        std::set<const Facade::Cluster*> m_cov_fit_scope;

        // doc pdvd/28 T1 (PDVD PR perf round 1): per-grouping coverage index
        // for is_cell_covered_by_foreign_blobs.  For every (apa,face), the
        // grouping's clusters owning at least one blob there, in
        // grouping->children() order, with the bounding box of those blobs
        // in (slice index, u/v/w wire index).  The box is a NECESSARY
        // condition of the per-blob predicate is_cell_covered_by_own_blobs,
        // so walking the index, skipping by box and re-applying the
        // unchanged predicate visits a superset of the covering clusters in
        // the original order: same first claimant, same OR.  Validity: a
        // (cluster pointer, nblobs) signature of children() is re-checked at
        // every rebuild_cov_fit_scope() (top of each fit pass; clusters are
        // never mutated inside a pass) and grouping + child count at every
        // use; reset_for_new_event() drops it.  Never iterated by pointer
        // (children() order only), so determinism is unaffected.
        struct CovBox {
            const Facade::Cluster* cluster = nullptr;
            int smin = 0, smax = 0;        // min slice_index_min .. max slice_index_max
            int wmin[3] = {0, 0, 0};       // per-plane min wire_index_min
            int wmax[3] = {0, 0, 0};       // per-plane max wire_index_max
        };
        struct CovIndex {
            const Facade::Grouping* grouping = nullptr;
            std::vector<std::pair<const Facade::Cluster*, size_t>> signature;
            std::map<std::pair<int, int>, std::vector<CovBox>> by_face;  // (apa,face) -> children order
        };
        mutable CovIndex m_cov_index;
        void ensure_cov_index(const Facade::Grouping* grouping) const;
        // WCT_TF_COVERAGE_CENSUS=1: log-only counters (calls, index entries
        // visited, box survivors, covering clusters), printed at each index
        // rebuild and at destruction.  Unset => no counting, no log line.
        mutable size_t m_cov_census[4] = {0, 0, 0, 0};

        // doc pr/50 (172230-class near-vertex robustness): graph vertex
        // positions (fit point when valid, else wcpt) and degrees at the
        // last rebuild_cov_fit_scope() call, in ordered_nodes order.  Used
        // by the deweight sentinel to report each firing's distance to the
        // nearest pattern vertex (and that vertex's degree) -- the census
        // evidence for classifying near-vertex vs far-ghost deweights.
        // Refreshed only while the fit_blob_coverage knob is on; read-only
        // positional data, never iterated by pointer.
        std::vector<std::pair<WireCell::Point, int>> m_cov_vtx_info;

        // Option 1: per-cluster edge descriptor cache to avoid full graph traversal
        std::unordered_map<Facade::Cluster*, std::vector<PR::edge_descriptor>> m_cluster_edges;
        std::vector<PR::edge_descriptor> m_all_edges;  ///< All segment edges in m_graph
        std::vector<PR::node_descriptor> m_ordered_nodes_vec;  ///< Nodes sorted by index, cached by build_cluster_edges
        void build_cluster_edges();                    ///< Rebuild m_cluster_edges, m_all_edges, and m_ordered_nodes_vec from m_graph
        const std::vector<PR::edge_descriptor>& get_segment_edges() const; ///< Return edges for current filter

        // Option 2: per-cluster charge data cache to avoid iterating full m_charge_data
        std::unordered_map<Facade::Cluster*, std::unordered_map<CoordReadout, ChargeMeasurement, CoordReadoutHash>> m_cluster_charge_data;

        std::set<Facade::Blob*> m_blobs;

        // input segment
        std::set<std::shared_ptr<PR::Segment> > m_segments;

        // input graph
        std::shared_ptr<PR::Graph> m_graph{nullptr};

        // Neutrino pattern-recognition results (set by TaggerCheckNeutrino)
        PR::VertexPtr        m_main_vertex{nullptr};
        PR::IndexedShowerSet m_showers;
        // doc sbnd_xin/docs/pr/40 round 9 B2 -- see the accessor comment.
        std::set<int>        m_bridged_cluster_ids;
        // doc sbnd_xin/docs/pr/92 -- see the accessor comment.
        std::set<int>        m_dropped_satellite_shower_ids;

        // Pi0 identification results (set by TaggerCheckNeutrino via set_pi0_data)
        PR::IndexedShowerSet                      m_pi0_showers;
        PR::ShowerIntMap                          m_map_shower_pio_id;
        std::map<int, std::vector<PR::ShowerPtr>> m_map_pio_id_showers;
        std::map<int, std::pair<double, int>>     m_map_pio_id_mass;

        // Kinematics and tagger features (set by TaggerCheckNeutrino)
        PR::KineInfo   m_kine_info{};
        PR::TaggerInfo m_tagger_info{};
        // doc sbnd_xin/docs/pr/75 -- diagnostic only, empty unless the
        // vertex_scoreboard knob was on.
        PR::VertexScoreboard m_vertex_scoreboard{};

        // =====================================================================
        // HYBRID CACHE IMPLEMENTATION
        // =====================================================================
        
        // Key types for caching
        using PlaneKey = std::tuple<int, int, int>;    // (apa, face, plane)
        using WireKey = std::tuple<int, int, int, int>; // (apa, face, plane, wire)
        
        // Hot cache: frequently accessed plane mappings (full plane cached)
        mutable std::map<PlaneKey, std::vector<int>> m_hot_cache;
        
        // Cold cache: individual wire lookups
        mutable std::map<WireKey, int> m_cold_cache;
        
        // Access frequency tracking
        mutable std::map<PlaneKey, int> m_access_count;
        
        // Cache statistics
        mutable CacheStats m_cache_stats = {0, 0, 0, 0, 0};
        
        // Configuration
        static constexpr int HOT_THRESHOLD = 50; // Access count to promote to hot cache
        
        // Helper methods
        void cache_entire_plane(int apa, int face, int plane) const;
        int fetch_channel_from_anode(int apa, int face, int plane, int wire) const;
    

        // ----------------------------------------
        // Internal Storage
        // ----------------------------------------
        std::unordered_map<CoordReadout, ChargeMeasurement, CoordReadoutHash> m_charge_data;  ///< Internal charge data storage using ChargeMeasurement struct
        std::unordered_map<CoordReadout, ChargeMeasurement, CoordReadoutHash> m_orig_charge_data; // saved original charge measurement, if modified

        std::map<Coord2D, std::set<int>> m_2d_to_3d;  ///< Internal 2D→3D mapping
        std::map<int, Point3DInfo> m_3d_to_2d;               ///< Internal 3D→2D mapping
    
        // Global (apa, time, channel) to blobs
        std::unordered_map<CoordReadout, std::unordered_set<Facade::Blob*>, CoordReadoutHash> global_rb_map;

        // Fitted 2D charge organized by (apa, face, plane) -> (wire, time)
        std::map<APAFacePlane, std::map<WireTime, FittedCharge2D>> m_fitted_charge_2d;

        /// Per-cluster snapshots of m_fitted_charge_2d captured at the end of
        /// every fill_fitted_charge_2d() call when m_cluster_filter is set.
        /// Refitting the same cluster replaces its entry in place, so this is
        /// still "latest fit wins per cluster" across pattern recognition.
        /// Merged into m_fitted_charge_2d by assemble_fitted_charge_2d().
        ///
        /// doc pr/109 §8: this was a std::map keyed by Facade::Cluster* and
        /// ORDERED BY IDENT (PR::ClusterPtrCmp).  Ident order was needed --
        /// the merge is last-writer-wins on shared cells and a pointer-ordered
        /// walk made pred_charge run-dependent (10.2% of cells moved between
        /// two `setarch -R` runs of SBND evt 388, doc pr/28 §4.3, §14) -- but
        /// ident order also made two live clusters that share an ident compare
        /// EQUAL, so the earlier snapshot was silently discarded.  That fired
        /// on uBooNE 5384-6528 and cost the main cluster its entire prediction
        /// (T_proj_data cid 19: 2141 cells, Σpred = 0).  A vector in capture
        /// order keeps both; assemble_fitted_charge_2d() still walks them in
        /// (ident, capture index) order, so the merged map is unchanged
        /// wherever no collision occurs.
        std::vector<ClusterFitted2D> m_cluster_fitted_charge_2d;

        // global geometry

        void BuildGeometry();
        void sync_from_graph();

        std::map<WirePlaneId , std::tuple<WireCell::Point, double, double, double>> wpid_params;
        std::map<WirePlaneId, std::pair<WireCell::Point, double> > wpid_U_dir;
        std::map<WirePlaneId, std::pair<WireCell::Point, double> > wpid_V_dir;
        std::map<WirePlaneId, std::pair<WireCell::Point, double> > wpid_W_dir;
        std::set<int> apas;

        // Time_width, Pitch_u, pitch_v, pitch_w, for each apa/face
        std::map<WirePlaneId, std::tuple<double, double, double, double >> wpid_geoms;

        // geometry information T, U, V, W for each apa/face
        std::map<WirePlaneId, std::tuple<double, double, double, double >> wpid_offsets;
        // T, slope_yu slope_zu, slope_yv slope_zv, slope_yw slope_zw 
        std::map<WirePlaneId, std::tuple<double, std::pair<double, double>, std::pair<double, double>, std::pair<double, double> >> wpid_slopes;

        // result
        std::vector<std::pair<WireCell::Point, std::shared_ptr<PR::Segment>>> fine_tracking_path;
        std::vector<double> dQ;
        std::vector<double> dx;
        std::vector<double> pu;
        std::vector<double> pv;
        std::vector<double> pw;
        std::vector<double> pt;
        std::vector<std::pair<int, int>> paf;
        std::vector<double> reduced_chi2;
    };

} // namespace WireCell::Clus

#endif // WIRECELLCLUS_TRACKFITTING_H