#include "WireCellClus/PRGraph.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/PRShower.h"
#include "WireCellClus/NeutrinoTaggerInfo.h"
#include "WireCellClus/IClusGeomHelper.h"
#include "WireCellClus/PRSegmentFunctions.h"

#include <array>
#include <cmath>

namespace WireCell::Clus::PR {

    /** BDT input features for the best pi0 candidate found by id_pi0_with_vertex (flag==1)
     *  or id_pi0_without_vertex (flag==2).  All angles in radians, energies and distances
     *  in WireCell internal units (multiply by 1/units::MeV, 1/units::cm for output).
     *  Default-initialised to zero (flag==0 means no pi0 found).
     */
    /// Type aliases for the pre-collected 2D charge maps used by cal_kine_charge.
    /// Collect once per shower_clustering_with_nv call via TrackFitting::collect_2D_charge()
    /// and reuse across all showers/merges to avoid O(N_hits) re-collection per shower.
    using ChargeMap = std::map<TrackFitting::CoordReadout, TrackFitting::ChargeMeasurement>;
    using WireMap   = std::map<std::pair<int,int>, std::vector<std::tuple<int,int,int>>>;

    /// Calibration constants of the charge -> kinetic-energy conversion in
    /// NeutrinoEnergyReco.cxx (cal_kine_charge and calculate_shower_kinematics).
    /// Every default below is the uBooNE-tuned number that used to be a literal
    /// at the three duplicated fudge/recom blocks and inside
    /// kine_charge_from_maps -- absent config keys are therefore byte-identical.
    /// None of them has been re-derived for any other detector; they are exposed
    /// so that a calibration can move them without a code change (see
    /// wcp-porting-img sbnd_xin/docs/pr/2 sec. 2e(iii)).
    ///
    /// The conversion is
    ///     E = sum_plane(w_p * Q_p) / sum(w) / recom / fudge * w_value[eV] * 1e-6 MeV
    /// with the (med,max) plane-asymmetry switch described at plane_asym_switch.
    struct KineChargeOptions {
        /// Average recombination survival fraction and residual scale factor for
        /// a track-like object.  uBooNE values; field- and calibration-dependent.
        double fudge_factor{0.95};
        double recom_factor{0.7};
        /// Same pair for an object flagged shower-like (higher local dE/dx =>
        /// more recombination, hence the smaller survival fraction).
        double shower_fudge_factor{0.8};
        double shower_recom_factor{0.5};
        /// Recombination survival for |pdg| == 2212 (proton).  The fudge factor
        /// deliberately stays at fudge_factor for protons -- only recombination
        /// is specialised, as in the prototype.
        double proton_recom_factor{0.35};
        /// Per-plane weights {U, V, W} of the charge average.  The uBooNE
        /// 0.25/0.25/1.0 encodes "induction planes are a quarter as trustworthy
        /// as collection".  A zero weight is legitimate ("ignore this plane");
        /// only an all-zero triple is rejected.
        std::array<double, 3> plane_weights{{0.25, 0.25, 1.0}};
        /// If the two largest per-plane charges disagree by more than this
        /// relative asymmetry, the three-plane weighted average is replaced by
        /// the (median, minimum) pair -- i.e. the largest plane is dropped as
        /// likely contaminated.  Dimensionless, uBooNE-tuned.
        double plane_asym_switch{0.04};
        /// Argon W-value in eV per electron-ion pair.  Duplicated here rather
        /// than read from the recombination model's Wi (docs/pr/2 sec. 2e(iii),
        /// "low risk, note only") -- if the model's Wi is ever wired in, this
        /// knob is what it replaces.
        double w_value{23.6};
    };

    struct Pi0KineFeatures {
        int    flag{0};      ///< 0=none, 1=with_vertex, 2=without_vertex
        double mass{0};      ///< reconstructed pi0 invariant mass
        double vtx_dis{0};   ///< distance from pi0 decay vertex to main vertex
        double energy_1{0};  ///< kine_charge of first shower
        double theta_1{0};   ///< polar angle of first shower direction (from z-axis)
        double phi_1{0};     ///< azimuthal angle of first shower direction
        double dis_1{0};     ///< distance from pi0 vertex to first shower start point
        double energy_2{0};  ///< kine_charge of second shower
        double theta_2{0};   ///< polar angle of second shower direction
        double phi_2{0};     ///< azimuthal angle of second shower direction
        double dis_2{0};     ///< distance from pi0 vertex to second shower start point
        double angle{0};     ///< opening angle between the two shower directions
    };

    class PatternAlgorithms{
        public:
        bool m_perf{false};  // if true, print per-step timing to stdout

        // Direction-weakness read fidelity switch (default false = legacy port
        // behavior, byte-identical).  The prototype's ONLY public accessor is the
        // score-thresholded ProtoSegment::is_dir_weak() (the raw member is
        // private); the port read the raw flag at every one of the ~83 sites
        // instead.  When true, all reads route through segment_is_dir_weak()
        // (PRSegmentFunctions.cxx), restoring the prototype semantics: a mu/p
        // segment whose particle_score is poor (or still the sentinel 100) is
        // treated as weak.  See wcp-porting-img sbnd_xin/docs/pr/6.
        bool m_dir_weak_use_score{false};
        // All PatternAlgorithms code must read direction weakness through this,
        // never seg->dir_weak() directly (flag propagation in break_segment is
        // the one intentional raw read, and it lives in PRSegmentFunctions.cxx).
        bool seg_dir_weak(SegmentPtr seg) const;

        // MIP dQ/dx scales, INTERNAL units (charge per WCT length unit).
        // Two distinct uBooNE-legacy roles, previously hard-coded everywhere
        // (see wcp-porting-img sbnd_xin/docs/pr/7 sec. 5 and pr/8):
        //  - m_mip_dqdx        (legacy 50e3 e/cm): the flat-template amplitude in
        //    the do_track_comp KS comparison and the segment_cal_4mom /
        //    segment_is_shower_trajectory scale.
        //  - m_mip_dqdx_median (legacy 43e3 e/cm): the scale every median-dQ/dx
        //    ratio threshold (x1.2/1.3/1.4/1.75...) and BDT normalization is
        //    quoted against.
        // Defaults = the uBooNE numbers => absent config keys are byte-identical.
        // SBND production sets 56000 / 48000 e/cm (owner 2026-07-30; 48000 is a
        // placeholder scaled by the uBooNE 43/50 ratio, pending an SBND MIP
        // median measurement).
        double m_mip_dqdx{50000/units::cm};
        double m_mip_dqdx_median{43000/units::cm};

        // Long shower-topology demote length (doc sbnd_xin/docs/pr/25 sec 3).
        // Passed to every segment_is_shower_topology call site: a segment the
        // topology test would flag kShowerTopology whose geometric length
        // exceeds this is demoted to a track instead, so it receives real
        // track PID rather than the hard-coded pdg=11/score=100.
        // Motivation: the test's only measurement axis satisfies
        // dir_3.xhat = sin(angle-to-drift), so for a near-isochronous segment
        // it measures spread along the DRIFT axis, where points sit on a
        // 0.313 cm time-slice lattice -- against a 0.4 cm "large spread" cut.
        // 86 of 91 long firings across a 572-event manifest carry no evidence
        // above that noise floor, and an owner hand-scan of all 10 such
        // segments on a selected nu-candidate main cluster found 10/10 tracks.
        // C++ default 0 => the guard never fires => byte-identical.
        // 50*units::cm reproduces the scan-supported rule for 9 of the 10
        // scanned events; ~45*units::cm covers all 10 (evt 400504 measures
        // 49.1 cm by segment_track_length, the length this guard uses).
        double m_shower_topo_demote_len{0};

        // Isochronous first-segment endpoint finding (doc sbnd_xin/docs/pr/24
        // round 2, SBND evt 271851).  When ON and the cluster is sheet-like
        // (long, with a small quantile-trimmed drift-x extent), the first PR
        // segment's endpoints are picked as the quantile-trimmed extremes
        // along the sheet's principal axis (each snapped to a Steiner
        // terminal) instead of the wire-footprint boundary metric, and the
        // local-PCA endpoint refinement (degenerate on a filled 2-D sheet) is
        // skipped on that branch.  No prototype counterpart; the prototype's
        // boundary metric (PR3DCluster_path.h:530-536, wire-delta coefficients
        // *0.0) is preserved byte-identically when off.  Precedents:
        // adjust_wcpoints_parallel (data/src/PR3DCluster.cxx:428, separation
        // only) and search_for_connection_isochronous
        // (pid/src/PR3DCluster_graph.h:1445, call site disabled upstream).
        // C++ default false => legacy path byte-identical.
        bool   m_iso_endpoint{false};
        double m_iso_endpoint_min_length{40 * units::cm};
        double m_iso_endpoint_max_xext{25 * units::cm};
        double m_iso_endpoint_xext_frac{0.35};
        double m_iso_endpoint_xext_quantile{0.02};
        // Round 3 (doc pr/24 sec 15).  Round 2 picked each endpoint as the band
        // point nearest the centroid of a 3 cm band at the QUANTILE-TRIMMED
        // axial extreme -- inward-biased by construction (8-15 cm on a
        // 7000-point track), which left a tip stub for find_other_segments to
        // claim as its own segment and put a 0.9 deg "vertex" in a straight
        // track (SBND mcp1k evt 284794, 59899).  Instead: the endpoint is the
        // TRUE, untrimmed axial extreme over all qualified points (so it can
        // never move inward -- "leaves no stub" holds by construction), and
        // lateral centring is applied only WITHIN the 3 cm end band, which is
        // what still keeps the pick off a sheet's edge corner.
        // m_iso_endpoint_min_aspect additionally requires real 2-D sheet-ness
        // (trimmed transverse/axial extent ratio), so the branch is inert on
        // 1-D track-like clusters rather than merely harmless.
        // m_iso_endpoint_tube_radius is DIAGNOSTIC ONLY: the DEBUG line reports
        // whether the chosen endpoint lies within it of the axis line.  It was
        // a hard filter in a first draft; measured on the 38 gated clusters
        // that filter pulled endpoints up to 28.6 cm inward on long/curved
        // objects (doc pr/24 sec 15), the same failure it was meant to cure.
        double m_iso_endpoint_tube_radius{4 * units::cm};
        double m_iso_endpoint_min_aspect{0.12};

        // Cathode kink veto (doc sbnd_xin/docs/pr/20 Part II, B0).  Passed to
        // segment_search_kink from break_segments: a candidate fit point within
        // m_cathode_kink_xcut of the cathode plane at m_cathode_x is skipped, so
        // the ~2 cm transverse cathode mismatch cannot invent a vertex that breaks
        // a crossing cosmic into two particles.  C++ default 0 => no point is ever
        // skipped => byte-identical to the pre-pr/20 behavior.
        double m_cathode_x{0};
        double m_cathode_kink_xcut{0};

        // Proton-template direction vote (doc pr/8; default false = legacy).
        // Thresholds are initial values pending the pr/8 sec. 6 calibration.
        bool   m_proton_dir_vote{false};
        double m_proton_dir_score_max{0.25};
        double m_proton_dir_asym_min{1.3};
        // Endpoint-trim retry (doc sbnd_xin/docs/pr/9 sec. 6 F1).  C++ default
        // false => legacy abstention path, byte-identical.
        bool   m_endpoint_trim_retry{false};
        // Minimum WCPT-path length for a segment to count as a leg in the
        // fit_vertex position fit (doc sbnd_xin/docs/pr/9 sec. 11 F3c, owner
        // 2026-07-30, deliberate prototype divergence): very short vertex-
        // activity stubs otherwise drag the fitted vertex (evt 172230: a
        // 0.62 cm drift-blur stub moved the true vertex 2.4 cm).  Measured in
        // wcpt space (a fresh stub's FIT cloud spreads past 1 cm).  If >=3
        // legs survive, the fit runs on the survivors; if the exclusion
        // leaves <=2, the vertex fit is skipped (an effectively two-leg
        // vertex was already fit by the plain two-leg pass, and refitting on
        // stub-contaminated re-tracked clouds reproduces the drag).  Excluded
        // segments stay in the graph and the particle flow.
        // C++ default 0 => no filtering, legacy byte-identical.
        double m_fit_vertex_min_seg_length{0};

        // ---- Detector-extent literals (doc sbnd_xin/docs/pr/2 sec. 2e(iv)) ----
        // The uBooNE active volume the prototype was written against is
        // y in [-116, +117] cm, z in [0, 1037] cm, x in [0, 256] cm
        // (prototype_base/pid/apps/wire-cell-prod-nue.cxx:417).  Every default
        // below is the prototype literal, so absent config keys are
        // byte-identical.  Values are INTERNAL units.
        //
        // cosmic_tagger() "reaches the top of the detector" tests, part D
        // (prototype NeutrinoID_cosmic_tagger.h:594-796).  Each is a fixed
        // DISTANCE BELOW THE TOP FACE in uBooNE terms (117 - cut), which is
        // what a downward cosmic entering the top satisfies regardless of
        // detector height -- so a taller detector moves them up, it does not
        // scale them.  SBND (top y = +200 cm) uses 183 / 185 / 163 / 133.
        double m_cosmic_y_top_main{100 * units::cm};    // :1175, 17 cm below top:
                                                        // the MAIN cluster's own highest
                                                        // point, relaxing the vertical-angle
                                                        // cut from 20 to 30 degrees.
        double m_cosmic_y_top_strict{102 * units::cm};  // :1190, 15 cm below top: event
                                                        // highest point, single-cosmic branch.
        double m_cosmic_y_top_loose{80 * units::cm};    // :1191, 37 cm below top: event
                                                        // highest point, global gate on the
                                                        // whole flagp_cosmic decision.
        double m_cosmic_y_small_piece{50 * units::cm};  // :1073, 67 cm below top: a <3 cm
                                                        // cluster counts as cosmic debris
                                                        // (acc_small_length) only if its PCA
                                                        // centre sits in the upper region.
        // Upstream-z preference in compare_main_vertices (:875) and
        // compare_main_vertices_global (:3001): each candidate is penalised by
        // (z - min_z) / m_vertex_z_prior_scale, competing with the +0.25-per-
        // track bonuses.  The uBooNE 200 cm gives ~5.2 penalty units across the
        // 1037 cm detector.  On a shorter detector the same 200 cm shrinks the
        // prior's dynamic range, which matters most in the _global comparison
        // (candidates from DIFFERENT clusters of the beam bundle, separations up
        // to the full detector length).  SBND uses 100 cm ~ 200 x 501/1037.
        double m_vertex_z_prior_scale{200 * units::cm};

        // Bundle the do_track_pid-related knobs for the free segment functions.
        TrackPidOptions track_pid_options() const {
            TrackPidOptions o;
            o.mip_dqdx = m_mip_dqdx;
            o.proton_dir_vote = m_proton_dir_vote;
            o.proton_dir_score_max = m_proton_dir_score_max;
            o.proton_dir_asym_min = m_proton_dir_asym_min;
            o.endpoint_trim_retry = m_endpoint_trim_retry;
            return o;
        }

        // Beam-line reference directions used by ssm_tagger, in the DETECTOR
        // frame.  They feed the 8 ssm_*_angle_{target,absorber} BDT features
        // via safe_acos(dir.dot(ref)) at NeutrinoTaggerSSM.cxx:908-909 (the
        // 10 cm initial-direction pair) and :1124-1125 (the nu/con_nu/prim_nu/
        // track momentum pairs).
        // Defaults = the prototype's hard-coded uBooNE numbers
        // (NeutrinoID_ssm_tagger.h): BNB target and NuMI absorber as seen from
        // MicroBooNE.  Absent config keys => byte-identical.  SBND needs its
        // own BNB-target vector and the NuMI-absorber features have no obvious
        // SBND meaning -- see wcp-porting-img sbnd_xin/docs/pr/2 sec. 2e(i);
        // the value itself is still an open gap, only its location moved.
        //
        // NOTE these are NOT unit vectors: |target| = 0.99866,
        // |absorber| = 1.00970.  safe_acos() clamps the dot product, so the
        // absorber angle saturates at exactly 0 for true angles out to ~7.9
        // deg.  Kept verbatim for prototype parity -- whoever supplies a
        // properly normalized SBND vector should expect these feature
        // distributions to shift even at the same physical direction.
        WireCell::Vector m_ssm_target_dir{0.46, 0.05, 0.885};
        WireCell::Vector m_ssm_absorber_dir{0.33, 0.75, -0.59};

        // Charge -> kinetic-energy calibration constants (NeutrinoEnergyReco.cxx).
        // Defaults = the uBooNE literals they replaced => byte-identical when the
        // config keys are absent.  See KineChargeOptions above.
        KineChargeOptions m_kine_charge{};

        // Muon median-dQ/dx-vs-length envelope {c0, c1, pivot, power}:
        //     dQ_dx_cut = c0 + c1 * pow(pivot / L, power)
        // a dimensionless multiple of m_mip_dqdx_median.  The default is the
        // prototype's empirical uBooNE stopping-muon refit (of the commented-out
        // 0.85 + 0.95*sqrt(25cm/L) family): short tracks are Bragg-dominated so
        // their MEDIAN dQ/dx sits well above MIP (2.5x at 5 cm), relaxing toward
        // the 0.8866 plateau asymptote.  Appears at NINE sites (numu x2,
        // vertex-finder, nue x4, ssm, cosmic), all as literals before this knob
        // (doc pr/2 sec 2e(iv) undercounted 3).  pivot is INTERNAL units.
        // Defaults = the uBooNE literals => absent config keys are
        // byte-identical.
        std::array<double, 4> m_muon_dqdx_curve{{0.8866, 0.9533, 18 * units::cm, 0.4234}};

        // The envelope, length in INTERNAL units.  Inline and shaped exactly
        // like the literal expression it replaced (c0 + c1*pow(180.0/L, c3)),
        // so the defaults are bit-identical.
        double muon_dqdx_cut(double length) const {
            return m_muon_dqdx_curve[0] +
                   m_muon_dqdx_curve[1] * std::pow(m_muon_dqdx_curve[2] / length, m_muon_dqdx_curve[3]);
        }
        // Same envelope for a length already in cm (the ssm_tagger convention,
        // whose sg_length values are divided by units::cm at the source).
        // pivot/units::cm = 18.0 exactly, so this too is bit-identical.
        double muon_dqdx_cut_cm(double length_cm) const {
            return m_muon_dqdx_curve[0] +
                   m_muon_dqdx_curve[1] * std::pow((m_muon_dqdx_curve[2] / units::cm) / length_cm, m_muon_dqdx_curve[3]);
        }

        // Single-photon stem dE/dx conversion (NeutrinoTaggerSinglePhoton.cxx).
        // Default false = the inline float inverse Modified-Box frozen at
        // uBooNE's field (A=1.0, B=0.255, rho=1.38, E=0.273 kV/cm), kept for
        // byte-identity.  When true, shw_sp_vec_{median,mean}_dedx go through
        // m_recomb_model->dE() instead -- the same configured model the rest of
        // the chain uses (docs/pr/2 sec 2e(i) third correctness item).
        bool m_sp_dedx_use_recomb_model{false};
        // The hard mip_quality cut on shw_sp_vec_mean_dedx (MeV/cm).  Stored as
        // float so the default compares bit-identically to the legacy 2.3f
        // literal.  uBooNE-tuned against the inline (0.273 kV/cm) dE/dx scale;
        // coupled to the knob above -- retune when switching models.
        float m_sp_mean_dedx_cut{2.3f};
        // The configured recombination model, set by TaggerCheckNeutrino::visit()
        // alongside the kine transfer.  Null when the component has none.
        IRecombinationModel::pointer m_recomb_model{};

        // 2D charge maps cached for the duration of shower_clustering_with_nv.
        // Populated once by collect_charge_maps(); reused by calculate_shower_kinematics
        // and all cal_kine_charge call sites to avoid O(N_hits) re-collection per shower.
        ChargeMap m_charge_2d_u, m_charge_2d_v, m_charge_2d_w;
        WireMap   m_map_apa_ch_plane_wires;

        // Populate the cached charge maps from track_fitter.
        // Call once at the start of shower_clustering_with_nv; the maps are then
        // valid for the entire call tree beneath it.
        void collect_charge_maps(TrackFitting& track_fitter);
        std::vector<VertexPtr> find_cluster_vertices(Graph& graph, const Facade::Cluster& cluster);
        std::vector<SegmentPtr> find_cluster_segments(Graph& graph, const Facade::Cluster& cluster);
        bool clean_up_graph(Graph& graph, const Facade::Cluster& cluster);

        SegmentPtr init_first_segment(Graph& graph, Facade::Cluster& cluster, Facade::Cluster* main_cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, bool flag_back_search = true);

        // Isochronous first-segment endpoints (m_iso_endpoint; doc
        // sbnd_xin/docs/pr/24 round 2).  Returns true and fills p1/p2 (both
        // snapped to Steiner terminals of "steiner_pc") when the cluster
        // passes the sheet gate; false => caller uses the legacy boundary
        // search unchanged.
        bool find_iso_first_segment_endpoints(const Facade::Cluster& cluster, Facade::geo_point_t& p1, Facade::geo_point_t& p2) const;

        // find the shortest path using steiner graph
        std::vector<Facade::geo_point_t> do_rough_path(const Facade::Cluster& cluster,Facade::geo_point_t& first_point, Facade::geo_point_t& last_point);
        std::vector<Facade::geo_point_t> do_rough_path_reg_pc(const Facade::Cluster& cluster, Facade::geo_point_t& first_point, Facade::geo_point_t& last_point,  std::string graph_name = "relaxed_pid");
        // create a segment given a path
        SegmentPtr create_segment_for_cluster(WireCell::Clus::Facade::Cluster& cluster, IDetectorVolumes::pointer dv, const std::vector<Facade::geo_point_t>& path_points, int dir = 0);
        // create a segment given two vertices, null, if failed
        SegmentPtr create_segment_from_vertices(Graph& graph, Facade::Cluster& cluster, VertexPtr v1, VertexPtr v2, IDetectorVolumes::pointer dv);
        // replace a segment and vertex with another segment and vertex, assuming the original vertex only connect to this segment
        bool replace_segment_and_vertex(Graph& graph, SegmentPtr& seg, VertexPtr& vtx, std::list<Facade::geo_point_t>& path_point_list, Facade::geo_point_t& break_point, IDetectorVolumes::pointer dv);
        bool replace_segment_and_vertex(Graph& graph, SegmentPtr& seg, VertexPtr old_vertex, VertexPtr new_vertex, IDetectorVolumes::pointer dv);
        bool break_segment_into_two(Graph& graph, VertexPtr vtx1, SegmentPtr seg, VertexPtr vtx2, std::list<Facade::geo_point_t>& path_point_list1, Facade::geo_point_t& break_point, std::list<Facade::geo_point_t>& path_point_list2, IDetectorVolumes::pointer dv, SegmentPtr& out_seg2);


        // return the point and its index in the steiner tree as a pair
        std::pair<Facade::geo_point_t, size_t> proto_extend_point(const Facade::Cluster& cluster, Facade::geo_point_t& p, Facade::geo_vector_t& dir, Facade::geo_vector_t& dir_other, bool flag_continue, std::vector<Facade::geo_point_t>* walk_history = nullptr);
        // return Steiner Graph path in wcps_list1 and wcps_list2
        bool proto_break_tracks(const Facade::Cluster& cluster, const Facade::geo_point_t& first_wcp, Facade::geo_point_t& curr_wcp, const Facade::geo_point_t& last_wcp, std::list<Facade::geo_point_t>& wcps_list1, std::list<Facade::geo_point_t>& wcps_list2, bool flag_pass_check);
        // breaking segments ...
        bool break_segments(Graph& graph, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, std::vector<SegmentPtr>& remaining_segments, float dis_cut = 0);
        // merge vertices within 0.1 cm after break_segments, then refit if changed
        bool merge_nearby_vertices(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        // merge two segments to one
        bool merge_two_segments_into_one(Graph& graph, SegmentPtr& seg1, VertexPtr& vtx, SegmentPtr& seg2, IDetectorVolumes::pointer dv);
        // merge vertex into another
        bool merge_vertex_into_another(Graph& graph, VertexPtr& vtx_from, VertexPtr& vtx_to, IDetectorVolumes::pointer dv);

        // get direction with  distance cut ... 
        Facade::geo_vector_t vertex_get_dir(VertexPtr& vertex, Graph& graph, double dis_cut = 5*units::cm);
        Facade::geo_vector_t vertex_segment_get_dir(VertexPtr& vertex, SegmentPtr& segment, Graph& graph, double dis_cut = 2*units::cm);


        // Structure examination
        void examine_structure(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv); // call examine_structure_1 and examine_structure_2      
        bool examine_structure_1(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool examine_structure_2(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);

        bool examine_structure_3(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool examine_structure_4(VertexPtr vertex, bool flag_final_vertex, Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);

        // identify other segments giving the graph ...
        void find_other_segments(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv , bool flag_break_track =true, double search_range = 1.5*units::cm, double scaling_2d = 0.8);
        VertexPtr find_vertex_other_segment(Graph& graph, Facade::Cluster& cluster, SegmentPtr seg, bool flag_forwrard, Facade::geo_point_t& wcp, TrackFitting& track_fitter, IDetectorVolumes::pointer dv );
        std::tuple<VertexPtr, SegmentPtr, Facade::geo_point_t> check_end_point(Graph& graph, Facade::Cluster& cluster, std::vector<Facade::geo_point_t>& tracking_path, bool flag_front = true, double vtx_cut1 = 0.9 * units::cm, double vtx_cut2 = 2.0 * units::cm, double sg_cut1  = 2.0 * units::cm, double sg_cut2  = 1.2 * units::cm);
        bool modify_vertex_isochronous(Graph& graph, Facade::Cluster& cluster, VertexPtr vtx, VertexPtr v1, SegmentPtr sg, VertexPtr v2, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool modify_segment_isochronous(Graph& graph, Facade::Cluster& cluster, SegmentPtr sg1, VertexPtr v1, SegmentPtr sg, VertexPtr v2, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, double dis_cut = 6*units::cm, double angle_cut = 15, double extend_cut = 15*units::cm);

        // examine segment
        void examine_segment(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool crawl_segment(Graph& graph, Facade::Cluster& cluster, SegmentPtr seg, VertexPtr vertex, TrackFitting& track_fitter, IDetectorVolumes::pointer dv );
        void examine_partial_identical_segments(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);

        //examine vertices
        void examine_vertices(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, VertexPtr main_vertex = nullptr);
        bool examine_vertices_1(Graph&graph, Facade::Cluster&cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, VertexPtr main_vertex = nullptr);
        bool examine_vertices_1p(Graph&graph, VertexPtr v1, VertexPtr v2, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool examine_vertices_2(Graph&graph, Facade::Cluster&cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, VertexPtr main_vertex = nullptr);
        bool examine_vertices_4(Graph&graph, Facade::Cluster&cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, VertexPtr main_vertex = nullptr);
        bool examine_vertices_4p(Graph&graph, VertexPtr v1, VertexPtr v2, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        Facade::geo_point_t get_local_extension(Facade::Cluster& cluster, const Facade::geo_point_t& wcp);
        void examine_vertices_3(Graph& graph, Facade::Cluster& main_cluster, std::pair<VertexPtr, VertexPtr> main_cluster_initial_pair_vertices, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);

        // master pattern recognition function
        bool find_proto_vertex(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, bool flag_break_track = true, int nrounds_find_other_tracks = 2, bool flag_back_search = true);
        
        void init_point_segment(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);

        // examine the structure of the patterns ... 
        bool examine_structure_final_1(Graph& graph, VertexPtr main_vertex, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool examine_structure_final_1p(Graph& graph, VertexPtr main_vertex, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool examine_structure_final_2(Graph& graph, VertexPtr main_vertex, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool examine_structure_final_3(Graph& graph, VertexPtr main_vertex, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool examine_structure_final(Graph& graph, VertexPtr main_vertex, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);

        // EM shower related
        void clustering_points(Graph& graph, Facade::Cluster& cluster, const IDetectorVolumes::pointer& dv, const std::string& cloud_name = "associate_points", double search_range = 1.2*units::cm, double scaling_2d = 0.7);
        void separate_track_shower(Graph&graph, Facade::Cluster& cluster);
        // Direction
        void determine_direction(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        std::pair<int, double> calculate_num_daughter_showers(Graph& graph, VertexPtr vertex, SegmentPtr segment, bool flag_count_shower = true);
        std::pair<int, double> calculate_num_daughter_tracks(Graph& graph, VertexPtr vtx, SegmentPtr sg, bool flag_count_shower = false, double length_cut = 0);
        // find_cont_muon_segment_nue: find adjacent segment continuing in same direction (MIP-like)
        std::pair<SegmentPtr, VertexPtr> find_cont_muon_segment_nue(Graph& graph, SegmentPtr sg, VertexPtr vtx, bool flag_ignore_dQ_dx = false);
        void examine_good_tracks(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data);
        // about fix maps
        void fix_maps_multiple_tracks_in(Graph& graph, Facade::Cluster& cluster);
        void fix_maps_shower_in_track_out(Graph& graph, Facade::Cluster& cluster);
        void improve_maps_one_in(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_strong_check = true);
        void improve_maps_shower_in_track_out(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_strong_check = true);
        void improve_maps_no_dir_tracks(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void improve_maps_multiple_tracks_in(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void judge_no_dir_tracks_close_to_showers(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, IDetectorVolumes::pointer dv);
        bool examine_maps(Graph&graph, Facade::Cluster& cluster);
        void examine_all_showers(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data);
        void shower_determining_in_main_cluster(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, IDetectorVolumes::pointer dv);
        void set_default_shower_particle_info(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);

        // PCA calculation
        std::pair<Facade::geo_point_t, Facade::geo_vector_t> calc_PCA_main_axis(std::vector<Facade::geo_point_t>& points);

        // vertex related functions 
        bool search_for_vertex_activities(Graph& graph, VertexPtr vertex, std::vector<SegmentPtr>& segments_set, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, double search_range = 1.5*units::cm);
        // existing_segments stays a plain std::set<SegmentPtr> (pointer-keyed)
        // ON PURPOSE.  It is iterated once, in the Case-5 block, where the body
        // only takes a running min over three distances -- order-insensitive.
        // Its find()/count() must stay POINTER identity: an index-ordered set
        // would compare by Segment::get_graph_index(), and that value is NOT
        // unique across live Segment objects.  PR::add_segment on a vertex pair
        // that already carries an edge takes the "edge already existed" path
        // (PRGraph.cxx:86-89): it overwrites g[desc].segment with the new
        // segment and copies the existing edge index into it.  The displaced
        // segment keeps that same index, so any SegmentPtr still held to it --
        // which is exactly what this set holds -- now compares EQUAL to a
        // different live segment.  Measured: making that swap moved
        // kine_reco_Enu on SBND evt 239794 from 2930 to 1687 MeV.
        bool eliminate_short_vertex_activities(Graph& graph, Facade::Cluster& cluster, VertexPtr main_vertex, std::set<SegmentPtr>& existing_segments, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        std::tuple<bool, int, int> examine_main_vertex_candidate(Graph& graph, VertexPtr vertex);
        VertexPtr compare_main_vertices_all_showers(Graph& graph, Facade::Cluster& cluster, std::vector<VertexPtr>& vertex_candidates, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        float calc_conflict_maps(Graph& graph, VertexPtr vertex);
        VertexPtr compare_main_vertices(Graph& graph, Facade::Cluster& cluster, std::vector<VertexPtr>& vertex_candidates);
        std::pair<SegmentPtr, VertexPtr> find_cont_muon_segment(Graph &graph, SegmentPtr sg, VertexPtr vtx, bool flag_ignore_dQ_dx = false);
        bool examine_direction(Graph& graph, VertexPtr vertex, VertexPtr main_vertex, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_final = false);

        bool fit_vertex(Facade::Cluster& cluster, VertexPtr vertex, VertexPtr main_vertex, std::vector<SegmentPtr>& sg_set, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        void improve_vertex(Graph& graph, Facade::Cluster& cluster, VertexPtr& main_vertex, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_search_vertex_activity = true , bool flag_final_vertex = false);
        void determine_main_vertex(Graph& graph, Facade::Cluster& cluster, VertexPtr& main_vertex, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void change_daughter_type(Graph& graph, VertexPtr vertex, SegmentPtr segment, int particle_type, double mass, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void examine_main_vertices_local(Graph& graph, std::vector<VertexPtr>& vertices, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);

        // cluster functions ...
        Facade::geo_vector_t calc_dir_cluster(Graph& graph, Facade::Cluster& cluster, const Facade::geo_point_t& orig_p, double dis_cut);
        Facade::Cluster* swap_main_cluster(Facade::Cluster& new_main_cluster, Facade::Cluster& old_main_cluster, std::vector<Facade::Cluster*>& other_clusters);
        void examine_main_vertices(Graph& graph, ClusterVertexMap& map_cluster_main_vertices, Facade::Cluster*& main_cluster, std::vector<Facade::Cluster*>& other_clusters);

        VertexPtr compare_main_vertices_global(Graph& graph, std::vector<VertexPtr>& vertex_candidates, Facade::Cluster& main_cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        Facade::Cluster* check_switch_main_cluster(Graph& graph, ClusterVertexMap map_cluster_main_vertices, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        Facade::Cluster* check_switch_main_cluster_2(Graph& graph, VertexPtr temp_main_vertex, Facade::Cluster* max_length_cluster, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters);
        VertexPtr determine_overall_main_vertex(Graph& graph, ClusterVertexMap map_cluster_main_vertices, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_dev_chain = true);

        // Deep-learning vertex refinement.  Returns true if the DL network changed
        // the selected vertex (in which case the traditional determine_overall_main_vertex
        // should NOT be called).  Returns false if the DL network was unavailable, failed,
        // or did not improve on the candidate vertices (fall back to traditional).
        bool determine_overall_main_vertex_DL(Graph& graph, ClusterVertexMap& map_cluster_main_vertices, Facade::Cluster*& main_cluster, std::vector<Facade::Cluster*>& other_clusters, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, const std::string& dl_weights, double dl_vtx_cut, double dQdx_scale = 0.1, double dQdx_offset = -1000.0, bool flag_rerank = false, int dl_vtx_top_k = 5, double dl_vtx_min_accept_score = 4.0, double dl_vtx_score_scale = 1000.0);

        // global information transfer
        void transfer_info_from_segment_to_cluster(Graph& graph, Facade::Cluster& cluster, const std::string& cloud_name = "associated_points");

        // print information
        void print_segs_info(Graph& graph, Facade::Cluster& cluster, VertexPtr vertex= nullptr);

        // shower related functions
        void update_shower_maps(IndexedShowerSet& showers, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters);
        void shower_clustering_with_nv_in_main_cluster(Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void shower_clustering_connecting_to_main_vertex(Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters);
        void shower_clustering_with_nv_from_main_cluster(Graph& graph, VertexPtr main_vertex, Facade::Cluster* main_cluster, IndexedShowerSet& showers, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters);
        void shower_clustering_with_nv_from_vertices(Graph& graph, VertexPtr main_vertex, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, IndexedShowerSet& showers, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void calculate_shower_kinematics(IndexedShowerSet& showers, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, Graph& graph, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void examine_merge_showers(IndexedShowerSet& showers, VertexPtr main_vertex, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, Graph& graph, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void shower_clustering_in_other_clusters(Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, ClusterVertexMap map_cluster_main_vertices, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_save = true);
        void examine_shower_1(Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, ClusterVertexMap map_cluster_main_vertices, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void examine_showers(Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, ClusterVertexMap map_cluster_main_vertices, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void id_pi0_with_vertex(int acc_segment_id, IndexedShowerSet& pi0_showers, ShowerIntMap& map_shower_pio_id, std::map<int, std::vector<ShowerPtr > >& map_pio_id_showers, std::map<int, std::pair<double, int> >& map_pio_id_mass, std::map<int, std::pair<int, int> >& map_pio_id_saved_pair, Pi0KineFeatures& pio_kine, Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, ClusterVertexMap map_cluster_main_vertices, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void id_pi0_without_vertex(int acc_segment_id, IndexedShowerSet& pi0_showers, ShowerIntMap& map_shower_pio_id, std::map<int, std::vector<ShowerPtr > >& map_pio_id_showers, std::map<int, std::pair<double, int> >& map_pio_id_mass, std::map<int, std::pair<int, int> >& map_pio_id_saved_pair, Pi0KineFeatures& pio_kine, Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, ClusterVertexMap map_cluster_main_vertices, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, IndexedSegmentSet& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void shower_clustering_with_nv(int acc_segment_id, IndexedShowerSet& pi0_showers, ShowerIntMap& map_shower_pio_id, std::map<int, std::vector<ShowerPtr > >& map_pio_id_showers, std::map<int, std::pair<double, int> >& map_pio_id_mass, std::map<int, std::pair<int, int> >& map_pio_id_saved_pair, Pi0KineFeatures& pio_kine, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, ClusterVertexMap map_cluster_main_vertices, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);


        // deghost related functions ...
        void order_clusters(Graph& graph, std::vector<Facade::Cluster*>& ordered_clusters, std::map<Facade::Cluster*, std::vector<SegmentPtr> >& map_cluster_to_segments, std::map<Facade::Cluster*, double>& map_cluster_total_length);
        void deghost_clusters(Graph& graph, std::vector<Facade::Cluster*>& all_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        void order_segments(std::vector<SegmentPtr>& ordered_segments, std::vector<SegmentPtr>& segments);
        void deghost_segments(Graph& graph, ClusterVertexMap map_cluster_main_vertices, std::vector<Facade::Cluster*>& all_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        void deghosting(Graph& graph, ClusterVertexMap& map_cluster_main_vertices, std::vector<Facade::Cluster*>& all_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv );

        // energy calculation ...
        double cal_corr_factor(WireCell::Point& pt, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);

        // kinematics output ...
        void init_tagger_info(TaggerInfo& ti);

        // cosmic tagger
        bool bad_reconstruction(Graph& graph, VertexPtr main_vertex, ShowerPtr shower,
                                bool flag_fill = false, TaggerInfo* ti = nullptr);
        bool cosmic_tagger(Graph& graph, VertexPtr main_vertex,
                           IndexedShowerSet& showers,
                           ShowerSegmentMap& map_segment_in_shower,
                           VertexShowerSetMap& map_vertex_to_shower,
                           IndexedSegmentSet& segments_in_long_muon,
                           Facade::Cluster* main_cluster,
                           std::vector<Facade::Cluster*>& all_clusters,
                           IDetectorVolumes::pointer dv,
                           TaggerInfo& ti);

        // numu tagger
        // count_daughters: counts track branches and total branches at the far end of a muon.
        // For a segment: BFS from the vertex closer to main_vertex through sg, count what is beyond.
        // For a shower:  finds the last segment of the long-muon chain, then counts from its near end.
        std::pair<int, int> count_daughters(Graph& graph, SegmentPtr sg, VertexPtr main_vertex);
        std::pair<int, int> count_daughters(Graph& graph, ShowerPtr shower, VertexPtr main_vertex,
                                            IndexedSegmentSet& segments_in_long_muon);
        // numu_tagger: fills TaggerInfo numu_cc_* BDT features and returns
        //   {flag_long_muon, max_muon_length}.
        // Prototype: WCPPID::NeutrinoID::numu_tagger() in NeutrinoID_numu_tagger.h.
        std::pair<bool, double> numu_tagger(Graph& graph,
                                            VertexPtr main_vertex,
                                            IndexedShowerSet& showers,
                                            IndexedSegmentSet& segments_in_long_muon,
                                            Facade::Cluster* main_cluster,
                                            TaggerInfo& ti);

        // nue_tagger: fills TaggerInfo NuE BDT features and returns flag_nue.
        // Prototype: WCPPID::NeutrinoID::nue_tagger() in NeutrinoID_nue_tagger.h.
        bool nue_tagger(Graph& graph,
                        Facade::Cluster* main_cluster,
                        VertexPtr main_vertex,
                        int apa, int face,
                        IndexedShowerSet& showers,
                        VertexShowerSetMap& map_vertex_to_shower,
                        IndexedShowerSet& pi0_showers,
                        ShowerIntMap& map_shower_pio_id,
                        std::map<int, std::vector<ShowerPtr>>& map_pio_id_showers,
                        std::map<int, std::pair<double,int>>& map_pio_id_mass,
                        IDetectorVolumes::pointer dv,
                        ParticleDataSet::pointer particle_data,
                        double muon_length,
                        TaggerInfo& ti);

        // singlephoton_tagger: fills TaggerInfo shw_sp_* BDT features and returns flag_sp.
        // Prototype: WCPPID::NeutrinoID::singlephoton_tagger() in NeutrinoID_singlephoton_tagger.h.
        // apa/face are derived internally from dv->contained_by(main_vertex_pt) so callers
        // do not need to know detector geometry details.
        bool singlephoton_tagger(Graph& graph,
                                 Facade::Cluster* main_cluster,
                                 VertexPtr main_vertex,
                                 IndexedShowerSet& showers,
                                 VertexShowerSetMap& map_vertex_to_shower,
                                 ShowerIntMap& map_shower_pio_id,
                                 std::map<int, std::vector<ShowerPtr>>& map_pio_id_showers,
                                 std::map<int, std::pair<double,int>>& map_pio_id_mass,
                                 IDetectorVolumes::pointer dv,
                                 TaggerInfo& ti);

        // ssm_tagger: fills TaggerInfo ssm_* and ssmsp_* features and returns flag_ssm.
        // Prototype: WCPPID::NeutrinoID::ssm_tagger() in NeutrinoID_ssm_tagger.h.
        bool ssm_tagger(Graph& graph,
                        VertexPtr main_vertex,
                        IndexedShowerSet& showers,
                        ShowerVertexMap& map_vertex_in_shower,
                        ShowerSegmentMap& map_segment_in_shower,
                        const Pi0KineFeatures& pio_kine,
                        int flag_ssmsp,
                        int& acc_segment_id,
                        const ParticleDataSet::pointer& particle_data,
                        const IRecombinationModel::pointer& recomb_model,
                        TaggerInfo& ti);

        KineInfo fill_kine_tree(VertexPtr main_vertex, IndexedShowerSet& showers,
                                const Pi0KineFeatures& pio_kine,
                                Graph& graph, TrackFitting& track_fitter,
                                IDetectorVolumes::pointer dv,
                                WireCell::IClusGeomHelper::pointer geom_helper,
                                const Clus::ParticleDataSet::pointer& particle_data,
                                const IRecombinationModel::pointer& recomb_model);
        // Convenience overloads: collect 2D charge maps internally (safe for isolated calls).
        double cal_kine_charge(ShowerPtr Shower, Graph& graph, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        double cal_kine_charge(SegmentPtr segment, Graph& graph, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        // Fast overload: reuse pre-collected 2D charge maps (avoids O(N_hits) collection per call).
        // Collect maps once with track_fitter.collect_2D_charge() and pass here when calling in a loop.
        double cal_kine_charge(ShowerPtr shower,
            const ChargeMap& charge_2d_u, const ChargeMap& charge_2d_v, const ChargeMap& charge_2d_w,
            const WireMap& map_apa_ch_plane_wires,
            TrackFitting& track_fitter, IDetectorVolumes::pointer dv);

    };
}
