#ifndef WIRECELL_CLUS_PRSEGMENTFUNCTIONS
#define WIRECELL_CLUS_PRSEGMENTFUNCTIONS

#include "WireCellClus/PRGraph.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/D4Vector.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellIface/IRecombinationModel.h"
#include "WireCellClus/ParticleDataSet.h"

namespace WireCell::Clus::PR {

    using geo_point_t = WireCell::Point;

    /// Replace the segment in the graph with two new segments that meet at a
    /// new vertex nearest to the point.
    ///
    /// The input segment is removed from the graph.
    ///
    /// The point must be withing max_dist of the segment.
    ///
    /// Returns true if the graph was modified.
    std::tuple<bool, std::pair<SegmentPtr, SegmentPtr>, VertexPtr> break_segment(Graph& graph, SegmentPtr seg, Point point, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, const IDetectorVolumes::pointer& dv,
                       double max_dist=1e9*units::cm);
    // patter recognition
    // cathode_kink_xcut > 0 makes the kink search skip candidate fit points within
    // that distance of the cathode plane at cathode_x, where the transverse
    // cathode mismatch fakes a turn and the para_angle guard (which exists to
    // protect isochronous imaging artifacts) is wide open because the crossing is
    // drift-x dominated -- see sbnd_xin/docs/pr/20 Part II.  Default 0 = no point
    // is ever skipped = byte-identical to the pre-pr/20 behavior.
    std::tuple<WireCell::Point, WireCell::Vector, WireCell::Vector, bool> segment_search_kink(SegmentPtr seg, WireCell::Point& start_p, const std::string& cloud_name = "fit", double dQ_dx_threshold = 43000/units::cm, double cathode_x = 0, double cathode_kink_xcut = 0 );

    /// Calculate track length from segment
    ///
    /// If flag == 1 and segment has fitted dx values, sum the dx values from fits.
    /// If flag == 0, calculate geometric length from wcpts.
    ///
    /// @param seg The segment to calculate length for
    /// @param flag Calculation method: 0=geometric from points, 1=from fit dx values
    /// @return Track length
    double segment_track_length(SegmentPtr seg, int flag = 0, int n1 = -1, int n2 = -1, WireCell::Vector dir_perp = WireCell::Vector(0,0,0));
    double segment_track_direct_length(SegmentPtr seg, int n1 = -1, int n2 = -1, WireCell::Vector dir_perp = WireCell::Vector(0,0,0));
    double segment_track_max_deviation(SegmentPtr seg, int n1 = -1, int n2 = -1);
    /// Calculate track length above dQ/dx threshold
    ///
    /// Extracts dQ and dx from segment's fits and calculates length above threshold.
    ///
    /// @param seg The segment containing fit data
    /// @param threshold dQ/dx threshold value
    /// @return Length of track segments above threshold
    double segment_track_length_threshold(SegmentPtr seg, double threshold = 75000./units::cm);
    /// Calculate track length from segment using geometric distance between points
    ///
    /// This is a convenience function that always uses geometric calculation
    /// regardless of available dx data.
    ///
    /// @param seg The segment to calculate length for
    /// @return Geometric track length
    double segment_geometric_length(SegmentPtr seg, int n1 = -1, int n2 = -1, WireCell::Vector dir_perp = WireCell::Vector(0,0,0));


    /// Calculate median dQ/dx for a segment
    ///
    /// Extracts dQ and dx from segment's fits and calculates median dQ/dx.
    ///
    /// @param seg The segment containing fit data
    /// @return Median dQ/dx value (0 if no valid fits)
    double segment_median_dQ_dx(SegmentPtr seg, int n1 = -1, int n2 = -1);
    double segment_rms_dQ_dx(SegmentPtr seg);

    /// doc sbnd_xin/docs/pr/40 -- shared evidence check for the F2/F3 knobs
    /// (shower_reclass_dqdx_guard, shower_topo_dqdx_guard).  Reuses the SAME
    /// median-dQ/dx thresholds segment_determine_dir_track's short-track
    /// fallback already trusts (PRSegmentFunctions.cxx, medium_dQ_dx >
    /// MIP_dQdx*1.75 => proton, < MIP_dQdx*1.2 => muon) to decide whether a
    /// segment's OWN charge profile is decisively non-electron-like, i.e. it
    /// should be spared from an unconditional track-to-electron
    /// reclassification.  A zero/absent median (no valid dQ/dx samples) never
    /// spares -- that is "no evidence", not "MIP-like evidence".
    bool segment_dqdx_spares_electron_reclass(SegmentPtr seg, double MIP_dQdx);

    /// doc sbnd_xin/docs/pr/40 round 2 F5 -- an electron cannot father a proton.
    ///
    /// True iff `seg` emanates from `main_vertex` (graph identity at either
    /// endpoint, not a distance cut -- every measured case sits at exactly
    /// d=0.00 cm) and its FAR endpoint's out-edges include a segment that is
    /// (a) already PID'd proton (pdg 2212, i.e. has_particle_info() and
    /// pdg()==2212) and (b) independently charge-confirmed by its own
    /// median dQ/dx (> 1.75x MIP_dQdx, the same threshold
    /// segment_dqdx_spares_electron_reclass and segment_determine_dir_track's
    /// own short-track fallback already trust).  `main_vertex` may be null
    /// (e.g. before stage 4 determines one), in which case this always
    /// returns false -- there is no vertex to require the segment emanate
    /// from.  Callers own the `flag_shower` check; this only judges the
    /// proton-daughter topology.
    bool segment_has_proton_daughter(Graph& graph, SegmentPtr seg, VertexPtr main_vertex, double MIP_dQdx);


    /// Create and associate a DynamicPointCloud with a segment from path points
    ///
    /// @param segment The segment to associate the DynamicPointCloud with
    /// @param path_points Vector of 3D points to process
    /// @param dv Detector volume for wire plane ID determination
    /// @param cloud_name Name for the DynamicPointCloud (default: "main")
    void create_segment_point_cloud(SegmentPtr segment,
                                    const std::vector<geo_point_t>& path_points,
                                    const IDetectorVolumes::pointer& dv,
                                    const std::string& cloud_name = "main",
                                    const std::vector<size_t>& global_indices = {});

    void create_segment_fit_point_cloud(SegmentPtr segment,
                                    const IDetectorVolumes::pointer& dv,
                                    const std::string& cloud_name = "fit");

    std::pair<double, WireCell::Point> segment_get_closest_point(SegmentPtr seg, const WireCell::Point& point, const std::string& cloud_name = "fit", const std::string& base_cloud_name = "main");
    std::tuple<double, double, double> segment_get_closest_2d_distances(SegmentPtr seg, const WireCell::Point& point, int apa, int face, const std::string& cloud_name = "fit");
    double segment_get_closest_2d_distance(SegmentPtr seg, const WireCell::Point& point, int apa, int face, int plane, const std::string& cloud_name = "fit");


    // PID related
    bool eval_ks_ratio(double ks1, double ks2, double ratio1, double ratio2);
    // skip_stop_samples: exclude that many samples nearest the hypothesized
    // stopping end (largest L, smallest residual range) from the comparison,
    // WITHOUT moving the template anchor (end_L).  0 = legacy, byte-identical.
    // Used by the endpoint-trim retry (doc sbnd_xin/docs/pr/9 sec. 6 F1): the
    // stopping tip's dQ/dx is unreliable when the endpoint is ill-defined
    // (diluted OR piled-up), and it is compared against the template's Bragg
    // maximum, so one bad tip sample can veto an otherwise clean decision.
    // empty_abstain (doc sbnd_xin/docs/pr/31 §12, F6 was P7): when the
    // comparison window holds no samples (or a dEdx template is missing), the
    // degenerate return's element 0 becomes 0.0 ("abstain") instead of 1.0
    // ("this orientation passed the direction gate").  0.0 is the prototype's
    // degenerate answer, verified by execution: zero-bin TH1F KolmogorovTest
    // returns 0 for both templates, so eval_ks_ratio's first gate
    // (ks1-ks2 >= 0.0) fails and results[0] is false.  false = legacy 1.0 =
    // byte-identical.
    std::vector<double> do_track_comp(std::vector<double>& L , std::vector<double>& dQ_dx, double compare_range, double offset_length, const Clus::ParticleDataSet::pointer& particle_data, double MIP_dQdx = 50000/units::cm, int skip_stop_samples = 0, bool empty_abstain = false);

    // Options for segment_do_track_pid / segment_determine_dir_track.
    // Defaults reproduce legacy behavior byte-for-byte (uBooNE values, vote off).
    // - mip_dqdx: the flat-template MIP amplitude handed to do_track_comp
    //   (uBooNE legacy 50e3 e/cm; distinct from the 43e3 median-threshold
    //   scale carried by the functions' MIP_dQdx parameter).
    // - proton_dir_vote: doc sbnd_xin/docs/pr/8 improvement (beyond prototype):
    //   when the muon-vs-flat direction gate abstains in both orientations,
    //   let the proton template declare the direction, guarded by
    //   score/asymmetry/free-end conditions.  Initial thresholds are
    //   placeholders pending the pr/8 sec. 6 calibration.
    // - start_n/end_n: vertex connectivity, filled by
    //   segment_determine_dir_track (guard G4: stop end must be free).
    // - endpoint_trim_retry: doc sbnd_xin/docs/pr/9 sec. 6 F1 (beyond
    //   prototype, owner-approved 2026-07-30): when the legacy decision
    //   abstains in both orientations, retry ONCE with exactly one sample
    //   excluded at each orientation's hypothesized stopping end (template
    //   anchor unchanged).  Dynamic by construction: a well-found endpoint
    //   means the first attempt already decides and nothing is trimmed;
    //   never trims more than 1 sample; value-agnostic (robust to diluted
    //   AND piled-up tips).  Runs BEFORE the proton_dir_vote fallback.
    // - track_comp_empty_abstain: doc sbnd_xin/docs/pr/31 §12 (F6, was P7).
    //   Forwarded to do_track_comp's empty_abstain (see its comment above).
    //   false = legacy "confirmed" filler = byte-identical.
    // - dir_track_median_local: doc sbnd_xin/docs/pr/31 §12 (F4, was P8).
    //   segment_determine_dir_track takes its median dQ/dx over the SAME local
    //   dQ_dx vector it hands to the PID (prototype ProtoSegment.cxx:1574-1576,
    //   :1602-1611: nth_element over a copy of dQ_dx), instead of
    //   segment_median_dQ_dx's filtered rebuild from fits() -- which drops
    //   invalid/dx<=0/dQ<0 samples the PID vector keeps as zeros, so the two
    //   disagree about the same pathological point in opposite directions and
    //   nth_element selects a different order statistic.  Same internal-unit
    //   scale either way (charge per internal length; see the unit-convention
    //   comment in segment_determine_dir_track).  false = filtered helper =
    //   byte-identical.
    // - track_pid_persist_dqdx: doc sbnd_xin/docs/pr/40 (F1, = doc pr/7 sec 5 /
    //   pr/31 P14/F8).  segment_determine_dir_track's final store
    //   (`if (pdg_code != 0 && ((dirsign==1 && end_n==1) || (dirsign==-1 &&
    //   start_n==1)))`) gates BOTH the type/mass persistence and the 4-mom
    //   recompute on the direction pointing at a free end -- so a dQ/dx-only
    //   recovery (medium_dQ_dx > 1.75x/< 1.2x MIP with dirsign left at 0, or a
    //   confident template PID whose stop end is not topologically free) is
    //   computed and then silently discarded; the segment exits with no
    //   particle_info at all.  The prototype (ProtoSegment.cxx:1637-1639)
    //   persists type+mass unconditionally and gates ONLY cal_4mom() on
    //   get_particle_4mom(3)>0 ("an energy was already computed").  true =
    //   restore that shape: store type+mass whenever pdg_code != 0, gate only
    //   the 4-momentum recompute on the existing free-end test.  false =
    //   legacy = byte-identical.  Measured (doc pr/40 Part 0): of 9 owner-
    //   reported track-to-electron cases, this persistence gate is why 4 of
    //   them (evt 74544, 174637, 267597, 269774) had NO particle_info at all
    //   by the time a wholesale shower-reclassification site looked at them.
    // - track_pid_persist_4mom: doc sbnd_xin/docs/pr/40 round 2 (F4).  When
    //   track_pid_persist_dqdx rescues a non-free-end store, the legacy
    //   4-momentum stub is D4Vector(mass,0,0,0) -- E = mass exactly, so
    //   Aux::ParticleInfo::kinetic_energy() (= E - mass) is ZERO for that
    //   segment, byte-for-byte, forever (SBND evt 174637 seg 9050 measured
    //   as "mu- 0 MeV" in the Bee PF tree once track_pid_persist_dqdx went
    //   SBND-default-ON).  segment_cal_4mom has NO dependence on the
    //   direction being a free end -- its only direction coupling is
    //   segment_cal_dir_3vector, which already returns a zero 3-vector when
    //   dirsign()==0 -- so the free-end gate on it was purely external and
    //   unnecessarily strict.  true = call segment_cal_4mom unconditionally
    //   whenever pdg_code != 0 (correct KE, momentum direction only as good
    //   as dirsign allows).  false = legacy rest-mass-only stub =
    //   byte-identical.  Independent of track_pid_persist_dqdx.
    struct TrackPidOptions {
        double mip_dqdx{50000/units::cm};
        bool   proton_dir_vote{false};
        double proton_dir_score_max{0.25};
        double proton_dir_asym_min{1.3};
        bool   endpoint_trim_retry{false};
        int    start_n{1};
        int    end_n{1};
        bool   track_comp_empty_abstain{false};
        bool   dir_track_median_local{false};
        bool   track_pid_persist_dqdx{false};
        bool   track_pid_persist_4mom{false};
    };

    // success, flag_dir, pdg_code, particle_score
    std::tuple<bool, int, int, double> segment_do_track_pid(SegmentPtr segment, std::vector<double>& L , std::vector<double>& dQ_dx, const Clus::ParticleDataSet::pointer& particle_data, double compare_range=35*units::cm, double offset_length = 0*units::cm, bool flag_force = false,  const TrackPidOptions& pid_opts = {});

    // direction calculation ...
    /// Returns true if the segment's direction is weakly determined.
    ///
    /// Mirrors prototype ProtoSegment::is_dir_weak(): checks particle-type-based
    /// score thresholds (muon, proton) and the static dir_weak flag.
    bool segment_is_dir_weak(SegmentPtr seg);
    WireCell::Vector segment_cal_dir_3vector(SegmentPtr seg);
    WireCell::Vector segment_cal_dir_3vector(SegmentPtr seg, WireCell::Point& p, double dis_cut);
    WireCell::Vector segment_cal_dir_3vector(SegmentPtr seg, int direction, int num_points, int start);
    void segment_determine_dir_track(SegmentPtr segment, int start_n, int end_n, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, double MIP_dQdx = 43000/units::cm, bool flag_print = false, const TrackPidOptions& pid_opts = {});

    // kinemiatics calculations ...
    double segment_cal_kine_dQdx(SegmentPtr seg, const IRecombinationModel::pointer& recomb_model);
    double cal_kine_dQdx(std::vector<double>& vec_dQ, std::vector<double>& vec_dx, const IRecombinationModel::pointer& recomb_model);
    double cal_kine_range(double L, int pdg_code, const Clus::ParticleDataSet::pointer& particle_data);
    // 4-momentum: E, px, py, pz
    WireCell::D4Vector<double> segment_cal_4mom(SegmentPtr segment, int pdg_code, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, double MIP_dQdx = 50000/units::cm);

    // EMshower PID
    void clustering_points_segments(std::vector<SegmentPtr> segments, const IDetectorVolumes::pointer& dv, const std::string& cloud_name = "associate_points", double search_range = 1.2*units::cm, double scaling_2d = 0.7);

    /// doc sbnd_xin/docs/pr/32 §11 F2 -- knob transport for
    /// segment_is_shower_trajectory's flag semantics.
    ///
    /// The prototype's ProtoSegment::is_shower_trajectory() opens with
    /// `flag_shower_trajectory = false` (ProtoSegment.cxx:544) and sets it true
    /// only if the test fires (:608), so EVERY call re-caches the label and a
    /// segment that no longer qualifies is DEMOTED.  The toolkit's port only
    /// ever calls set_flags (PRSegmentFunctions.cxx), so kShowerTrajectory is
    /// monotone: once set it survives every later negative test.  That is why
    /// the two improve_vertex call sites recompute instead of reading the flag
    /// (doc pr/32 §10.3) -- reading it would return an OR over history.
    ///
    /// true => mirror the prototype: clear the bit on entry, set it on a
    /// positive result.  false => today's set-only behaviour, byte-identical.
    ///
    /// This is a free function reached from three files with no component
    /// configuration in scope, so the value travels through a process-wide
    /// flag written once per event by TaggerCheckNeutrino before any graph is
    /// built -- the same transport as PR::g_graph_endpoint_policy (doc pr/30
    /// §11 P8).  Read-mostly; never written concurrently with clustering.
    extern bool g_shower_traj_refresh_flag;

    bool segment_is_shower_trajectory(SegmentPtr seg, double step_size = 10*units::cm, double mip_dQ_dx = 50000 / units::cm);
    void segment_determine_shower_direction_trajectory(SegmentPtr segment, int start_n, int end_n, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, double MIP_dQdx = 43000/units::cm, bool flag_print = false, const TrackPidOptions& pid_opts = {});
    
    // median_local (doc sbnd_xin/docs/pr/31 §12, F4): forwarded to the interior
    // segment_determine_dir_track call on the short-segment branch as
    // TrackPidOptions::dir_track_median_local ONLY -- that interior call
    // otherwise passes a default-constructed TrackPidOptions, and forwarding a
    // caller's full option set there would unconditionally change proton_dir_vote
    // et al. at that site.  false = legacy = byte-identical.
    bool segment_determine_shower_direction(SegmentPtr segment, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, const std::string& cloud_name = "associate_points", double MIP_dQdx = 43000/units::cm, double rms_cut= 0.4*units::cm, double mip_dqdx = 50000/units::cm, bool median_local = false);
    // sbnd_xin doc pr/25 sec 3: `demote_len` > 0 makes the existing long-track
    // guard unconditional -- any segment this function would flag
    // kShowerTopology whose geometric length (segment_track_length(seg,0),
    // the same measure the legacy >50cm guard uses) exceeds `demote_len` is
    // demoted to a track instead.  Default 0 == off == byte-identical.
    // Motivation: an owner hand-scan of every long shower-topology segment on
    // a selected nu-candidate main cluster in the 572-event valfast manifest
    // (2026-08-03, 10/10 events) found NO showers -- all tracks.  See sec 3.8.
    // reset (doc sbnd_xin/docs/pr/31 §12, F3 was P13): the prototype's
    // is_shower_topology opens with two assignments BEFORE its early returns
    // (ProtoSegment.cxx:319-321): flag_shower_topology = tmp_val (all callers
    // pass false => the flag is CLEARED on every entry) and flag_dir = 0.  The
    // toolkit's port assigns a local and never clears the segment flag, and its
    // four early returns skip dirsign() -- so a segment that qualified on an
    // earlier pass keeps a stale kShowerTopology flag and a stale direction.
    // true => clear the flag bit (unset_flags, not clear_flags -- other bits
    // survive) and zero dirsign at entry, mirroring the prototype; the tail
    // re-sets both when the answer is still yes.  false = set-only legacy =
    // byte-identical.
    // dqdx_guard (doc sbnd_xin/docs/pr/40, F3): the local vec_dQ_dx this
    // function already builds (fits[i].dQ/dx, normalized by MIP_dQ_dx) is
    // otherwise DEAD -- read only as .size() (pr/31 GOTCHA 5) -- so the
    // 5-branch geometric spread test that decides flag_shower_topology never
    // once consults the segment's own charge.  true = after that test (and
    // the existing length-based demotions) decide flag_shower_topology=true,
    // override it back to false when segment_dqdx_spares_electron_reclass
    // says the segment's median dQ/dx is decisively proton- or muon-like.
    // This does NOT touch flag_dir or the geometric test itself -- a
    // genuinely isochronous/EM-like spread pattern with ambiguous charge
    // still sets the flag.  false = legacy = byte-identical.  Measured (doc
    // pr/40 Part 0): rescues both owner-reported cases where the topology
    // test itself, not a pdg-11 write elsewhere, was the mechanism (evt
    // 256587 seg 11079 at 1.26x MIP, evt 489330 seg 4018 at 2.73x MIP -- both
    // shorter than the existing shower_topo_demote_len's reach).
    bool segment_is_shower_topology(SegmentPtr seg, bool tmp_val=false, double MIP_dQ_dx = 43000/units::cm,
                                    double demote_len = 0, bool reset = false, bool dqdx_guard = false);
}

#endif
