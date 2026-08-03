#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/ClusteringFuncsMixins.h"
#include "WireCellClus/ParticleDataSet.h"
#include "WireCellClus/FiducialUtils.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellAux/Logger.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellClus/TrackFitting.h"  
#include "WireCellClus/TrackFittingPresets.h"
#include "WireCellClus/PRSegmentFunctions.h"

#include "WireCellIface/IScalarFunction.h"
#include "WireCellUtil/KSTest.h"

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

class TaggerCheckNeutrino : public Aux::Logger, public IConfigurable, public Clus::IEnsembleVisitor, private Clus::NeedDV, private Clus::NeedPCTS, private Clus::NeedRecombModel, private Clus::NeedParticleData, private Clus::NeedClusGeomHelper {
public:
    TaggerCheckNeutrino() : Aux::Logger("TaggerCheckNeutrino", "clus") {
        // Initialize with default preset
        m_track_fitter = std::make_shared<TrackFitting>(TrackFittingPresets::create_with_current_values());
    }
    virtual ~TaggerCheckNeutrino() {}
    virtual void configure(const WireCell::Configuration& config) ;

    virtual Configuration default_configuration() const ;

    virtual void visit(Ensemble& ensemble) const;


    private:
        std::string m_grouping_name{"live"};
        std::string m_trackfitting_config_file;  // Path to TrackFitting config file
        bool m_perf{false};  // if true, print per-step timing to stdout
        // Route direction-weakness reads through segment_is_dir_weak() (the
        // faithful port of prototype ProtoSegment::is_dir_weak(), score
        // thresholds included).  C++ default false = legacy raw-flag reads,
        // byte-identical.  See PatternAlgorithms::m_dir_weak_use_score.
        bool m_dir_weak_use_score{false};
        // MIP dQ/dx scales in e/cm (uBooNE legacy defaults => absent keys are
        // byte-identical).  Two roles: mip_dqdx = flat-template amplitude
        // (legacy 50e3); mip_dqdx_median = median-ratio threshold scale
        // (legacy 43e3).  See PatternAlgorithms::m_mip_dqdx{,_median}.
        double m_mip_dqdx{50000.0};
        double m_mip_dqdx_median{43000.0};
        // Proton-template direction vote (doc sbnd_xin/docs/pr/8; default
        // false = legacy abstention).  Thresholds pending pr/8 calibration.
        bool   m_proton_dir_vote{false};
        double m_proton_dir_score_max{0.25};
        double m_proton_dir_asym_min{1.3};
        // Endpoint-trim retry (doc sbnd_xin/docs/pr/9 sec. 6 F1; default false
        // = legacy).  On abstention only, retry the dQ/dx PID once with 1
        // sample excluded at the hypothesized stopping end (ill-defined
        // endpoints dilute or inflate the tip dQ/dx).  Runs before the vote.
        bool   m_endpoint_trim_retry{false};
        // Minimum wcpt-path length (cm) for a segment to count as a leg in the
        // fit_vertex position fit (doc sbnd_xin/docs/pr/9 sec. 11 F3c; default
        // 0 = legacy include-all, byte-identical).  Stops vertex-activity
        // stubs from dragging the vertex: >=3 surviving legs => fit on the
        // survivors; <=2 => skip the fit (two-leg position already fit).
        double m_fit_vertex_min_seg_length{0};
        // Cathode kink veto (doc sbnd_xin/docs/pr/20 Part II B0), both in cm.
        // cathode_kink_xcut = 0 => segment_search_kink sees every fit point =>
        // byte-identical to the pre-pr/20 behavior.  cathode_x is the cathode
        // plane in the T0-corrected frame (SBND 0; the same convention
        // ClusteringCathodeConnect's cathode_x uses).
        double m_cathode_x{0};
        double m_cathode_kink_xcut{0};
        // Long shower-topology demote length, cm (doc sbnd_xin/docs/pr/25
        // sec 3).  0 => the guard never fires => byte-identical.  50 is the
        // scan-supported operating point (9/10 owner-scanned events; ~45
        // covers all 10).  Threaded to PatternAlgorithms and thence to every
        // segment_is_shower_topology call site -- both NeutrinoTrackShowerSep
        // and NeutrinoVertexFinder, so one segment cannot be classified two
        // ways within one event.
        double m_shower_topo_demote_len{0};
        // Isochronous first-segment endpoint finding (doc sbnd_xin/docs/pr/24
        // round 2, SBND evt 271851); lengths in cm.  false => legacy
        // wire-footprint boundary endpoints, byte-identical.  When on, a
        // long cluster whose quantile-trimmed drift-x extent is small (a
        // filled 2-D sheet) gets its first-segment endpoints from the sheet's
        // principal axis instead, and the local-PCA endpoint refinement is
        // skipped on that branch.
        bool   m_iso_endpoint{false};
        double m_iso_endpoint_min_length{40};   // cm
        double m_iso_endpoint_max_xext{25};     // cm
        double m_iso_endpoint_xext_frac{0.35};
        double m_iso_endpoint_xext_quantile{0.02};
        // Detector-extent literals (docs/pr/2 sec. 2e(iv)), all in cm.
        // Defaults = the uBooNE prototype values (active volume y in
        // [-116,+117], z in [0,1037]) => absent keys are byte-identical.
        // cosmic_tagger()'s four "reaches the top" cuts, quoted as a distance
        // below the top face (117 - value): 17 / 15 / 37 / 67 cm.
        double m_cosmic_y_top_main{100};
        double m_cosmic_y_top_strict{102};
        double m_cosmic_y_top_loose{80};
        double m_cosmic_y_small_piece{50};
        // Denominator of the upstream-z vertex-score penalty (z-min_z)/scale in
        // compare_main_vertices{,_global}.  See PatternAlgorithms.
        double m_vertex_z_prior_scale{200};
        // SSM beam-line reference directions in the detector frame, {x,y,z}.
        // Defaults = the prototype's uBooNE BNB-target / NuMI-absorber vectors
        // (absent keys => byte-identical).  SBND has no value for either yet
        // (docs/pr/2 sec. 2e(i)).  See PatternAlgorithms::m_ssm_target_dir.
        std::vector<double> m_ssm_target_dir{0.46, 0.05, 0.885};
        std::vector<double> m_ssm_absorber_dir{0.33, 0.75, -0.59};
        // Charge -> kinetic-energy calibration constants of NeutrinoEnergyReco
        // (docs/pr/2 sec. 2e(iii)).  Defaults = the uBooNE literals they
        // replaced => absent keys are byte-identical.  Dimensionless except
        // m_kine_w_value (eV).  See PR::KineChargeOptions.
        double m_kine_fudge_factor{0.95};
        double m_kine_recom_factor{0.7};
        double m_kine_shower_fudge_factor{0.8};
        double m_kine_shower_recom_factor{0.5};
        double m_kine_proton_recom_factor{0.35};
        std::vector<double> m_kine_plane_weights{0.25, 0.25, 1.0};  // {U,V,W}
        double m_kine_plane_asym_switch{0.04};
        double m_kine_w_value{23.6};  // eV per electron-ion pair
        // Muon median-dQ/dx-vs-length envelope {c0, c1, pivot_cm, power}:
        // cut = c0 + c1*pow(pivot/L, power), a multiple of mip_dqdx_median,
        // used at nine tagger sites (docs/pr/2 sec. 2e(iv)).  Defaults = the
        // prototype's empirical uBooNE stopping-muon refit => absent key is
        // byte-identical.  See PatternAlgorithms::m_muon_dqdx_curve.
        std::vector<double> m_muon_dqdx_curve{0.8866, 0.9533, 18.0, 0.4234};
        // Single-photon stem dE/dx conversion (docs/pr/2 sec. 2e(i) third
        // correctness item).  Default false = the inline uBooNE-field inverse
        // Modified-Box, byte-identical; true routes shw_sp_vec_*_dedx through
        // the configured recombination model.  The mean-dedx threshold
        // (MeV/cm) is coupled to that choice; default = the legacy 2.3.
        bool   m_sp_dedx_use_recomb_model{false};
        double m_sp_mean_dedx_cut{2.3};
        std::string m_dl_weights;              // path to SCN vertex .pth weights file (empty = DL disabled)
        double m_dl_vtx_cut{25.0};             // max distance (mm) from DL prediction to accept candidate vertex (default 2.5 cm)
        double m_dQdx_scale{0.1};              // scale factor applied to dQ before passing to SCN network
        double m_dQdx_offset{-1000.0};         // offset applied after scaling: q_in = dQ * scale + offset
        bool   m_dl_vtx_rerank{true};              // if true, use top-K + soft re-rank; if false, use legacy single-best argmax
        int    m_dl_vtx_top_k{5};                // number of top voxels from DL inference to re-rank (only when rerank enabled)
        double m_dl_vtx_min_accept_score{4.0};    // minimum composite re-rank score to accept DL vertex (only when rerank enabled; matches the 4.0 advertised in default_configuration -- docs/pr/2 sec 8.4)
        double m_dl_vtx_score_scale{1000.0};      // scale factor applied to the raw DL score term in the composite re-rank score
        double m_beam_window_low{0};   // beam window [low, high) on cluster_t0 (matched flash time, WCT units).
        double m_beam_window_high{0};  // low >= high (default) disables the gate: uBooNE single-main behavior.
        bool m_nu_skip_cosmic{false};  // if true (beam-gate only), skip in-window mains already tagged
                                       // cosmic upstream: flag_TGM, flag_STM, or lm_flag > 0.
        bool m_nu_skip_cosmic_bundle{false};  // if true, that verdict vetoes the whole flash bundle
                                              // (every in-window main sharing matched_flash_gid with a
                                              // cosmic-tagged main), not just the tagged main itself.
                                              // Inert unless m_nu_skip_cosmic.
        double m_nu_skip_cosmic_bundle_min_length{0};  // cm.  > 0 spares an UNTAGGED in-window main at
                                                       // least this long from the bundle veto (its
                                                       // cosmic-tagged bundle-mate stays excluded --
                                                       // companions are associated-only).  0 = legacy:
                                                       // the veto removes every bundle-mate.
                                                       // Inert unless m_nu_skip_cosmic_bundle.
                                                       // docs/pr/16 design A (SBND 15).
        bool m_skip_cosmic_companions{false};  // doc pr/20 Part I P4.  If true, a COMPANION
                                               // (other_clusters member) that is TGM- or
                                               // STM-tagged and at least
                                               // m_cosmic_companion_min_length long is dropped
                                               // from the neutrino's companion list.  Nothing
                                               // tags a companion today unless the taggers run
                                               // with evaluate_demoted_mains (P3), so this is
                                               // inert without it.
        double m_cosmic_companion_min_length{0};  // cm.  A tagged companion SHORTER than this
                                                  // stays in regardless of verdict, so a
                                                  // mis-tagged short neutrino daughter can never
                                                  // be silently dropped.  0 = no floor: drop on
                                                  // the verdict alone.  Own tuning, deliberately
                                                  // NOT inheriting nu_skip_cosmic_bundle_min_length
                                                  // -- a different question.
        mutable std::shared_ptr<TrackFitting> m_track_fitter;

        void load_trackfitting_config(const std::string& config_file);

};