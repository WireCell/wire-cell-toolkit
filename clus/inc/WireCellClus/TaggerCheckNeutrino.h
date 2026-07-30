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
        // SSM beam-line reference directions in the detector frame, {x,y,z}.
        // Defaults = the prototype's uBooNE BNB-target / NuMI-absorber vectors
        // (absent keys => byte-identical).  SBND has no value for either yet
        // (docs/pr/2 sec. 2e(i)).  See PatternAlgorithms::m_ssm_target_dir.
        std::vector<double> m_ssm_target_dir{0.46, 0.05, 0.885};
        std::vector<double> m_ssm_absorber_dir{0.33, 0.75, -0.59};
        std::string m_dl_weights;              // path to SCN vertex .pth weights file (empty = DL disabled)
        double m_dl_vtx_cut{25.0};             // max distance (mm) from DL prediction to accept candidate vertex (default 2.5 cm)
        double m_dQdx_scale{0.1};              // scale factor applied to dQ before passing to SCN network
        double m_dQdx_offset{-1000.0};         // offset applied after scaling: q_in = dQ * scale + offset
        bool   m_dl_vtx_rerank{true};              // if true, use top-K + soft re-rank; if false, use legacy single-best argmax
        int    m_dl_vtx_top_k{5};                // number of top voxels from DL inference to re-rank (only when rerank enabled)
        double m_dl_vtx_min_accept_score{0.0};    // minimum composite re-rank score to accept DL vertex (only when rerank enabled)
        double m_dl_vtx_score_scale{1000.0};      // scale factor applied to the raw DL score term in the composite re-rank score
        double m_beam_window_low{0};   // beam window [low, high) on cluster_t0 (matched flash time, WCT units).
        double m_beam_window_high{0};  // low >= high (default) disables the gate: uBooNE single-main behavior.
        bool m_nu_skip_cosmic{false};  // if true (beam-gate only), skip in-window mains already tagged
                                       // cosmic upstream: flag_TGM, flag_STM, or lm_flag > 0.
        mutable std::shared_ptr<TrackFitting> m_track_fitter;

        void load_trackfitting_config(const std::string& config_file);

};