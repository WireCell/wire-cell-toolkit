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

class TaggerCheckNeutrino : public Aux::Logger, public IConfigurable, public Clus::IEnsembleVisitor, private Clus::NeedDV, private Clus::NeedPCTS, private Clus::NeedRecombModel, private Clus::NeedParticleData {
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
        std::string m_dl_weights;      // path to SCN vertex .pth weights file (empty = DL disabled)
        double m_dl_vtx_cut{20.0};     // max distance (mm) from DL prediction to accept candidate vertex (default 2 cm)
        double m_dQdx_scale{0.1};      // scale factor applied to dQ before passing to SCN network
        double m_dQdx_offset{-1000.0}; // offset applied after scaling: q_in = dQ * scale + offset
        mutable std::shared_ptr<TrackFitting> m_track_fitter;

        void load_trackfitting_config(const std::string& config_file);

};