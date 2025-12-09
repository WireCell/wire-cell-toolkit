#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/ClusteringFuncsMixins.h"
#include "WireCellClus/ParticleDataSet.h"
#include "WireCellClus/FiducialUtils.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Logging.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellClus/TrackFitting.h"  
#include "WireCellClus/TrackFittingPresets.h"
#include "WireCellClus/PRSegmentFunctions.h"

#include "WireCellIface/IScalarFunction.h"
#include "WireCellUtil/KSTest.h"

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

class TaggerCheckNeutrino : public IConfigurable, public Clus::IEnsembleVisitor, private Clus::NeedDV, private Clus::NeedPCTS, private Clus::NeedRecombModel, private Clus::NeedParticleData {
public:
    TaggerCheckNeutrino() {
        // Initialize with default preset
        m_track_fitter = TrackFittingPresets::create_with_current_values();
    }
    virtual ~TaggerCheckNeutrino() {}
    virtual void configure(const WireCell::Configuration& config) ;
    
    virtual Configuration default_configuration() const ;

    virtual void visit(Ensemble& ensemble) const;
    

    private:
        std::string m_grouping_name{"live"};
        std::string m_trackfitting_config_file;  // Path to TrackFitting config file
        mutable TrackFitting m_track_fitter; 
        
        void load_trackfitting_config(const std::string& config_file);

};