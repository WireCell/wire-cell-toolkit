/// This API provides some mixin classes for "clustering classes" to handle
/// common behavior.

#ifndef WIRECELLCLUS_CLUSTERINGFUNCSMIXINS
#define WIRECELLCLUS_CLUSTERINGFUNCSMIXINS

#include "WireCellClus/IPCTransform.h"
#include "WireCellClus/ParticleDataSet.h"

#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellIface/IRecombinationModel.h"
#include "WireCellIface/IFiducial.h"
#include "WireCellUtil/Configuration.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"

#include "WireCellUtil/PointTree.h"

namespace WireCell::Clus {

    // A mixin to get an IDetectorVolumes
    // 
    // Configuration:
    //
    // "detector_volumes" : string type/name of an IDetectorVolumes
    struct NeedDV {
        IDetectorVolumes::pointer m_dv;
        void configure(const WireCell::Configuration &cfg);
    };

    // A mixin to get a IRecombinationModel
    struct NeedRecombModel {
        IRecombinationModel::pointer m_recomb_model;
        void configure(const WireCell::Configuration &cfg);
    };

    // A mixin to get particle_data_set
    class NeedParticleData {
        ParticleDataSet::pointer m_particle_data;
        std::string m_particle_data_name{"ParticleDataSet"};
        
    public:
        void configure(const WireCell::Configuration& config) {
            m_particle_data_name = get<std::string>(config, "particle_dataset", m_particle_data_name);
            auto configurable = Factory::find_tn<IConfigurable>(m_particle_data_name);
            m_particle_data = std::dynamic_pointer_cast<ParticleDataSet>(configurable);
            if (!m_particle_data) {
                THROW(ValueError() << errmsg{"Failed to find or cast ParticleDataSet: " + m_particle_data_name});
            }
        }
        
    protected:
        ParticleDataSet::pointer particle_data() { return m_particle_data; }
        const ParticleDataSet::pointer particle_data() const { return m_particle_data; }
    };

    // A mixin to get an IFiducial
    // 
    // Configuration:
    //
    // "fiducial" : string type/name of an IFiducial
    struct NeedFiducial {
        IFiducial::pointer m_fiducial;
        void configure(const WireCell::Configuration &cfg);
    };

    // A mixin to get an IPCTransformSet
    // 
    // Configuration:
    //
    // "pc_transforms" : string type/name of an IPCTransformSet
    struct NeedPCTS {
        Clus::IPCTransformSet::pointer m_pcts;
        void configure(const WireCell::Configuration &cfg);
    };

    // A mixin for things that need to be configured for a scope (pc name and coord names)
    // 
    // Configuration:
    //
    // "pc_name" : string name for a "node-local" PointCloud::Dataset
    // "coords" : array of string giving names of arrays in the PC to use as coordinates
    struct NeedScope {
        PointCloud::Tree::Scope m_scope;
        NeedScope(const std::string &pcname = "3d", 
                  const std::vector<std::string> &coords = {"x", "y", "z"},
                  size_t depth = 0);
        void configure(const WireCell::Configuration &cfg);
    };

}
#endif
