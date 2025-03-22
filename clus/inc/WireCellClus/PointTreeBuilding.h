/** Sample blobs to make point cloud tree and output as tensors.
 * Same as PointTreeBuilding but use ICluster as input.
*/
#ifndef WIRECELL__CLUS_POINTTREEBUILDING
#define WIRECELL__CLUS_POINTTREEBUILDING

#include "WireCellIface/IClusterFaninTensorSet.h"
#include "WireCellIface/IBlobSampler.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellAux/Logger.h"
#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/Units.h"
#include "WireCellClus/IClusGeomHelper.h"


namespace WireCell::Clus {

    class PointTreeBuilding : public Aux::Logger, public IClusterFaninTensorSet, public IConfigurable
    {
      public:
        PointTreeBuilding();
        virtual ~PointTreeBuilding();

        // INode, override because we get multiplicity at run time.
        virtual std::vector<std::string> input_types();

        // IConfigurable
        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        virtual bool operator()(const input_vector& invec, output_pointer& tensorset);

      private:
        // sampling for live/dead
        using node_ptr = WireCell::PointCloud::Tree::Points::node_ptr;
        node_ptr sample_live(const WireCell::ICluster::pointer cluster, const double tick, const double angle_u, const double angle_v, const double angle_w) const;
        node_ptr sample_dead(const WireCell::ICluster::pointer cluster, const double tick) const;
        // add CT point cloud to the root/Grouping
        void add_ctpc(node_ptr& root, const WireCell::ICluster::pointer cluster) const;
        // wind -> xbeg, xend
        void add_dead_winds(node_ptr& root, const WireCell::ICluster::pointer cluster) const;

        size_t m_multiplicity {2};
        std::vector<std::string> m_tags;
        size_t m_count{0};


        // double m_tick {0.5*units::us};
        // double m_drift_speed {1.101*units::millimeter/units::us};
        // double m_time_offset {-1600 * units::us};
        double m_dead_threshold {1e10};
        // double m_angle_u {1.0472}; // 60 degrees
        // double m_angle_v {-1.0472}; // -60 degrees
        // double m_angle_w {0}; // 0 degrees

        // the anode to be processed
        IAnodePlane::pointer m_anode;

        IDetectorVolumes::pointer m_dv;

        // the face to be processed
        int m_face{0};

        // the geometry helper
        // IClusGeomHelper::pointer m_geomhelper;
        
        /** Configuration: "samplers"

            An object with attributes providing names of
            IBlobSamplers.  The attribute names will be used to name
            the point cloud produced by the samplers.
        */
        std::map<std::string, IBlobSampler::pointer> m_samplers;

        /** Config: "datapath"

            Set the datapath for the tensor representing the point
            cloud tree.  If a %d format code is found it wil be
            interpolated with the IBlobSet::ident() value.
         */
        std::string m_datapath = "pointtrees/%d";

    };
}

#endif
