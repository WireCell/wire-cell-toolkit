/** A set of point cloud transforms for clustering.
    
    These interfaces are not at all general purpose so it is buried inside clus. 

 */
#ifndef WIRECELLCLUS_IPCTRANSFORM
#define WIRECELLCLUS_IPCTRANSFORM

#include "WireCellUtil/IComponent.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/PointCloudDataset.h"
#include "WireCellUtil/PointTree.h"
#include <vector>

namespace WireCell::Clus {

    class IPCTransform : public IComponent<IPCTransform> {
    public:
        using Dataset = WireCell::PointCloud::Dataset;
        using Array = WireCell::PointCloud::Array;

        virtual ~IPCTransform() {}

        virtual Point forward(const Point& pos_raw, double clustser_t0, int face, int apa) const = 0;
        virtual Point backward(const Point& pos_cor, double clustser_t0, int face, int apa) const = 0;
        virtual bool filter(const Point& pos_cor, double clustser_t0, int face, int apa) const = 0;

        virtual Dataset forward(const Dataset& pc_raw, const std::vector<std::string>& arr_raw_names, const std::vector<std::string>& arr_cor_names, double clustser_t0, int face, int apa) const = 0;
        virtual Dataset backward(const Dataset& pc_cor, const std::vector<std::string>& arr_cor_names, const std::vector<std::string>& arr_raw_names, double clustser_t0, int face, int apa) const = 0;
        virtual Dataset filter(const Dataset& pc_cor, const std::vector<std::string>& arr_cor_names, double clustser_t0, int face, int apa) const = 0;

        /// Return the point cloud scope that this transform produces.
        /// This scope is registered as Cluster::m_scopes[correction_name] by
        /// Cluster::add_corrected_points().  Downstream components (e.g.
        /// ClusteringRetile) must be configured with a matching scope.
        virtual PointCloud::Tree::Scope output_scope() const = 0;

        /// Return the subset of corrected array names that must be persisted in
        /// each blob's local "3d" point cloud.  Only coordinates that actually
        /// change need to be stored; unchanged coordinates already exist in the
        /// original blob PC under their raw names.
        ///
        /// Example — T0Correction shifts x only, so {"x_t0cor"} is sufficient.
        /// A full 3D SCE correction would return {"x_sce", "y_sce", "z_sce"}.
        virtual std::vector<std::string> stored_array_names() const = 0;
    };

    class IPCTransformSet : public IComponent<IPCTransformSet> {
       public:
        virtual ~IPCTransformSet() {}
        virtual IPCTransform::pointer pc_transform(const std::string &name) const = 0;
    };
}

#endif
