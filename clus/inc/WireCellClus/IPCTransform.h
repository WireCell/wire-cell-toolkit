/** A set of point cloud transforms for clustering.
    
    These interfaces are not at all general purpose so it is buried inside clus. 

 */
#ifndef WIRECELLCLUS_IPCTRANSFORM
#define WIRECELLCLUS_IPCTRANSFORM

#include "WireCellUtil/IComponent.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/PointCloudDataset.h"
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

    };

    class IPCTransformSet : public IComponent<IPCTransformSet> {
       public:
        virtual ~IPCTransformSet() {}
        virtual IPCTransform::pointer pc_transform(const std::string &name) const = 0;
    };
}

#endif
