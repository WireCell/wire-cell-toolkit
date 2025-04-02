#ifndef WIRECELL_POINTCLOUDTRANSFORM
#define WIRECELL_POINTCLOUDTRANSFORM

#include "WireCellUtil/PointCloudDataset.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/Configuration.h"
#include <functional>

namespace WireCell::PointCloud {

    /**
     * Transform provides operations for transforming and filtering 
     * point cloud datasets.
     */
    class Transform {
      public:
        using pointer = std::shared_ptr<Transform>;

        // Default constructor
        Transform() = default;
        
        // Virtual destructor for proper inheritance
        virtual ~Transform() = default;
        
        virtual Point forward(const Point& pos, double clustser_t0, int face, int apa) const = 0;
        virtual Point backward(const Point& pos, double clustser_t0, int face, int apa) const = 0;
        virtual bool filter(const Point& pos, double clustser_t0, int face, int apa) const = 0;

        virtual Dataset forward(const Dataset& pc, const std::vector<std::string>& arr_names, double clustser_t0, int face, int apa) const = 0;
        virtual Dataset backward(const Dataset& pc, const std::vector<std::string>& arr_names, double clustser_t0, int face, int apa) const = 0;
        virtual Dataset filter(const Dataset& pc, const std::vector<std::string>& arr_names, double clustser_t0, int face, int apa) const = 0;
    };

}  // namespace WireCell::PointCloud

#endif // WIRECELL_POINTCLOUDTRANSFORM
