#include "WireCellUtil/PointCloudDataset.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/Units.h"
#include "WireCellIface/IBlob.h"

namespace WireCell::Aux {


    // Some crufty stuff used in PointTreeBuilding and UbooneClusterSource.
    PointCloud::Dataset
    make_scaler_dataset(const IBlob::pointer iblob, const Point& center,
                        const int npoints = 0, const double tick_span = 0.5*units::us);

    // Calculate the average position of a point cloud tree.
    WireCell::Point calc_blob_center(const PointCloud::Dataset& ds);

}
