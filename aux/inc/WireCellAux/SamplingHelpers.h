#include "WireCellUtil/PointCloudDataset.h"
#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/Units.h"
#include "WireCellIface/IBlob.h"
#include "WireCellIface/IBlobSet.h"

namespace WireCell::Aux {


    // Some crufty stuff used in PointTreeBuilding and UbooneClusterSource.
    PointCloud::Dataset
    make_scalar_dataset(const IBlob::pointer iblob, const Point& center,
                        const int npoints = 0, const double tick_span = 0.5*units::us);

    // Calculate the average position of a point cloud tree.
    WireCell::Point calc_blob_center(const PointCloud::Dataset& ds);

    // Calculate a dataset of blob corners
    PointCloud::Dataset make_corner_dataset(const IBlob::pointer iblob);

    double time2drift(const IAnodeFace::pointer anodeface, const double time_offset, const double drift_speed,
                      double time);

    void add_ctpc(PointCloud::Tree::Points::node_t& root, const IBlobSet::vector ibsv,
                  const IAnodeFace::pointer iface, const int face = 0,
                  const double time_offset = -1600 * units::us + 6 * units::mm / (1.101 * units::mm / units::us),
                  const double drift_speed = 1.101 * units::mm / units::us,
                  const double tick = 0.5 * units::us,
                  const double dead_threshold = 1e10);
}
