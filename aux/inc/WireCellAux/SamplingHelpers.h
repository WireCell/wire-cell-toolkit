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

    PointCloud::Dataset make2dds(const PointCloud::Dataset& ds3d, const double angle);

    // Calculate the average position of a point cloud tree.
    WireCell::Point calc_blob_center(const PointCloud::Dataset& ds);

    // Calculate a dataset of blob corners
    // if drift is true, the corner x would be drifted x insetead of the wire plane x
    PointCloud::Dataset make_corner_dataset(const IBlob::pointer iblob, const bool drift = false, const double time_offset = 0,
                            const double drift_speed = 1.6 * units::mm / units::us);



    double time2drift(const IAnodeFace::pointer anodeface, const double time_offset, const double drift_speed,
                      double time);

    void add_ctpc(PointCloud::Tree::Points::node_t& root, const IBlobSet::vector ibsv,
                  const IAnodeFace::pointer iface, const int face = 0,
                  const double time_offset = -1600 * units::us + 6 * units::mm / (1.101 * units::mm / units::us),
                  const double drift_speed = 1.101 * units::mm / units::us,
                  const double tick = 0.5 * units::us,
                  const double dead_threshold = 1e10);

    void add_dead_winds(PointCloud::Tree::Points::node_t& root, const IBlobSet::vector ibsv,
                    const IAnodeFace::pointer iface, const int face = 0,
                    const double time_offset = -1600 * units::us + 6 * units::mm / (1.101 * units::mm / units::us),
                    const double drift_speed = 1.101 * units::mm / units::us,
                    const double tick = 0.5 * units::us,
                    const double dead_threshold = 1e10
    );
}
