#include "WireCellUtil/PointCloudDataset.h"
#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/Units.h"

#include "WireCellIface/IBlob.h"
#include "WireCellIface/IBlobSampler.h"
#include "WireCellIface/IBlobSet.h"

#include <vector>

namespace WireCell::Aux {

    /** Return a "sampling" of a live blob.

        These roll up the relevant fill_*() methods for "live" blobs.

        If ident<0 will use iblob.ident().

        An empty set of PCs is returned if there are no sample points.
     */
    PointCloud::Tree::named_pointclouds_t
    sample_live(const IBlobSampler::pointer& sampler, const IBlob::pointer& blob,
                const std::vector<double>& angles,
                const double tick=500*units::ns, int ident=-1);

    /** Return a "sampling" of a dead blob.

        These roll up the relevant fill_*() methods for "dead" blobs.
    */
    PointCloud::Tree::named_pointclouds_t
    sample_dead(const IBlob::pointer& blob, const double tick=500*units::ns);


    /** Add per-plane 2D point cloud arrays to the given pc that represent the
     * plane coordinates derived from the "x", "y" and "z" arrays in the given
     * "pc".  The "pattern" argument is used to name the arrays of the 2D PCs.
     * The default is likely best kept as-is unless there is some custom needs.
     * The first argument to the string interpolation on pattern is the index
     * that runs over the angles (ie, plane index).  The second argument is the
     * 2D coordinate character ('x' or 'y').
     */
    void fill_2dpcs(PointCloud::Dataset& pc, 
                    const std::vector<double>& angles,
                    const std::string& pattern="2dp%1%_%2%");
    

    /// Fill various types of per-blob "scalar" info.
    void fill_scalar_blob(PointCloud::Dataset& scalar, const IBlob& iblob, const double tick=500*units::ns);
    /// Fill "center" related items with values
    void fill_scalar_center(PointCloud::Dataset& scalar, const PointCloud::Dataset& pc3d);
    /// As above but fill with zeros.
    void fill_scalar_center(PointCloud::Dataset& scalar);
    /// Transfer "aux" values.  See #426.
    void fill_scalar_aux(PointCloud::Dataset& scalar, const PointCloud::Dataset& aux);
    /// As above but fill with zeros
    void fill_scalar_aux(PointCloud::Dataset& scalar);

    // Calculate the average position of a point cloud tree.
    WireCell::Point calc_blob_center(const PointCloud::Dataset& ds);

    // Calculate a dataset of blob corners
    PointCloud::Dataset make_corner_dataset(const IBlob& iblob);



    double time2drift(IAnodeFace::pointer anodeface, const double time_offset, const double drift_speed,
                      double time);

    void add_ctpc(PointCloud::Tree::Points::node_t& root,
                  const ISlice::vector& slices,
                  IAnodeFace::pointer iface, const int face = 0,
                  const double time_offset = -1600 * units::us + 6 * units::mm / (1.101 * units::mm / units::us),
                  const double drift_speed = 1.101 * units::mm / units::us,
                  const double tick = 0.5 * units::us,
                  const double dead_threshold = 1e10);

    void add_dead_winds(PointCloud::Tree::Points::node_t& root,
                        const ISlice::vector& slices,
                        IAnodeFace::pointer iface, const int face = 0,
                        const double time_offset = -1600 * units::us + 6 * units::mm / (1.101 * units::mm / units::us),
                        const double drift_speed = 1.101 * units::mm / units::us,
                        const double tick = 0.5 * units::us,
                        const double dead_threshold = 1e10
    );
}
