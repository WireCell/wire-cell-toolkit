#include "WireCellAux/SamplingHelpers.h"

using namespace WireCell;
using WireCell::PointCloud::Dataset;
using WireCell::PointCloud::Array;
WireCell::PointCloud::Dataset
Aux::make_scalar_dataset(const IBlob::pointer iblob, const Point& center,
                         const int npoints, const double tick)
{
    Dataset ds;
    // Warning, these types must match consumers.  In particular, PointTreeBuilding.
    ds.add("charge", Array({(double)iblob->value()}));
    ds.add("center_x", Array({(double)center.x()}));
    ds.add("center_y", Array({(double)center.y()}));
    ds.add("center_z", Array({(double)center.z()}));
    ds.add("npoints", Array({(int)npoints}));
    const auto& islice = iblob->slice();
    // fixme: possible risk of roundoff error + truncation makes _min == _max?
    ds.add("slice_index_min", Array({(int)(islice->start()/tick)})); // unit: tick
    ds.add("slice_index_max", Array({(int)((islice->start()+islice->span())/tick)}));
    const auto& shape = iblob->shape();
    const auto& strips = shape.strips();
    /// ASSUMPTION: is this always true?
    std::unordered_map<RayGrid::layer_index_t, std::string> layer_names = {
        {2, "u"},
        {3, "v"},
        {4, "w"}
    };
    for (const auto& strip : strips) {
        if(layer_names.find(strip.layer) == layer_names.end()) {
            continue;
        }
        ds.add(layer_names[strip.layer]+"_wire_index_min", Array({(int)strip.bounds.first}));
        ds.add(layer_names[strip.layer]+"_wire_index_max", Array({(int)strip.bounds.second}));
    }
    return ds;
}

// Calculate the average position of a point cloud tree.
Point Aux::calc_blob_center(const Dataset& ds)
{
    const auto& arr_x = ds.get("x")->elements<Point::coordinate_t>();
    const auto& arr_y = ds.get("y")->elements<Point::coordinate_t>();
    const auto& arr_z = ds.get("z")->elements<Point::coordinate_t>();
    const size_t len = arr_x.size();
    if(len == 0) {
        raise<ValueError>("empty point cloud");
    }
    Point ret(0,0,0);
    for (size_t ind=0; ind<len; ++ind) {
        ret += Point(arr_x[ind], arr_y[ind], arr_z[ind]);
    }
    ret = ret / len;
    return ret;
}

