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



/// extract corners
Dataset Aux::make_corner_dataset(const IBlob::pointer iblob)
{
    using float_t = double;

    Dataset ds;
    const auto& shape = iblob->shape();
    const auto& crossings = shape.corners();
    const auto& anodeface = iblob->face();
    const auto& coords = anodeface->raygrid();

    // ray center
    WireCell::Point center;
    for (const auto& crossing : crossings) {
        const auto& [one, two] = crossing;
        auto pt = coords.ray_crossing(one, two);
        center += pt;
    }
    center = center / crossings.size();

    std::vector<float_t> corner_x;
    std::vector<float_t> corner_y;
    std::vector<float_t> corner_z;
        
    for (const auto& crossing : crossings) {
        const auto& [one, two] = crossing;
        RayGrid::coordinate_t o1 = one;
        RayGrid::coordinate_t o2 = two;
        // {
        //     auto pt = coords.ray_crossing(one, two);
        //     auto is_higher = [&](const RayGrid::coordinate_t& c) {
        //         if (c.layer < 2) {
        //             return false;
        //         }
        //         const double diff = (pt - center).dot(coords.pitch_dirs()[c.layer]);
        //         return diff > 0;
        //     };
        //     if (is_higher(one)) o1.grid += 1;
        //     if (is_higher(two)) o2.grid += 1;
        // }
        auto opt = coords.ray_crossing(o1, o2);
        corner_x.push_back(opt.x());
        corner_y.push_back(opt.y());
        corner_z.push_back(opt.z());
        // std::cout << "orig: " << one.layer << " " << one.grid << ", " << two.layer << " " << two.grid
        //           << " new: " << o1.grid << " " << o2.grid << " corner: " << opt.x() << " " << opt.y() << " "
        //           << opt.z() << std::endl;
    }

    ds.add("x", Array(corner_x));
    ds.add("y", Array(corner_y));
    ds.add("z", Array(corner_z));
        
    return ds;
}
