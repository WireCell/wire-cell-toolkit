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

double Aux::time2drift(const IAnodeFace::pointer anodeface, const double time_offset, const double drift_speed, double time) {
    // std::cout << "time2drift: " << time << " " << time_offset << " " << drift_speed << std::endl;
    const Pimpos* colpimpos = anodeface->planes()[2]->pimpos();
    double xsign = colpimpos->axis(0)[0];
    double xorig = anodeface->planes()[2]->wires().front()->center().x();
    const double drift = (time + time_offset)*drift_speed;
    /// TODO: how to determine xsign?
    // std::cout << "drift: " << xorig + xsign*drift << std::endl;
    return xorig + xsign*drift;
}

template <typename T>
using mapfp_t = std::unordered_map<int, std::unordered_map<int, T>>;
void Aux::add_ctpc(PointCloud::Tree::Points::node_t& root, const WireCell::IBlobSet::vector ibsv,
                   const IAnodeFace::pointer iface, const int face, const double time_offset, const double drift_speed,
                   const double tick, const double dead_threshold)
{
    mapfp_t<std::vector<double>> ds_x, ds_y, ds_charge, ds_charge_err;
    mapfp_t<std::vector<int>> ds_cident, ds_wind, ds_slice_index;
    mapfp_t<double> proj_centers;
    mapfp_t<double> pitch_mags;

    const auto& coords = iface->raygrid();
    const int ndummy_layers = 2;
    // skip dummy layers so the vector matches 0, 1, 2 plane order
    for (int layer=ndummy_layers; layer<coords.nlayers(); ++layer) {
        const auto& pitch_dir = coords.pitch_dirs()[layer];
        const auto& center = coords.centers()[layer];
        double proj_center = center.dot(pitch_dir);
        proj_centers[iface->which()][layer-ndummy_layers] = proj_center;
        pitch_mags[iface->which()][layer-ndummy_layers] = coords.pitch_mags()[layer];
    }

    size_t nslices = 0;
    for (const auto& ibs : ibsv) {
        const auto& slice = ibs->slice();
        {
            // auto& slice = std::get<slice_t>(cgnode.ptr);
            ++nslices;
            const auto slice_index = slice->start()/tick;
            const auto& activity = slice->activity();
            for (const auto& [ichan, charge] : activity) {
                if(charge.uncertainty() > dead_threshold) {
                    // if (charge.value() >0)
                    // std::cout << "Test: dead_threshold " << dead_threshold << " charge.uncertainty() " << charge.uncertainty() << " " << charge.value() << " " << ichan << " " << slice_index << std::endl;
                    continue;
                } 
                const auto& cident = ichan->ident();
                const auto& wires = ichan->wires();
                for (const auto& wire : wires) {
                    const auto& wind = wire->index();
                    const auto& plane = wire->planeid().index();
                    // log->debug("slice {} chan {} charge {} wind {} plane {} face {}", slice_index, cident, charge, wind, plane, wire->planeid().face());
                    // const auto& face = wire->planeid().face();
                    // const auto& face = m_face;
                    /// FIXME: is this the way to get face?

//                    std::cout << "Test: " << slice->start() <<  " " << slice_index << " " << tp.time_offset << " " << tp.drift_speed << std::endl;

                    const auto& x = time2drift(iface, time_offset, drift_speed, slice->start());
                    const double y = pitch_mags.at(face).at(plane)* (wind +0.5) + proj_centers.at(face).at(plane); // the additon of 0.5 is to match with the convetion of WCP (X. Q.)

                    // if (abs(wind-815) < 2 or abs(wind-1235) < 2 or abs(wind-1378) < 2) {
                    //     log->debug("slice {} chan {} charge {} wind {} plane {} face {} x {} y {}", slice_index, cident, charge,
                    //                wind, plane, face, x, y);
                    // }
                    ds_x[face][plane].push_back(x);
                    ds_y[face][plane].push_back(y);
                    ds_charge[face][plane].push_back(charge.value());
                    ds_charge_err[face][plane].push_back(charge.uncertainty());
                    ds_cident[face][plane].push_back(cident);
                    ds_wind[face][plane].push_back(wind);
                    ds_slice_index[face][plane].push_back(slice_index);
                }
            }
            // log->debug("ds_x.size() {}", ds_x.size());
        }
    }
    // log->debug("got {} slices", nslices);

    for (const auto& [face, planes] : ds_x) {
        for (const auto& [plane, x] : planes) {
            // log->debug("ds_x {} ds_y {} ds_charge {} ds_charge_err {} ds_cident {} ds_wind {} ds_slice_index {}",
            //            x.size(), ds_y[face][plane].size(), ds_charge[face][plane].size(),
            //            ds_charge_err[face][plane].size(), ds_cident[face][plane].size(), ds_wind[face][plane].size(),
            //            ds_slice_index[face][plane].size());
            // std::cout << "[yuhw] " << " face " << face << " plane " << plane << " x " << x.size() << " y "
            //           << ds_y[face][plane].size() << " charge " << ds_charge[face][plane].size() << " charge_err "
            //           << ds_charge_err[face][plane].size() << " cident " << ds_cident[face][plane].size() << " wind "
            //           << ds_wind[face][plane].size() << " slice_index " << ds_slice_index[face][plane].size()
            //           << std::endl;
            Dataset ds;
            ds.add("x", Array(x));
            ds.add("y", Array(ds_y[face][plane]));
            ds.add("charge", Array(ds_charge[face][plane]));
            ds.add("charge_err", Array(ds_charge_err[face][plane]));
            ds.add("cident", Array(ds_cident[face][plane]));
            ds.add("wind", Array(ds_wind[face][plane]));
            ds.add("slice_index", Array(ds_slice_index[face][plane]));
            const std::string ds_name = String::format("ctpc_f%dp%d", face, plane);
            // root->insert(Points(named_pointclouds_t{{ds_name, std::move(ds)}}));
            root.value.local_pcs().emplace(ds_name, ds);
            // log->debug("added point cloud {} with {} points", ds_name, x.size());
        }
    }
    // for (const auto& [name, pc] : root->value.local_pcs()) {
    //     log->debug("contains point cloud {} with {} points", name, pc.get("x")->size_major());
    // }
}

void Aux::add_dead_winds(PointCloud::Tree::Points::node_t& root, const IBlobSet::vector ibsv,
    const IAnodeFace::pointer iface, const int face ,
    const double time_offset ,
    const double drift_speed ,
    const double tick, const double dead_threshold){

    std::set<int> faces;
    std::set<int> planes;

    std::map<std::pair<int,int>, std::unordered_map<int, std::pair<double, double>> > map_dead_winds;
    mapfp_t<std::vector<double>> xbegs, xends;
    mapfp_t<std::vector<int>> winds;

    for (const auto& ibs : ibsv) {
        const auto& slice = ibs->slice();
        {
            // const auto& slice_index = slice->start()/tick;
            const auto& activity = slice->activity();
            for (const auto& [ichan, charge] : activity) {
                if(charge.uncertainty() < dead_threshold) continue;
                const auto& wires = ichan->wires();
                for (const auto& wire : wires) {
                    const auto& wind = wire->index();
                    const auto& plane = wire->planeid().index();
                    //                     const auto& x = time2drift(iface, time_offset, drift_speed, slice->start());
                    const auto& xbeg = time2drift(iface, time_offset, drift_speed, slice->start());
                    const auto& xend = time2drift(iface, time_offset, drift_speed, slice->start() + slice->span());

                    auto& dead_winds = map_dead_winds[std::make_pair(face, plane)];
                    if (dead_winds.find(wind) == dead_winds.end()) {
                        dead_winds[wind] = {std::min(xbeg,xend)-0.1*units::cm, std::max(xbeg,xend) + 0.1*units::cm};
                    } else {
                        const auto& [xbeg_now, xend_now] = dead_winds[wind];
                        dead_winds[wind] = {std::min(std::min(xbeg,xend)-0.1*units::cm, xbeg_now), std::max(std::max(xbeg,xend) + 0.1*units::cm, xend_now)};
                    }
                    faces.insert(face);
                    planes.insert(plane);

                }
            }
        }
    }
    
    for (const auto& face : faces) {
        for (const auto& plane : planes) {
            for (const auto& [wind, xbeg_xend] : map_dead_winds[std::make_pair(face, plane)]) {
                xbegs[face][plane].push_back(xbeg_xend.first);
                xends[face][plane].push_back(xbeg_xend.second);
                winds[face][plane].push_back(wind);
            }
        }
    }
    for (const auto& face : faces) {
        for (const auto& plane : planes) {
            Dataset ds;
            ds.add("xbeg", Array(xbegs[face][plane]));
            ds.add("xend", Array(xends[face][plane]));
            ds.add("wind", Array(winds[face][plane]));
            const std::string ds_name = String::format("dead_winds_f%dp%d", face, plane);
            root.value.local_pcs().emplace(ds_name, ds);
            // log->debug("added point cloud {} with {} points", ds_name, xbeg.size());
        }
    }

}
