#include "WireCellClus/PointTreeBuilding.h"
#include "WireCellImg/Projection2D.h"
#include "WireCellClus/Facade.h"
#include "WireCellClus/Facade_Util.h"
#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/RayTiling.h"
#include "WireCellUtil/GraphTools.h"
#include "WireCellUtil/NamedFactory.h"

#include "WireCellAux/ClusterHelpers.h"
#include "WireCellAux/TensorDMpointtree.h"
#include "WireCellAux/TensorDMcommon.h"

WIRECELL_FACTORY(PointTreeBuilding, WireCell::Clus::PointTreeBuilding,
                 WireCell::INamed,
                 WireCell::IClusterFaninTensorSet,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::GraphTools;
using namespace WireCell::Clus;
/// FIXME : move this to aux some day
using WireCell::Img::Projection2D::get_geom_clusters;
using namespace WireCell::Aux;
using namespace WireCell::Aux::TensorDM;
using namespace WireCell::PointCloud;
using namespace WireCell::PointCloud::Tree;
using WireCell::PointCloud::Dataset;
using WireCell::PointCloud::Array;

PointTreeBuilding::PointTreeBuilding()
    : Aux::Logger("PointTreeBuilding", "clus")
{
}


PointTreeBuilding::~PointTreeBuilding()
{
}

std::vector<std::string> Clus::PointTreeBuilding::input_types()
{
    const std::string tname = std::string(typeid(input_type).name());
    std::vector<std::string> ret(m_multiplicity, tname);
    return ret;
}

void PointTreeBuilding::configure(const WireCell::Configuration& cfg)
{
    int m = get<int>(cfg, "multiplicity", (int) m_multiplicity);
    if (m != 1 and m != 2) {
        raise<ValueError>("multiplicity must be 1 or 2");
    }
    m_multiplicity = m;

    m_tags.resize(m);

    // Tag entire input frame worth of traces in the output frame.
    auto jtags = cfg["tags"];
    for (int ind = 0; ind < m; ++ind) {
        m_tags[ind] = convert<std::string>(jtags[ind], "");
    }

    m_datapath = get(cfg, "datapath", m_datapath);

    log->debug("using anode plane: {}", cfg["anode"].asString());
    m_anode = Factory::find_tn<IAnodePlane>(cfg["anode"].asString());
    if (!m_anode) {
        raise<ValueError>("failed to get anode plane");
    }

    m_face = get<int>(cfg, "face", 0);
    log->debug("using face: {}", m_face);
    if (m_anode->face(m_face) == nullptr) {
        raise<ValueError>("failed to get face %d", m_face);
    }

    m_geomhelper = Factory::find_tn<IClusGeomHelper>(cfg["geom_helper"].asString());

    auto samplers = cfg["samplers"];
    if (samplers.isNull()) {
        raise<ValueError>("add at least one entry to the \"samplers\" configuration parameter");
    }

    for (auto name : samplers.getMemberNames()) {
        auto tn = samplers[name].asString();
        if (tn.empty()) {
            raise<ValueError>("empty type/name for sampler \"%s\"", name);
        }
        log->debug("point cloud \"{}\" will be made by sampler \"{}\"",
                   name, tn);
        m_samplers[name] = Factory::find_tn<IBlobSampler>(tn); 
    }
    if (m_samplers.find("3d") == m_samplers.end()) {
        raise<ValueError>("m_samplers must have \"3d\" sampler");
    }

}


WireCell::Configuration PointTreeBuilding::default_configuration() const
{
    Configuration cfg;
    // eg:
    //    cfg["samplers"]["samples"] = "BlobSampler";
    cfg["datapath"] = m_datapath;
    return cfg;
}

namespace {

// unused....
#if 0
    std::string dump_ds(const WireCell::PointCloud::Dataset& ds) {
        std::stringstream ss;
        for (const auto& key : ds.keys()) {;
            const auto& arr = ds.get(key);
            ss << " {" << key << ":" << arr->dtype() << ":" << arr->shape()[0] << "} ";
            // const auto& arr = ds.get(key)->elements<float>();
            // for(auto elem : arr) {
            //     ss << elem << " ";
            // }
        }
        return ss.str();
    }
    // dump a NaryTree node
    std::string dump_node(const WireCell::PointCloud::Tree::Points::node_t* node)
    {
        std::stringstream ss;
        ss << "node: " << node;
        if (node) {
            const auto& lpcs = node->value.local_pcs();
            ss << " with " << lpcs.size() << " local pcs";
            for (const auto& [name, pc] : lpcs) {
                ss << " " << name << ": " << dump_ds(pc);
            }
        } else {
            ss << " null";
        }
        return ss.str();
    }
    // dump childs of a NaryTree node
    std::string dump_children(const WireCell::PointCloud::Tree::Points::node_t* root)
    {
        std::stringstream ss;
        ss << "NaryTree: " << root->nchildren() << " children";
        const auto first = root->children().front();
        ss << dump_node(first);
        return ss.str();
    }
#endif

    // Calculate the average position of a point cloud tree.
    Point calc_blob_center(const Dataset& ds)
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
    /// TODO: add more info to the dataset
    Dataset make_scaler_dataset(const IBlob::pointer iblob, const Point& center, const int npoints = 0, const double tick_span = 0.5*units::us)
    {
        using float_t = Facade::float_t;
        using int_t = Facade::int_t;
        Dataset ds;
        ds.add("charge", Array({(float_t)iblob->value()}));
        ds.add("center_x", Array({(float_t)center.x()}));
        ds.add("center_y", Array({(float_t)center.y()}));
        ds.add("center_z", Array({(float_t)center.z()}));
	ds.add("npoints", Array({(int_t)npoints}));
        const auto& islice = iblob->slice();
        // fixme: possible risk of roundoff error + truncation makes _min == _max?
        ds.add("slice_index_min", Array({(int_t)(islice->start()/tick_span)})); // unit: tick
        ds.add("slice_index_max", Array({(int_t)((islice->start()+islice->span())/tick_span)}));
        const auto& shape = iblob->shape();
        const auto& strips = shape.strips();
        /// ASSUMPTION: is this always true?
        std::unordered_map<RayGrid::layer_index_t, std::string> layer_names = {
            {2, "u"},
            {3, "v"},
            {4, "w"}
        };
        for (const auto& strip : strips) {
            // std::cout << "layer " << strip.layer << " bounds " << strip.bounds.first << " " << strip.bounds.second << std::endl;
            if(layer_names.find(strip.layer) == layer_names.end()) {
                continue;
            }
            ds.add(layer_names[strip.layer]+"_wire_index_min", Array({(int_t)strip.bounds.first}));
            ds.add(layer_names[strip.layer]+"_wire_index_max", Array({(int_t)strip.bounds.second}));
        }
        return ds;
    }

    /// extract corners
    /// if drift is true, the corner x would be drifted x insetead of the wire plane x
    Dataset make_corner_dataset(const IBlob::pointer iblob, const bool drift = false, const double time_offset = 0,
                                const double drift_speed = 1.6 * units::mm / units::us)
    {
        using float_t = double;

        Dataset ds;
        const auto& shape = iblob->shape();
        const auto& crossings = shape.corners();
        const auto& anodeface = iblob->face();
        const auto& coords = anodeface->raygrid();

        // ray center
        Facade::geo_point_t center;
        for (const auto& crossing : crossings) {
            const auto& [one, two] = crossing;
            auto pt = coords.ray_crossing(one, two);
            center += pt;
        }
        center = center / crossings.size();

        std::vector<float_t> corner_x;
        std::vector<float_t> corner_y;
        std::vector<float_t> corner_z;

        std::vector<float_t> corner_x_start;
        std::vector<float_t> corner_x_end;
        
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
            if (drift) {
                double x_start = Facade::time2drift(iblob->face(), time_offset, drift_speed, iblob->slice()->start());
                double x_end = Facade::time2drift(iblob->face(), time_offset, drift_speed, iblob->slice()->start() + iblob->slice()->span());
                corner_x_start.push_back(x_start);
                corner_x_end.push_back(x_end);
                corner_y.push_back(opt.y());
                corner_z.push_back(opt.z());
                // std::cout << "[delete] start" << iblob->slice()->start() << " time offset: " << time_offset << " drift speed: " << drift_speed << " start: " << x_start << " end: " << x_end << std::endl;

            } else {
                corner_x.push_back(opt.x());
                corner_y.push_back(opt.y());
                corner_z.push_back(opt.z());
            }
            // std::cout << "orig: " << one.layer << " " << one.grid << ", " << two.layer << " " << two.grid
            //           << " new: " << o1.grid << " " << o2.grid << " corner: " << opt.x() << " " << opt.y() << " "
            //           << opt.z() << std::endl;
        }
        if (drift) {
            corner_x.reserve(corner_x_start.size() + corner_x_end.size());
            corner_x.insert(corner_x.end(), corner_x_start.begin(), corner_x_start.end());
            corner_x.insert(corner_x.end(), corner_x_end.begin(), corner_x_end.end());
            corner_y.reserve(corner_y.size() * 2);
            corner_y.insert(corner_y.end(), corner_y.begin(), corner_y.end());
            corner_z.reserve(corner_z.size() * 2);
            corner_z.insert(corner_z.end(), corner_z.begin(), corner_z.end());
        }

        ds.add("x", Array(corner_x));
        ds.add("y", Array(corner_y));
        ds.add("z", Array(corner_z));
        
        return ds;
    }
}

static Dataset make2dds (const Dataset& ds3d, const double angle) {
    Dataset ds;
    const auto& x = ds3d.get("x")->elements<Facade::float_t>();
    const auto& y = ds3d.get("y")->elements<Facade::float_t>();
    const auto& z = ds3d.get("z")->elements<Facade::float_t>();
    std::vector<Facade::float_t> x2d(x.size());
    std::vector<Facade::float_t> y2d(y.size());
    for (size_t ind=0; ind<x.size(); ++ind) {
        const auto& xx = x[ind];
        const auto& yy = y[ind];
        const auto& zz = z[ind];
        x2d[ind] = xx;
        y2d[ind] = cos(angle) * zz - sin(angle) * yy;
    }
    ds.add("x", Array(x2d));
    ds.add("y", Array(y2d));
    return ds;
}

// Points::node_ptr PointTreeBuilding::sample_live(const WireCell::ICluster::pointer icluster, const double tick, const double angle_u, const double angle_v, const double angle_w) const {
void PointTreeBuilding::sample_live(Points::node_ptr& root, const WireCell::ICluster::pointer icluster) const {
    auto grouping = root->value.facade<Facade::Grouping>();
    const auto& tp = grouping->get_params();
    using int_t = Facade::int_t;
    const auto& gr = icluster->graph();
    log->debug("load cluster {} at call={}: {}", icluster->ident(), m_count, dumps(gr));

    auto clusters = get_geom_clusters(gr);
    log->debug("got {} clusters", clusters.size());
    size_t nblobs = 0;
    // Points::node_ptr root = std::make_unique<Points::node_t>();
    auto& sampler = m_samplers.at("3d");
    for (auto& [cluster_id, vdescs] : clusters) {
        auto cnode = root->insert();
        for (const auto& vdesc : vdescs) {
            const char code = gr[vdesc].code();
            if (code != 'b') {
                continue;
            }
            const IBlob::pointer iblob = std::get<IBlob::pointer>(gr[vdesc].ptr);
            named_pointclouds_t pcs;
            /// TODO: use nblobs or iblob->ident()?  A: Index.  The sampler takes blob->ident() as well.
            auto [pc3d, aux] = sampler->sample_blob(iblob, nblobs);
            pcs.emplace("3d", pc3d);
            pcs.emplace("2dp0", make2dds(pc3d, tp.angle_u));
            pcs.emplace("2dp1", make2dds(pc3d, tp.angle_v));
            pcs.emplace("2dp2", make2dds(pc3d, tp.angle_w));
            pcs.emplace("corner", make_corner_dataset(iblob, true, tp.time_offset, tp.drift_speed));
            const Point center = calc_blob_center(pcs["3d"]);
            // std::cout << "[delete] center: " << center.x() << " " << center.y() << " " << center.z() << std::endl;
            auto scaler_ds = make_scaler_dataset(iblob, center, pcs["3d"].get("x")->size_major(), tp.tick);
            int_t max_wire_interval = aux.get("max_wire_interval")->elements<int_t>()[0];
            int_t min_wire_interval = aux.get("min_wire_interval")->elements<int_t>()[0];
            int_t max_wire_type = aux.get("max_wire_type")->elements<int_t>()[0];
            int_t min_wire_type = aux.get("min_wire_type")->elements<int_t>()[0];
            scaler_ds.add("max_wire_interval", Array({(int_t)max_wire_interval}));
            scaler_ds.add("min_wire_interval", Array({(int_t)min_wire_interval}));
            scaler_ds.add("max_wire_type", Array({(int_t)max_wire_type}));
            scaler_ds.add("min_wire_type", Array({(int_t)min_wire_type}));
            pcs.emplace("scalar", std::move(scaler_ds));
            cnode->insert(Points(std::move(pcs)));

            ++nblobs;
        }
    }
    
    log->debug("sampled {} live blobs to tree with {} children", nblobs, root->nchildren());
    // return root;
}

Points::node_ptr PointTreeBuilding::sample_dead(const WireCell::ICluster::pointer icluster, const double tick) const {
    using int_t = Facade::int_t;
    const auto& gr = icluster->graph();
    log->debug("load cluster {} at call={}: {}", icluster->ident(), m_count, dumps(gr));

    auto clusters = get_geom_clusters(gr);
    log->debug("got {} clusters", clusters.size());
    size_t nblobs = 0;
    Points::node_ptr root = std::make_unique<Points::node_t>();
    // if (m_samplers.find("dead") == m_samplers.end()) {
    //     raise<ValueError>("m_samplers must have \"dead\" sampler");
    // }
    // auto& sampler = m_samplers.at("dead");
    for (auto& [cluster_id, vdescs] : clusters) {
        auto cnode = root->insert(std::make_unique<Points::node_t>());
        for (const auto& vdesc : vdescs) {
            const char code = gr[vdesc].code();
            if (code != 'b') {
                continue;
            }
            auto iblob = std::get<IBlob::pointer>(gr[vdesc].ptr);
            named_pointclouds_t pcs;
            auto scaler_ds = make_scaler_dataset(iblob, {0,0,0}, 0, tick);
            scaler_ds.add("max_wire_interval", Array({(int_t)-1}));
            scaler_ds.add("min_wire_interval", Array({(int_t)-1}));
            scaler_ds.add("max_wire_type", Array({(int_t)-1}));
            scaler_ds.add("min_wire_type", Array({(int_t)-1}));
            pcs.emplace("scalar", scaler_ds);
            pcs.emplace("corner", make_corner_dataset(iblob));
            // for (const auto& [name, pc] : pcs) {
            //     log->debug("{} -> keys {} size_major {}", name, pc.keys().size(), pc.size_major());
            // }
            cnode->insert(Points(std::move(pcs)));
            ++nblobs;
        }
        /// DEBUGONLY
        // if (nblobs > 1) {
        //     break;
        // }
    }
    
    log->debug("sampled {} dead blobs to tree with {} children", nblobs, root->nchildren());
    return root;
}

void PointTreeBuilding::add_ctpc(Points::node_ptr& root, const WireCell::ICluster::pointer icluster) const {
    using slice_t = WireCell::cluster_node_t::slice_t;
    using float_t = Facade::float_t;
    using int_t = Facade::int_t;

    const auto& cg = icluster->graph();
    // log->debug("add_ctpc load cluster {} at call={}: {}", icluster->ident(), m_count, dumps(cg));

    auto grouping = root->value.facade<Facade::Grouping>();
    const auto& tp = grouping->get_params();
    const auto& proj_centers = grouping->proj_centers();
    const auto& pitch_mags = grouping->pitch_mags();
    /// DEBUGONLY: remove these prints after debugging
    // for(const auto& [face, mags] : pitch_mags) {
    //     for(const auto& [pind, mag] : mags) {
    //         log->debug("face {} pind {} pitch_mag {}", face, pind, mag);
    //     }
    // }
    // for (const auto& [face, centers] : proj_centers) {
    //     for(const auto& [pind, center] : centers) {
    //         log->debug("face {} pind {} center {}", face, pind, center);
    //     }
    // }

    Facade::mapfp_t<std::vector<float_t>> ds_x, ds_y, ds_charge, ds_charge_err;
    Facade::mapfp_t<std::vector<int_t>> ds_cident, ds_wind, ds_slice_index;

    size_t nslices = 0;
    for (const auto& vdesc : GraphTools::mir(boost::vertices(cg))) {
        const auto& cgnode = cg[vdesc];
        if (cgnode.code() == 's') {
            auto& slice = std::get<slice_t>(cgnode.ptr);
            ++nslices;
            const auto& slice_index = slice->start()/tp.tick;
            const auto& activity = slice->activity();
            for (const auto& [ichan, charge] : activity) {
                if(charge.uncertainty() > m_dead_threshold) {
                    // if (charge.value() >0)
                    // std::cout << "Test: m_dead_threshold " << m_dead_threshold << " charge.uncertainty() " << charge.uncertainty() << " " << charge.value() << " " << ichan << " " << slice_index << std::endl;
                    continue;
                } 
                const auto& cident = ichan->ident();
                const auto& wires = ichan->wires();
                for (const auto& wire : wires) {
                    const auto& wind = wire->index();
                    const auto& plane = wire->planeid().index();
                    // log->debug("slice {} chan {} charge {} wind {} plane {} face {}", slice_index, cident, charge, wind, plane, wire->planeid().face());
                    // const auto& face = wire->planeid().face();
                    const auto& face = m_face;
                    /// FIXME: is this the way to get face?

//                    std::cout << "Test: " << slice->start() <<  " " << slice_index << " " << tp.time_offset << " " << tp.drift_speed << std::endl;

                    const auto& x = Facade::time2drift(m_anode->face(face), tp.time_offset, tp.drift_speed, slice->start());
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
            //            x.size(), ds_y[face][plane].size(), ds_charge[face][plane].size(), ds_charge_err[face][plane].size(),
            //            ds_cident[face][plane].size(), ds_wind[face][plane].size(), ds_slice_index[face][plane].size());
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
            root->value.local_pcs().emplace(ds_name, ds);
            // log->debug("added point cloud {} with {} points", ds_name, x.size());
        }
    }
    // for (const auto& [name, pc] : root->value.local_pcs()) {
    //     log->debug("contains point cloud {} with {} points", name, pc.get("x")->size_major());
    // }
}

void PointTreeBuilding::add_dead_winds(Points::node_ptr& root, const WireCell::ICluster::pointer icluster) const {
    using slice_t = WireCell::cluster_node_t::slice_t;
    using float_t = Facade::float_t;
    using int_t = Facade::int_t;
    const auto& cg = icluster->graph();
    auto grouping = root->value.facade<Facade::Grouping>();
    const auto& tp = grouping->get_params();
    std::set<int> faces;
    std::set<int> planes;
    for (const auto& vdesc : GraphTools::mir(boost::vertices(cg))) {
        const auto& cgnode = cg[vdesc];
        if (cgnode.code() != 's') continue;
        auto& slice = std::get<slice_t>(cgnode.ptr);
        // const auto& slice_index = slice->start()/m_tick;
        const auto& activity = slice->activity();
        for (const auto& [ichan, charge] : activity) {
            if(charge.uncertainty() < m_dead_threshold) continue;
            // log->debug("m_dead_threshold {} charge.uncertainty() {}", m_dead_threshold, charge.uncertainty());
            // const auto& cident = ichan->ident();
            const auto& wires = ichan->wires();
            for (const auto& wire : wires) {
                const auto& wind = wire->index();
                const auto& plane = wire->planeid().index();
                // const auto& face = wire->planeid().face();
                // log->debug("dead chan {} charge {} wind {} plane {} face {}", ichan->ident(), charge, wind, plane, wire->planeid().face());
                const auto& face = m_face;
                /// FIXME: is this the way to get face?
                const auto& xbeg = Facade::time2drift(m_anode->face(face), tp.time_offset, tp.drift_speed, slice->start());
                const auto& xend = Facade::time2drift(m_anode->face(face), tp.time_offset, tp.drift_speed, slice->start() + slice->span());
                // if (true) {
                //     log->debug("dead chan {} slice_index_min {} slice_index_max {} charge {} xbeg {} xend {}", ichan->ident(),
                //                slice_index, (slice->start() + slice->span()) / m_tick, charge, xbeg, xend);
                // }
                faces.insert(face);
                planes.insert(plane);

                auto & dead_winds = grouping->get_dead_winds(face, plane);
                // fix a bug how do we know the smaller or bigger value of xbeg and xend?
                // if (dead_winds.find(wind) == dead_winds.end()) {
                //     dead_winds[wind] = {xbeg, xend};
                // } else {
                //     const auto& [xbeg_now, xend_now] = dead_winds[wind];
                //     dead_winds[wind] = {std::min(xbeg, xbeg_now), std::max(xend, xend_now)};
                // }
                if (dead_winds.find(wind) == dead_winds.end()) {
                    dead_winds[wind] = {std::min(xbeg,xend)-0.1*units::cm, std::max(xbeg,xend) + 0.1*units::cm};
                } else {
                    const auto& [xbeg_now, xend_now] = dead_winds[wind];
                    dead_winds[wind] = {std::min(std::min(xbeg,xend)-0.1*units::cm, xbeg_now), std::max(std::max(xbeg,xend) + 0.1*units::cm, xend_now)};
                }
                


                // if (cident == 903) {
                //     log->debug("wind {} xbeg {} xend {}", wind, dead_winds[wind].first, dead_winds[wind].second);
                // }
            }
        }
    }
    /// DEBUGONLY:
    // for (int pind = 0; pind < 2; ++pind) {
    //     const auto& dead_winds = grouping->get_dead_winds(0, pind);
    //     for (const auto& [wind, xbeg_xend] : dead_winds) {
    //         log->debug("dead wind {} xbeg {} xend {}", wind, xbeg_xend.first, xbeg_xend.second);
    //     }
    // }
    log->debug("got dead winds {} {} {} ", grouping->get_dead_winds(0, 0).size(), grouping->get_dead_winds(0, 1).size(),
               grouping->get_dead_winds(0, 2).size());

    Facade::mapfp_t<std::vector<float_t>> xbegs, xends;
    Facade::mapfp_t<std::vector<int_t>> winds;
    for (const auto& face : faces) {
        for (const auto& plane : planes) {
            for (const auto& [wind, xbeg_xend] : grouping->get_dead_winds(face, plane)) {
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
            // root->insert(Points(named_pointclouds_t{{ds_name, std::move(ds)}}));
            root->value.local_pcs().emplace(ds_name, ds);
            // log->debug("added point cloud {} with {} points", ds_name, xbeg.size());
        }
    }
    for (const auto& [name, pc] : root->value.local_pcs()) {
        if (name.find("dead_winds") != std::string::npos) {
            log->debug("contains point cloud {} with {} points", name, pc.get("xbeg")->size_major());
        }
    }

}

bool PointTreeBuilding::operator()(const input_vector& invec, output_pointer& tensorset)
{
    tensorset = nullptr;

    size_t neos = 0;
    for (const auto& in : invec) {
        if (!in) {
            ++neos;
        }
    }
    if (neos == invec.size()) {
        // all inputs are EOS, good.
        log->debug("EOS at call {}", m_count++);
        return true;
    }
    if (neos) {
        raise<ValueError>("missing %d input tensors ", neos);
    }

    if (invec.size() != m_multiplicity) {
        raise<ValueError>("unexpected multiplicity got %d want %d",
                          invec.size(), m_multiplicity);
        return true;
    }


    const auto& iclus_live = invec[0];
    const int ident = iclus_live->ident();
    std::string datapath = m_datapath;
    if (datapath.find("%") != std::string::npos) {
        datapath = String::format(datapath, ident);
    }

    const auto& tp_json = m_geomhelper->get_params(m_anode->ident(), m_face);
    Points::node_ptr root_live = std::make_unique<Points::node_t>();
    auto grouping = root_live->value.facade<Facade::Grouping>();
    grouping->set_anode(m_anode);
    grouping->set_params(tp_json);
    // {
    //     const auto& tp = grouping->get_params();
    //     std::cout << "[delete] sample_live: tp.time_offset: " << tp.time_offset << " tp.drift_speed: " << tp.drift_speed << std::endl;
    // }
    sample_live(root_live, iclus_live);
    // Points::node_ptr root_live = sample_live(iclus_live, tp_json["tick"].asDouble(), tp_json["angle_u"].asDouble(),
    //                                          tp_json["angle_v"].asDouble(), tp_json["angle_w"].asDouble());
    add_ctpc(root_live, iclus_live);
    add_dead_winds(root_live, iclus_live);
    /// TODO: remove after debugging
    // {
    //     for (const auto& [name, pc] : root_live->value.local_pcs()) {
    //         log->debug("contains point cloud {} with {} points", name, pc.get("x")->size_major());
    //     }
    //     /// test ctpc_f0p0 exists
    //     grouping->kd2d(0,0);

    //     /// find test point on ctpc
    //     const auto ctest = grouping->children().front();
    //     const auto p3ds = ctest->points();
    //     log->debug("p3ds.size() {}", p3ds[0].size());
    //     {
    //         const auto winds = ctest->wire_indices();
    //         log->debug("winds.size() {}", winds[0].size());
    //         log->debug("ctest point x {} y {} z {}", p3ds[0][0], p3ds[1][0], p3ds[2][0]);
    //         log->debug("ctest winds {} {} {}", winds[0][0], winds[1][0], winds[2][0]);
    //         const double radius = 0.6 * units::cm;
    //         auto ret0 = grouping->get_closest_points({p3ds[0][0], p3ds[1][0], p3ds[2][0]}, radius, 0, 0);
    //         auto ret1 = grouping->get_closest_points({p3ds[0][0], p3ds[1][0], p3ds[2][0]}, radius, 0, 1);
    //         auto ret2 = grouping->get_closest_points({p3ds[0][0], p3ds[1][0], p3ds[2][0]}, radius, 0, 2);
    //         log->debug("closest points u {} v {} w {}", ret0.size(), ret1.size(), ret2.size());
    //         const auto& ctpc = root_live->value.local_pcs().at("ctpc_f0p0");
    //         const auto& x = ctpc.get("x")->elements<Facade::float_t>();
    //         const auto& y = ctpc.get("y")->elements<Facade::float_t>();
    //         const auto& slice_index = ctpc.get("slice_index")->elements<Facade::int_t>();
    //         const auto& wind = ctpc.get("wind")->elements<Facade::int_t>();
    //         for (const auto& [ind, dist] : ret0) {
    //             log->debug("ind {} dist {} x {} y {} slice_index {} wind {}", ind, dist, x[ind], y[ind], slice_index[ind], wind[ind]);
    //         }
    //     }

    //     {
    //         const auto tw0 = grouping->convert_3Dpoint_time_ch({p3ds[0][0], p3ds[1][0], p3ds[2][0]}, 0, 0);
    //         const auto tw1 = grouping->convert_3Dpoint_time_ch({p3ds[0][0], p3ds[1][0], p3ds[2][0]}, 0, 1);
    //         const auto tw2 = grouping->convert_3Dpoint_time_ch({p3ds[0][0], p3ds[1][0], p3ds[2][0]}, 0, 2);
    //         log->debug("tind {} wind {}", std::get<0>(tw0), std::get<1>(tw0));
    //         log->debug("tind {} wind {}", std::get<0>(tw1), std::get<1>(tw1));
    //         log->debug("tind {} wind {}", std::get<0>(tw2), std::get<1>(tw2));
    //     }
    //     {
    //         bool d0 = grouping->get_closest_dead_chs({p3ds[0][0], p3ds[1][0], p3ds[2][0]}, 1, 0, 0);
    //         bool d1 = grouping->get_closest_dead_chs({p3ds[0][0], p3ds[1][0], p3ds[2][0]}, 1, 0, 1);
    //         bool d2 = grouping->get_closest_dead_chs({p3ds[0][0], p3ds[1][0], p3ds[2][0]}, 1, 0, 2);
    //         log->debug("dead chs {} {} {}", d0, d1, d2);

    //         bool is_good = grouping->is_good_point({p3ds[0][0], p3ds[1][0], p3ds[2][0]}, 0);
    //         log->debug("is_good_point {}", is_good);
    //     }
    //     // exit(0);
    // }
    // {
    //     auto grouping = root_live->value.facade<Facade::Grouping>();
    //     auto children = grouping->children(); // copy
    //     sort_clusters(children);
    //     size_t count=0;
    //     for(const auto* cluster : children) {
    //         bool sane = cluster->sanity(log);
    //         log->debug("live cluster {} {} sane:{}", count++, *cluster, sane);
    //     }
    // }
    auto tens_live = as_tensors(*root_live.get(), datapath+"/live");
    log->debug("Made {} live tensors", tens_live.size());
    for(const auto& ten : tens_live) {
        log->debug("tensor {} {}", ten->metadata()["datapath"].asString(), ten->size());
        break;
    }

    if (m_multiplicity == 2) {
        const auto& iclus_dead = invec[1];
        /// FIXME: what do we expect?
        if(ident != iclus_dead->ident()) {
            raise<ValueError>("ident mismatch between live and dead clusters");
        }
        Points::node_ptr root_dead = sample_dead(iclus_dead, tp_json["tick"].asDouble());
        /// DEBUGONLY:
        // {
        //     Facade::Grouping& dead_grouping = *root_dead->value.facade<Facade::Grouping>();
        //     // std::cout<< "dumping\n";
        //     // for (const auto cluster : dead_grouping.children()) {
        //     //     std::cout << cluster->dump() << std::endl;
        //     // }
        // }
        auto tens_dead = as_tensors(*root_dead.get(), datapath+"/dead");
        log->debug("Made {} dead tensors", tens_dead.size());
        for(const auto& ten : tens_dead) {
            log->debug("tensor {} {}", ten->metadata()["datapath"].asString(), ten->size());
            break;
        }
        /// TODO: is make_move_iterator faster?
        tens_live.insert(tens_live.end(), tens_dead.begin(), tens_dead.end());
    }

    log->debug("Total outtens {} tensors", tens_live.size());
    tensorset = as_tensorset(tens_live, ident);

    m_count++;
    return true;
}

        

