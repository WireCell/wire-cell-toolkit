#include "WireCellImg/PointTreeBuilding.h"
#include "WireCellImg/Projection2D.h"
#include "WireCellImg/PointCloudFacade.h"
#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/RayTiling.h"
#include "WireCellUtil/GraphTools.h"
#include "WireCellUtil/NamedFactory.h"

#include "WireCellAux/ClusterHelpers.h"
#include "WireCellAux/TensorDMpointtree.h"
#include "WireCellAux/TensorDMcommon.h"

WIRECELL_FACTORY(PointTreeBuilding, WireCell::Img::PointTreeBuilding,
                 WireCell::INamed,
                 WireCell::IClusterFaninTensorSet,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::GraphTools;
using namespace WireCell::Img;
using WireCell::Img::Projection2D::get_geom_clusters;
using namespace WireCell::Aux;
using namespace WireCell::Aux::TensorDM;
using namespace WireCell::PointCloud;
using namespace WireCell::PointCloud::Tree;
using WireCell::PointCloud::Dataset;
using WireCell::PointCloud::Array;

PointTreeBuilding::PointTreeBuilding()
    : Aux::Logger("PointTreeBuilding", "img")
{
}


PointTreeBuilding::~PointTreeBuilding()
{
}

std::vector<std::string> Img::PointTreeBuilding::input_types()
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
        using float_t = double;
        using int_t = int;
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
    Dataset make_corner_dataset(const IBlob::pointer iblob)
    {
        using float_t = double;

        Dataset ds;
        const auto& shape = iblob->shape();
        const auto& crossings = shape.corners();
        const auto& anodeface = iblob->face();
        std::vector<float_t> corner_x;
        std::vector<float_t> corner_y;
        std::vector<float_t> corner_z;
        
        for (const auto& crossing : crossings) {
            const auto& coords = anodeface->raygrid();
            const auto& [one, two] = crossing;
            auto pt = coords.ray_crossing(one, two);
            corner_x.push_back(pt.x());
            corner_y.push_back(pt.y());
            corner_z.push_back(pt.z());
        }
        
        ds.add("x", Array(corner_x));
        ds.add("y", Array(corner_y));
        ds.add("z", Array(corner_z));
        
        return ds;
    }
}

Points::node_ptr PointTreeBuilding::sample_live(const WireCell::ICluster::pointer icluster) const {
    const auto& gr = icluster->graph();
    log->debug("load cluster {} at call={}: {}", icluster->ident(), m_count, dumps(gr));

    auto clusters = get_geom_clusters(gr);
    log->debug("got {} clusters", clusters.size());
    size_t nblobs = 0;
    Points::node_ptr root = std::make_unique<Points::node_t>();
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
            pcs.emplace("3d", sampler->sample_blob(iblob, nblobs));
            const Point center = calc_blob_center(pcs["3d"]);
            const auto scaler_ds = make_scaler_dataset(iblob, center, pcs["3d"].get("x")->size_major(), m_tick);
            pcs.emplace("scalar", std::move(scaler_ds));
            cnode->insert(Points(std::move(pcs)));

            ++nblobs;
        }
    }
    
    log->debug("sampled {} live blobs to tree with {} children", nblobs, root->nchildren());
    return root;
}

Points::node_ptr PointTreeBuilding::sample_dead(const WireCell::ICluster::pointer icluster) const {
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
            // pcs.emplace("dead", sampler->sample_blob(iblob, nblobs));
            pcs.emplace("scalar", make_scaler_dataset(iblob, {0,0,0}, 0, m_tick));
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

    Points::node_ptr root_live = sample_live(iclus_live);
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
        Points::node_ptr root_dead = sample_dead(iclus_dead);
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

        

