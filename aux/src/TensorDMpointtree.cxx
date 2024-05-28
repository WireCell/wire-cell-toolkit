#include "WireCellAux/TensorDMpointtree.h"
#include "WireCellAux/TensorDMcommon.h"
#include "WireCellAux/SimpleTensor.h"
#include <iostream>


using namespace WireCell;
using namespace WireCell::Aux;
using namespace WireCell::Aux::TensorDM;

WireCell::Aux::TensorDM::named_pointclouds_t
WireCell::Aux::TensorDM::as_pcnamedset(const ITensor::vector& tens, const std::string& datapath)
{
    TensorIndex ti(tens);
    return as_pcnamedset(ti, datapath);
}

WireCell::Aux::TensorDM::named_pointclouds_t
WireCell::Aux::TensorDM::as_pcnamedset(const TensorIndex& ti, const std::string& datapath)
{
    named_pointclouds_t ret;
    auto top = ti.at(datapath, "pcnamedset");
    const auto& md = top->metadata();
    auto items = md["items"];
    for (const auto& name : items.getMemberNames()) {
        const auto path = items[name].asString();

        auto ds = as_dataset(ti, path);
        ret.emplace(name, ds);
    }
    return ret;
}

// Type for storing the parentage indices.
using parentage_index_type = int;
// Type for major axis size of local PCs.  2 gigapoints is enough for anyone.
using lpcmaps_index_type = int;


ITensor::vector WireCell::Aux::TensorDM::as_tensors(
    const WireCell::PointCloud::Tree::Points::node_t& root,
    const std::string& datapath)
{
    auto descent = root.depth();

    // Lookup from node to index.
    std::unordered_map<const WireCell::PointCloud::Tree::Points::node_t*, size_t> nodesindex;

    // Store the flattened tree structure.  Each element represents a node in
    // the DFS order and holds index representing parent or self-index if root.
    std::vector<parentage_index_type> parentage;

    // Concatenated local point clouds in DFS order.  Key is PC name.
    PointCloud::Tree::named_pointclouds_t pointclouds;

    // Sizes of local PC contribution to the pointclouds.  Key is PC name.  Will
    // become a pcdataset.
    std::map<std::string, std::vector<lpcmaps_index_type>> lpcmaps;

    // Visit each node in tree in depth-first descent order.
    for (const auto& noderef : descent) {
        const auto* node = &noderef;
        const size_t index = nodesindex.size();
        nodesindex[node] = index;

        // Extend parentage map
        parentage.resize(index + 1);
        if (node->parent) {
            // DFS: we've already seen this parent so no need to find().
            parentage[index] = nodesindex[node->parent];
        }
        else {
            // No parent, root nodes have their own self as "parent".
            parentage[index] = index;
        }

        // Extend pointclouds and lpcmaps
        for (const auto& [pcname,pcds] : node->value.local_pcs()) {
            auto& lpcmap = lpcmaps[pcname];
            lpcmap.resize(index + 1, 0);
            lpcmap[index] = pcds.size_major();

            auto& cpc = pointclouds[pcname];
            cpc.append(pcds);
        }
    }

    // Metadata for the pctree.
    Configuration md;
    md["datapath"] = datapath;
    md["datatype"] = "pctree";
    // datapath to pcnamedset holding the concatenated local PC
    md["pointclouds"] = "";
    // datapath to pcdataset holding the local PC size maps
    md["lpcmaps"] = "";

    ITensor::vector ret;
    ret.push_back(nullptr);     // fill in below

    if (pointclouds.size()) {
        auto tens = pcnamedset_as_tensors(pointclouds.begin(), pointclouds.end(),
                                          datapath + "/pointclouds");
        md["pointclouds"] = tens[0]->metadata()["datapath"];
        ret.insert(ret.end(), tens.begin(), tens.end());

        PointCloud::Dataset ds;
        for (auto& [pcname, lpcmap] : lpcmaps) {
            ds.add(pcname, PointCloud::Array(lpcmap));
        }
        tens = as_tensors(ds, datapath + "/lpcmaps");
        md["lpcmaps"] = tens[0]->metadata()["datapath"];
        ret.insert(ret.end(), tens.begin(), tens.end());
    }

    ITensor::shape_t shape = { parentage.size() };
    ret[0] = std::make_shared<SimpleTensor>(shape, parentage.data(), md);
    return ret;
}

std::unique_ptr<WireCell::PointCloud::Tree::Points::node_t>
WireCell::Aux::TensorDM::as_pctree(const ITensor::vector& tens,
                                   const std::string& datapath)
{
    TensorIndex ti(tens);
    return as_pctree(ti, datapath);
}

std::unique_ptr<WireCell::PointCloud::Tree::Points::node_t>
WireCell::Aux::TensorDM::as_pctree(const TensorIndex& ti,
                                   const std::string& datapath)
{
    using WireCell::PointCloud::Tree::Points;

    auto top = ti.at(datapath, "pctree");
    const auto parentage = reinterpret_cast<const parentage_index_type*>(top->data());
    const size_t nnodes = top->size() / sizeof(parentage_index_type);

    std::vector<Points::node_t*> nodes;

    // Build the tree
    Points::node_t::owned_ptr root = std::make_unique<Points::node_t>();
    nodes.push_back(root.get());
    for (size_t index=1; index<nnodes; ++index) {
        size_t parent_index = parentage[index];
        if (parent_index == index) {
            // while multi-root trees are supported in general, this function
            // only returns the first.
            break;
        }
        auto* parent = nodes[parent_index];
        auto* node = parent->insert();
        nodes.push_back(node);
    }

    auto const& md = top->metadata();        
    auto pointclouds = as_pcnamedset(ti, md["pointclouds"].asString());
    auto lpcmaps_ds = as_dataset(ti, md["lpcmaps"].asString());

    // Loop cross product of (PC name,node)
    for (const auto& [pcname, pcds] : pointclouds) {

        // Get local PC map vector as a span on DS array.
        auto lpcmap = lpcmaps_ds.get(pcname)->elements<lpcmaps_index_type>();

        // Count number of points seen as we walk through the node list.
        int offset = 0;
        for (size_t index=0; index < nnodes; ++index) {
            const int npoints = lpcmap[index];
            if (!npoints) {
                continue;       // omit an empty local PC
            }
            
            auto node = nodes[index];

            // make the local PC from slice of the aggregate.
            auto& lpcs = node->value.local_pcs();
            lpcs[pcname] = pcds.slice(offset, npoints);

            offset += npoints;
        }
    }

    return root;
}
