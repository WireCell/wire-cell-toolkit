#include "WireCellAux/TensorDMcluster.h"
#include "WireCellAux/ClusterArrays.h"
#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleCluster.h"
#include "WireCellUtil/Exceptions.h"

using namespace WireCell;
using namespace WireCell::Aux;
using namespace WireCell::Aux::ClusterArrays;

template<typename ElementType>
ITensor::pointer make_array_tensor(const boost::multi_array<ElementType, 2>& arr,
                                   const std::string& datatype,
                                   const std::string& datapath)
{
    Json::Value md = Json::objectValue;
    md["datatype"] = datatype;
    md["datapath"] = datapath;
    std::vector<size_t> shape = {arr.shape()[0], arr.shape()[1]};
    return std::make_shared<SimpleTensor>(shape, arr.data(), md);
}

ITensor::vector TensorDM::as_tensors(ICluster::pointer cluster,
                                     const std::string& datapath)
{
    ITensor::vector ret;

    // cluster
    const int ident = cluster->ident();
    Json::Value md = Json::objectValue;
    md["datatype"] = "cluster";
    md["datapath"] = datapath;
    md["ident"] = ident;
    md["nodes"] = Json::objectValue;
    md["edges"] = Json::objectValue;
    ret.push_back(nullptr);     // fill in below

    node_array_set_t nas;
    edge_array_set_t eas;
    to_arrays(cluster->graph(), nas, eas);

    {                           // translate node arrays to tensors
        const std::string root = datapath + "/nodes/";
        for (const auto& [nc, na] : nas) {
            const std::string ncs(1, nc);
            const std::string path = root + ncs;
            md["nodes"][ncs] = path;
            ret.push_back(make_array_tensor<double>(na, "clnodeset", path));
        }
    }

    {                           // translate edge arrays to tensors
        const std::string root = datapath + "/edges/";
        for (const auto& [ec, ea] : eas) {
            const std::string ecs = to_string(ec);
            const std::string path = root + ecs;
            md["edges"][ecs] = path;
            ret.push_back(make_array_tensor<int>(ea, "cledgeset", path));
        }
    }

    // Build main "cluster" tensor
    ret[0] = std::make_shared<SimpleTensor>(md);

    return ret;
}

ICluster::pointer TensorDM::as_cluster(const ITensor::vector& tens,
                                       const IAnodePlane::vector& anodes,
                                       const std::string& datapath)
{
    TensorIndex ti(tens);
    return as_cluster(ti, anodes, datapath);
}
ICluster::pointer TensorDM::as_cluster(const TensorIndex& ti,
                                       const IAnodePlane::vector& anodes,
                                       const std::string& datapath)
{
    ITensor::pointer top = ti.at(datapath, "cluster");

    node_array_set_t nas;
    edge_array_set_t eas;

    auto topmd = top->metadata();
    const int ident = topmd["ident"].asInt();

    {                           // copy over node arrays
        auto paths = topmd["nodes"];
        for (const auto& code : paths.getMemberNames()) {
            const std::string dpath = paths[code].asString();
            auto ten = ti.at(dpath); // fixme: datatype
            const char nc = code[0];
            const double* data = reinterpret_cast<const double*>(ten->data());
            nas.emplace(nc, boost::const_multi_array_ref<double,2>(data, ten->shape()));
        }
    }

    {                           // copy over edge arrays
        auto paths = topmd["edges"];
        for (const auto& code : paths.getMemberNames()) {
            const std::string dpath = paths[code].asString();
            auto ten = ti.at(dpath); // fixme: datatype
            auto ec = edge_code(code[0], code[1]);
            const int* data = reinterpret_cast<const int*>(ten->data());
            eas.emplace(ec, boost::const_multi_array_ref<int,2>(data, ten->shape()));
        }

    }

    auto cgraph = to_cluster(nas, eas, anodes);
    return std::make_shared<SimpleCluster>(std::move(cgraph), ident);
    
}
