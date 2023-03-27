#include "WireCellAux/TensorDM.h"
#include "WireCellAux/SimpleTensor.h"


using namespace WireCell;
using namespace WireCell::Aux;
using namespace WireCell::Aux::TensorDM;

WireCell::PointGraph
WireCell::Aux::TensorDM::as_pointgraph(const ITensor::vector& tens,
                                       const std::string& datapath)
{
    ITensor::pointer top = nullptr;
    std::unordered_map<std::string, ITensor::pointer> located;
    for (const auto& iten : tens) {
        const auto& tenmd = iten->metadata();
        const auto dtype = tenmd["datatype"].asString();
        const auto dpath = tenmd["datapath"].asString();

        if (!top and dtype == "pcgraph") {
            if (datapath.empty() or datapath == dpath) {
                top = iten;
            }
            continue;
        }
        if (dtype == "pcdataset") {
            located[dpath] = iten;
        }
        continue;
    }

    if (!top) {
        THROW(ValueError() << errmsg{"no array of datatype \"pcgraph\" at datapath \"" + datapath + "\"" });
    }

    const auto& topmd = top->metadata();
    auto nodes = as_dataset(tens, topmd["nodes"].asString());
    auto edges = as_dataset(tens, topmd["edges"].asString());
    return PointGraph(nodes, edges);
}

ITensor::vector
WireCell::Aux::TensorDM::as_tensors(const PointGraph& pcgraph,
                                    const std::string datapath)
{
    Configuration md;
    md["datatype"] = "pcgraph";
    md["datapath"] = datapath;
    md["nodes"] = datapath + "/nodes";
    md["edges"] = datapath + "/edges";

    ITensor::vector tens;
    tens.push_back(std::make_shared<SimpleTensor>(md));

    auto nodes = as_tensors(pcgraph.nodes(), datapath + "/nodes");
    tens.insert(tens.end(), nodes.begin(), nodes.end());

    auto edges = as_tensors(pcgraph.edges(), datapath + "/edges");
    tens.insert(tens.end(), edges.begin(), edges.end());

    return tens;
}

