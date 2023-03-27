#include "WireCellImg/PointGraphProcessor.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellAux/TensorDM.h"
#include "WireCellAux/SimpleTensorSet.h"

WIRECELL_FACTORY(PointGraphProcessor, WireCell::Img::PointGraphProcessor,
                 WireCell::INamed,
                 WireCell::ITensorSetFilter,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Img;
using namespace WireCell::Aux;
using namespace WireCell::Aux::TensorDM;

PointGraphProcessor::PointGraphProcessor()
    : Aux::Logger("PointGraphProcessor", "img")
{
}

void PointGraphProcessor::configure(const WireCell::Configuration& cfg)
{
    m_inpath = get(cfg, "inpath", m_inpath);
    m_outpath = get(cfg, "outpath", m_outpath);
    m_visitors.clear();
    for (const auto& jvis : cfg["visitors"]) {
        m_visitors.push_back(Factory::find_tn<IPointGraphVisitor>(jvis.asString()));
    }
}

WireCell::Configuration PointGraphProcessor::default_configuration() const
{
    Configuration cfg;
    cfg["inpath"] = m_inpath;
    cfg["outpath"] = m_outpath;
    cfg["visitors"] = Json::arrayValue;
    return cfg;
}


bool PointGraphProcessor::operator()(const input_pointer& ints, output_pointer& outts)
{
    outts = nullptr;
    if (!ints) {
        log->debug("EOS at call {}", m_count++);
        return true;
    }

    const auto& tens = *ints->tensors();
    const auto pcgs = match_at(tens, m_inpath);
    if (pcgs.size() != 1) {
        outts = std::make_shared<SimpleTensorSet>(ints->ident());
        log->warn("unexpected number of tensors ({}) matching {} at call={}",
                  pcgs.size(), m_inpath, m_count++);
        return true;
    }

    // Convert tensor set to point graph
    PointGraph pg;

    const auto md0 = pcgs[0]->metadata();
    const std::string datatype = md0["datatype"].asString();
    const std::string datapath = md0["datapath"].asString();

    if (datatype == "pcgraph") {
        pg = as_pointgraph(tens, datapath);
    }
    else if (datatype == "pcdataset") {
        auto nodes = as_dataset(tens, datapath);
        pg = PointGraph(nodes);
    }
    else {
        outts = std::make_shared<SimpleTensorSet>(ints->ident());
        log->warn("no pcgraph at datapath \"{}\" at call={}",
                  datapath, m_count++);
        return true;
    }


    for (auto& visitor : m_visitors) {
        visitor->visit_pointgraph(pg);
    }

    std::string outpath = m_outpath;
    if (outpath.find("%") != std::string::npos) {
        outpath = String::format(outpath, ints->ident());
    }
    auto out_tens = as_tensors(pg, outpath);

    // Convert point graph to tensor set
    auto stv = std::make_shared<ITensor::vector>(out_tens.begin(), out_tens.end());
    outts = std::make_shared<SimpleTensorSet>(ints->ident(), ints->metadata(), stv);

    log->debug("processed point graph at call={}", m_count++);

    return true;
}
