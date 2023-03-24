#include "WireCellImg/ClusterFanin.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellAux/SimpleCluster.h"
#include "WireCellUtil/NamedFactory.h"
WIRECELL_FACTORY(ClusterFanin, WireCell::Img::ClusterFanin,
                 WireCell::INamed,
                 WireCell::IClusterFanin, WireCell::IConfigurable)

using namespace WireCell;

using WireCell::Aux::SimpleCluster;

Img::ClusterFanin::ClusterFanin()
    : Aux::Logger("ClusterFanin", "glue")
{
}
Img::ClusterFanin::~ClusterFanin()
{
}


// IConfigurable
WireCell::Configuration Img::ClusterFanin::default_configuration() const
{
    Configuration cfg;
    cfg["multiplicity"] = (int) m_multiplicity;
    cfg["ident_source"] = (int) m_ident_source;
    return cfg;
}

void Img::ClusterFanin::configure(const WireCell::Configuration& cfg)
{
    const int mult = get<int>(cfg, "multiplicity", (int) m_multiplicity);
    if (mult <= 0) {
        log->critical("illegal multiplicity={}", mult);
        THROW(ValueError() << errmsg{"multiplicity must be positive"});
    }
    m_multiplicity = mult;

    int isrc = get<int>(cfg, "ident_source", (int) m_ident_source);
    if (isrc < 0 or isrc >= mult) {
        log->critical("given illegal ident_source={} with multiplicity={}", isrc, mult);
        THROW(ValueError() << errmsg{"illegal ident source"});
    }
    m_ident_source = isrc;
}

// IFanin
bool Img::ClusterFanin::operator()(const input_vector& invec, output_pointer& out)
{
    out = nullptr;
    size_t neos = 0;
    for (const auto& fr : invec) {
        if (!fr) {
            ++neos;
        }
    }
    if (neos) {
        log->debug("EOS at call={} with {}", m_count, neos);
        ++m_count;
        return true;
    }
    if (invec.size() != m_multiplicity) {
        log->critical("input vector size={} my multiplicity={}", invec.size(), m_multiplicity);
        THROW(ValueError() << errmsg{"input vector size mismatch"});
    }

    cluster_graph_t res;
    for (const auto& ic : invec) {
        const auto& og = ic->graph();

        boost::copy_graph(og, res); // need some converter function?
    }

    const int ident = invec[m_ident_source]->ident();

    out = std::make_shared<SimpleCluster>(res, ident);
        
    const size_t nnodes = boost::num_vertices(res);
    const size_t nedges = boost::num_edges(res);
    log->debug("fan {} clusters to ident={} with Nv={}, Ne={} at call={}",
               m_multiplicity, ident, nnodes, nedges, m_count);

    ++m_count;
    return true;
}


// INode, override because we get multiplicity at configure time.
std::vector<std::string> Img::ClusterFanin::input_types()
{
    const std::string tname = std::string(typeid(input_type).name());
    std::vector<std::string> ret(m_multiplicity, tname);
    return ret;
}
