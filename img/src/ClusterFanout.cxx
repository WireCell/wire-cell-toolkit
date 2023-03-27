#include "WireCellImg/ClusterFanout.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellAux/SimpleCluster.h"

WIRECELL_FACTORY(ClusterFanout, WireCell::Img::ClusterFanout,
                 WireCell::INamed,
                 WireCell::IClusterFanout, WireCell::IConfigurable)

using namespace WireCell;
using WireCell::Aux::SimpleCluster;

Img::ClusterFanout::ClusterFanout(size_t multiplicity)
    : Aux::Logger("ClusterFanout", "glue")
    , m_multiplicity(multiplicity)
{
}
Img::ClusterFanout::~ClusterFanout() {}

WireCell::Configuration Img::ClusterFanout::default_configuration() const
{
    Configuration cfg;
    // How many output ports
    cfg["multiplicity"] = (int) m_multiplicity;
    return cfg;
}
void Img::ClusterFanout::configure(const WireCell::Configuration& cfg)
{
    int m = get<int>(cfg, "multiplicity", (int) m_multiplicity);
    if (m <= 0) {
        THROW(ValueError() << errmsg{"multiplicity must be positive"});
    }
    m_multiplicity = m;
}

std::vector<std::string> Img::ClusterFanout::output_types()
{
    const std::string tname = std::string(typeid(output_type).name());
    std::vector<std::string> ret(m_multiplicity, tname);
    return ret;
}

bool Img::ClusterFanout::operator()(const input_pointer& in, output_vector& outv)
{
    outv.resize(m_multiplicity);
    if (!in) {  //  pass on EOS
        log->debug("EOS at call={}", m_count);
    }
    for (size_t ind = 0; ind < m_multiplicity; ++ind) {
        outv[ind] = in;
    }
    ++m_count;
    return true;
}
