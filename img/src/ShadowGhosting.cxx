#include "WireCellImg/ShadowGhosting.h"
#include "WireCellAux/ClusterShadow.h"

#include "WireCellUtil/NamedFactory.h"


WIRECELL_FACTORY(ShadowGhosting, WireCell::Img::ShadowGhosting,
                 WireCell::INamed,
                 WireCell::IClusterFilter, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Img;
using namespace WireCell::Aux;


ShadowGhosting::ShadowGhosting()
    : Aux::Logger("ShadowGhosting", "img")
{
}

ShadowGhosting::~ShadowGhosting()
{
}

void ShadowGhosting::configure(const WireCell::Configuration& cfg)
{
    m_shadow_type = get<std::string>(cfg, "shadow_type", m_shadow_type);
}

WireCell::Configuration ShadowGhosting::default_configuration() const
{
    WireCell::Configuration cfg;
    cfg["shadow_type"] = m_shadow_type;
    return cfg;
}

bool ShadowGhosting::operator()(const input_pointer& in, output_pointer& out)
{
    // out = nullptr;
    out = in;                   // pass-through for now
    if (!in) {
        log->debug("EOS at call={}", m_count);
        ++m_count;
        return true;
    }
    
    const auto& cgraph = in->graph();

    auto bsg = BlobShadow::shadow(cgraph, m_shadow_type[0]);
    log->debug("nblobs={}", boost::num_vertices(bsg));

    ClusterShadow::blob_cluster_map_t clusters;
    auto csg = ClusterShadow::shadow(cgraph, bsg, clusters);

    log->debug("nclusters={}", boost::num_vertices(csg));
    ++m_count;
    return true;
}

