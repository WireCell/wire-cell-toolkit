/// ClusterFlashDump implementation
#include "WireCellAux/ClusterFlashDump.h"
#include "WireCellAux/TensorDMpointtree.h"

#include "WireCellUtil/NamedFactory.h"



WIRECELL_FACTORY(ClusterFlashDump,
                 WireCell::Aux::ClusterFlashDump,
                 WireCell::INamed,
                 WireCell::ITensorSetSink)

using namespace WireCell;
using WireCell::Aux::TensorDM::as_pctree;

Aux::ClusterFlashDump::ClusterFlashDump()
    : Aux::Logger("ClusterFlashDump", "aux")
{
}

Aux::ClusterFlashDump::~ClusterFlashDump()
{
}

void Aux::ClusterFlashDump::configure(const WireCell::Configuration& cfg)
{
    m_datapath = get(cfg, "datapath", m_datapath);
    if (m_datapath.empty()) {
        raise<ValueError>("datapath must be provided");
    }
}

WireCell::Configuration Aux::ClusterFlashDump::default_configuration() const
{
    Configuration cfg;
    cfg["datapath"] = m_datapath;    // required
    return cfg;
}

bool Aux::ClusterFlashDump::operator()(const ITensorSet::pointer& in)
{
    if (!in) {
        return true;            // eos
    }

    const int ident = in->ident();

    std::string datapath = m_datapath;
    if (datapath.find("%") != std::string::npos) {
        datapath = String::format(datapath, ident);
    }

    const auto& tens = *in->tensors();
    auto root = as_pctree(tens, datapath);

    log->debug("got pc-tree ident {} with {} L2 children (clusters) in {} tensors",
               ident, root->nchildren(), tens.size());

    return true;
}
