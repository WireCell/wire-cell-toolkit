#include "WireCellClus/PointTreeMerging.h"
#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/NamedFactory.h"

#include "WireCellAux/TensorDMpointtree.h"
#include "WireCellAux/TensorDMcommon.h"

WIRECELL_FACTORY(PointTreeMerging, WireCell::Clus::PointTreeMerging,
                 WireCell::INamed,
                 WireCell::ITensorSetFanin,
                 WireCell::IConfigurable)

using namespace WireCell;


Clus::PointTreeMerging::PointTreeMerging()
    : Aux::Logger("PointTreeMerging", "clus")
{
}

std::vector<std::string> Clus::PointTreeMerging::input_types()
{
    const std::string tname = std::string(typeid(input_type).name());
    std::vector<std::string> ret(m_multiplicity, tname);
    return ret;
}

void Clus::PointTreeMerging::configure(const WireCell::Configuration& cfg)
{
}

WireCell::Configuration Clus::PointTreeMerging::default_configuration() const
{
    Configuration cfg;
    return cfg;
}

void Clus::PointTreeMerging::finalize()
{
}

bool Clus::PointTreeMerging::operator()(const input_vector& invec, output_pointer& out)
{
    out = nullptr;
    // check input size
    if (invec.size() != m_multiplicity) {
        raise<ValueError>("unexpected multiplicity got %d want %d", invec.size(), m_multiplicity);
        return true;
    }

    // boilerplate for EOS handling
    size_t neos = 0;
    for (const auto& in : invec) {
        if (!in) { ++neos; }
    }
    if (neos == invec.size()) {
        // all inputs are EOS, good.
        log->debug("EOS at call {}", m_count++);
        return true;
    }
    if (neos) { raise<ValueError>("missing %d input tensors ", neos); }

    // TODO: actual impl.
    out = invec[0];
    return true;
}