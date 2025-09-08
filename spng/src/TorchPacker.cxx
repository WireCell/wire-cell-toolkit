#include "WireCellSpng/TorchPacker.h"
// #include "WireCellAux/TensUtil.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"

#include <algorithm>

WIRECELL_FACTORY(TorchPacker, WireCell::SPNG::TorchPacker, WireCell::ITorchPacker, WireCell::IConfigurable)

using namespace WireCell;

SPNG::TorchPacker::TorchPacker(size_t multiplicity)
    : Aux::Logger("TorchPacker", "spng")
    , m_multiplicity(multiplicity)
{
}
SPNG::TorchPacker::~TorchPacker() {}

WireCell::Configuration SPNG::TorchPacker::default_configuration() const
{
    Configuration cfg;
    // How many input ports
    cfg["multiplicity"] = (int) m_multiplicity;
    return cfg;
}
void SPNG::TorchPacker::configure(const WireCell::Configuration& cfg)
{
    auto m = get<int>(cfg, "multiplicity", m_multiplicity);
    if (m <= 0) {
        raise<ValueError>("TorchPacker multiplicity must be positive definite");
    }
    m_multiplicity = m;
    m_count = get<int>(cfg, "count", m_count);
}

std::vector<std::string> SPNG::TorchPacker::input_types()
{
    const std::string tname = std::string(typeid(input_type).name());
    std::vector<std::string> ret(m_multiplicity, tname);
    return ret;
}

bool SPNG::TorchPacker::operator()(const input_vector& invec, output_pointer& out)
{
    if (invec.size() != m_multiplicity) {
        log->error("unexpected multiplicity, got:{} want:{}", invec.size(), m_multiplicity);
        raise<ValueError>("unexpected multiplicity, got:%d want:%d", invec.size(), m_multiplicity);
    }

    out = nullptr;
    size_t neos = std::count(invec.begin(), invec.end(), nullptr);
    if (neos) {
        log->debug("EOS in {} of {} at call={}",
                   neos, m_multiplicity, m_count);
        ++m_count;
        return true;
    }

    auto itv = std::make_shared<ITorchTensor::vector>();
    for (auto iten : invec) {
        itv->push_back(iten);
        log->trace("tag: {}, type: {}", iten->metadata()["tag"], iten->metadata()["type"]);
    }

    // FIXME: is there a meaningful way to set md?
    Configuration md;
    out = std::make_shared<SimpleTorchTensorSet>(m_count, md, itv);
    ++m_count;
    return true;
}
