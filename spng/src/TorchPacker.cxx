#include "WireCellSpng/TorchPacker.h"
// #include "WireCellAux/TensUtil.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"

WIRECELL_FACTORY(TorchPacker, WireCell::SPNG::TorchPacker, WireCell::ITorchPacker, WireCell::IConfigurable)

using namespace WireCell;

SPNG::TorchPacker::TorchPacker(size_t multiplicity)
  : m_multiplicity(multiplicity)
  , log(Log::logger("sig"))
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
    m_cfg = cfg;
    auto m = get<int>(cfg, "multiplicity", m_multiplicity);
    if (m <= 0) {
        THROW(ValueError() << errmsg{"TorchPacker multiplicity must be positive"});
    }
    m_multiplicity = m;
}

std::vector<std::string> SPNG::TorchPacker::input_types()
{
    const std::string tname = std::string(typeid(input_type).name());
    std::vector<std::string> ret(m_multiplicity, tname);
    return ret;
}

bool SPNG::TorchPacker::operator()(const input_vector& invec, output_pointer& out)
{
    out = nullptr;
    size_t neos = 0;
    for (const auto& fr : invec) {
        if (!fr) {
            ++neos;
        }
    }
    if (neos == invec.size()) {
        return true;
    }
    if (neos) {
        std::cerr << "SPNG::TorchPacker: " << neos << " input tensors missing\n";
    }

    if (invec.size() != m_multiplicity) {
        std::cerr << "SPNG::TorchPacker: got unexpected multiplicity, got:" << invec.size()
                  << " want:" << m_multiplicity << std::endl;
        THROW(ValueError() << errmsg{"unexpected multiplicity"});
    }

    ITorchTensor::vector* itv = new ITorchTensor::vector;
    for (auto iten : invec) {
        itv->push_back(iten);
        log->trace("tag: {}, type: {}", iten->metadata()["tag"], iten->metadata()["type"]);
    }

    // TODO: set md and ident
    Configuration set_md;
    out = std::make_shared<SimpleTorchTensorSet>(0, set_md, ITorchTensor::shared_vector(itv));

    return true;
}
