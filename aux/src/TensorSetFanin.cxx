#include "WireCellAux/TensorSetFanin.h"
#include "WireCellAux/TensorDMcommon.h"
#include "WireCellAux/SimpleTensorSet.h"

#include "WireCellUtil/NamedFactory.h"


WIRECELL_FACTORY(TensorSetFanin, WireCell::Aux::TensorSetFanin,
                 WireCell::INamed,
                 WireCell::IConfigurable,
                 WireCell::ITensorSetFanin)


using namespace WireCell;


Aux::TensorSetFanin::TensorSetFanin()
  : Aux::Logger("TensorSetFanin", "aux")
{
}


Aux::TensorSetFanin::~TensorSetFanin()
{
}


WireCell::Configuration Aux::TensorSetFanin::default_configuration() const
{
    Configuration cfg;
    cfg["ident_port"] = m_ident_port;
    // also: "tensor_order" and "multiplicity"
    return cfg;
}


void Aux::TensorSetFanin::configure(const WireCell::Configuration& cfg)
{
    m_multiplicity = get(cfg, "multiplicity", -1);
    if (m_multiplicity < 0) {
        raise<ValueError>("illegal 'multiplicity' value: %d", m_multiplicity);
    }

    m_ident_port = get(cfg, "ident_port", m_ident_port);
    if (m_ident_port < 0 || m_ident_port >= m_multiplicity) {
        raise<ValueError>("illegal 'ident_port' value: %d", m_ident_port);
    }

    m_tensor_order = get(cfg, "tensor_order", m_tensor_order);
    if (m_tensor_order.size() != (size_t)m_multiplicity) {
        log->warn("configured to ignore some ports, I hope this is what you want");
    }

}


std::vector<std::string> Aux::TensorSetFanin::input_types() 
{
    const std::string tname = std::string(typeid(input_type).name());
    std::vector<std::string> ret(m_multiplicity, tname);
    return ret;
}


bool Aux::TensorSetFanin::operator()(const input_vector& in_tss, output_pointer& out_ts)
{
    out_ts = nullptr;
    if (in_tss.empty()) {
        return true;  // eos
    }
    for (const auto& tsptr : in_tss) {
        if (!tsptr) {           // we assume input is synced so any EOS is EOS
            return true;
        }
    }

    const int ident = in_tss[m_ident_port]->ident();

    // FIXME: strategies to fill metadata from the input are needed.

    ITensor::vector tensors;

    for (int iport : m_tensor_order) {
        ITensorSet::pointer ts = in_tss[iport];
        log->debug("port={} ident={} md={}", iport, ts->ident(), ts->metadata());
        const auto& tens = ts->tensors();
        tensors.insert(tensors.end(), tens->begin(), tens->end());
    }
    out_ts = Aux::TensorDM::as_tensorset(tensors, ident); // no md
    return true;
}

