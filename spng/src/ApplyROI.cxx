#include "WireCellSpng/ApplyROI.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"
#include "WireCellSpng/ITorchSpectrum.h"
#include "WireCellSpng/Util.h"

WIRECELL_FACTORY(SPNGApplyROI, WireCell::SPNG::ApplyROI,
                 WireCell::INamed,
                 WireCell::SPNG::ITorchTensorSetFanin, WireCell::IConfigurable)

WireCell::SPNG::ApplyROI::ApplyROI()
  : Aux::Logger("SPNGApplyROI", "spng") {

}

WireCell::SPNG::ApplyROI::~ApplyROI() {};


void WireCell::SPNG::ApplyROI::configure(const WireCell::Configuration& config) {

    m_ROI_tensor_index = get(config, "ROI_tensor_index", m_ROI_tensor_index);
    m_value_tensor_index = get(config, "value_tensor_index", m_value_tensor_index);

    log->debug("Will apply ROI_tensor {} to value_tensor {}", m_ROI_tensor_index, m_value_tensor_index);

    m_output_tensor_tag = get(config, "output_tensor_tag", m_output_tensor_tag);

    if (!config.isMember("output_set_tag")) {
        THROW(ValueError()
            << errmsg{"Must provide output_set_tag in ApplyROI configuration"});
    }
    m_output_set_tag = config["output_set_tag"];
}
std::vector<std::string> WireCell::SPNG::ApplyROI::input_types() {
    if (!m_multiplicity) {
        return std::vector<std::string>();
    }
    const std::string tname = std::string(typeid(input_type).name());
    std::vector<std::string> ret(m_multiplicity, tname);
    return ret;
}

bool WireCell::SPNG::ApplyROI::operator()(const input_vector& inv, output_pointer& out) {
    out = nullptr;
    size_t neos = 0;
    for (const auto& in : inv) {
        if (!in) {
            ++neos;
        }
    }
    if (neos) {
        log->debug("EOS with {}", neos);
        return true;
    }


    if (m_multiplicity > 0 && (inv.size() != (size_t)m_multiplicity)) {
        THROW(ValueError()
            << errmsg{"Expected input multiplicity ({" + std::to_string(inv.size()) +
            "}) to match configuration ({" + std::to_string(m_multiplicity) + "})"});
    }
    log->debug("Running ApplyROI");

    auto ROI_tensor_clone = (*inv[m_ROI_tensor_index]->tensors())[0]->tensor().clone();
    auto value_tensor_clone = (*inv[m_value_tensor_index]->tensors())[0]->tensor().clone();

    //Mask out the value
    value_tensor_clone = value_tensor_clone * ROI_tensor_clone;

    //FFT on requested dim
    // TODO: set md?
    Configuration set_md, tensor_md;
    set_md["tag"] = m_output_set_tag;
    tensor_md["tag"] = m_output_tensor_tag;
    std::vector<ITorchTensor::pointer> itv{
        std::make_shared<SimpleTorchTensor>(value_tensor_clone, tensor_md)
    };
    out = std::make_shared<SimpleTorchTensorSet>(
        0, set_md,
        std::make_shared<std::vector<ITorchTensor::pointer>>(itv)
    );

    return true;
}