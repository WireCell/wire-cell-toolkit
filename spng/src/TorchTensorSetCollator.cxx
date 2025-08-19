#include "WireCellSpng/TorchTensorSetCollator.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellAux/FrameTools.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellIface/INamed.h"
#include <cuda.h>

WIRECELL_FACTORY(TorchTensorSetCollator, WireCell::SPNG::TorchTensorSetCollator,
                 WireCell::INamed,
                 WireCell::SPNG::ITorchTensorSetFanin)


using namespace WireCell;

WireCell::Configuration SPNG::TorchTensorSetCollator::default_configuration() const
{
    Configuration cfg;
    cfg["multiplicity"] = m_multiplicity;
    cfg["output_set_tag"] = m_output_set_tag;


    return cfg;
}

void SPNG::TorchTensorSetCollator::configure(const WireCell::Configuration& config)
{
    //Get the anode to make a channel map for output
    m_multiplicity = get(config, "multiplicity", m_multiplicity);
    m_output_set_tag = get(config, "output_set_tag", m_output_set_tag.asString());
    log->debug("Will tag output set with {}", m_output_set_tag);
}

std::vector<std::string> SPNG::TorchTensorSetCollator::input_types()
{
    if (!m_multiplicity) {
        return std::vector<std::string>();
    }
    const std::string tname = std::string(typeid(input_type).name());
    std::vector<std::string> ret(m_multiplicity, tname);
    return ret;
}


SPNG::TorchTensorSetCollator::TorchTensorSetCollator()
    : Aux::Logger("TorchTensorSetCollator", "spng") {}

bool SPNG::TorchTensorSetCollator::operator()(const input_vector& inv, output_pointer& out) {
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


    if (m_multiplicity > 0 && (inv.size() != m_multiplicity)) {
        THROW(ValueError()
            << errmsg{"Expected input multiplicity ({" + std::to_string(inv.size()) +
            "}) to match configuration ({" + std::to_string(m_multiplicity) + "})"});
    }

    std::vector<ITorchTensor::pointer> itv;
    for (auto in : inv) {
        auto this_set_tag = in->metadata()["tag"];
        for (auto tensor : (*in->tensors())) {
            auto input_metadata = tensor->metadata();
            auto this_tensor_tag = input_metadata["tag"];
            std::string combined_tag = this_set_tag.asString() + ":" + this_tensor_tag.asString();
            log->debug("Adding tensor with tag {}", combined_tag);
            
            Configuration tensor_md;
            if (input_metadata.isMember("channel_map")) {
                tensor_md["channel_map"] = input_metadata["channel_map"];
                std::cout << "Storing channel map" << std::endl;
            }

            tensor_md["tag"] = combined_tag;
            itv.push_back(
                std::make_shared<SimpleTorchTensor>(tensor->tensor(), tensor_md)
            );
        }
    }

    Configuration output_md;
    output_md["tag"] = m_output_set_tag;
    log->debug("Writing out Tensor Set with output set tag {} {}", output_md["tag"], m_output_set_tag);
    out = std::make_shared<SimpleTorchTensorSet>(
        inv.at(0)->ident(), output_md,
        std::make_shared<std::vector<ITorchTensor::pointer>>(itv)
    );
    


    return true;
}
