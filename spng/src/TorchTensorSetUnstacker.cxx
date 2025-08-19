#include "WireCellSpng/TorchTensorSetUnstacker.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellAux/FrameTools.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellIface/INamed.h"
#include <cuda.h>

WIRECELL_FACTORY(TorchTensorSetUnstacker, WireCell::SPNG::TorchTensorSetUnstacker,
                 WireCell::INamed,
                 WireCell::SPNG::ITorchTensorSetFanout)


using namespace WireCell;

WireCell::Configuration SPNG::TorchTensorSetUnstacker::default_configuration() const
{
    Configuration cfg;
    cfg["multiplicity"] = m_multiplicity;

    return cfg;
}

void SPNG::TorchTensorSetUnstacker::configure(const WireCell::Configuration& config)
{
    m_multiplicity = get(config, "multiplicity", m_multiplicity);

}

std::vector<std::string> SPNG::TorchTensorSetUnstacker::output_types()
{
    const std::string tname = std::string(typeid(ITorchTensorSet).name());
    log->debug("Got {}", m_multiplicity);
    std::vector<std::string> ret(m_multiplicity, tname);
    return ret;
}


SPNG::TorchTensorSetUnstacker::TorchTensorSetUnstacker()
    : Aux::Logger("TorchTensorSetUnstacker", "spng") {}

bool SPNG::TorchTensorSetUnstacker::operator()(const input_pointer& in, output_vector& outv) {
    outv.resize(m_multiplicity);
    //Default null ptrs
    for (size_t ind = 0; ind < m_multiplicity; ++ind) {
        outv[ind] = nullptr;
    }

    if (!in) {
        // log->debug("EOS ");
        return true;
    }

    //TODO -- implement more sophisticated access?
    std::vector<torch::Tensor> output_tensors = torch::split(
        (*in->tensors())[0]->tensor(), 1, 0
    );
    //Make each output tensor a clone since torch::split returns a view
    for (auto & tensor : output_tensors) {
        tensor = tensor.clone();
    }
    // if (m_multiplicity != output_tensors.size()) {
    // THROW
    // }

    Configuration tensor_md, set_md;
    
    // std::cout << tensor_md["channel_map"] << std::endl;

    for (int output_index = 0; output_index < outv.size(); ++output_index) {
        //Clone the tensor to take ownership of the memory and put into 
        //output 
        tensor_md["channel_map"] = in->tensors()->at(0)->metadata()["channel_map"][output_index];
        std::vector<ITorchTensor::pointer> itv{
            std::make_shared<SimpleTorchTensor>(
                output_tensors[output_index],
                tensor_md
            ) //.clone())
        };
        outv[output_index] = std::make_shared<SimpleTorchTensorSet>(
            in->ident(), set_md,
            std::make_shared<std::vector<ITorchTensor::pointer>>(itv)
        );
    }
    return true;
}
