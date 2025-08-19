#include "WireCellSpng/TorchTensorSetTagger.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"
#include "WireCellSpng/ITorchSpectrum.h"

WIRECELL_FACTORY(SPNGTorchTensorSetTagger, WireCell::SPNG::TorchTensorSetTagger,
                 WireCell::INamed,
                 WireCell::ITorchTensorSetFilter, WireCell::IConfigurable)

WireCell::SPNG::TorchTensorSetTagger::TorchTensorSetTagger()
  : Aux::Logger("SPNGTorchTensorSetTagger", "spng") {

}

WireCell::SPNG::TorchTensorSetTagger::~TorchTensorSetTagger() {};


void WireCell::SPNG::TorchTensorSetTagger::configure(const WireCell::Configuration& config) {

    if (config.isMember("tag_list")) {
        // m_tag_list = 
        auto tag_names = config["tag_list"].getMemberNames();
        for (auto name : tag_names) {
            m_tag_list[name] = config["tag_list"][name].asString();
            log->debug("Will tag TorchTensorSet with {}: {}", name, m_tag_list[name]);
            std::cout << name << " " << m_tag_list[name] << std::endl;
        }
    }

    m_allow_retagging = get(config, "allow_retagging", m_allow_retagging);
    if (m_allow_retagging) {
        log->debug("Will retag");
    }
}

bool WireCell::SPNG::TorchTensorSetTagger::operator()(const input_pointer& in, output_pointer& out) {
    out = nullptr;
    if (!in) {
        log->debug("EOS ");
        return true;
    }
    log->debug("Running TorchTensorSetTagger");

    // out = in;
    auto metadata = in->metadata();
    for (const auto & [key, val] : m_tag_list) {
        if (metadata.isMember(key) && !m_allow_retagging) {
            THROW(ValueError()
                  << errmsg{"Attempting to overwrite tag " + key});
        }
        metadata[key] = val;
    }
    //Clone the tensor to take ownership of the memory and put into 
    //output 
    // std::vector<ITorchTensor::pointer> itv;
    // std::cout << "Making vec" << std::endl;
    // for (auto torch_tensor : *(in->tensors())) {
    //     std::cout << "Pushing back " << std::endl;
    //     itv.push_back(
    //         std::make_shared<SimpleTorchTensor>(torch_tensor->tensor().clone())
    //     );
    // };
    out = std::make_shared<SimpleTorchTensorSet>(
        in->ident(), metadata,
        in->tensors()
        // std::make_shared<std::vector<ITorchTensor::pointer>>(itv)
    );

    return true;
}