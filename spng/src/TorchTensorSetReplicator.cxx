#include "WireCellSpng/TorchTensorSetReplicator.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellAux/FrameTools.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellIface/INamed.h"
#include <cuda.h>

WIRECELL_FACTORY(TorchTensorSetReplicator, WireCell::SPNG::TorchTensorSetReplicator,
                 WireCell::INamed,
                 WireCell::SPNG::ITorchTensorSetFanout)


using namespace WireCell;

WireCell::Configuration SPNG::TorchTensorSetReplicator::default_configuration() const
{
    Configuration cfg;
    cfg["multiplicity"] = m_multiplicity;

    return cfg;
}

void SPNG::TorchTensorSetReplicator::configure(const WireCell::Configuration& config)
{
    m_multiplicity = get(config, "multiplicity", m_multiplicity);
}

std::vector<std::string> SPNG::TorchTensorSetReplicator::output_types()
{
    const std::string tname = std::string(typeid(ITorchTensorSet).name());
    log->debug("Got {}", m_multiplicity);
    std::vector<std::string> ret(m_multiplicity, tname);
    return ret;
}


SPNG::TorchTensorSetReplicator::TorchTensorSetReplicator()
    : Aux::Logger("TorchTensorSetReplicator", "spng") {}

bool SPNG::TorchTensorSetReplicator::operator()(const input_pointer& in, output_vector& outv) {
    outv.resize(m_multiplicity);
    //Default null ptrs
    for (size_t ind = 0; ind < m_multiplicity; ++ind) {
        outv[ind] = in;
    }
    return true;

    // //Nothing in, nothing out
    // if (!in) {  //  pass on EOS
    //     log->debug("Exiting");
    //     return true;
    // }

    // for (size_t ind = 0; ind < m_multiplicity; ++ind) {
    //     outv[ind] = in;
    // }

    // return true;
}
