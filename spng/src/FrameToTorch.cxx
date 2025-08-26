#include "WireCellSpng/FrameToTorch.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellAux/FrameTools.h"
#include "WireCellAux/SimpleFrame.h"

WIRECELL_FACTORY(FrameToTorch, WireCell::SPNG::FrameToTorch,
                 WireCell::INamed,
                 WireCell::SPNG::IFrameToTorch)


using namespace WireCell;

WireCell::Configuration SPNG::FrameToTorch::default_configuration() const
{
    Configuration cfg;
    cfg["Test"] = true;
    return cfg;
}

bool SPNG::FrameToTorch::operator()(const input_pointer& in, output_pointer& out) {
    out = nullptr;
    if (!in) return true;
    const long ntraces = in->traces()->size();
    
    log->debug("Ntraces: {}", ntraces);

    if (ntraces == 0) {
        log->debug("No traces, exiting");
        return true;
    }

    
    const long nticks = (*in->traces())[0]->charge().size();
    log->debug("Making tensor of shape: {} {}", ntraces, nticks);
    torch::Tensor output_tensor = torch::zeros({ntraces, nticks});
    log->debug("Putting");
    for (size_t i = 0; i < ntraces; ++i) {
        for (size_t j = 0; j < nticks; ++j) {
            output_tensor.index_put_({(int64_t)i, (int64_t)j}, (*in->traces())[i]->charge()[j]);
        }    
    }
    log->debug("Put");


    out = SimpleTorchTensor::pointer(
        new SimpleTorchTensor(output_tensor.clone()));

    return true;
}
