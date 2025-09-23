#include "WireCellSpng/TorchService.h"
#include "WireCellSpng/Util.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Persist.h"
#include "WireCellSpng/TorchScript.h"

//do we rely on openMP...probably not..
//#include <omp.h>
#include <ATen/Parallel.h>

WIRECELL_FACTORY(SPNGTorchService, 
                 WireCell::SPNG::TorchService,
                 WireCell::SPNG::ITorchForward, // interface 1
                 WireCell::IConfigurable)

using namespace WireCell;



SPNG::TorchService::TorchService()
    : Aux::Logger("SPNGTorchService", "spng")
{

}

Configuration SPNG::TorchService::default_configuration() const
{
    Configuration cfg;

    // TorchScript model
    cfg["model"] = "model.ts";

    // one of: {cpu, gpu, gpuN} where "N" is a GPU number.  "gpu"
    // alone will use GPU 0.
    cfg["device"] = "cpu";
    
    return cfg;
}

void SPNG::TorchService::configure(const WireCell::Configuration& cfg)
{
    // Device setup
    auto devname = get<std::string>(cfg, "device", "cpu");
    log->debug("TorchService device configuration: {}", devname);
    
    if(devname == "cpu"){
        m_device = torch::Device(torch::kCPU);
    }
    else if(devname == "gpu" || devname == "gpu0") {
        m_device = torch::Device(torch::kCUDA, 0);
    }
    else if(devname.substr(0, 3) == "gpu" && devname.size() > 3) {
        int devnum = 0;
        try {
            devnum = std::stoi(devname.substr(3));
        } catch (const std::exception& e) {
            log->warn("Failed to parse GPU device number from '{}', using GPU 0", devname);
            devnum = 0;
        }
        m_device = torch::Device(torch::kCUDA, devnum);
    }
    else {
        log->warn("Unknown device type '{}', defaulting to CPU", devname);
        m_device = torch::Device(torch::kCPU);
    }

    log->debug("TorchService configured for device: {}", m_device.str());
    
    const std::string model_path = Persist::resolve(cfg["model"].asString());
    
    if (model_path.empty()) {
        log->critical("no TorchScript model file provided");
        THROW(ValueError() << errmsg{"no TorchScript model file provided"});
    }

    // Check CUDA availability
    if(!torch::cuda::is_available() && m_device.is_cuda()) {
        log->critical("CUDA is not available but GPU is requested. Please check your configuration.");
        THROW(ValueError() << errmsg{"CUDA is not available but GPU is requested."});
    }

    // Use NoGradGuard for inference mode
    torch::NoGradGuard no_grad;
    
    try {
        log->debug("Loading TorchScript model from: {} to device {}", model_path, m_device.str());
        

        m_module = torch::jit::load(model_path, m_device);
        
        log->debug("Loaded TorchScript model from: {} to device: {}", model_path, m_device.str()); 
    }
    catch (const c10::Error& e) {
        log->critical("error loading model: \"{}\" to device \"{}\": {}",
                      model_path, m_device.str(), e.what());
        throw;              
    }
    catch (const std::exception& e) {
        log->critical("error loading model: \"{}\" to device \"{}\": {}",
                      model_path, m_device.str(), e.what());
        THROW(ValueError() << errmsg{"error loading TorchScript model"} <<
                            errmsg{" " + std::string(e.what())});
    }
}


ITorchTensorSet::pointer SPNG::TorchService::forward(const ITorchTensorSet::pointer& in) const
{

    log->debug("TorchService::forward function entered");


    if (!in) {
        log->critical("TorchService::forward received a null input pointer");
        THROW(ValueError() << errmsg{"TorchService::forward received a null input pointer"});
    }

    try {
        log->debug("TorchService::forward called with input: {}", in->ident());
    }
    catch (const std::exception& e) {
        log->critical("Exception while accessing in->ident(): {}", e.what());
        THROW(ValueError() << errmsg{"Exception while accessing in->ident()"} <<
                            errmsg{" " + std::string(e.what())});
    }


    
    std::vector<torch::IValue> iival = SPNG::from_itensor(in, m_device.is_cuda());

    torch::IValue oival;

    // Use NoGradGuard for inference (saves memory and improves performance)
    torch::NoGradGuard no_grad;

    try {
        log->debug("TorchService::forward running model with {} inputs", iival.size());
        oival = m_module.forward(iival);
        log->debug("TorchService::forward model execution completed successfully");
    }
    catch (const c10::Error& err) {
        log->critical("PyTorch C10 error running model on device {}: {}", m_device.str(), err.what());
        THROW(ValueError() << errmsg{"PyTorch C10 error running model on device"} <<
                            errmsg{" " + std::string(err.what())});
    }

    ITorchTensorSet::pointer ret = SPNG::to_itensor({oival});
    return ret;
}
