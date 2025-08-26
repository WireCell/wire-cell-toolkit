#include "WireCellUtil/NamedFactory.h"
#include "WireCellIface/IConfigurable.h"
//#include "WireCellAux/Logger.h"
#include "WireCellSpng/Torch.h"  // One-stop header.
#include "WireCellSpng/SimpleTorchTensorSet.h"
#include "WireCellSpng/ITorchTensorSet.h"
#include "WireCellSpng/ITorchForward.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellSpng/Util.h"
#include <chrono>
#include <torch/script.h>
#include <iostream>
#include <filesystem>
#include <signal.h>
#include <fstream>


// Signal handler to help debug crashes
void crash_handler(int sig) {
    std::cerr << "Segmentation fault caught! Signal: " << sig << std::endl;
    exit(1);
}

//TODO: Add a dummy abstract class so that TestInterface Inherits from it.

class TestBase {
public:
    virtual ~TestBase() = default;
    virtual void dummyMethod() = 0;
};

// Just some test interface...
//If the TestInterface does not inherit anything from ITorchForward or IConfigurable, no error during torch::jit::load...
class TestInterface:public WireCell::SPNG::ITorchForward,
                    public WireCell::IConfigurable,
                    public TestBase
{
public:
    TestInterface();
    virtual ~TestInterface();
    //~TestInterface();
    virtual WireCell::ITorchTensorSet::pointer forward(const WireCell::ITorchTensorSet::pointer& input) const;
    virtual void configure(const WireCell::Configuration& cfg);
    virtual void dummyMethod() override;

private:
    mutable std::unique_ptr<torch::jit::script::Module> m_module;
    //mutable torch::jit::script::Module m_module;
    torch::Device m_device;
    bool m_configured = false;

};


TestInterface::TestInterface():
    m_device(torch::kCPU)
{}


TestInterface::~TestInterface()
{
}

void TestInterface::dummyMethod() {
    std::cout << "Dummy method called" << std::endl;
}

WireCell::ITorchTensorSet::pointer TestInterface::forward(const WireCell::ITorchTensorSet::pointer& input) const
{
    // Perform the forward pass using the TorchScript model
    auto ret = std::make_shared<WireCell::SPNG::SimpleTorchTensorSet>(0);
    return ret;
}


void TestInterface::configure(const WireCell::Configuration& cfg)
{
    std::cout << "TestInterface::configure called" << std::endl;
    
    try {
        if (!cfg.isMember("model")) {
            throw std::runtime_error("No 'model' key in configuration");
        }
        
        std::string model_path = cfg["model"].asString();
        std::cout << "Model path: " << model_path << std::endl;
        
        // Check file exists and is readable
        if (!std::filesystem::exists(model_path)) {
            throw std::runtime_error("Model file does not exist: " + model_path);
        }
        
        std::cout << "File exists, checking size..." << std::endl;
        auto file_size = std::filesystem::file_size(model_path);
        std::cout << "File size: " << file_size << " bytes" << std::endl;
        
        if (file_size == 0) {
            throw std::runtime_error("Model file is empty: " + model_path);
        }
        
        std::cout << "About to call torch::jit::load..." << std::endl;
        std::cout.flush();  // Force output before potential crash
        
        // Load model with explicit error handling
            // Create module on heap
    m_module = std::make_unique<torch::jit::script::Module>(
        torch::jit::load(model_path, m_device)
    );

    std::cout << "torch::jit::load completed successfully" << std::endl;

    // Set to evaluation mode
    //m_module.eval();
        std::cout << "Model set to evaluation mode" << std::endl;
        
        m_configured = true;
        std::cout << "Configuration completed successfully" << std::endl;
        
    } catch (const c10::Error& e) {
        std::cerr << "PyTorch C10 error in configure: " << e.what() << std::endl;
        throw;
    } catch (const std::exception& e) {
        std::cerr << "Standard exception in configure: " << e.what() << std::endl;
        throw;
    } catch (...) {
        std::cerr << "Unknown exception in configure" << std::endl;
        throw;
    }
}

int main(int argc, const char* argv[])
{
    signal(SIGSEGV, crash_handler);
    // Create an instance of TestInterface using WireCell::Factory
    auto test_interface = std::make_shared<TestInterface>();
    WireCell::Configuration cfg;
    const std::string model_path = argv[1];
    //const std::string model_path="/nfs/data/1/abashyal/spng/spng_dev_050525/Pytorch-UNet/ts-model-2.3/unet-l23-cosmic500-e50.ts";
    cfg["model"] = model_path;
    //check if the path exists
    
    if (!std::filesystem::exists(cfg["model"].asString())) {
        std::cerr << "Model path does not exist: " << cfg["model"].asString() << std::endl;
        return 1;
    }
    
    std::cout<<"TestInterface: Configuring with model path: "<<model_path<<std::endl;
    std::cout.flush(); //output before the crash..
    //also test torch::jit::stand alone...
    //torch::jit::script::Module module;
    /*
    if(torch::cuda::is_available()) {
        try {
            torch::Device device = torch::Device(torch::kCUDA, 0);
            auto module = torch::jit::load(model_path, device);
            std::cout << "Standalone torch::jit::load succeeded." << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Error loading model: " << e.what() << std::endl;
            return 1;
        }
    }
    */
    //now test with the interface...
     test_interface->configure(cfg);

return 0;

}