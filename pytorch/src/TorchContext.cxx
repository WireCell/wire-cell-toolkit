#include "WireCellPytorch/TorchContext.h"
#include "WireCellUtil/NamedFactory.h"

using namespace WireCell;
using namespace WireCell::Pytorch;

TorchContext::TorchContext() {}
TorchContext::~TorchContext() { }
TorchContext::TorchContext(const std::string& devname,
                           const std::string& semname)
{
    connect(devname, semname);
}
void TorchContext::connect(const std::string& devname,
                           const std::string& semname)
{
    // Use almost 1/2 the memory and 3/4 the time.
    // fixme: this does effectively nothing as it asserts only inside this code scope.
    torch::NoGradGuard no_grad;

    if (devname == "cpu") {
        m_dev = torch::Device(torch::kCPU);
    }
    else {
        int devnum = 0;
        if (devname.size() > 3) {
            devnum = atoi(devname.substr(3).c_str());
        }
        m_dev = torch::Device(torch::kCUDA, devnum);
    }
    m_devname = devname;

    std::string s_tn = "Semaphore:torch-" + devname;
    if (not semname.empty()) {
        s_tn = semname;
    }

    m_sem = Factory::lookup_tn<ISemaphore>(s_tn);
}

