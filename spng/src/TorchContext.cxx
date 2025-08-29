#include "WireCellSpng/TorchContext.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"

namespace WireCell::SPNG {


    TorchContext::TorchContext()
    {
        connect();
    }

    void TorchContext::connect(const std::string& devname, const std::string& semname)
    {
        if (devname == "cpu") {
            m_dev = torch::Device(torch::kCPU);
        } // all rest are cuda
        else if (torch::cuda::is_available()) {

            if (devname == "gpu" || devname == "cuda") {
                m_dev = torch::Device(torch::kCUDA);
            }
            else if (String::startswith(devname, "gpu")) {
                int devnum = atoi(devname.substr(3).c_str());
                m_dev = torch::Device(torch::kCUDA, devnum);
            }
            else {              // assume integer
                int devnum = atoi(devname.c_str());
                m_dev = torch::Device(torch::kCUDA, devnum);
            }
        }
        else {
            raise<RuntimeError>("No torch device matches: \"%s\"", devname);
        }
        m_devname = devname;
            

        auto m_semname = semname;
        if (m_semname.empty()) {
            m_semname = "Semaphore:torch-" + m_devname;
        }
        m_sem = Factory::lookup_tn<ISemaphore>(m_semname);
    }

    void TorchContext::enter() const { 
        if (m_sem) m_sem->acquire(); 

        m_agm = torch::GradMode::is_enabled();
        torch::GradMode::set_enabled(false);
    }
    void TorchContext::exit() const {
        torch::GradMode::set_enabled(m_agm);

        if (m_sem) m_sem->release();
    }

}
