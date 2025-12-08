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
            
        
        m_semname = semname;
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

    // A per-C++ thread "stream" is needed to submit concurrent tasks to a GPU.
    // If no stream is used, the default is used.
    // Tasks are placed in the default stream queue but GPU may still execute those tasks concurrently.
    // Using streams should be faster.
    //
    // torch::cuda::Stream stream;
    // torch::cuda::StreamGuard guard(stream); 
    // torch::cuda::Stream original_stream = torch::cuda::current_stream();
    // if (use_custom) {
    //         torch::cuda::set_stream(my_custom_stream);
    // }
    // if (use_custom) { // restore
    //         torch::cuda::set_stream(original_stream);
    // }
    //
    // But, even better would be to batch across C++ threads into a single
    // kernel.  This would require thread safe queue in C++ to collect batches
    // from inputs.  It would require to track association of batches to
    // requester for return.  std::promise could help this.  It would require
    // some kind of Nagle's theorem to wait some short time to collect batches.
    // De-batching must take care that the output tensors are actually done.
    // Finding the optimum Nagle delay is ... empirical.  I don't think it's
    // worth it.

}
