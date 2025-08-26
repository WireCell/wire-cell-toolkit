/** A mixin class to provide a torch context
    
 */

//#include "WireCellIface/ISemaphore.h"
#include "WireCellSpng/Torch.h"
//#include <torch/torch.h>

namespace WireCell::SPNG {

    class TorchContext {
      public:

        // The "devname" is "cpu" or "gpu" or "gpuN" where N is a GPU
        // number.  If "semname" is given, use it for semaphore,
        // otherwise use canonically tn=Semaphore:torch-<devname>.
        TorchContext(const std::string& devname);
        TorchContext();
        ~TorchContext();

        // Default constructor makes context with no device nor
        // semaphore.  This will make the "connection" to them.
        void connect(const std::string& devname);

        torch::Device device() const { return m_dev; }
        std::string devname() const { return m_devname; }

        bool is_gpu() const { return m_dev.is_cuda();}

        // Context manager methods.  Caller should prefer using a
        // TorchSemaphore class but if called directly, caller MUST
        // balance an enter() with an exit().  These can and should be
        // used in multi-thread run stage.

      private:

        torch::Device m_dev{torch::kCUDA}; //kCUDA or kCPU
        std::string m_devname{""};
     };

    /// Use like:
    ///
    /// void mymeth() {
    ///   TorchSemaphore sem(m_ctx);
    ///   ... more code may return/throw
    /// } // end of scope
    /*
    class TorchSemaphore {
        const TorchContext& m_th;
      public:
        TorchSemaphore(const TorchContext& th) : m_th(th) {
            m_th.enter();
        }
        ~TorchSemaphore() {
            m_th.exit();
        }
    };            
     */
}
