/** Provides TorchContext and TorchSemaphore as ways to manage how SPNG code
 * uses Torch devices.
 *
 * It provides
 *
 * - A configurable "device" object to use when creating tensors and no input
 * tensors are available to supply it.
 *
 * - A semaphore to cooperatively limit the number of tasks on a torch device.
 *
 * - A context manager to handle semaphore and to turn off auto grad mode. 
 *
 *
 * Note: this is a revamped version of classes of the same name from wct/pytorch/.
 */

#ifndef WIRECELL_SPNG_TORCHCONTEXT
#define WIRECELL_SPNG_TORCHCONTEXT

#include "WireCellIface/ISemaphore.h"
#include "WireCellSpng/Torch.h"

namespace WireCell::SPNG {

    class TorchContext {
      public:

        // Assume default device (cpu).
        TorchContext();
        ~TorchContext() = default;

        /// Set or reset device and/or semaphore.
        ///
        /// The device name can be "cpu" (default) or "cuda" (use default CUDA
        /// device) or "gpu" (alias for "cuda") or "gpuN" or "N" where N is a
        /// GPU index number.
        ///
        /// The semnam can be empty or a WCT component "type:name".  If empty,
        /// the "Semaphore" type with an instance "torch-<devname>" is
        /// constructed.
        void connect(const std::string& devname="cpu", const std::string& semname="");

        torch::Device device() const { return m_dev; }
        std::string devname() const { return m_devname; }
        std::string semname() const { return m_semname; }

        bool is_gpu() const { return m_dev.is_cuda();}

        /// Manual context manager methods.
        ///
        /// If used, these should bracket any code that will access the device.
        ///
        /// For a scope-based context manager, use a TorchSemaphore.
        void enter() const;
        void exit() const;

      private:

        torch::Device m_dev{torch::kCPU};
        std::string m_devname{"cpu"}, m_semname{""};
        ISemaphore::pointer m_sem;
        mutable bool m_agm{false};   // to store autogradmode through the context

     };

    /// Use like:
    ///
    /// void mymeth() {
    ///   TorchSemaphore sem(m_ctx);
    ///   ... more code may return/throw
    /// } // end of scope
    class TorchSemaphore {
        const TorchContext& m_th;

      public:
        TorchSemaphore() = delete;

        TorchSemaphore(const TorchContext& th) : m_th(th) {
            m_th.enter();
        }
        ~TorchSemaphore() {
            m_th.exit();
        }
    };            
}
#endif
