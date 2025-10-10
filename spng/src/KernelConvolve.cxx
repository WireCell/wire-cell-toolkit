#include "WireCellSpng/KernelConvolve.h"
#include "WireCellSpng/Convo.h"
#include "WireCellSpng/SimpleTorchTensor.h"

#include "WireCellUtil/HanaJsonCPP.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellIface/IFieldResponse.h"
#include "WireCellIface/IWaveform.h"

#include "boost/container_hash/hash.hpp"


WIRECELL_FACTORY(SPNGKernelConvolve,
                 WireCell::SPNG::KernelConvolve,
                 WireCell::ITorchSpectrum,
                 WireCell::IConfigurable)

using namespace WireCell::HanaJsonCPP;                         

namespace WireCell::SPNG {

    KernelConvolve::KernelConvolve()
        : Logger("KernelConvolve", "spng")
    {}

    KernelConvolve::KernelConvolve(const KernelConvolveConfig& cfg)
        : Logger("KernelConvolve", "spng")
        , m_cfg(cfg)
    {
        configme();
    }

    WireCell::Configuration KernelConvolve::default_configuration() const
    {
        auto cfg = this->ContextBase::default_configuration();
        auto cfg2 = this->Logger::default_configuration();
        update(cfg, cfg2);
        cfg2 = to_json(m_cfg);
        update(cfg, cfg2);
        return cfg;
    }

    void KernelConvolve::configure(const WireCell::Configuration& config)
    {
        this->ContextBase::configure(config);
        this->Logger::configure(config);
        from_json(m_cfg, config);
        configme();
    }

    void KernelConvolve::configme()
    {
        if (m_cfg.kernel.empty()) {
            raise<ValueError>("the 'kernel' parameter is required");
        }
        m_kernel = Factory::find_tn<ITorchSpectrum>(m_cfg.kernel);
        if (m_cfg.axis.size() != 2) {
            log->warn("number of configured axis: {} is not 2, using defaults for missing or ignoring extra",
                m_cfg.axis.size());
        }
        m_cfg.axis.resize(2);

        auto kshape = m_kernel->shape();
        m_roll.resize(2);
        for (size_t dim=0; dim<2; ++dim) {
            const auto& acfg = m_cfg.axis[dim];

            if (acfg.crop < -2) {
                raise<ValueError>("illegal crop configured: %d", acfg.crop);
            }

            m_roll[dim] = acfg.roll;
            if (acfg.roll_mode == "decon") {
                m_roll[dim] += kshape[dim];
            }
        }
    }

    bool KernelConvolve::operator()(const input_pointer& in, output_pointer& out)
    {
        out = nullptr;
        if (!in) {
            logit("EOS");
            ++m_count;
            return true;
        }

        bool batched = true;

        // Assure the tensor is batched.  Everything that touches "tensor" must
        // take care to consider indices {1,2} to be dimensions {0,1}!
        auto tensor = in->tensor();
        if (tensor.dim() == 2) {
            tensor = tensor.unsqueeze(0);
            batched = false;    // squeeze on output 
        }
        auto tensor_shape = tensor.sizes().vec(); // size 3

        // Consider non batch dimensions!
        std::vector<int64_t> basic_shape(2);
        std::vector<int64_t> convolve_shape(2);
        auto kernel_shape = m_kernel->shape();
        for (size_t dim=0; dim<2; ++dim) {
            const auto& acfg = m_cfg.axis[dim];
            const auto tshape = tensor_shape[dim+1]; // skip batch dim;

            if (acfg.cyclic) {
                convolve_shape[dim] = basic_shape[dim] = tshape;
                continue;
            }
            // linear
            convolve_shape[dim] = basic_shape[dim] = linear_shape(tshape, kernel_shape[dim]);
            if (!m_cfg.faster) {
                continue;
            }
            convolve_shape[dim] = m_faster(convolve_shape[dim]);
        }

        /// Do padding of input tensor
        for (size_t dim=0; dim<2; ++dim) {
            const auto tshape = tensor_shape[dim+1]; // skip batch dim;
            if (tshape == convolve_shape[dim]) {
                // avoid a copy inside resize()
                continue;
            }
            log->debug("resize: {} dim={} size={}",
                       to_string(tensor), dim, convolve_shape[dim]);
            tensor = LMN::resize(tensor, convolve_shape[dim], dim+1);
        }

        // applies to last 2 dimensions by default
        tensor = torch::fft::fft2(tensor);

        auto kernel = m_kernel->spectrum(convolve_shape);

        // This is supposed to not change data shared by other shallow copies.
        kernel = kernel.unsqueeze(0);


        // The reason of our existence, I give you, the CONVOLUTION:
        tensor = tensor * kernel;


        // applies to last 2 dimensions by default
        tensor = torch::real(torch::fft::ifft2(tensor));

        /// Do crop and/or roll
        for (size_t dim=0; dim<2; ++dim) {        
            const int crop = m_cfg.axis[dim].crop;
            if (crop > 0) {  // absolute crop
                tensor = LMN::resize(tensor, crop, dim+1);
            }
            else if (crop == -1) {
                tensor == LMN::resize(tensor, basic_shape[dim], dim+1);
            }
            else if (crop == -2) {
                tensor = LMN::resize(tensor, tensor_shape[dim], dim+1);
            }
            // else, crop==0 and no crop

            if (m_roll[dim]) {
                tensor = torch::roll(tensor, m_roll[dim], dim+1);
            }
        }

        if (! batched) {
            tensor = tensor.squeeze(0);
        }

        out = std::make_shared<SimpleTorchTensor>(tensor); // fixme: md?

        ++m_count;
        return true;
    }

}
