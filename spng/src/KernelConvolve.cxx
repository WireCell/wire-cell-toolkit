#include "WireCellSpng/KernelConvolve.h"
#include "WireCellSpng/Convo.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/TdmTools.h"

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
        auto kernel_tensor = m_kernel->spectrum(kshape);
        if (has_nan(kernel_tensor)) {
            raise<ValueError>("kernel has NaN");
        }


        m_roll.resize(2);
        for (size_t dim=0; dim<2; ++dim) {
            const auto& acfg = m_cfg.axis[dim];

            if (acfg.crop < -2) {
                raise<ValueError>("illegal crop configured: %d", acfg.crop);
            }
            if (kshape[dim] == 0) {
                log->warn("kernel dimension {} has zero size, is this really what you want?", dim);
            }

            m_roll[dim] = acfg.roll;
            if (acfg.roll_mode == "decon") {
                m_roll[dim] += kshape[dim];
            }
        }
        log->debug("using tag: {} and datapath format: {}", m_cfg.tag, m_cfg.datapath_format);


    }

    bool KernelConvolve::operator()(const input_pointer& in, output_pointer& out)
    {
        out = nullptr;
        if (!in) {
            logit("EOS");
            ++m_count;
            return true;
        }

        logit(in, "input");

        // Fixme: uplift tag/datapath_format config and application to a base
        // class.
        Configuration md;
        if (m_cfg.tag.size()) {
            md["tag"] = m_cfg.tag;
        }
        md = TDM::derive_metadata(md, in->metadata(), m_cfg.datapath_format);

        auto tensor = in->tensor();
        if (has_nan(tensor)) {
            log->critical("input tensor has NaNs {}", to_string(tensor));
            raise<ValueError>("input tensor has NaNs");
        }

        // An upstream source may provide an empty tensor due to no activity
        // (specifically sim which makes neither signal nor noise).  Better that
        // a branch/merge protect us from this case but we try to act in good
        // faith and just pass it along to the next sucker^W node.
        if (tensor.size(-1) <= 0 || tensor.size(-2) <= 0) {
            log->warn("empty tensor dimensions at call={}, passing it along.  Are we sparse processing?", m_count);
            out = std::make_shared<SimpleTorchTensor>(tensor, md);
            logit(out, "empty");
            ++m_count;
            // fixme: should still do output derivation
            return true;
        }


        // Assure the tensor is batched.  Everything that touches "tensor" must
        // take care to consider indices {1,2} to be dimensions {0,1}!
        bool batched = true;     // start with assuming ndim==3
        if (tensor.dim() == 2) { // unless shown otherwise
            tensor = tensor.unsqueeze(0);
            batched = false;    // squeeze on output 
        }
        auto tensor_shape = tensor.sizes().vec(); // size 3

        if (tensor_shape.size() != 3) {
            log->critical("illegal number of input tensor dimensions at call=%d: %d",
                          m_count, tensor_shape.size());
            raise<ValueError>("illegal number of input tensor dimensions");
        }

        // Consider non batch dimensions!
        std::vector<int64_t> basic_shape(2);
        std::vector<int64_t> convolve_shape(2);
        auto kernel_shape = m_kernel->shape();
        for (size_t dim=0; dim<2; ++dim) {
            const auto& acfg = m_cfg.axis[dim];
            const auto dim_size = tensor_shape[dim+1]; // skip batch dim;

            if (kernel_shape[dim] == 0) {
                log->warn("shape: natural kernel size of dimension {} is zero, using input size {}", dim, dim_size);
                convolve_shape[dim] = basic_shape[dim] = dim_size;
                continue;
            }

            if (acfg.cyclic) {
                convolve_shape[dim] = basic_shape[dim] = dim_size;
                log->debug("shape: dim={} cyclic size: {}", dim, dim_size);
                continue;
            }
            // linear
            convolve_shape[dim] = basic_shape[dim] = linear_shape(dim_size, kernel_shape[dim]);
            if (!m_cfg.faster) {
                log->debug("shape: dim={} linear size: {}->{}", dim, dim_size, convolve_shape[dim]);
                continue;
            }
            convolve_shape[dim] = m_faster(convolve_shape[dim]);
            log->debug("shape: dim={} faster size: {}->{}->{}",
                       dim, dim_size, basic_shape[dim], convolve_shape[dim]);
        }

        /// Do padding of input tensor
        for (size_t dim=0; dim<2; ++dim) {
            const auto dim_size = tensor_shape[dim+1]; // skip batch dim;
            if (dim_size == convolve_shape[dim]) {
                // avoid a copy inside resize()
                continue;
            }
            log->debug("resize: dim={} {} -> {}",
                       dim, dim_size, convolve_shape[dim]);
            tensor = LMN::resize(tensor, convolve_shape[dim], dim+1);
        }

        // applies to last 2 dimensions by default
        tensor = torch::fft::fft2(tensor);

        auto kernel = m_kernel->spectrum(convolve_shape);
        if (has_nan(kernel)) {
            log->critical("kernel has NaNs {}", to_string(kernel));
            raise<ValueError>("kernel has NaNs");
        }

        // This is supposed to not change data shared by other shallow copies.
        kernel = kernel.unsqueeze(0);


        // The reason of our existence, I give you, the CONVOLUTION:
        tensor = tensor * kernel;


        // applies to last 2 dimensions by default
        tensor = torch::real(torch::fft::ifft2(tensor));

        /// Do crop and/or roll
        for (size_t dim=0; dim<2; ++dim) {        
            const int crop = m_cfg.axis[dim].crop;
            log->debug("before: crop={} dim={} tensor={}", crop, dim, to_string(tensor));

            if (crop > 0) {  // absolute crop
                tensor = LMN::resize(tensor, crop, dim+1);
            }
            else if (crop == -1) {
                // crop away the "faster" padding keep the basic convoulutional size
                tensor == LMN::resize(tensor, basic_shape[dim], dim+1);
            }
            else if (crop == -2) {
                // crop faster and convolutional padding
                const auto dim_size = tensor_shape[dim+1];
                tensor = LMN::resize(tensor, dim_size, dim+1);
            }
            // else, crop==0 and no crop
            log->debug("after: crop={} dim={} tensor={}", crop, dim, to_string(tensor));

            if (m_roll[dim]) {
                tensor = torch::roll(tensor, m_roll[dim], dim+1);
            }
        }

        if (! batched) {
            tensor = tensor.squeeze(0);
        }

        out = std::make_shared<SimpleTorchTensor>(tensor, md);

        logit(out, "output");
        ++m_count;
        return true;
    }

}
