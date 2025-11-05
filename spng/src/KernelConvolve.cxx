#include "WireCellSpng/KernelConvolve.h"
#include "WireCellSpng/Convo.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/TdmTools.h"
#include "WireCellSpng/HanaConfigurable.h"

#include "WireCellUtil/HanaJsonCPP.h"
#include "WireCellUtil/String.h"

#include "WireCellUtil/NamedFactory.h"

#include "boost/container_hash/hash.hpp"


WIRECELL_FACTORY(SPNGKernelConvolve,
                 WireCell::SPNG::KernelConvolve,
                 WireCell::ITorchTensorFilter,
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
        auto cfg = WireCell::default_configuration_bases<KernelConvolve, ContextBase, Logger>(this);
        auto cfg2 = to_json(m_cfg);
        update(cfg, cfg2);
        return cfg;
    }

    void KernelConvolve::configure(const WireCell::Configuration& config)
    {
        WireCell::configure_bases<KernelConvolve, ContextBase, Logger>(this, config);
        from_json(m_cfg, config);
        configme();
    }

    torch::Tensor KernelConvolve::convolve(torch::Tensor tensor, torch::Tensor kernel)
    {
        // 2D
        if (m_cfg.axis[0].dft && m_cfg.axis[1].dft) {
            return torch::real(torch::fft::ifft2(torch::fft::fft2(tensor) * kernel));
        }

        int64_t dim = 0;
        if (m_cfg.axis[1].dft) {
            dim = 1;
        };
        log->debug("convolve 1d along axis dim {}", dim);

        // 1D
        auto spec = torch::fft::fft(tensor, std::nullopt, dim+1);
        auto conv = spec * kernel;
        return torch::real(torch::fft::ifft(conv, std::nullopt, dim+1));
    }

    void KernelConvolve::configme()
    {
        if (m_cfg.kernel.empty()) {
            raise<ValueError>("the 'kernel' parameter is required");
        }
        m_kernel = Factory::find_tn<ITorchSpectrum>(m_cfg.kernel);

        // sanity check the kernel
        auto kshape = m_kernel->shape();
        auto kernel_tensor = m_kernel->spectrum(kshape);
        if (has_nan(kernel_tensor)) {
            raise<ValueError>("kernel has NaN");
        }
        const size_t ndims = m_kernel->shape().size();
        if (ndims != 2) {
            raise<ValueError>("convolution requires 2D kernel, got (%d)", ndims);
        }

        const size_t naxes = m_cfg.axis.size();
        if (naxes != 2) {
            raise<ValueError>("convolution requires 2D axes, got (%d)", naxes);
        }


        // Collect and check per input dimension values.
        m_roll.clear();
        m_crop.clear();
        m_roll.resize(2, 0);
        m_crop.resize(2, 0);

        for (size_t dim=0; dim<naxes; ++dim) {
            auto& acfg = m_cfg.axis[dim];

            if (!String::has({"head","tail","none"}, acfg.padding)) {
                raise<ValueError>("unsupported padding mode: \"%s\" for axis %d", acfg.padding, dim);
            }

            if (acfg.crop < -2) {
                raise<ValueError>("illegal crop configured: %d for axis %d", acfg.crop, dim);
            }
            m_crop[dim] = acfg.crop;

            m_roll[dim] = acfg.roll;
            if (acfg.roll_mode == "decon") {
                m_roll[dim] += kshape[dim];
            }
            log->debug("convolution axis dim={} roll={} crop={} natural size={} dft={}",
                       dim, m_roll[dim], m_crop[dim], kshape[dim], acfg.dft);

        }
        log->debug("using tag: {} and datapath format: {}", m_cfg.tag, m_cfg.datapath_format);
    }

    bool KernelConvolve::operator()(const input_pointer& in, output_pointer& out)
    {
        out = nullptr;
        if (!in) {
            logit("EOS");
            next_count();
            return true;
        }

        logit(in, "input");

        TorchSemaphore sem(context());

        // Fixme: uplift tag/datapath_format config and application to a base
        // class.
        Configuration md;
        if (m_cfg.tag.size()) {
            md["tag"] = m_cfg.tag;
        }
        md = TDM::derive_metadata(md, in->metadata(), m_cfg.datapath_format);
        
        // We will update time as we head-pad or do a roll.
        double time = get<double>(md, "time", 0.0);
        const double sample_period = get<double>(md, "period", 1.0); // 1.0 is certainly wrong
        // fixme: I'm ignoring tbin.

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
            log->warn("empty tensor dimensions at call={}, passing it along.  Are we sparse processing?",
                      get_count());
            out = std::make_shared<SimpleTorchTensor>(tensor, md);
            logit(out, "empty");
            next_count();
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
                          get_count(), tensor_shape.size());
            raise<ValueError>("illegal number of input tensor dimensions");
        }

        // For saving out debug tensors
        using tensor_map = torch::Dict<std::string, torch::Tensor>;
        tensor_map to_save;
        auto maybe_save = [&](torch::Tensor ten, std::string name) {
            if (m_cfg.debug_filename.size()) {
                if (!batched) {
                    ten = ten.squeeze(0);
                }
                to_save.insert(name, ten);
            }
        };

        // Consider non batch dimensions!
        std::vector<int64_t> basic_shape(2);
        std::vector<int64_t> convolve_shape(2);
        auto kernel_shape = m_kernel->shape();
        for (size_t dim=0; dim<2; ++dim) {
            const auto& acfg = m_cfg.axis[dim];
            const size_t batched_dim = dim+1;
            const auto dim_size = tensor_shape[batched_dim];

            if (acfg.baseline) {
                auto medians = std::get<0>(torch::median(tensor, batched_dim, false));
                medians = medians.unsqueeze(batched_dim);
                tensor = tensor - medians;
                log->debug("median baseline subtraction for dimension {} of size {}", dim, dim_size);

                maybe_save(tensor, fmt::format("median_subtracted_dim{}", dim));
            }

            // Find the padding sizes.
            convolve_shape[dim] = basic_shape[dim] = dim_size;

            // Ways we may not actually pad:
            if (acfg.padding == "none") {
                log->debug("shape: dim={} \"{}\" padding, size: {}, dft={}",
                           dim, acfg.padding, dim_size, acfg.dft);
                continue;
            }
            if (kernel_shape[dim] == 0) {
                log->warn("shape: natural kernel size of dimension {} is zero, using input size {}",
                          dim, dim_size);
                continue;
            }
            if (acfg.cyclic) {
                log->debug("shape: dim={} cyclic size: {}", dim, dim_size);
                continue;
            }

            // Achieve linear convolution.
            convolve_shape[dim] = basic_shape[dim] = linear_shape(dim_size, kernel_shape[dim]);
            if (!m_cfg.faster) {
                log->debug("shape: dim={} linear size: {}->{}", dim, dim_size, convolve_shape[dim]);
                continue;
            }

            // Use a somewhat larger, faster DFT size.
            convolve_shape[dim] = m_faster(convolve_shape[dim]);
            log->debug("shape: dim={} faster size: {}->{}->{}",
                       dim, dim_size, basic_shape[dim], convolve_shape[dim]);
        }

        /// Do the actual padding of input tensor
        for (size_t dim=0; dim<2; ++dim) {
            const auto& acfg = m_cfg.axis[dim];
            const size_t batched_dim = dim+1;
            const auto dim_size = tensor_shape[batched_dim];
            if (dim_size == convolve_shape[dim]) {
                continue;
            }
            log->debug("resize: dim={} {} -> {}",
                       dim, dim_size, convolve_shape[dim]);
            // axis, index, size
            if (acfg.padding == "tail") {
                tensor = resize_tensor_tail(tensor, batched_dim, convolve_shape[dim]);
            }
            else if (acfg.padding == "head") {
                // "head" padding 
                tensor = resize_tensor_head(tensor, batched_dim, convolve_shape[dim]);
                if (dim == 1) {
                    time -= convolve_shape[dim] * sample_period;
                }
            }
            /// else padding is "none"

            maybe_save(tensor, fmt::format("resized_dim{}", dim));
        }

        auto kernel = m_kernel->spectrum(convolve_shape);
        if (has_nan(kernel)) {
            log->critical("kernel has NaNs {}", to_string(kernel));
            raise<ValueError>("kernel has NaNs");
        }

        // This is supposed to not change data shared by other shallow copies.
        kernel = kernel.unsqueeze(0);

        // Our main event of the evening.
        tensor = convolve(tensor, kernel);

        maybe_save(tensor, "raw_convo");

        /// Do crop and/or roll
        for (size_t dim=0; dim<2; ++dim) {        
            const size_t batched_dim = dim + 1; // the dimension in the batched tensor

            const int crop = m_crop[dim];
            log->debug("axis before: dim={} crop={} tensor={} dft={}", dim, crop, to_string(tensor), m_cfg.axis[dim].dft);

            if (m_roll[dim]) {
                tensor = torch::roll(tensor, m_roll[dim], batched_dim);
                if (dim == 1) {
                    time -= m_roll[dim] * sample_period;
                }
            }

            if (crop > 0) {  // absolute crop to user provided size
                tensor = resize_tensor_tail(tensor, batched_dim, crop);

            }
            else if (crop == -1) {
                // crop away the "faster" padding keep the basic convoulutional size
                tensor == resize_tensor_tail(tensor, batched_dim, basic_shape[dim]);
            }
            else if (crop == -2) {
                // crop both the "faster" and the convolutional padding
                const auto dim_size = tensor_shape[batched_dim];
                tensor = resize_tensor_tail(tensor, batched_dim, dim_size);
            }
            // else, crop==0 and no crop
            log->debug("axis after: dim={} crop={} tensor={} dft={}", dim, crop, to_string(tensor), m_cfg.axis[dim].dft);
        }

        if (! batched) {
            tensor = tensor.squeeze(0);
        }

        if (m_cfg.debug_filename.size()) {
            std::string filename = fmt::format(m_cfg.debug_filename, fmt::arg("ident", get_count()));
            auto data = torch::pickle_save(to_save);
            std::ofstream output_file(filename, std::ios::binary);
            output_file.write(data.data(), data.size());
        }

        md["time"] = time;
        out = std::make_shared<SimpleTorchTensor>(tensor, md);

        logit(out, "output");
        next_count();
        return true;
    }

}
