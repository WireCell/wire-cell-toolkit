#include "WireCellSpng/FilterKernel.h"
#include "WireCellSpng/TorchLMN.h"
#include "WireCellSpng/Util.h"

#include "WireCellUtil/Response.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(SPNGFilterKernel,
                 WireCell::SPNG::FilterKernel,
                 WireCell::ITorchSpectrum,
                 WireCell::IConfigurable)

namespace WireCell::SPNG {

    using WireCell::SPNG::nhalf;

    FilterKernel::FilterKernel()
        : Logger("FilterKernel", "spng")
    {
    }
    FilterKernel::FilterKernel(const FilterKernelConfig& cfg)
        : Logger("FilterKernel", "spng")
        , m_cfg(cfg)
    {
        configme();
    }

    FilterKernel::shape_t FilterKernel::shape() const
    {
        shape_t shape;
        const size_t ndims = m_cfg.axis.size();
        for (size_t ind=0; ind<ndims; ++ind) {
            shape.push_back(m_cfg.axis[ind].size);
        }
        return shape;
    }


    void FilterKernel::configme()
    {
        const size_t ndim = m_cfg.axis.size();

        if (ndim == 0) {
            raise<ValueError>("FilterKernel has empty 'axis' parameter, check config?");
        }
        if (ndim > 3) {
            raise <ValueError>("a filter kernel only supports up to 3D (currently)");
        }

        m_funcs.clear();

        for (size_t ind=0; ind<ndim; ++ind) {
            const auto& axis = m_cfg.axis[ind];

            axis_functions funcs;

            for (auto filter : axis.filters) {

                log->debug("add to filter axis {} kind=\"{}\" "
                           "size={} period={} scale={} power={} ignore_baseline={}",
                           ind, filter.kind, axis.size, axis.period, filter.scale,
                           filter.power, axis.ignore_baseline);
                
                if (filter.kind == "flat") {
                    funcs.emplace_back([filter](double freq) -> float {
                        return filter.scale;
                    });
                    continue;
                }
            
                // the "pass" filters have parameter constraints.
                if (axis.period <= 0) {
                    raise<ValueError>("period must be positive");
                }
                if (filter.scale <= 0) {
                    raise<ValueError>("scale must be positive");
                }

                if (filter.kind == "lowpass") { // aka "hf"
                    funcs.emplace_back([filter](double freq) -> float {
                        return exp(-0.5 * pow(freq / filter.scale, filter.power));
                    });
                }
                else if (filter.kind == "highpass") { // aka "lf"
                    funcs.emplace_back([filter](double freq) -> float {
                        return 1 - exp(-pow(freq / filter.scale, filter.power));
                    });
                }
                else {
                    raise<ValueError>("unsupported filter kind: %s", filter.kind);
                }
            }
            m_funcs.emplace_back([funcs](double freq) -> float {
                double val = 1.0;
                for (const auto& func : funcs) {
                    val *= func(freq);
                }
                return val;
            });
        }
            
        if (m_cfg.debug_filename.size()) {
            log->debug("writing debug file: {}", m_cfg.debug_filename);
            write_debug(m_cfg.debug_filename);
        }

    }

    void FilterKernel::write_debug(const std::string& filename) const
    {
        /// There is no "native" size for the filter since it is a sampled
        /// analytical function.  Ultimately it is requested to be the padded
        /// size required for M*F/R.  Here we will kludge that size without
        /// going overboard..
        std::vector<int64_t> shape = {100,2000};
        
        using tensor_map = torch::Dict<std::string, torch::Tensor>;
        tensor_map to_save;

        const size_t ndims = m_funcs.size();
        for (size_t idim=0; idim<ndims; ++idim) {
            auto s = spectrum_axis(shape, idim);
            to_save.insert(fmt::format("filter1d{}", idim), s);
        }
        to_save.insert("filter2d", spectrum(shape));
        auto data = torch::pickle_save(to_save);
        std::ofstream output_file(filename, std::ios::binary);
        output_file.write(data.data(), data.size());
    }


    void FilterKernel::configure(const WireCell::Configuration& cfg)
    {
        this->ContextBase::configure(cfg);
        HanaJsonCPP::from_json(m_cfg, cfg);
        configme();
    }

    WireCell::Configuration FilterKernel::default_configuration() const
    {
        auto cfg = HanaJsonCPP::to_json(m_cfg);
        auto cfg2 = this->ContextBase::default_configuration();
        update(cfg, cfg2);
        return cfg;
    }


    torch::Tensor FilterKernel::spectrum_axis(const shape_t & shape, int64_t axis) const
    {
        const auto& cfg = m_cfg.axis[axis];

        const int64_t size = shape[axis];            

        if (cfg.filters.empty()) {
            // Note, this "sampling" does not get into a final kernel.
            return torch::ones(size, torch::kFloat);
        }

        const auto& sample = m_funcs[axis];

        // Build on CPU.  We will send to device just before return.
        auto s = torch::zeros(size, torch::kFloat);

        if (! cfg.ignore_baseline) {
            s[0] = sample(0.0);
        }

        // From 0 to sample frequency
        Binning bins(size, 0, 1.0/cfg.period);

        // Number of non-zero/non-nyquist "positive frequencies".
        const size_t half = nhalf(size);
        for (size_t ind=1; ind <= half; ++ind) {
            float freq = bins.edge(ind);
            s[ind] = s[size-ind] = sample(freq);
        }
        if (size%2 == 0) {      // a nyquist bin exists.
            int64_t nyquist_sample = size/2;
            float nyquist_frequency = bins.edge(nyquist_sample);
            s[half+1] = sample(nyquist_frequency);
        }
        return s;
    }

    torch::Tensor FilterKernel::spectrum(const shape_t & shape) const
    {
        const size_t ndims = m_funcs.size();

        if (shape.size() != ndims) {
            raise<ValueError>("dimension mismatch, want: %d got: %d.  Check config or code?",
                              ndims, shape.size());
        }
        if (shape.size() > 3) {
            raise<ValueError>("no support for %d dimension, current max is 3.  Extend the code.",
                              shape.size());
        }

        
        std::vector<torch::Tensor> full_axis, none_axis;
        std::vector<int64_t> none_dims;

        for (size_t idim=0; idim<ndims; ++idim) {
            auto s = spectrum_axis(shape, idim);

            log->debug("filter spectrum dim={} size={} nfuncs=\"{}\"",
                       idim, shape[idim], m_cfg.axis[idim].filters.size());

            if (m_cfg.axis[idim].filters.empty()) {
                none_axis.push_back(s);
                none_dims.push_back(idim);
            }
            else {
                full_axis.push_back(s);
            }
        }

        // To handle unsqueezed axes, we must have at least one sampled axis.
        // If all are "none" we use the last one to provide a "real" sampled
        // tensor and will unsqueeze to get the rest.
        if (full_axis.empty()) {
            full_axis.push_back(none_axis.back());
            none_axis.pop_back();
            none_dims.pop_back();
        }

        // 1D special case must end up with full_axis as the sole holder.
        if (ndims == 1) {
            return to(full_axis[0]);
        }

        torch::Tensor result;

        if (full_axis.size() == 1) {
            // No outer product needed.
            result = full_axis[0];
        }
        else {
            // Craft the einsum "i,j,k->ijk" directive to do outer product;

            // Technically, this supports up to 6D of "real" axes.  We could
            // relax the 3D clamp in configme().
            const std::string index_letters = "ijklmn";
            std::string from, to, comma="";
            for (size_t ind=0; ind<full_axis.size(); ++ind) {
                from += comma;
                comma = ",";
                from.push_back(index_letters[ind]);
                to.push_back(index_letters[ind]);
            }                
            std::string ein = from + "->" + to;
            log->debug("filter outer product \"{}\"", ein);
            result = torch::einsum(ein, full_axis);
        }

        // This relies dimensions sorted in increasing order.
        for (auto dim : none_dims) {
            log->debug("unsqueeze dim={}", dim);
            result = result.unsqueeze(dim);
        }

        return to(result);
    }
    


}
