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

    void FilterKernel::configme()
    {
        if (m_cfg.axis.empty()) {
            raise<ValueError>("FilterKernel has empty 'axis' parameter, check config?");
        }

        m_funcs.clear();
        size_t ndim = m_cfg.axis.size();
        for (size_t ind=0; ind<ndim; ++ind) {
            const auto& axis = m_cfg.axis[ind];
            
            if (axis.kind == "identity") {
                m_funcs.emplace_back([axis](double freq) -> float {
                    return 1.0;
                });
                continue;
            }

            // the "pass" filters require legal parameters.

            if (axis.period <= 0) {
                raise<ValueError>("period must be positive");
            }
            if (axis.scale <= 0) {
                raise<ValueError>("scale must be positive");
            }

            if (axis.kind == "lowpass") { // aka "hf"
                m_funcs.emplace_back([axis](double freq) -> float {
                    return exp(-0.5 * pow(freq / axis.scale, axis.power));
                });
            }
            else if (axis.kind == "highpass") { // aka "lf"
                m_funcs.emplace_back([axis](double freq) -> float {
                    return 1 - exp(-pow(freq / axis.scale, axis.power));
                });
            }
            else {
                raise<ValueError>("unsupported axis kind: %s", axis.kind);
            }
            log->debug("filter axis {} with {} period={} scale={} power={}",
                       ind, axis.kind, axis.period, axis.scale, axis.power);
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

        size_t ndims = m_funcs.size();
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
        const auto& sample = m_funcs[axis];

        int64_t size = shape[axis];

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
        size_t ndims = m_funcs.size();

        if (shape.size() != ndims) {
            raise<ValueError>("dimension mismatch, want: %d got: %d.  Check config or code?",
                              ndims, shape.size());
        }
        if (shape.size() > 3) {
            raise<ValueError>("no support for %d dimension, current max is 3.  Extend the code.",
                              shape.size());
        }

        
        std::vector<torch::Tensor> axis;

        for (size_t idim=0; idim<ndims; ++idim) {
            auto s = spectrum_axis(shape, idim);
            axis.push_back(s);
        }

        // 1D
        if (ndims == 1) {
            return to(axis[0]);
        }
        
        // 2D
        if (ndims == 2) {
            return to(torch::outer(axis[0], axis[1]));
        }

        // 3D
        return to(torch::einsum("i,j,k->ijk", axis));

        // Future extension for 4D etc can be added easily by giving the
        // appropriate, larger einsum string.  Actually, all dimensions can be
        // symmetric with a call to einsum("...", oneds) with a per-dimension
        // string....
    }
    


}
