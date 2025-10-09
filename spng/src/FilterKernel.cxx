#include "WireCellSpng/FilterKernel.h"

#include "WireCellUtil/Response.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(SPNGFilterKernel,
                 WireCell::SPNG::FilterKernel,
                 WireCell::ITorchSpectrum,
                 WireCell::IConfigurable)

namespace WireCell::SPNG {

    FilterKernel::FilterKernel() {
    }
    FilterKernel::FilterKernel(const FilterKernelConfig& cfg)
        : m_cfg(cfg)
    {
        configme();
    }

    void FilterKernel::configme()
    {
        if (m_cfg.axis.empty()) {
            raise<ValueError>("FilterKernel has empty 'axis' parameter, check config?");
        }

        m_funcs.clear();
        for (const auto& axis : m_cfg.axis) { 
            
            if (axis.kind != "lowpass" && axis.kind != "highpass") {
                raise<ValueError>("unsupported kind of filter: %s", axis.kind);
            }
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
            else {
                m_funcs.emplace_back([axis](double freq) -> float {
                    return 1 - exp(-pow(freq / axis.scale, axis.power));
                });
            }
        }
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
            const auto& cfg = m_cfg.axis[idim];
            const auto& sample = m_funcs[idim];

            int64_t size = shape[idim];

            // Build on CPU.  We will send to device just before return.
            auto s = torch::zeros(size, torch::kFloat);

            if (! cfg.ignore_baseline) {
                s[0] = sample(0.0);
            }

            // From 0 to sample frequency
            Binning bins(size, 0, 1.0/cfg.period);

            const size_t half = (size+1)/2;
            for (size_t ind=1; ind<half; ++ind) {
                float freq = bins.edge(ind);
                s[ind] = s[size-ind] = sample(freq);
            }
            if (size%2 == 0) {
                float nyquist = bins.edge(half);
                s[half] = bins.edge(nyquist);
            }

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
