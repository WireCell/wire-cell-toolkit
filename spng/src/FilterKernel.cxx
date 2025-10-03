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
        // call me after setting my m_cfg by constructor or configure().

        if (m_cfg.kind != "lowpass" && m_cfg.kind != "highpass") {
            raise<ValueError>("unsupported kind of filter: %s", m_cfg.kind);
        }
        if (m_cfg.period <= 0) {
            raise<ValueError>("period must be positive");
        }
        if (m_cfg.scale <= 0) {
            raise<ValueError>("scale must be positive");
        }

        if (m_cfg.kind == "lowpass") { // aka "hf"
            m_func = [=](double freq) {
                return exp(-0.5 * pow(freq / m_cfg.scale, m_cfg.power));
            };
        }
        else {
            m_func = [=](double freq) {
                return 1 - exp(-pow(freq / m_cfg.scale, m_cfg.power));
            };
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
        if (shape.size() != 1) {
            raise<ValueError>("only 1D supported by FilterKernel, got %d", shape.size());
        }

        // Build on CPU, will send to device just before return.
        auto s = torch::zeros(shape, torch::kFloat);

        if (! m_cfg.ignore_baseline) {
            s[0] = m_func(0.0);
        }

        // From 0 to sample frequency
        const size_t nbins = shape[0];
        Binning bins(nbins, 0, 1/m_cfg.period);

        const size_t half = (nbins+1)/2;
        for (size_t ind=1; ind<half; ++ind) {
            float freq = bins.edge(ind);
            s[ind] = s[nbins-ind] = m_func(freq);
        }
        if (nbins%2 == 0) {
            float nyquist = bins.edge(half);
            s[half] = bins.edge(nyquist);
        }

        return to(s);
    }
    


}
