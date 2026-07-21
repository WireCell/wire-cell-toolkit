#include "WireCellSpng/RatioKernel.h"
#include "WireCellSpng/TorchLMN.h"
#include "WireCellSpng/Util.h"

#include "WireCellUtil/Response.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(SPNGRatioKernel,
                 WireCell::SPNG::RatioKernel,
                 WireCell::ITorchSpectrum,
                 WireCell::IConfigurable)

namespace WireCell::SPNG {

    using WireCell::SPNG::nhalf;

    RatioKernel::RatioKernel()
        : Logger("RatioKernel", "spng")
    {
    }

    RatioKernel::shape_t RatioKernel::shape() const
    {
        auto shape_num = m_num->shape();
        auto shape_den = m_den->shape();

        return {std::min(shape_num[0], shape_den[0]),
                std::max(shape_num[1], shape_den[1])};
    }

    WireCell::Configuration RatioKernel::default_configuration() const
    {
        auto cfg = this->ContextBase::default_configuration();
        cfg["numerator"] = "";
        cfg["denominator"] = "";
        return cfg;
    }

    void RatioKernel::configure(const WireCell::Configuration& cfg)
    {
        this->ContextBase::configure(cfg);
        m_num = Factory::find_tn<ITorchSpectrum>(cfg["numerator"].asString());
        m_den = Factory::find_tn<ITorchSpectrum>(cfg["denominator"].asString());
    }

    torch::Tensor RatioKernel::spectrum(const shape_t & shape) const
    {
        auto num = to(m_num->spectrum(shape));
        auto den = to(m_den->spectrum(shape));
        return num/den;
    }
    


}
