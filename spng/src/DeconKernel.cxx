#include "WireCellSpng/DeconKernel.h"

#include "WireCellUtil/HanaJsonCPP.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellIface/IFieldResponse.h"
#include "WireCellIface/IWaveform.h"

#include "boost/container_hash/hash.hpp"


// First use of new configuration struct pattern!
//
// Here we register the structure of the struct with hana.  Thanks, hana.
BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::DeconKernelConfig,
                        field_response,
                        elec_response,
                        time_filter,
                        channel_filter,
                        plane_id,
                        capacity);
                         
using namespace WireCell::HanaJsonCPP;                         

namespace WireCell::SPNG {

    DeconKernel::DeconKernel()
        : Logger("DeconKernel")
        , m_cache(1)
    {}

    void DeconKernel::configure(const WireCell::Configuration& config)
    {
        this->Logger::configure(config);
        this->ContextBase::configure(config);

        // use the hana serializer
        from_json(m_cfg, config);

        if (m_cfg.plane_id < 0 || m_cfg.plane_id>3) {
            raise<ValueError>("illegal plane ID: %d", m_cfg.plane_id);
        }
        m_cache.reset_capacity(m_cfg.capacity);


        m_time_filter = Factory::find_tn<IFilterWaveform>(m_cfg.time_filter);
        m_channel_filter = Factory::find_tn<IFilterWaveform>(m_cfg.channel_filter);

        auto ifr = Factory::find_tn<IFieldResponse>(m_cfg.field_response);
        auto ier = Factory::find_tn<IWaveform>(m_cfg.elec_response);


        
        torch::Tensor m_natural_spectrum;
        torch::Tensor m_response_waveform;


    }

    WireCell::Configuration DeconKernel::default_configuration() const
    {
        auto cfg = this->Logger::default_configuration();
        auto cfg2 = this->ContextBase::default_configuration();
        update(cfg, cfg2);

        // Use the hana serializer
        cfg2 = to_json(m_cfg);

        update(cfg, cfg2);
        return cfg;
    }

    size_t DeconKernel::make_cache_key(const shape_t& shape) const {
        size_t h = 0;
        for (const auto& s : shape) {
            boost::hash_combine(h, s);
        }
        return h;
    }

    torch::Tensor DeconKernel::spectrum() const
    {
        return m_natural_spectrum; // hard-cache
    }

    torch::Tensor DeconKernel::spectrum(const std::vector<int64_t> & shape) const
    {
        torch::Tensor ten;
        return ten;
    }

    std::vector<int64_t> DeconKernel::shape() const {
        std::vector<int64_t> ret;
        return ret;
    }

    std::vector<int64_t> DeconKernel::shifts() const {
        std::vector<int64_t> ret;
        return ret;
    }
}
