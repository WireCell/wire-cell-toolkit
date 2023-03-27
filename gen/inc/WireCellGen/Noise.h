#ifndef WIRECELL_GEN_NOISE
#define WIRECELL_GEN_NOISE

#include "WireCellAux/Logger.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IRandom.h"
#include "WireCellIface/IDFT.h"

#include "WireCellUtil/NamedFactory.h"

#include <map>
#include <string>
#include <unordered_set>

namespace WireCell::Gen {


    /// Add*Noise and *NoiseSource all share some commonalities.
    class NoiseBase : public Aux::Logger,
                      public IConfigurable
    {
      public:
        NoiseBase(const std::string& tname);
        virtual ~NoiseBase();        
        /// IConfigurable
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;

      protected:
        IRandom::pointer m_rng;
        IDFT::pointer m_dft;

        std::unordered_set<std::string> m_model_tns;
        size_t m_nsamples{9600};
        double m_rep_percent{0.02};
        size_t m_count{0};

        double m_bug202{0.0};   // non-zero emulates issue #202
    };

    /// Incoherent and coherent noise each have their own type of
    /// model
    template<class IModel>
    class NoiseBaseT : public NoiseBase {
      public:
        NoiseBaseT(const std::string& tname) : NoiseBase(tname) {}
        virtual ~NoiseBaseT() {}
        virtual void configure(const WireCell::Configuration& config) {
            NoiseBase::configure(config);
            for (const auto& mtn : m_model_tns) {
                m_models[mtn] = Factory::find_tn<IModel>(mtn);
            }
        }

      protected:
        std::map<std::string, typename IModel::pointer> m_models;

    };

}

#endif
