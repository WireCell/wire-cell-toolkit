#pragma once

#include "WireCellSpng/Logger.h"
#include "WireCellSpng/ContextBase.h"
#include "WireCellSpng/HanaConfigurable.h"
#include "WireCellSpng/TensorTools.h"

#include "WireCellIface/IFaninNode.h"
#include "WireCellIface/IFanoutNode.h"

namespace WireCell::SPNG {

    /** @brief Configuration common to fans.
     */
    struct FanConfig {
        /// The size of the fan.
        int multiplicity = 2;
    };

}

BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::FanConfig, multiplicity);

namespace WireCell::SPNG {


    /** @brief Base class providing a fanin DFP node.

        The sub class should implement at least fanin_combine() to apply a fan
        operation.

        The data type of objects in the fan is IFansType and output is
        IScalarType which defaults to IFansType.

    */
    template<typename IFansType, typename IScalarType = IFansType>  
    struct FaninBase : public Logger,
                       public ContextBase,
                       public virtual IConfigurable,
                       public virtual IFaninNode<IFansType, IScalarType, 0>
    {
        using fan_type = IFaninNode<IFansType, IScalarType, 0>;
        using input_vector = typename fan_type::input_vector;
        using output_pointer = typename fan_type::output_pointer;


        // If you must call this default constructor, best to call Logger::init().
        FaninBase() = default;
        FaninBase(const std::string& type_name)
            : Logger(type_name, "spng") {
        }

        FaninBase(const std::string& type_name,
                  const std::string& group_name)
            : Logger(type_name, group_name) {
        }

        FanConfig m_config;

        // IConfigurable API.
        virtual void configure(const WireCell::Configuration& jconfig) {
            WireCell::configure_bases<FaninBase<IFansType, IScalarType>, Logger, ContextBase>(this, jconfig);
            HanaJsonCPP::from_json(m_config, jconfig);

            this->debug("fanin multiplicity: {}", 
                        this->m_config.multiplicity);
        }
                
        virtual WireCell::Configuration default_configuration() const {
            auto cfg = WireCell::default_configuration_bases<FaninBase<IFansType, IScalarType>, Logger, ContextBase>(this);
            update(cfg, HanaJsonCPP::to_json(m_config));
            return cfg;
        }

        /// Subclass must implement to set out which starts as nullptr.
        virtual void fanin_combine(const input_vector& inv, output_pointer& out) = 0;

        /// Let subclass do EOS processing or override default logging.
        virtual void fanin_eos() { this->logit("EOS"); };

        virtual bool operator()(const input_vector& inv, output_pointer& out) {
            TorchSemaphore sem(this->context());

            size_t multiplicity = m_config.multiplicity;

            if (inv.empty()) {
                this->critical("empty vector input to fanin, something dreadfully wrong");
                raise<ValueError>("fanin empty vector input to fanin, something dreadfully wrong");
            }
            if (inv.size() != multiplicity) {
                raise<ValueError>("fanin unexpected multiplicity, got:%d want:%d",
                                  inv.size(), multiplicity);
            }

            out=nullptr;
            size_t neos = std::count(inv.begin(), inv.end(), nullptr);
            if (neos) {
                fanin_eos();
                this->next_count();
                return true;
            }

            auto device_inv = to_device(inv, device());

            for (const auto& one : device_inv) {
                logit(one, "fanin");
            }

            fanin_combine(device_inv, out);

            this->next_count();
            return true;
        }

        /// Supply the IFaninNode 
        virtual std::vector<std::string> input_types() {
            const std::string tname = std::string(typeid(IFansType).name());
            std::vector<std::string> ret(m_config.multiplicity, tname);
            return ret;
        }
        virtual std::string signature() { return typeid(fan_type).name(); }

    };

    /** @brief Base class providing a fanout DFP node.

        The sub class should implement at least fanout_separate();

        The type in the fan is IFansType and output is IScalarType which defaults to IFansType.
    */
    template<typename IFansType, typename IScalarType = IFansType>  
    struct FanoutBase : public Logger,
                        public ContextBase,
                        public virtual IConfigurable,
                        public virtual IFanoutNode<IScalarType, IFansType, 0>
    {
        using fan_type = IFanoutNode<IScalarType, IFansType, 0>;
        using input_pointer = typename fan_type::input_pointer;
        using output_vector = typename fan_type::output_vector;

        // If you must call this default constructor, best to call Logger::init().
        FanoutBase() = default;
        FanoutBase(const std::string& type_name)
            : Logger(type_name, "spng") {
        }

        FanoutBase(const std::string& type_name,
                   const std::string& group_name)
            : Logger(type_name, group_name) {
        }

        FanConfig m_config;

        // IConfigurable API.
        virtual void configure(const WireCell::Configuration& jconfig) {
            configure_bases<FanoutBase<IFansType, IScalarType>, Logger, ContextBase>(this, jconfig);
            HanaJsonCPP::from_json(m_config, jconfig);

            this->debug("fanin multiplicity: {}", 
                        this->m_config.multiplicity);
        }
                
        virtual WireCell::Configuration default_configuration() const {
            auto cfg = default_configuration_bases<FanoutBase<IFansType, IScalarType>, Logger, ContextBase>(this);
            update(cfg, HanaJsonCPP::to_json(m_config));
            return cfg;
        }

        /// Subclass must implement.  outv is pre-sized and holds nullptrs.
        virtual void fanout_separate(const input_pointer& in, output_vector& outv) = 0;

        /// Let subclass do EOS processing or override default logging.
        virtual void fanout_eos() { this->logit("EOS"); }

        virtual bool operator()(const input_pointer& in, output_vector& outv) {
            TorchSemaphore sem(this->context());

            size_t multiplicity = m_config.multiplicity;
            if (outv.size() > 0 && outv.size() != multiplicity) {
                raise<ValueError>("fanout got unexpected multiplicity, got:%d want:%d",
                                  outv.size(), multiplicity);
            }
            // Start with consistent EOS marker.
            outv.clear();
            outv.resize(multiplicity, nullptr);
            if (! in) {
                fanout_eos();
            this->next_count();
            return true;
            }

            auto device_in = to_device(in, device());
            fanout_separate(device_in, outv);

            this->next_count();
            return true;
        }

        /// Supply IFanoutNode interface.
        virtual std::vector<std::string> output_types() {
            const std::string tname = std::string(typeid(IFansType).name());
            std::vector<std::string> ret(m_config.multiplicity, tname);
            return ret;
        }
        virtual std::string signature() { return typeid(fan_type).name(); }
    };





    // deprecated 
    class FanBase : public Logger {
    public:
        FanBase(const std::string& logname, const std::string& pkgnam="spng");
        virtual ~FanBase() = default;

        // IConfigurable
        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

    protected:

        /// Configuration
        ///
        /// - multiplicity :: the fan-out multiplicity.
        size_t m_multiplicity{2};
    };
}

