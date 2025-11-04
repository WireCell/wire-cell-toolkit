/**
   These fans keep the type unchanged between input and output.
 */

#ifndef WIRECELL_SPNG_FANBASE
#define WIRECELL_SPNG_FANBASE

#include "WireCellSpng/ITorchTensorSet.h"
#include "WireCellSpng/Logger.h"
#include "WireCellSpng/ContextBase.h"
#include "WireCellSpng/HanaConfigurable.h"
#include "WireCellIface/IConfigurable.h"
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

        The sub class should implement at least fanin_combine().
    */
    template<typename IDataType>  
    struct FaninBase : public HanaConfigurable<FaninBase<IDataType>, FanConfig, Logger, ContextBase>,
                       public virtual IFaninNode<IDataType, IDataType, 0>
    {
        using fan_type = IFaninNode<IDataType, IDataType, 0>;
        using input_vector = typename fan_type::input_vector;
        using output_pointer = typename fan_type::output_pointer;


        FaninBase(const std::string& type_name="SPNGFaninTensorSets",
                  const std::string& group_name="spng") {
            auto& cfg = this->Logger::hana_config();
            cfg.type_name = type_name;
            cfg.group_name = group_name;
        }


        void configured() {
            this->debug("fanin multiplicity: {}", 
                        this->m_hana_config.multiplicity);
        }

        /// Subclass must implement to set out which starts as nullptr.
        virtual void fanin_combine(const input_vector& inv, output_pointer& out) = 0;

        /// Let subclass do EOS processing or logging.
        virtual void fanin_eos() {
            this->logit("EOS");
        };

        virtual bool operator()(const input_vector& inv, output_pointer& out) {
            TorchSemaphore sem(this->context());

            size_t multiplicity = this->hana_config().multiplicity;

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
            }
            else {
                fanin_combine(inv, out);
            }
            this->next_count();
            return true;
        }

        /// Supply the IFaninNode 
        virtual std::vector<std::string> input_types() {
            int multiplicity = this->hana_config().multiplicity;
            const std::string tname = std::string(typeid(IDataType).name());
            std::vector<std::string> ret(multiplicity, tname);
            return ret;
        }
        virtual std::string signature() { return typeid(FaninBase<IDataType>).name(); }

    };

    /** @brief Base class providing a fanout DFP node.

        The sub class should implement at least fanout_separate();
    */
    template<typename IDataType>  
    struct FanoutBase : public HanaConfigurable<FanoutBase<IDataType>, FanConfig, Logger, ContextBase>,
                        public virtual IFanoutNode<IDataType, IDataType, 0>
    {
        using fan_type = IFanoutNode<IDataType, IDataType, 0>;
        using input_pointer = typename fan_type::input_pointer;
        using output_vector = typename fan_type::output_vector;

        FanoutBase(const std::string& type_name="SPNGFaninTensorSets",
                   const std::string& group_name="spng") {
            auto& cfg = this->Logger::hana_config();
            cfg.type_name = type_name;
            cfg.group_name = group_name;
        }

        void configured() {
            this->debug("fanout multiplicity: {}", 
                        this->m_hana_config.multiplicity);
        }

        /// Subclass must implement.  outv is pre-sized and holds nullptrs.
        virtual void fanout_separate(const input_pointer& in, output_vector& outv) = 0;

        /// Let subclass do EOS processing or logging.
        virtual void fanout_eos() {
            this->logit("EOS");
        };

        virtual bool operator()(const input_pointer& in, output_vector& outv) {
            TorchSemaphore sem(this->context());

            size_t multiplicity = this->hana_config().multiplicity;
            if (outv.size() > 0 && outv.size() != multiplicity) {
                raise<ValueError>("fanout got unexpected multiplicity, got:%d want:%d",
                                  outv.size(), multiplicity);
            }
            // Start with consistent EOS marker.
            outv.clear();
            outv.resize(multiplicity, nullptr);
            if (! in) {
                fanout_eos();
            }
            else {
                fanout_separate(in, outv);
            }
            this->next_count();
            return true;
        }

        /// Supply IFanoutNode interface.
        virtual std::vector<std::string> output_types() {
            int multiplicity = this->hana_config().multiplicity;
            const std::string tname = std::string(typeid(IDataType).name());
            std::vector<std::string> ret(multiplicity, tname);
            return ret;
        }
        virtual std::string signature() { return typeid(FanoutBase<IDataType>).name(); }
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
#endif
