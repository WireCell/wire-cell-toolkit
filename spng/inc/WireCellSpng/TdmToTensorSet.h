#pragma once

#include "WireCellSpng/ITorchToTensorSet.h"
#include "WireCellSpng/Logger.h"
#include "WireCellSpng/ContextBase.h" // may not be needed, we must safely call .to(cpu).
#include "WireCellUtil/HanaJsonCPP.h"


namespace WireCell::SPNG {

    struct TdmToTensorSetConfig {
        /// An array of "match" objects to determine which tensors in the input
        /// set are included for consideration.  If empty (the default) then all
        /// are input tensors included, otherwise, only those that match are
        /// included.  Included tensors are de-duplicated in first-seen order. 
        std::vector<Configuration> include_rules = {};
        
        /// An array of "match" objects to determine which included tensors are
        /// ultimately excluded.
        std::vector<Configuration> exclude_rules = {};

        // something to muck with the tensor set metadata

    };
}

BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::TdmToTensorSetConfig, include_rules, exclude_rules);

namespace WireCell::SPNG {

    ///
    struct TdmToTensorSet: public Logger,
                              public ContextBase,
                              public ITorchToTensorSet {

        TdmToTensorSet();
        virtual ~TdmToTensorSet();

        // IFunction
        virtual bool operator()(const input_pointer& in, output_pointer& out);
        
        // IConfigurable
        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

    private:
        TdmToTensorSetConfig m_config;
    };

}
