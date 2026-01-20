#pragma once

#include "WireCellSpng/Logger.h"
#include "WireCellSpng/ContextBase.h"
#include "WireCellSpng/ITorchTensorSetFilter.h"

#include "WireCellUtil/HanaJsonCPP.h"

namespace WireCell::SPNG {

    /// A TensorSetFilter provides a base class for tensor set filters.
    ///
    /// It provides "standard" logger and torch context functionality.
    struct TensorSetFilter: public ContextBase,
                            public Logger,
                            virtual public ITorchTensorSetFilter,
                            virtual public IConfigurable {

        TensorSetFilter(const std::string& type_name, const std::string& group_name="spng");
        virtual ~TensorSetFilter() = default;

        /// ITorchTensorSetFilter
        virtual bool operator()(const input_pointer& in, output_pointer& out);

        // IConfigurable
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;
        
        /// Subclass must implement.
        ///
        /// This base class takes care of EOS and this method is only called if
        /// "in" is not nullptr.  This class also assures the context base
        /// semaphore.
        virtual ITorchTensor::pointer filter_tensor(const ITorchTensor::pointer& in) = 0;


    };

}
