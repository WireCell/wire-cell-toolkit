#pragma once

#include "WireCellSpng/Logger.h"
#include "WireCellSpng/ContextBase.h"
#include "WireCellSpng/ITorchTensorFilter.h"

#include "WireCellUtil/HanaJsonCPP.h"

namespace WireCell::SPNG {

    struct TensorFilterConfig {

        /// Set the "tag" metadata attribute on the produced tensor.
        std::string tag="";

        /// Apply a format string to the metadata to produce a datapath for the
        /// output tensor.  If not provided (default) the produced tensor will
        /// retain the same datapath as the input tensor.  If provided, the
        /// input datapath is stored as the "derived_from" attribute.
        std::string datapath_format="";

        /// If set, save output tensor to the named file.  File is produced by
        /// torch::pickle_save().
        std::string debug_filename = "";
        
    };

}
BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::TensorFilterConfig,
                        tag, datapath_format, debug_filename);

namespace WireCell::SPNG {

    struct TensorFilter: public ContextBase,
                         public Logger,
                         virtual public ITorchTensorFilter,
                         virtual public IConfigurable {

        TensorFilter(const std::string& type_name, const std::string& group_name="spng");
        virtual ~TensorFilter() = default;

        /// ITorchTensorFilter
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

        /// Add a tensor to be saved in the debug file, if one was named.  Subclass may call.
        void maybe_save(torch::Tensor tensor, std::string name);


    private:
        TensorFilterConfig m_config;

        using tensor_map = torch::Dict<std::string, torch::Tensor>;
        tensor_map m_to_save;

    };

}
