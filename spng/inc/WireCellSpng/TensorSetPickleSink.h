#pragma once

#include "WireCellSpng/Logger.h"
#include "WireCellSpng/ITorchTensorSetSink.h"
#include "WireCellUtil/HanaJsonCPP.h"

// fixme: should factor out meaning from names.
namespace WireCell::SPNG {
    struct FilenameConfig {
        std::string filename;
    };
}
BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::FilenameConfig, filename);

namespace WireCell::SPNG {

    struct TensorSetPickleSink : public Logger,
                                 public virtual IConfigurable,
                                 public virtual ITorchTensorSetSink {
        TensorSetPickleSink();
        virtual ~TensorSetPickleSink() = default;

        // IConfigurable API.
        virtual void configure(const WireCell::Configuration& jconfig);
        virtual WireCell::Configuration default_configuration() const;

        // ITorchTensorSetSink API.
        virtual bool operator()(const ITorchTensorSet::pointer& in);

    private:
        FilenameConfig m_config;
        std::ofstream m_output;

    };
}
