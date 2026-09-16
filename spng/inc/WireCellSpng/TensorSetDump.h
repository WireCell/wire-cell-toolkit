#pragma once

#include "WireCellSpng/Logger.h"
#include "WireCellSpng/ITorchTensorSetSink.h"
#include "WireCellUtil/HanaJsonCPP.h"

// fixme: should factor out meaning from names.
namespace WireCell::SPNG {

    struct TensorSetDump : public Logger,
                                 public virtual IConfigurable,
                                 public virtual ITorchTensorSetSink {
        TensorSetDump();
        virtual ~TensorSetDump() = default;

        // IConfigurable API.
        virtual void configure(const WireCell::Configuration& jconfig);
        virtual WireCell::Configuration default_configuration() const;

        // ITorchTensorSetSink API.
        virtual bool operator()(const ITorchTensorSet::pointer& in);

    };
}
