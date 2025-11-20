#pragma once

#include "WireCellSpng/Logger.h"
#include "WireCellSpng/ContextBase.h"
#include "WireCellSpng/TorchLMN.h"
#include "WireCellSpng/ITorchTensorFilter.h"

#include "WireCellUtil/HanaJsonCPP.h"
#include "WireCellIface/IConfigurable.h"

namespace WireCell::SPNG {
    struct ResampleConfig {
        
    };
}

namespace WireCell::SPNG {

    class Resample : public ContextBase,
                     public Logger,
                     public ITorchTensorFilter,
                     virtual public IConfigurable {
    };
}
