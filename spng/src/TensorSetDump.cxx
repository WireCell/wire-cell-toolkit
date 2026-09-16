#include "WireCellSpng/TensorSetDump.h"
#include "WireCellSpng/Torch.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(SPNGTensorSetDump,
                 WireCell::SPNG::TensorSetDump,
                 WireCell::ITorchTensorSetSink,
                 WireCell::IConfigurable,
                 WireCell::INamed);


namespace WireCell::SPNG {

    TensorSetDump::TensorSetDump()
        : Logger("TensorSetDump", "spng")
    {
    }

    void TensorSetDump::configure(const WireCell::Configuration& jconfig)
    {
        this->Logger::configure(jconfig);
    }

    WireCell::Configuration TensorSetDump::default_configuration() const
    {
        return this->Logger::default_configuration();
    }

    bool TensorSetDump::operator()(const ITorchTensorSet::pointer& in)
    {
        if (!in) {
            logit("EOS");
            next_count();
            return true;
        }
        logit(in, "sinking");

        next_count();
        return true;
    }


}
