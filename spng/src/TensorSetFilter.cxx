#include "WireCellSpng/TensorSetFilter.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"
#include "WireCellSpng/HanaConfigurable.h"

using namespace WireCell::HanaJsonCPP;

namespace WireCell::SPNG {

    TensorSetFilter::TensorSetFilter(const std::string& type_name, const std::string& group_name)
        : Logger(type_name, group_name)
    {
    }

    bool TensorSetFilter::operator()(const input_pointer& in, output_pointer& out)
    {
        out = nullptr;
        if (!in) {
            logit("EOS");
            next_count();
            return true;
        }

        logit(in, "input");

        TorchSemaphore sem(context());

        out = filter_tensor(in);

        logit(out, "output");
        next_count();
        return true;
    }

    void TensorSetFilter::configure(const WireCell::Configuration& config)
    {
        log->debug("{}", config);
        WireCell::configure_bases<TensorSetFilter, ContextBase, Logger>(this, config);
    }

    WireCell::Configuration TensorSetFilter::default_configuration() const
    {
        return WireCell::default_configuration_bases<TensorSetFilter, ContextBase, Logger>(this);
    }
}
