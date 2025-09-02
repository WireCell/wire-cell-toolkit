/** An SPNG::FunctionNode provides a base class to SPNG ITorchTensorFilter nodes.
 *
 * It overrides:
 *
 *   operator()(ITorchTensorSet::pointer in, ITorchTensorSet::pointer& out)
 *
 * and calls:
 *
 *   operator()(TensorIndex& in, ITensorSet::pointer& out)
 *
 * The index will only contain parent and their children tensors for which the
 * parent tensor passed the configured tensor selector.
 *
 * This indexed operator is called in a semaphore and nograd protected context.
 *
 * A subclass must assure IConfigurable and ITorchTensorSet are in the list of
 * interfaces given to WIRECELL_FACTORY().
 *
 * If the subclass itself is configurable it must marshal the config object from
 * the base default_configuration() and provide the config object to the base
 * configure().  See this class's methods if you need an example of how to do
 * this forwarding of the config object.
 */

#ifndef WIRECELL_SPNG_FUNCTIONNODE
#define WIRECELL_SPNG_FUNCTIONNODE

#include "WireCellSpng/ITorchTensorSetFilter.h"
#include "WireCellSpng/ContextBase.h"
#include "WireCellSpng/TensorSelector.h"
#include "WireCellSpng/TensorIndex.h"
#include "WireCellIface/IConfigurable.h"

namespace WireCell::SPNG {

    class FunctionNode : public ContextBase, public TensorSelector, public WireCell::ITorchTensorSetFilter { 
    public:
        FunctionNode() = default;
        virtual ~FunctionNode() = default;

        // IConfigurable
        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        // ITorchTensorSetFilter
        virtual bool operator()(const input_pointer& in, output_pointer& out);

        virtual bool operator()(const TensorIndex& ti, output_pointer& out) = 0;

        /// Configuration: see TorchContext and TensorSelector.
        ///
        /// Any tensor selection is applied only to parent tensors (see SPNG
        /// torch data model).  Parent and children live and die together.
        ///
        /// If a parent fails to be explicitly accepted or rejected, it is
        /// selected (accepted).
        ///
        /// Currently, renaming is not supported.
        ///
        /// The input tensors be placed in the context's device() and device()
        /// may be used when a subclass needs to make a new tensor on the same
        /// device.  Generally, a subclass should NOT make tensors on any other
        /// device.

    };
}

#endif
