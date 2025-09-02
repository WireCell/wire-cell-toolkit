/** The SPNG::BaseNode class bundles some common behavior for SPNG data flow
 * graph nodes.
 *
 * This includes
 * - A TorchContext for torch device and semaphore.
 * - Facilities to work with ITorchTensorSet
 *
 * See SPNG::FunctionNode and SPNG::{Fanout,Fanin}Node and perhaps other classes
 * that include SPNG::BaseNode and provide additional functionality as base
 * classes.
 *
 * The SPNG::BaseNode is an IConfigurable.  When used as a base for your own
 * component class, be sure to include IConfigurable in the list of interfaces
 * passed to the `WIRECELL_FACTORY()` macro all.  If your subclass requires
 * IConfigurable itself, be sure to forward configuration objects to/from the
 * base.
*/

#ifndef WIRECELL_SPNG_BASENODE
#define WIRECELL_SPNG_BASENODE


#include "WireCellUtil/Configuration.h"
#include "WireCellSpng/TorchContext.h"

namespace WireCell::SPNG {

    class BaseNode {
    public:

        BaseNode() = default;
        virtual ~BaseNode() = default;

        /// IConfigurable.

        /// Configurable subclasses should call this method from their own
        /// default_configuration() and use update(cfgA, cfgB) to merge.
        virtual WireCell::Configuration default_configuration() const;

        /// Configurable subclasses should call this method from their own
        /// configurure() to forward the cfg object.
        virtual void configure(const WireCell::Configuration& cfg);

        /// Configuration:
        ///
        /// - device :: optional string giving the torch device to use.
        ///             Defaults to "cpu".  See TorchContext.
        ///
        /// - semaphore :: optional string giving the semaphore name.
        ///                Defaults to none, derive from "device".  See TorchContext.

        const TorchContext& context() const { return m_ctx; }


    protected:

        /// Subclass may use the torch context like:
        ///
        /// void mymeth() {
        ///   TorchSemaphore sem(m_ctx);
        ///   ... code using a torch device
        /// } // end of scope
        TorchContext m_ctx;

    };
}

#endif

