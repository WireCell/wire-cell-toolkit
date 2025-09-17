/** The SPNG::ContextBase provides configurable torch device and semaphore.
 *
 * It is expected that every SPNG data flow graph node that performs torch
 * operations will inherit from ContexBase.
 *
 * The ContexBase is an IConfigurable.  Any subclass that is itself a component
 * must include IConfigurable interface in the subclass's WIRECELL_FACTORY()
 * macro call.
 *
 * A subclass that is itself configurable must combine the result of the
 * ContextBase's default_configuration() return to its own and must forward the
 * configuration object received by its own configure() method to that method on
 * the ContextBase.
*/

#ifndef WIRECELL_SPNG_CONTEXTBASE
#define WIRECELL_SPNG_CONTEXTBASE


#include "WireCellIface/IConfigurable.h"
#include "WireCellSpng/TorchContext.h"

namespace WireCell::SPNG {

    class ContextBase : public virtual IConfigurable {
    public:

        ContextBase() = default;
        virtual ~ContextBase() = default;

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

        /// Access the device and the semaphore.
        ///
        /// The subclass likely should use the semaphore by via TorchSemaphore
        /// wrapper that lives for the duration of a code-scope.
        ///
        /// void mymeth() {
        ///   TorchSemaphore sem(context());
        ///   ... code using a torch device
        /// } // end of scope
        const TorchContext& context() const { return m_ctx; }

        /// Shorthand method to get the context's device.
        torch::Device device() const { return m_ctx.device(); }


    protected:

        TorchContext m_ctx;

    };
}

#endif

