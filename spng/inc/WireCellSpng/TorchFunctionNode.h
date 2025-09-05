#ifndef WIRECELL_SPNG_TORCHFUNCTIONNODE
#define WIRECELL_SPNG_TORCHFUNCTIONNODE

#include "WireCellSpng/FunctionNode.h"
#include "WireCellSpng/ContextBase.h"

namespace WireCell::SPNG {

    /** A TorchFunctionNode adds torch device, semaphore and autograd guards.

        
        This is intended to be used as a base class for a subclass to implement
        transform_tensors() with code that calls Torch operators.  If your
        subclass does not require any Torch operators, inherit from FunctionNode
        instead.

        This class provides a device() method so that the subclass may know the
        user-configured torch device.  This is required when the subclass has no
        way to discern the device from input tensors when constructing new
        tensors.

        The transform_tensors() method is called in a contex that is guarded by
        a semaphore to limit the number of simultaneous tasks on a GPU and to
        assure torch autograd is disabled.  It also assures that all tensors
        provided to the method are on the device().


        This class is an IConfigurable and an ITorchTensorSetFilter.  These two
        types must be included in the list of interfaces passed to the
        subclass's call to the WIRECELL_FACTORY() CPP macro.
       
        If the subclass is itself configurable it must marshal the config object
        from this base class's default_configuration() and forward the config
        object to this base class's configure().  This classes methods provide
        examples of how this marshalling must be done.  See this class's
        IConfigurable implementation for an example of this marshaling of the
        configuration object.
     
     */
    class TorchFunctionNode : public FunctionNode, public ContextBase
    {
    public:

        /// Call this from subclass constructor to pass unique logname.
        TorchFunctionNode(const std::string& logname, const std::string& pkgnam="spng");
        virtual ~TorchFunctionNode() = default;

        // IConfigurable to marshal configuration to/from base classes.
        virtual WireCell::Configuration default_configuration() const;
        virtual void configure(const WireCell::Configuration& cfg);

        /// Override this "system" method from FunctionNode in order to call
        /// transform_tensors() in a context guarded by the torch semaphore and
        /// no-auto-grad-mode.  
        virtual TensorIndex sys_transform_tensors(TensorIndex ti) const;
        
        /// Subclass should override:
        // virtual TensorIndex transform_tensors(TensorIndex ti) const;


    };
}

#endif
