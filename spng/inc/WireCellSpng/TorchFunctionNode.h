#ifndef WIRECELL_SPNG_TORCHFUNCTIONNODE
#define WIRECELL_SPNG_TORCHFUNCTIONNODE

#include "WireCellSpng/FunctionNode.h"
#include "WireCellSpng/ContextBase.h"

namespace WireCell::SPNG {

    /** A TorchFunctionNode adds torch device, semaphore and autograd guards.
        
        This class is not particularly useful stand-alone (use FunctionNode
        instead) and is intended to be used as a base class in order to
        implement transform_tensors() with code that calls Torch operators.

        If implementing transform_tensors() without any Torch operators inherit
        from FunctionNode instead.

        This class provides a device() method so that the subclass may know the
        user-configured torch device.  This is required when the subclass has no
        way to discern the device from input tensors when constructing new
        tensors.

        The transform_tensors() method is called in a contex that is guarded by
        a semaphore to limit the number of simultaneous tasks on a GPU and to
        assure torch autograd is disabled.  It also assures that all tensors
        provided to the method are on the device().

        Your subclass must assure the following:
        
        1) Include this list of interfaces in the subclass's WIRECELL_FACTORY().
     
            INamed, IConfigurable, ITorchTensorSetFilter
     
        2) Have the subclass constructor call `TorchFunctionNode("myname")`
        constructor to set a custom logging name.
     
        3) If the subclass is itself an IConfigurable, it must marshal the
        configuration object to/from FunctionNode::IConfigurable methods.  See
        implementation of these methods in FunctionNode for examples how to do
        this.


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
