/** This provides FunctionNode operating on SPNG TDM tensor sets.

    See the datamodel.org document for more information on the SPNG Tensor Data
    Model and the C++ class support.

 */

#ifndef WIRECELL_SPNG_FUNCTIONNODE
#define WIRECELL_SPNG_FUNCTIONNODE

#include "WireCellSpng/ITorchTensorSetFilter.h"
#include "WireCellSpng/TensorSelector.h"
#include "WireCellSpng/TensorRenaming.h"
#include "WireCellSpng/TensorIndex.h"
#include "WireCellIface/IConfigurable.h"

namespace WireCell::SPNG {

    /** Provide standard selection and datapath renaming for tensor sets and a
     * base class for more rich operations.
     *
     * In particular the transform_tensors() method may be overridden in a base
     * class to operate on selected tensors.  This should be done only in the
     * case that the base class will NOT use any Torch operations.  To safely
     * apply Torch operations, use the TorchFunctionNode as your base class.
     * 
     * This class is an IConfigurable and an ITorchTensorSetFilter.  If used as
     * a base class for your own data flow graph node, these two types must be
     * included in the list of interfaces passed to your WIRECELL_FACTORY() CPP
     *
     * If the subclass is itself configurable it must marshal the config object
     * from this base class's default_configuration() and forward the config
     * object to this base class's configure().  This classes methods provide
     * examples of how this marshalling must be done.
     */       
    class FunctionNode : public IConfigurable, public WireCell::ITorchTensorSetFilter { 
    public:
        FunctionNode() = default;
        virtual ~FunctionNode() = default;

        // IConfigurable
        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        // ITorchTensorSetFilter
        virtual bool operator()(const ITorchTensorSet::pointer& in, output_pointer& out) const;

        /// The FunctionNode API.  These are called by the operator() in the
        /// order given here.  Any may be overridden by a subclass.  Such
        /// overrides should likely still call these base class methods in order
        /// to retain standard behaviors.

        /// Index the input tensors.
        virtual TensorIndex index_tensors(const ITorchTensorSet::pointer& in) const;

        /// Apply standard and configurable selection rules.
        virtual TensorIndex select_tensors(TensorIndex ti) const;
        
        /// Apply a transform.  Here, this is a no-op.  It may be overridden by
        /// subclass.
        virtual TensorIndex transform_tensors(TensorIndex ti) const;

        // This method is NOT intended for user override but is overriden by
        // TorchFunctionNode to establish the torch context prior to calling
        // transform_tensors().
        virtual TensorIndex sys_transform_tensors(TensorIndex ti) const;

        /// Apply standard and configurable rename rules.
        virtual TensorIndex rename_tensors(TensorIndex ti) const;

        /// Copy the indexed tensors into a tensor set.  
        ITorchTensorSet::pointer pack_tensors(TensorIndex ti) const;


    protected:

        /// Configuration: "tensor_selector".  See TensorSelector class for details.
        TensorSelector m_selector;
        /// Configuration: "tensor_renaming".  See TensorRenaming class for details.
        TensorRenaming m_renaming;

    };
}

#endif
