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
#include "WireCellAux/Logger.h"
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
     * When FunctionNode is used as a base class, assure the following:
     *
     * 1) Include this list of interfaces in the subclass's WIRECELL_FACTORY().
     *
     *     INamed, IConfigurable, ITorchTensorSetFilter
     *
     * 2) Have the subclass constructor call FunctionNode("myname") constructor
     * to set a custom logging name.
     *
     * 3) If the subclass is itself an IConfigurable, it must marshal the
     * configuration object to/from FunctionNode::IConfigurable methods.  See
     * implementation of these methods in FunctionNode for examples how to do
     * this.
     */       
    class FunctionNode : public Aux::Logger, public IConfigurable, public WireCell::ITorchTensorSetFilter { 
    public:

        // Subclass SHOULD call this to provide a subclass-specific log name
        FunctionNode(const std::string& logname="function", const std::string& pkgnam="spng");
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

        /// Apply standard and configurable rename rules.
        virtual TensorIndex rename_tensors(TensorIndex ti) const;

        /// Copy the indexed tensors into a tensor set.  
        virtual ITorchTensorSet::pointer pack_tensors(TensorIndex ti) const;

        // The "sys_*" methods call their like-named siblings.  These should NOT
        // be overridden by subclasses (with a special exception for
        // TorchFunctionNode providing sys_transform_tensors).  They are used to
        // "inject" code to for regardless of subclass overrides of non-sys_*
        // methods.
        virtual TensorIndex sys_index_tensors(const ITorchTensorSet::pointer& in) const;
        virtual TensorIndex sys_transform_tensors(TensorIndex ti) const;
        virtual ITorchTensorSet::pointer sys_pack_tensors(TensorIndex ti) const;


    protected:

        /// Configuration: "tensor_selector".  See TensorSelector class for details.
        TensorSelector m_selector;
        /// Configuration: "tensor_renaming".  See TensorRenaming class for details.
        TensorRenaming m_renaming;

        /// Emit standard log line for the state of the tensor index.
        virtual void maybe_log(const TensorIndex& ti, const std::string& context="") const;

        /// Configuration: quiet.  Set true to not call any logging.  Default is false.
        bool m_quiet{false};

        /// Keep track of calls
        mutable size_t m_count{0};

    };
}

#endif
