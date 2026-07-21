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
#include "WireCellSpng/Logger.h"
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
    class FunctionNode : public Logger,
                         public virtual IConfigurable,
                         public ITorchTensorSetFilter { 
    public:

        // Subclass SHOULD call this to provide a subclass-specific log name
        FunctionNode(const std::string& logname="SPNGFunctionNode", const std::string& pkgnam="spng");
        virtual ~FunctionNode() = default;

        // IConfigurable
        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        /// Subclass should call this method to register each "kind" of tensor
        /// it expects to find.  The "kind" string is defined by the subclass
        /// and will be reflected to an attribute of an object that may be
        /// provided by the user configuration parameter "input_tensors".
        ///
        /// The subclass should issue all of its register_input() calls prior to
        /// calling the base class default_configuration() in order that any
        /// default_datapaths can be added to the default configuration object.
        void register_input(const std::string& kind, const std::string& default_datapath="");

        /// After the subclass calls this base class's configure() it may call
        /// input_datapath() to learn the user's configured choice.
        ///
        /// The datapath may be used with the TensorIndex to get the tensor.
        ///
        /// Any error results in an empty string being returned.
        std::string input_datapath(const std::string& kind) const;

        /// Symmetric with register_input() / "output_tensors" configuration
        /// parameter.  See its comments above for details.
        void register_output(const std::string& kind, const std::string& default_datapath="");

        /// Symmetric with input_datapath().
        ///
        /// The datapath may be used to add a new tensor to the output
        /// TensorIndex.
        std::string output_datapath(const std::string& kind) const;


        // ITorchTensorSetFilter
        virtual bool operator()(const ITorchTensorSet::pointer& in, output_pointer& out);

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
        virtual TensorIndex sys_select_tensors(TensorIndex ti) const;
        virtual TensorIndex sys_transform_tensors(TensorIndex ti) const;
        virtual TensorIndex sys_rename_tensors(TensorIndex ti) const;
        virtual ITorchTensorSet::pointer sys_pack_tensors(TensorIndex ti) const;


    protected:

        /// Configuration: "tensor_selector".  See TensorSelector class for details.
        TensorSelector m_selector;
        /// Configuration: "tensor_renaming".  See TensorRenaming class for details.
        TensorRenaming m_renaming;

        /// Configuration: "keep_unselected" (default true)
        ///
        /// IF true, any tensors not explicitly selected or rejected are
        /// selected.  If false, they are rejected.
        bool m_keep_unselected{true};

        /// Configuration: "select_parents" (default true)
        ///
        /// If true, the selector is applied only to ultimate parents (children
        /// of the tree root node).  Any tree descendants will follow the
        /// parent.  If false, the selector is applied to all tensors.
        bool m_select_parents{true};

        /// Configuration: "combine_policy" (string, default="union_replace")
        ///
        /// Control how tensor index that is produced by the "transform" call
        /// is combined with the input tensor index.  Supported policies are:
        ///
        /// - union_replace :: a transformed tensor will replace an input tensor that has the same datapath.
        ///
        /// - union_keep :: an input tensor of the same datapath as a transformed tensor is kept.
        ///
        /// - input_only :: the transform is pure side-effect and the input tensor index is not modified.
        ///
        /// - produced_only :: keep only the tensors that survive selection and transformation.
        std::string m_combine_policy{"union_replace"};

        std::map<std::string, std::string> m_input_datapath, m_output_datapath;

    };
}

#endif
