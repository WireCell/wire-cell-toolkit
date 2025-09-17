/** This provides FaninNode operating on SPNG TDM tensor sets.

    See the datamodel.org document for more information on the SPNG Tensor Data
    Model and the C++ class support.

 */
#ifndef WIRECELL_SPNG_FANOUTNODE
#define WIRECELL_SPNG_FANOUTNODE

#include "WireCellSpng/ITorchTensorSetFanin.h"
#include "WireCellSpng/FanBase.h"

namespace WireCell::SPNG {

    /// Provide a general ITorchTensorSet FaninNode.
    ///
    /// The is operates on TDM-compliant and non-TDM ITorchTensorSets.
    ///
    /// This class is intended to be used as-is to provide general purpose
    /// fanin semantics.
    ///
    /// See FanBase for configuration parameters.
    ///
    /// It is possible but discouraged to use this class as a base and override
    /// the separate_tensors() method to provide novel behavior.  Developers are
    /// urged instead to provide the novel behavior in the form of a
    /// FunctionNode or a TorchFunctionNode to operate on pre/post fanned tensor
    /// sets.  If used as a base class, there are various requirements that must
    /// be satisfied by the subclass.
    class FaninNode : public FanBase, public ITorchTensorSetFanin {
    public:
        FaninNode(const std::string& logname="SPNGFaninNode", const std::string& pkgnam="spng");
        virtual ~FaninNode() = default;

        // INode, override because we get multiplicity at run time.
        // FIXME: sigh, all the way up into iface, this method should be const!
        virtual std::vector<std::string> input_types();

        // IFanin
        virtual bool operator()(const ITorchTensorSet::vector& inv, ITorchTensorSet::pointer& out);

        /// Combine "multiplicity" number of sets into one.
        ///
        /// The default implementation is to form the ordered union.
        ///
        /// If more rich combination is required, attempt to use preceding or
        /// following FunctionNodes to perform filters.
        ///
        /// Last resort, this class may be inherited and this method overridden
        virtual ITorchTensorSet::pointer combine_tensors(const  ITorchTensorSet::vector& inv) const;
        virtual ITorchTensorSet::pointer sys_combine_tensors(const  ITorchTensorSet::vector& inv) const;


    };
}

#endif
