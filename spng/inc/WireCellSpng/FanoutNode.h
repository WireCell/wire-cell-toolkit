/** This provides FanoutNode operating on SPNG TDM tensor sets.

    See the datamodel.org document for more information on the SPNG Tensor Data
    Model and the C++ class support.

 */
#ifndef WIRECELL_SPNG_FANOUTNODE
#define WIRECELL_SPNG_FANOUTNODE

#include "WireCellSpng/ITorchTensorSetFanout.h"
#include "WireCellSpng/FanBase.h"

namespace WireCell::SPNG {

    /// Provide a general ITorchTensorSet FanoutNode.
    ///
    /// The is operates on TDM-compliant and non-TDM ITorchTensorSets.
    ///
    /// This class is intended to be used as-is to provide general purpose
    /// fanout semantics.
    ///
    /// See FanBase for configuration parameters.
    /// 
    /// It is possible but discouraged to use this class as a base and override
    /// the separate_tensors() method to provide novel behavior.  Developers are
    /// urged instead to provide the novel behavior in the form of a
    /// FunctionNode or a TorchFunctionNode to operate on pre/post fanned tensor
    /// sets.  If used as a base class, there are various requirements that must
    /// be satisfied by the subclass.
    class FanoutNode : public FanBase, public ITorchTensorSetFanout {
    public:
        FanoutNode(const std::string& logname="SPNGFanoutNode", const std::string& pkgnam="spng");
        virtual ~FanoutNode() = default;

        // INode, override because we get multiplicity at run time.
        // FIXME: sigh, all the way up into iface, this method should be const!
        virtual std::vector<std::string> output_types();

        // IFanout
        virtual bool operator()(const ITorchTensorSet::pointer& in, ITorchTensorSet::vector& outv);

        /// Separate the set into multiple a "multiplicity" number of sets.
        ///
        /// The default implementation copy the input pointer to each output port.
        ///
        /// If a more rich separation is required, attempt to use preceding or
        /// following FunctionNodes to perform filters.
        ///
        /// Last resort, this class may be inherited and this method overridden.
        virtual ITorchTensorSet::vector separate_tensors(const ITorchTensorSet::pointer& in) const;
        virtual ITorchTensorSet::vector sys_separate_tensors(const ITorchTensorSet::pointer& in) const;


    private:

        /// Configuration: See FanBase for more.

    };
}

#endif
