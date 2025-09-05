/** This provides FanoutNode operating on SPNG TDM tensor sets.

    See the datamodel.org document for more information on the SPNG Tensor Data
    Model and the C++ class support.

 */
#ifndef WIRECELL_SPNG_FANOUTNODE
#define WIRECELL_SPNG_FANOUTNODE

#include "WireCellIface/IConfigurable.h"
#include "WireCellSpng/ITorchTensorSetFanout.h"

namespace WireCell::SPNG {

    /// Provide a general ITorchTensorSet FanoutNode.
    ///
    /// The input and output ITorchTensorSets adhere to the SPNG tensor data
    /// model (TDM).
    ///
    /// This class is intended to be used as-is to provide general purpose
    /// fanout semantics.  It is possible to use this class as a base and
    /// override the separate_tensors() method to provide novel behavior.
    /// However, developers are urged to provide the novel behavior by using or
    /// developing an FunctionNode or a TorchFunctionNode to operate on pre/post
    /// fanned tensor sets.
    class FanoutNode : public IConfigurable, public ITorchTensorSetFanout {
    public:
        FanoutNode() = default;
        virtual ~FanoutNode() = default;

        // IConfigurable
        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        // INode, override because we get multiplicity at run time.
        // FIXME: sigh, all the way up into iface, this method should be const!
        virtual std::vector<std::string> output_types();

        // IFanout
        virtual bool operator()(const ITorchTensorSet::pointer& in, ITorchTensorSet::vector& outv) const;

        /// Separate the set into multiple a "multiplicity" number of sets.  The
        /// default implementation is to simply copy the input pointer to all
        /// outputs.  If more rich separation is required, this method may be
        /// overridden in a base class.  However, developers are urged to
        /// instead follow each fan out output port with standard FunctionNode
        /// to select and/or rename desired subset of the tensor set.
        virtual ITorchTensorSet::vector separate_tensors(const ITorchTensorSet::pointer& in) const;


        /// Configuration
        ///
        /// - multiplicity :: the fan-out multiplicity.

    private:
        int m_multiplicity{2};

    };
}

#endif
