#ifndef WIRECELL_SPNG_ITENSORTOTORCHSET
#define WIRECELL_SPNG_ITENSORTOTORCHSET

#include "WireCellIface/IFunctionNode.h"
#include "WireCellIface/ITensorSet.h"
#include "WireCellSpng/ITorchTensorSet.h"

namespace WireCell {
namespace SPNG {
    /** A component takes 1 input ITensorSet and produces one
     * ITorchTensorSet.

     */
    class ITensorToTorchSet : public IFunctionNode<ITensorSet, ITorchTensorSet> {
       public:
        virtual ~ITensorToTorchSet() {};

        virtual std::string signature() { return typeid(ITensorToTorchSet).name(); }

        // Subclass must implement:
        // virtual std::vector<std::string> output_types() = 0;
        // and the already abstract:
        // virtual bool operator()(const input_pointer& in, output_vector& outv);
    };
}  // namespace SPNG
}  // namespace WireCell

#endif