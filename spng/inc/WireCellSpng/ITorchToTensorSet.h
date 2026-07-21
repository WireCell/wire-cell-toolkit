#ifndef WIRECELL_SPNG_ITORCHTOTENSORSET
#define WIRECELL_SPNG_ITORCHTOTENSORSET

#include "WireCellIface/IFunctionNode.h"
#include "WireCellIface/ITensorSet.h"
#include "WireCellSpng/ITorchTensorSet.h"

namespace WireCell {
namespace SPNG {
    /** A component takes 1 input ITensorSet and produces one
     * ITorchTensorSet.

     */
    class ITorchToTensorSet : public IFunctionNode<ITorchTensorSet, ITensorSet> {
       public:
        virtual ~ITorchToTensorSet() {};

        virtual std::string signature() { return typeid(ITorchToTensorSet).name(); }

        // Subclass must implement:
        // virtual std::vector<std::string> output_types() = 0;
        // and the already abstract:
        // virtual bool operator()(const input_pointer& in, output_vector& outv);
    };
}  // namespace SPNG
}  // namespace WireCell

#endif