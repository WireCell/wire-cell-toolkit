#ifndef WIRECELL_ITORCHTENSORFILTER
#define WIRECELL_ITORCHTENSORFILTER

#include "WireCellUtil/IComponent.h"
#include "WireCellIface/IFunctionNode.h"
#include "WireCellSpng/ITorchTensor.h"

namespace WireCell {

    /** A torch tensor filter is something that applies some transformation
     * on its input torch tensor to produce and output torch tennsor.
     */
    class ITorchTensorFilter : public IFunctionNode<ITorchTensor, ITorchTensor> {
       public:
        typedef std::shared_ptr<ITorchTensorFilter> pointer;

        virtual ~ITorchTensorFilter() {};

        virtual std::string signature() { return typeid(ITorchTensorFilter).name(); }

        // supply:
        // virtual bool operator()(const input_pointer& in, output_pointer& out);
    };

}  // namespace WireCell

#endif