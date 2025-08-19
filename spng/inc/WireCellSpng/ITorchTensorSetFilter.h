#ifndef WIRECELL_ITORCHTENSORSETFILTER
#define WIRECELL_ITORCHTENSORSETFILTER

#include "WireCellUtil/IComponent.h"
#include "WireCellIface/IFunctionNode.h"
#include "WireCellSpng/ITorchTensorSet.h"

namespace WireCell {

    /** A torch tensor filter is something that applies some transformation
     * on its input torch tensor to produce and output torch tennsor.
     */
    class ITorchTensorSetFilter : public IFunctionNode<ITorchTensorSet, ITorchTensorSet> {
       public:
        typedef std::shared_ptr<ITorchTensorSetFilter> pointer;

        virtual ~ITorchTensorSetFilter() {};

        virtual std::string signature() { return typeid(ITorchTensorSetFilter).name(); }

        // supply:
        // virtual bool operator()(const input_pointer& in, output_pointer& out);
    };

}  // namespace WireCell

#endif