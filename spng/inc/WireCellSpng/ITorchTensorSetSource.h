#ifndef WIRECELL_ITORCHTENSORSETSOURCE
#define WIRECELL_ITORCHTENSORSETSOURCE

#include "ITorchTensorSet.h"
#include "WireCellIface/ISourceNode.h"
#include "WireCellUtil/IComponent.h"

namespace WireCell {

    /** A frame source is something that generates IFrames.
     */
    class ITorchTensorSetSource : public ISourceNode<ITorchTensorSet> {
       public:
        typedef std::shared_ptr<ITorchTensorSetSource> pointer;

        virtual ~ITorchTensorSetSource(){};

        // supply:
        // virtual bool operator()(ITorchTensorSet::pointer& frame);
    };

}  // namespace WireCell

#endif