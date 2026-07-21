#ifndef WIRECELL_ITORCHTENSORSETSINK
#define WIRECELL_ITORCHTENSORSETSINK

#include "WireCellSpng/ITorchTensorSet.h"
#include "WireCellIface/ISinkNode.h"
#include "WireCellUtil/IComponent.h"

namespace WireCell {

    /** A frame sink is something that generates IFrames.
     */
    class ITorchTensorSetSink : public ISinkNode<ITorchTensorSet> {
       public:
        typedef std::shared_ptr<ITorchTensorSetSink> pointer;

        virtual ~ITorchTensorSetSink(){};

        // subclass supply:
        // virtual bool operator()(const ITorchTensorSet::pointer& in);
    };

}  // namespace WireCell

#endif
