/** Consume IBlobSet until EOS and then emit an ITensorSet 
 */

#ifndef WIRECELL_ICLUSTERING
#define WIRECELL_ICLUSTERING

#include "WireCellIface/IQueuedoutNode.h"
#include "WireCellIface/IBlobSet.h"
#include "WireCellIface/ITensorSet.h"

namespace WireCell {

    class IBlobTensoring : public IQueuedoutNode<IBlobSet, ITensorSet> {
       public:
        typedef std::shared_ptr<IBlobTensoring> pointer;

        virtual ~IBlobTensoring() {}

        virtual std::string signature() { return typeid(IBlobTensoring).name(); }

        /// supply:
        // virtual bool operator()(const input_pointer& in, output_queue& outq) = 0;
    };
}  // namespace WireCell

#endif
