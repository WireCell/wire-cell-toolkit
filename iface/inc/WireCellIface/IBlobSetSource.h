#ifndef WIRECELL_IBLOBSETSOURCE
#define WIRECELL_IBLOBSETSOURCE

#include "WireCellIface/ISourceNode.h"
#include "WireCellIface/IBlobSet.h"

namespace WireCell {

    class IBlobSetSource : public ISourceNode<IBlobSet> {
       public:
        typedef std::shared_ptr<IBlobSetSource> pointer;

        virtual ~IBlobSetSource() {}

        virtual std::string signature() { return typeid(IBlobSetSource).name(); }

        // supply:
        // virtual bool operator()(IBlobSet::pointer& cluster);
    };
}  // namespace WireCell

#endif
