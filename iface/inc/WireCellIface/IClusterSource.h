#ifndef WIRECELL_ICLUSTERSOURCE
#define WIRECELL_ICLUSTERSOURCE

#include "WireCellIface/ISourceNode.h"
#include "WireCellIface/ICluster.h"

namespace WireCell {

    class IClusterSource : public ISourceNode<ICluster> {
       public:
        typedef std::shared_ptr<IClusterSource> pointer;

        virtual ~IClusterSource() {}

        virtual std::string signature() { return typeid(IClusterSource).name(); }

        // supply:
        // virtual bool operator()(ICluster::pointer& cluster);
    };
}  // namespace WireCell

#endif
