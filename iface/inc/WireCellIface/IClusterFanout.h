#ifndef WIRECELL_IFACE_ICLUSTERFANOUT
#define WIRECELL_IFACE_ICLUSTERFANOUT

#include "WireCellIface/IFanoutNode.h"
#include "WireCellIface/ICluster.h"

namespace WireCell {
    /** A cluster fan-out component takes 1 input cluster and produces one
     * cluster on each of its output ports.
     */
    class IClusterFanout : public IFanoutNode<ICluster, ICluster, 0> {
       public:
        virtual ~IClusterFanout();

        virtual std::string signature() { return typeid(IClusterFanout).name(); }

        // Subclass must implement:
        virtual std::vector<std::string> output_types() = 0;
        // and the already abstract:
        // virtual bool operator()(const input_pointer& in, output_vector& outv);
    };
}  // namespace WireCell

#endif
