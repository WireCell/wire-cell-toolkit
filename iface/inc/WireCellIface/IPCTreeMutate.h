#ifndef WIRECELL_IPCTREEMUTATE
#define WIRECELL_IPCTREEMUTATE

#include "WireCellUtil/IComponent.h"
#include "WireCellUtil/PointTree.h"

namespace WireCell {

    /** An IPCTreeMutate accepts a point cloud tree node, may mutate it and may
     * produces a new one.
     *
     * Caller takes ownership of the returned node via unique pointer.  
     *
     */
    class IPCTreeMutate : public IComponent<IPCTreeMutate> {
    public:
        virtual ~IPCTreeMutate() {}

        using node_t = PointCloud::Tree::Points::node_t;

        virtual std::unique_ptr<node_t> mutate(node_t& node) const = 0;

    };

}  // namespace WireCell

#endif
