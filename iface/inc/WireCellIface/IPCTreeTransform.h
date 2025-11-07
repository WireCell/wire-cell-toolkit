#ifndef WIRECELL_IPCTREETRANSFORM
#define WIRECELL_IPCTREETRANSFORM

#include "WireCellUtil/IComponent.h"
#include "WireCellUtil/PointTree.h"

namespace WireCell {

    /** An IPCTreeTransform accepts a point cloud tree node and produces a new
     * one.
     *
     * Caller takes ownership of the returned node via unique pointer.  
     *
     * Implementations should provide a contract as to what, if any, node facade
     * is assumed for both input and produced output nodes.
     */
    class IPCTreeTransform : public IComponent<IPCTreeTransform> {
    public:
        virtual ~IPCTreeTransform() {}

        using node_t = PointCloud::Tree::Points::node_t;

        virtual std::unique_ptr<node_t> transform(const node_t& node) const = 0;

    };

}  // namespace WireCell

#endif
