/** ITruthPointLabeler joins a point-cloud tree tensor set with the
 * cluster it was built from and the (drifted) depos that made the
 * cluster, and produces a tensor set.
 *
 * Used by truth labelling to attach per-point simulation truth to the
 * point tree (a local point cloud per blob node) without touching the
 * existing point clouds.
 */

#ifndef WIRECELL_ITRUTHPOINTLABELER
#define WIRECELL_ITRUTHPOINTLABELER

#include "WireCellUtil/IComponent.h"
#include "WireCellIface/IJoinNode.h"
#include "WireCellIface/ICluster.h"
#include "WireCellIface/IDepoSet.h"
#include "WireCellIface/ITensorSet.h"

namespace WireCell {

    class ITruthPointLabeler : public IJoinNode<std::tuple<ITensorSet, ICluster, IDepoSet>, ITensorSet> {
       public:
        typedef std::shared_ptr<ITruthPointLabeler> pointer;

        virtual ~ITruthPointLabeler();

        virtual std::string signature() { return typeid(ITruthPointLabeler).name(); }

        // subclass supply:
        // virtual bool operator()(const input_tuple_type& intup,
        //                         output_pointer& out);
    };

}  // namespace WireCell

#endif
