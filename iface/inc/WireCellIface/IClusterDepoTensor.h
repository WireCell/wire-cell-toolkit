/** IClusterDepoTensor inputs an ICluster and an IDepoSet to produce an
 * ITensorSet.
 *
 * The join of a cluster with the (drifted) depos that made it: used
 * by truth labelling to attach per-blob simulation truth (which depo,
 * which G4 track, how much charge) to the cluster as tensors.
 */

#ifndef WIRECELL_ICLUSTERDEPOTENSOR
#define WIRECELL_ICLUSTERDEPOTENSOR

#include "WireCellUtil/IComponent.h"
#include "WireCellIface/IJoinNode.h"
#include "WireCellIface/ICluster.h"
#include "WireCellIface/IDepoSet.h"
#include "WireCellIface/ITensorSet.h"

namespace WireCell {

    class IClusterDepoTensor : public IJoinNode<std::tuple<ICluster, IDepoSet>, ITensorSet> {
       public:
        typedef std::shared_ptr<IClusterDepoTensor> pointer;

        virtual ~IClusterDepoTensor();

        virtual std::string signature() { return typeid(IClusterDepoTensor).name(); }

        // subclass supply:
        // virtual bool operator()(const input_tuple_type& intup,
        //                         output_pointer& out);
    };

}  // namespace WireCell

#endif
