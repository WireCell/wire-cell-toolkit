/** IBlobDepoFill inputs an ICluster and an IDepoSet to produce an
 * ICluster.
 */

#ifndef WIRECELL_IBLOBDEPOFILL
#define WIRECELL_IBLOBDEPOFILL

#include "WireCellUtil/IComponent.h"
#include "WireCellIface/IJoinNode.h"
#include "WireCellIface/ICluster.h"
#include "WireCellIface/IDepoSet.h"

namespace WireCell {

    class IBlobDepoFill : public IJoinNode<std::tuple<ICluster, IDepoSet>, ICluster> {
       public:
        typedef std::shared_ptr<IBlobDepoFill> pointer;

        virtual ~IBlobDepoFill();

        virtual std::string signature() { return typeid(IBlobDepoFill).name(); }

        // subclass supply:
        // virtual bool operator()(const input_tuple_type& intup,
        //                         output_pointer& out);
    };

}  // namespace WireCell

#endif
