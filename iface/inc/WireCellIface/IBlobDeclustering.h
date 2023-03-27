/**  A blob declustering consumes an ICluster and produces an
 *  IBlobSet.
 */

#ifndef WIRECELL_IBLOBDECLUSTERING
#define WIRECELL_IBLOBDECLUSTERING

#include "WireCellIface/IFunctionNode.h"
#include "WireCellIface/ICluster.h"
#include "WireCellIface/IBlobSet.h"

namespace WireCell {

    class IBlobDeclustering : public IFunctionNode<ICluster, IBlobSet> {

    public:
        typedef std::shared_ptr<IBlobDeclustering> pointer;

        virtual ~IBlobDeclustering() {}

        virtual std::string signature() { return typeid(IBlobDeclustering).name(); }

        // supply:
        // virtual bool operator()(const input_pointer& in, output_pointer& out);
    };
}  // namespace WireCell

#endif
