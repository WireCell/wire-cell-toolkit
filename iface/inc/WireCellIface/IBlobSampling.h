/**  A blob sampling consumes a set of blobs and produces
 *  a tensor set of point clouds datasets sampling the blobs.
 */

#ifndef WIRECELL_IBLOBSAMPLING
#define WIRECELL_IBLOBSAMPLING

#include "WireCellIface/IFunctionNode.h"
#include "WireCellIface/ITensorSet.h"
#include "WireCellIface/IBlobSet.h"

namespace WireCell {

    class IBlobSampling : public IFunctionNode<IBlobSet, ITensorSet> {

    public:
        typedef std::shared_ptr<IBlobSampling> pointer;

        virtual ~IBlobSampling() {}

        virtual std::string signature() { return typeid(IBlobSampling).name(); }

        // supply:
        // virtual bool operator()(const input_pointer& in, output_pointer& out);
    };
}  // namespace WireCell

#endif
