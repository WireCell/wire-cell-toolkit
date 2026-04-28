/** BlobLabeling performs labeling operations on IBlobSet.
 * 
 * This component takes an IBlobSet as input and produces a labeled IBlobSet as output.
 * Currently, input and output are the same (pass-through behavior).
 */

#ifndef WIRECELLIMG_BLOBLABELING
#define WIRECELLIMG_BLOBLABELING

#include "WireCellIface/IFunctionNode.h"
#include "WireCellIface/IBlobSet.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"

namespace WireCell {
    namespace Img {

        class BlobLabeling : public Aux::Logger, public IFunctionNode<IBlobSet, IBlobSet>, public IConfigurable {
           public:
            BlobLabeling();
            virtual ~BlobLabeling();

            // IConfigurable
            virtual void configure(const WireCell::Configuration& cfg);
            virtual WireCell::Configuration default_configuration() const;

            // IFunctionNode
            virtual bool operator()(const input_pointer& in, output_pointer& out);

           private:
            // Configuration parameters can be added here as needed
        };

    }  // namespace Img
}  // namespace WireCell

#endif