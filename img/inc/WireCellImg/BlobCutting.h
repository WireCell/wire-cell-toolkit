/** BlobCutting performs cutting operations on IBlobSet.
 * 
 * This component takes an IBlobSet as input and produces a cut IBlobSet as output.
 * Currently, input and output are the same (pass-through behavior).
 */

#ifndef WIRECELLIMG_BLOBCUTTING
#define WIRECELLIMG_BLOBCUTTING

#include "WireCellIface/IFunctionNode.h"
#include "WireCellIface/IBlobSet.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"

namespace WireCell {
    namespace Img {

        class BlobCutting : public Aux::Logger, public IFunctionNode<IBlobSet, IBlobSet>, public IConfigurable {
           public:
            BlobCutting();
            virtual ~BlobCutting();

            // IConfigurable
            virtual void configure(const WireCell::Configuration& cfg);
            virtual WireCell::Configuration default_configuration() const;

            // IFunctionNode
            virtual bool operator()(const input_pointer& in, output_pointer& out);

           private:
            // Configuration parameters
            int m_length_threshold;   // Strip length threshold for cutting
            double m_nudge;          // Numerical precision parameter
            int m_max_depth;         // Maximum recursion depth to prevent stack overflow
        };

    }  // namespace Img
}  // namespace WireCell

#endif