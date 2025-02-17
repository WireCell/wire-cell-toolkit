/** Some tools that operate on anode plane related interfaces.
 */

#ifndef WIRECELL_PLANETOOLS
#define WIRECELL_PLANETOOLS

#include "WireCellIface/IAnodePlane.h"

namespace WireCell {
    namespace Aux {

        /// Return an ordered vector of channels for the given wire
        /// plane index (in [0,1,2]).  Order is "face major" in that
        /// a two-faced anode will have face 0 channels first
        /// followed by face 1 channels.  W/in a face the order is
        /// determined by the order of channels w/in the face.
        IChannel::vector plane_channels(IAnodePlane::pointer anode,
                                        int wire_plane_index);

        struct WirePlaneInfo {
            int start_index;
            int end_index;
            int total_wires;
        };
        // Get wire information for a specific plane from IAnodeFace
        WirePlaneInfo get_wire_plane_info(IAnodeFace::pointer face, WirePlaneLayer_t layer);


    }
}
#endif
