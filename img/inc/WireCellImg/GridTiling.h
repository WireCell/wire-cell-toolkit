/** This tiling algorithm makes use of RayGrid for the heavy lifting.
 *
 * It must be configured with a face on which to focus.  Only wires
 * (segments) in that face which are connected to channels with data
 * in a slice will be considered when tiling.
 *
 * It does not "know" about dead channels.  If your detector has them
 * you may place a component upstream which artifically inserts
 * non-zero signal on dead channels in slice input here.
 *
 * The resulting IBlobs have ident numbers which increment starting
 * from zero.  The ident is reset when EOS is received.
 */
#ifndef WIRECELLIM_GRIDTILING
#define WIRECELLIM_GRIDTILING

#include "WireCellIface/ITiling.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IAnodeFace.h"

#include "WireCellAux/Logger.h"

namespace WireCell {
    namespace Img {

        class GridTiling : public Aux::Logger, public ITiling, public IConfigurable {
           public:
            GridTiling();
            virtual ~GridTiling();
            virtual void configure(const WireCell::Configuration& cfg);
            virtual WireCell::Configuration default_configuration() const;

            virtual bool operator()(const input_pointer& slice, output_pointer& blobset);

           private:
            /** Config: threshold

                Channel activity in slice must have at least this
                value (charge) to be considered in the tiling.
            */
            double m_threshold{0.0};
            /** Config: nudge

                Effectively move any 2-layer crossing point toward the
                center of the existing blob by this fraction of a
                pitch in a 3rd layer prior to testing if that crossing
                point is inside the strip of that 3rd layer.  This
                fights floating-point imprecision in wire geometry.
            */
            double m_nudge{1e-3};

            // Count blobs in each contiguous stream to assign blob
            // ident numbers.
            size_t m_blobs_seen{0};

            IAnodePlane::pointer m_anode;
            IAnodeFace::pointer m_face;

            IBlobSet::pointer make_empty(const input_pointer& slice);

        };
    }  // namespace Img
}  // namespace WireCell

#endif
