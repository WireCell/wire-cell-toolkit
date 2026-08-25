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
            // ident numbers.  Reset at EOS and, for a multi-event process,
            // at each frame boundary (see m_last_frame_ident).
            size_t m_blobs_seen{0};

            // The frame ident of the last slice seen, or -1 before any.  A
            // process that streams MANY events (doc 76 round 2 group mode)
            // gets no EOS between them, so m_blobs_seen would keep counting
            // and event N's blobs would be identified from an offset instead
            // of from 0.  That is not cosmetic: blob idents key the
            // unordered_map<int,...>/unordered_set<int> containers in
            // InSliceDeghosting, ProjectionDeghosting, BlobGrouping and
            // LocalGeomClustering, whose ITERATION ORDER depends on the key
            // values, so an offset flips order-dependent deghosting decisions
            // and changes blob COUNTS and CHARGES on some events (doc 81
            // round 0).  Inert for a one-event process: the sequence already
            // starts at 0 and there is only ever one frame ident.
            int m_last_frame_ident{-1};

            IAnodePlane::pointer m_anode;
            IAnodeFace::pointer m_face;

            IBlobSet::pointer make_empty(const input_pointer& slice);

        };
    }  // namespace Img
}  // namespace WireCell

#endif
