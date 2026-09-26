/** BlobCutting -- split wide blobs into sub-blobs (the isochronous "sub-blob generator").

    IBlobSet -> IBlobSet function node.  Every blob whose U, V or W strip is
    wider than `length_threshold` wires is recursively bisected at the
    mid-point of its widest wire-plane strip until no strip is wider than the
    threshold (or `max_depth` bisections were done, or the widest strip is not
    wider than `min_length`).  Each half is rebuilt with RayGrid::Blob::add on
    the same strips except the halved one, then RayGrid::drop_invalid ->
    RayGrid::prune -> RayGrid::drop_invalid; halves that end up with fewer
    than three corners are dropped.  If neither half survives, the blob is
    kept whole.  Blobs that need no cut pass through as the same IBlob.

    Strips are half-open [lo, hi) so the two halves are disjoint in the
    halved plane and their union is the parent's strip.  A sub-blob is an
    ordinary Aux::SimpleBlob with value = parent value / n_children and the
    parent's uncertainty, slice and face; idents come from a per-frame counter
    that starts at `ident_base` (default 1<<20, above GridTiling's per-event
    count) and restarts at every frame boundary and at EOS, so idents stay
    unique within an event and identical between runs.

    What this node does NOT do: it does not re-tile activity (sub-blobs are
    purely geometric subdivisions of the parent's strips, which by tiling are
    fully active), it does not build slice/wire/channel or measure edges
    (BlobClustering does that downstream from shape().strips()), and it does
    not know about masked or dummy planes: a blob from a two-view slicing pass
    has its dummy-plane strip cut like any other, so place the node on the
    three-plane ("active") path only.

    Ported from Xuyang Ning's XN_apply-pointcloud branch (toolkit commit
    19a7f8f5, img/src/BlobCutting.cxx, used by pdhd/img_simple.jsonnet), with
    the split logic unchanged; see wcfm doc 03 for the differences (linkage,
    logging, ident scheme, `min_length`).
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
            // A wire-plane strip wider than this many wires triggers a cut.
            int m_length_threshold{20};
            // Tolerance passed to RayGrid::Blob::add / prune (same as GridTiling's nudge).
            double m_nudge{0.01};
            // Maximum number of bisections along any one parent->child chain.
            int m_max_depth{10};
            // A blob whose widest wire-plane strip is not wider than this is never cut.
            int m_min_length{2};
            // First sub-blob ident of each frame.
            int m_ident_base{1 << 20};

            int m_next_ident{0};
            int m_last_frame_ident{-1};
            bool m_have_frame{false};
            size_t m_count{0};
        };

    }  // namespace Img
}  // namespace WireCell

#endif
