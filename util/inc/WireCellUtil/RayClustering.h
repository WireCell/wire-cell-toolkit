#ifndef WIRECELL_RAYCLUSTERING
#define WIRECELL_RAYCLUSTERING

#include "WireCellUtil/RayTiling.h"

#include <functional>
#include <unordered_set>
#include <vector>

namespace WireCell {
    namespace RayGrid {

        typedef typename blobs_t::const_iterator blobref_t;
        typedef typename std::vector<blobref_t> blobvec_t;
        typedef typename std::vector<blobvec_t> blobproj_t;

        struct blobref_hash {
            size_t operator()(blobref_t blobref) const { return std::hash<int*>()((int*) &*blobref); }
        };
        typedef typename std::unordered_set<blobref_t, blobref_hash> blobset_t;

        // project blobs along rays in layer
        blobproj_t projection(const blobvec_t& blobs, layer_index_t layer);

        // Return blobs in projection bound by the half open [ibeg,iend) range.
        blobvec_t select(const blobproj_t& proj, grid_range_t range);

        // Return blobs in projection which overlap with blob along
        // ray lines of given layer and all lower layers.
        blobvec_t overlap(const blobref_t& blob, const blobproj_t& proj, layer_index_t layer, grid_index_t tolerance, const bool verbose);
        inline blobvec_t overlap(const blobref_t& blob, const blobproj_t& proj, layer_index_t layer) {
            return overlap(blob, proj, layer, 0, false);
        }

        // Return true if a's strips are all inside b's strips or vice versa.
        bool surrounding(const blobref_t& a, const blobref_t& b);

        // Helper to make a vector of references from the backing vector.
        blobvec_t references(const blobs_t& blobs);

        // High level API.

        typedef std::function<void(blobref_t& a, blobref_t& b)> associator_t;

        // overlap function signature
        using overlap_t = std::function<blobvec_t(const blobref_t&, const blobproj_t&, layer_index_t)>;

        // Associate blobs from set one with blobs from set two which
        // mutually overlap.  The associator function will be called
        // with every pair of two overlapping blobs.
        void associate(
            const blobs_t& one, const blobs_t& two, associator_t func,
            overlap_t ol = static_cast<blobvec_t (*)(const blobref_t&, const blobproj_t&, layer_index_t)>(&overlap));

    }  // namespace RayGrid
}  // namespace WireCell

#endif
