// RayTiling.h
#ifndef WIRECELL_SPNG_RAYGRID_RAYTILING_H
#define WIRECELL_SPNG_RAYGRID_RAYTILING_H

#include <torch/torch.h>
#include <vector>
#include <string>
#include <numeric> // For std::iota
#include <algorithm> // For std::all_of

// Forward declaration of Coordinates class
namespace WireCell {
    namespace Spng {
        namespace RayGrid {
            class Coordinates; // Forward declaration
        }
    }
}


namespace WireCell::Spng::RayGrid {

    // --- Helper functions ---

    /**
     * @brief Returns a blob strip pair view index tensor for a blob of the given number of views.
     * @param nviews The number of views.
     * @return A tensor of shape (npairs, 2) where each row is a pair of view indices.
     * @throws std::invalid_argument if nviews < 2.
     */
    torch::Tensor strip_pairs(long nviews);

    /**
     * @brief Replicates a tensor n-times along a new first dimension.
     * @param ten The input tensor.
     * @param n The number of times to replicate.
     * @return The expanded tensor.
     */
    torch::Tensor expand_first(const torch::Tensor& ten, long n);

    /**
     * @brief Returns a (4,2) tensor giving the canonical ordering of the 4 crossings of pairs of edges.
     * @return A tensor of shape (4, 2) with edge indices.
     */
    torch::Tensor strip_pair_edge_indices();

    /**
     * @brief Returns a blob crossings tensor.
     * @param blobs The (nblobs, nstrips, 2) blobs tensor.
     * @return A tensor of shape (nblobs, npairs, 4, 2) holding ray indices.
     */
    torch::Tensor blob_crossings(const torch::Tensor& blobs);

    /**
     * @brief Flattens the crossings tensor and returns view and ray indices.
     * @param crossings The (nblobs, npairs, 4, 2) crossings tensor.
     * @param nviews The number of views from which the crossings were built.
     * @return A tuple of four 1D tensors: (v1, r1, v2, r2), each of shape (nblobs * npairs * 4,).
     * @throws std::invalid_argument if crossings tensor has an unexpected shape.
     */
    std::tuple<torch::Tensor, torch::Tensor, torch::Tensor, torch::Tensor>
    flatten_crossings(const torch::Tensor& crossings, long nviews);

    /**
     * @brief Returns True if the crossing of (v1,r1) and (v2,r2) are inside the
     * inclusive range [rbegin3,rend3] in view v3.
     * @param coords The Coordinates object.
     * @param v1 View index of the first ray.
     * @param r1 Ray index of the first ray.
     * @param v2 View index of the second ray.
     * @param r2 Ray index of the second ray.
     * @param v3 View index of the third strip.
     * @param rbegin3 Beginning ray index of the third strip.
     * @param rend3 Ending ray index (half-open) of the third strip.
     * @param nudge A small value to nudge pitches for boundary checks.
     * @return A boolean tensor indicating if the crossing is inside.
     */
    torch::Tensor crossing_in_other(const Coordinates& coords,
                                    const torch::Tensor& v1, const torch::Tensor& r1,
                                    const torch::Tensor& v2, const torch::Tensor& r2,
                                    const torch::Tensor& v3,
                                    const torch::Tensor& rbegin3, const torch::Tensor& rend3,
                                    double nudge = 1e-3);

    /**
     * @brief Finds all consecutive regions values above threshold and returns their half-open ranges.
     * @param activity A 1D real or boolean torch::Tensor.  If real, threshold is applied to form bool.
     * @return A 2D torch::Tensor of shape (N_runs, 2), where N_runs is the number of consecutive True segments.
     * Each row is [start_index, end_index_half_open). Returns an empty tensor if no True runs are found.
     */
    torch::Tensor get_true_runs(torch::Tensor activity, double threshold=0.0);

    /**
     * @brief Expands existing blobs into new blobs based on intersections with
     * consecutive True regions in the new 'activity' tensor.
     * @param blobs The existing (nblobs_old, nviews_old, 2) tensor.
     * @param lo A (nblobs_old,) tensor of low bounds for the new view, per blob.
     * @param hi A (nblobs_old,) tensor of high (half-open) bounds for the new view, per blob.
     * @param activity A 1D boolean torch::Tensor representing the new view's activity.
     * @return A new blobs tensor of shape (nblobs_new, nviews_new, 2).
     */
    torch::Tensor expand_blobs_with_activity(const torch::Tensor& blobs,
                                             const torch::Tensor& lo,
                                             const torch::Tensor& hi,
                                             const torch::Tensor& activity);

    /**
     * @brief Marks blob crossings as being inside or not.
     * @param coords A Coordinates object
     * @param blobs A "blobs" tensor.  This determines the device for the return value.
     * @param crossings A "crossings" tensor.
     * @return An "insides" tensor.
     */
    torch::Tensor blob_insides(const Coordinates& coords, const torch::Tensor& blobs, const torch::Tensor& crossings);

    /**
     * @brief Calculates the half-open ray index range (lo, hi) for each blob in the new view,
     * considering only crossings where 'insides' is True.
     * @param coords The Coordinates object.
     * @param crossings The (nblobs, npairs, 4, 2) tensor reflecting the *current* blobs and views.
     * @param nviews The total number of views that *currently* exist in 'crossings'.
     * @param insides The (nblobs, npairs, 4) boolean tensor indicating which crossings are "inside".
     * @return A tuple (lo, hi) where lo and hi are 1D tensors of shape (nblobs,).
     */
    std::tuple<torch::Tensor, torch::Tensor>
    blob_bounds(const Coordinates& coords, const torch::Tensor& crossings,
                long nviews, const torch::Tensor& insides);

    /**
     * @brief Clamps per-blob bounds to be consistent with an activity of length nmeasures.
     * @param lo A 1D tensor of low bounds.
     * @param hi A 1D tensor of high bounds.
     * @param nmeasures The maximum allowed measure (exclusive upper bound).
     * @return A tuple (lo, hi) with clamped bounds.
     */
    std::tuple<torch::Tensor, torch::Tensor>
    bounds_clamp(torch::Tensor lo, torch::Tensor hi, long nmeasures);

    /**
     * @brief Returns the initial trivial blobs tensor for a 2-view case.
     * @return A tensor of shape (1, 2, 2) representing the initial blob.
     */
    torch::Tensor trivial_blobs();


    // --- Tiling class ---

    class Tiling {
    public:
        // The blobs tensor. (N blobs, N strips, 2 ray indices)
        torch::Tensor blobs;
        // The crossing tensor. (N blobs, N strip pairs, 4 strip edge pairs, 2 strips, 2 rays)
        torch::Tensor crossings;
        // The insides tensor. (N blobs, N strip pairs, 4 strip edge pairs)
        torch::Tensor insides;

        /**
         * @brief Constructs a Tiling object.
         * @param blobs The blobs tensor.
         * @param crossings The crossings tensor.
         * @param insides The insides tensor.
         */
        Tiling(const torch::Tensor& blobs, const torch::Tensor& crossings, const torch::Tensor& insides);

        /**
         * @brief Returns the number of view layers of this tiling.
         * @return The number of views.
         */
        long nviews() const;

        /**
         * @brief Returns the number of blobs in this tiling.
         * @return The number of blobs.
         */
        long nblobs() const;

        /**
         * @brief Return new tiling with a subset of blobs given by indices.
         * @param ind A 1D tensor of indices to select.
         * @return A new Tiling object containing the selected subset.
         */
        Tiling select(const torch::Tensor& ind) const;

        /**
         * @brief Returns a string representation of the Tiling object.
         * @param depth The level of detail for the string representation.
         * @return A string describing the Tiling.
         */
        std::string as_string(int depth = 1) const;
    };

    /**
     * @brief Returns tiling solution for trivial 2 view case.
     * @return A Tiling object representing the trivial solution.
     */
    Tiling trivial();

    /**
     * @brief Applies activity to blobs to make new blobs with one more view.
     * @param coords The Coordinates object.
     * @param blobs The existing blobs tensor.
     * @param activity A 1D boolean tensor representing the new view's activity.
     * @return The new blobs tensor
     */
    torch::Tensor
    apply_activity(const Coordinates& coords, const torch::Tensor& blobs,
                   const torch::Tensor& activity);


    /**
     * @brief Applies activity to blobs to make new blobs with one more view.
     * @param coords The Coordinates object.
     * @param blobs The existing blobs tensor.
     * @param activity A 1D boolean tensor representing the new view's activity.
     * @return The half-open bounds which correspond to the new blobs which would be 
     * returned -- this is used within apply_activity
     */
    std::tuple<torch::Tensor, torch::Tensor>
    retrieve_new_blob_bounds(const Coordinates& coords, const torch::Tensor& blobs,
                        const torch::Tensor& activity);

} // namespace WireCell::Spng::RayGrid

#endif // WIRECELL_SPNG_RAYGRID_RAYTILING_H
