// RayTiling.cpp
#include "WireCellSpng/RayTiling.h"
#include "WireCellSpng/RayGrid.h" // For WireCell::SPNG::RayGrid::Coordinates
#include <iostream>

using namespace torch::indexing;

namespace WireCell::SPNG::RayGrid {


    // --- Helper functions implementation ---

    torch::Tensor strip_pairs(long nviews) {
        if (nviews < 2) {
            throw std::invalid_argument("too few views to make pairs: " + std::to_string(nviews));
        }
        if (nviews == 2) {
            return torch::tensor({{0, 1}}, torch::kLong);
        }

        torch::Tensor prior = strip_pairs(nviews - 1);
        // Python: torch.tensor(list(product(range(nviews-1), [nviews-1])))
        // This creates pairs like (0, nviews-1), (1, nviews-1), ..., (nviews-2, nviews-1)
        std::vector<long> flat_more_pairs;
        long num_new_pairs = nviews - 1;
        for (long i = 0; i < num_new_pairs; ++i) {
            flat_more_pairs.push_back(i);
            flat_more_pairs.push_back(nviews - 1);
        }
        // Create the tensor with the correct shape
        torch::Tensor more = torch::tensor(flat_more_pairs, torch::kLong).reshape({num_new_pairs, 2});
        return torch::vstack({prior, more});
    }

    torch::Tensor expand_first(const torch::Tensor& ten, long n) {
        // Python: shape = [n] + [-1] * ten.ndim
        // return ten.unsqueeze(0).expand(*shape)
        std::vector<long> shape_vec;
        shape_vec.push_back(n);
        for (long i = 0; i < ten.dim(); ++i) {
            shape_vec.push_back(-1); // -1 means infer dimension size
        }
        return ten.unsqueeze(0).expand(shape_vec);
    }

    torch::Tensor strip_pair_edge_indices() {
        return torch::tensor({{0, 0}, {0, 1}, {1, 0}, {1, 1}}, torch::kLong);
    }

    torch::Tensor blob_crossings(const torch::Tensor& blobs) {
        long nstrips = blobs.size(1);
        // (npairs)

        auto sp = strip_pairs(nstrips);
        std::vector<torch::Tensor> views_unbound = sp.unbind(/*dim=*/1);
        torch::Tensor views_a = views_unbound[0];
        torch::Tensor views_b = views_unbound[1];
        // views shape (npairs,)

        std::vector<torch::Tensor> edges_unbound = strip_pair_edge_indices().unbind(/*dim=*/1);
        torch::Tensor edges_a = edges_unbound[0];
        torch::Tensor edges_b = edges_unbound[1];
        // edges shape (4, )

        // (nblobs, npairs, 4)
        // Python: blobs[:, views_a, :][..., edges_a]
        // C++: Use advanced indexing
        auto v_a = blobs.index({Slice(), views_a, Slice()});
        auto v_b = blobs.index({Slice(), views_b, Slice()});
        // v shape (nblobs, npairs, 2)

        auto rays_a = v_a.index({Ellipsis, edges_a});
        auto rays_b = v_b.index({Ellipsis, edges_b});
        // ray shape (nblobs, npairs, 4)

        auto bc = torch::stack({rays_a, rays_b}, /*dim=*/3);
        // bc shape (nblobs, npairs, 4, 2)
        return bc;
    }

    std::tuple<torch::Tensor, torch::Tensor, torch::Tensor, torch::Tensor>
    flatten_crossings(const torch::Tensor& crossings, long nviews) {
        long nblobs = crossings.size(0);
        long npairs = crossings.size(1);
        long num_crossing_types = crossings.size(2); // Should be 4

        // Validate crossings shape
        if (crossings.sizes() != torch::IntArrayRef({nblobs, npairs, num_crossing_types, 2})) {
            std::stringstream ss;
            ss << "crossings tensor has unexpected shape: " << crossings.sizes() << ". Expected {" << nblobs << ", " << npairs << ", 4, 2}.";
            throw std::invalid_argument( ss.str() );
        }

        // 2. Extract r1 and r2 directly from crossings
        // Flatten across nblobs, npairs, and 4
        torch::Tensor r1 = crossings.index({torch::indexing::Slice(), torch::indexing::Slice(), torch::indexing::Slice(), 0}).reshape(-1);
        torch::Tensor r2 = crossings.index({torch::indexing::Slice(), torch::indexing::Slice(), torch::indexing::Slice(), 1}).reshape(-1);

        // 3. Create v1 and v2, repeated for nblobs and 4 crossing bounds
        //    First, get the base v1 and v2 from strip_pairs
        torch::Tensor pairs_tensor = strip_pairs(nviews);
        torch::Tensor base_v1 = pairs_tensor.index({torch::indexing::Slice(), 0}); // Shape (npairs,)
        torch::Tensor base_v2 = pairs_tensor.index({torch::indexing::Slice(), 1}); // Shape (npairs,)

        //    Repeat these for each blob and each of the 4 crossing bounds
        // (npairs,) -> (1, npairs, 1, 1) -> (nblobs, npairs, 4, 1)
        torch::Tensor v1 = base_v1.unsqueeze(0).unsqueeze(2).unsqueeze(2).expand({nblobs, -1, 4, -1}).reshape(-1);
        torch::Tensor v2 = base_v2.unsqueeze(0).unsqueeze(2).unsqueeze(2).expand({nblobs, -1, 4, -1}).reshape(-1);

        return std::make_tuple(v1, r1, v2, r2);
    }


    torch::Tensor crossing_in_other(const Coordinates& coords,
                                    const torch::Tensor& v1, const torch::Tensor& r1,
                                    const torch::Tensor& v2, const torch::Tensor& r2,
                                    const torch::Tensor& v3,
                                    const torch::Tensor& rbegin3, const torch::Tensor& rend3,
                                    double nudge) {
        // pitches = coords.pitch_location((v1,r1), (v2,r2), v3)
        torch::Tensor pitches = coords.pitch_location(v1, r1, v2, r2, v3);

        // pinds_lo = coords.pitch_index(pitches+nudge, v3)
        torch::Tensor pinds_lo = coords.pitch_index(pitches + nudge, v3);
        // pinds_hi = coords.pitch_index(pitches-nudge, v3)
        torch::Tensor pinds_hi = coords.pitch_index(pitches - nudge, v3);

        torch::Tensor pinds = coords.pitch_index(pitches, v3);
        torch::Tensor in_other = (pinds >= rbegin3) & (pinds < rend3);

        return in_other;
    }

    torch::Tensor get_true_runs(torch::Tensor activity, double threshold) {
        auto device = activity.device();

        auto Loptions = torch::TensorOptions().dtype(torch::kLong).device(device);

        if (activity.numel() == 0) {
            return torch::empty({0, 2}, Loptions);
        }

        if (activity.dtype() != torch::kBool) {
            activity = activity > threshold;
        }


        // Convert boolean to int (True=1, False=0)
        torch::Tensor activity_int = activity.to(torch::kInt);
        auto Ioptions = torch::TensorOptions().dtype(torch::kInt).device(device);

        // Pad with zeros at both ends to ensure all runs are 'closed' by zeros.
        torch::Tensor padded_activity_int = torch::cat({
                torch::tensor({0}, Ioptions),
                activity_int,
                torch::tensor({0}, Ioptions)
            });

        // Calculate the difference to find transitions (0->1 for starts, 1->0 for ends)
        torch::Tensor diff = padded_activity_int.diff();

        // Find start indices of True runs: where diff changes from 0 to 1
        torch::Tensor start_indices = (diff == 1).nonzero().squeeze(-1);

        // Find end indices (exclusive) of True runs: where diff changes from 1 to 0
        torch::Tensor end_indices_half_open = (diff == -1).nonzero().squeeze(-1);

        if (start_indices.numel() == 0) {
            return torch::empty({0, 2}, Loptions);
        }

        return torch::stack({start_indices, end_indices_half_open}, /*dim=*/-1);
    }

    torch::Tensor expand_blobs_with_activity(
        const torch::Tensor& blobs,
        const torch::Tensor& lo,
        const torch::Tensor& hi,
        const torch::Tensor& activity
        ) {
        long nblobs_old = blobs.size(0);
        long nviews_old = blobs.size(1);
        long nviews_new = nviews_old + 1;
        torch::Device device = blobs.device();

        // Step 1: Get all consecutive True runs from the new activity tensor.
        torch::Tensor true_runs_all = get_true_runs(activity);

        if (true_runs_all.numel() == 0) {
            return torch::empty({0, nviews_new, 2}, blobs.dtype()).to(device);
        }

        long num_true_runs = true_runs_all.size(0);

        // Step 2: Calculate intersections for every old blob with every True run.
        torch::Tensor lo_expanded = lo.unsqueeze(1).expand({-1, num_true_runs});
        torch::Tensor hi_expanded = hi.unsqueeze(1).expand({-1, num_true_runs});

        torch::Tensor true_runs_all_lo = true_runs_all.index({torch::indexing::Slice(), 0}).unsqueeze(0).expand({nblobs_old, -1});
        torch::Tensor true_runs_all_hi = true_runs_all.index({torch::indexing::Slice(), 1}).unsqueeze(0).expand({nblobs_old, -1});

        torch::Tensor intersection_lo = torch::max(lo_expanded, true_runs_all_lo);
        torch::Tensor intersection_hi = torch::min(hi_expanded, true_runs_all_hi);

        torch::Tensor valid_intersections_mask = intersection_lo < intersection_hi;

        // Step 3: Identify the (old_blob_index, true_run_index) pairs that form new blobs.
        // torch::Tensor flat_blob_indices = torch::empty({0}, torch::kLong);
        // torch::Tensor flat_true_run_indices = torch::empty({0}, torch::kLong);
        // std::tie(flat_blob_indices, flat_true_run_indices) = torch::nonzero(valid_intersections_mask, /*as_tuple=*/true);
        // We need to extract the columns for blob and true_run indices.
        torch::Tensor nonzero_indices = torch::nonzero(valid_intersections_mask);
        torch::Tensor flat_blob_indices = nonzero_indices.index({torch::indexing::Slice(), 0});
        torch::Tensor flat_true_run_indices = nonzero_indices.index({torch::indexing::Slice(), 1});

        long total_new_blobs = flat_blob_indices.numel();

        if (total_new_blobs == 0) {
            return torch::empty({0, nviews_new, 2}, blobs.dtype()).to(device);
        }

        // Step 4: Construct the new_blobs tensor.
        torch::Tensor old_blob_data = blobs.index({flat_blob_indices});

        torch::Tensor new_view_ranges_flat = torch::stack({
                intersection_lo.index({flat_blob_indices, flat_true_run_indices}),
                intersection_hi.index({flat_blob_indices, flat_true_run_indices})
            }, /*dim=*/-1);

        torch::Tensor new_view_ranges_reshaped = new_view_ranges_flat.unsqueeze(1);

        torch::Tensor new_blobs = torch::cat({old_blob_data, new_view_ranges_reshaped}, /*dim=*/1);

        return new_blobs;
    }

    torch::Tensor blob_insides(const Coordinates& coords, const torch::Tensor& blobs, const torch::Tensor& crossings)
    {
        long nviews = blobs.size(1);
        auto device = blobs.device();

        long nblobs = crossings.size(0);
        long npairs = crossings.size(1);
        long num_crossing_types = crossings.size(2); // Should be 4

        if (nviews < 3) {
            // If nviews < 3, there are no 'other' strips to check against (nviews-2 is < 1).
            // In this case, all crossings are trivially "inside" all other strips (because there are none).
            return torch::ones({nblobs, npairs, num_crossing_types}, torch::kBool).to(device);
        }

        // 1. Prepare base flattened crossing data (v1, r1, v2, r2)
        torch::Tensor v1_base, r1_base, v2_base, r2_base;
        std::tie(v1_base, r1_base, v2_base, r2_base) = flatten_crossings(crossings, nviews);
    
        // long num_flattened_crossings = v1_base.numel();

        // 2. Determine 'other' view indices for each pair
        torch::Tensor current_pairs = strip_pairs(nviews).to(device);
    
        torch::Tensor all_view_indices = torch::arange(nviews, torch::kLong).to(device);

        // (npairs, nviews) boolean mask, True where view is one of the pair views
        torch::Tensor is_in_pair_mask =
            (all_view_indices.unsqueeze(0) == current_pairs.index({torch::indexing::Slice(), 0}).unsqueeze(1)) |
            (all_view_indices.unsqueeze(0) == current_pairs.index({torch::indexing::Slice(), 1}).unsqueeze(1));

        // (npairs, nviews) boolean mask, True where view is NOT one of the pair views
        //torch::Tensor is_other_view_mask = !is_in_pair_mask;
        torch::Tensor is_other_view_mask = torch::logical_not(is_in_pair_mask);

        // List of 1D tensors, each holding the 'other' view indices for a specific pair
        std::vector<torch::Tensor> other_views_list;
        for (long p_idx = 0; p_idx < npairs; ++p_idx) {
            other_views_list.push_back(all_view_indices.index({is_other_view_mask.index({p_idx})}));
        }

        long n_others_per_pair = nviews - 2;
    
        torch::Tensor other_views_tensor = torch::stack(other_views_list, /*dim=*/0); // (npairs, n_others_per_pair)

        // 3. Prepare full flattened inputs for crossing_in_other
        torch::Tensor v1_expanded = v1_base.unsqueeze(1);
        torch::Tensor r1_expanded = r1_base.unsqueeze(1);
        torch::Tensor v2_expanded = v2_base.unsqueeze(1);
        torch::Tensor r2_expanded = r2_base.unsqueeze(1);

        torch::Tensor v3_base_per_pair_and_crossing = other_views_tensor.unsqueeze(1).expand({-1, num_crossing_types, -1});
        torch::Tensor v3_base_per_crossing = v3_base_per_pair_and_crossing.reshape({npairs * num_crossing_types, n_others_per_pair});
        torch::Tensor v3_flat = v3_base_per_crossing.unsqueeze(0).expand({nblobs, -1, -1}).reshape(-1);

        torch::Tensor blob_indices_overall = torch::arange(nblobs, torch::kLong).to(device);
        torch::Tensor blob_indices_flat = blob_indices_overall.unsqueeze(1).expand({-1, npairs * num_crossing_types * n_others_per_pair}).reshape(-1);

        torch::Tensor rbegin3_flat = blobs.index({blob_indices_flat, v3_flat, 0});
        torch::Tensor rend3_flat = blobs.index({blob_indices_flat, v3_flat, 1});

        torch::Tensor v1_final = v1_expanded.expand({-1, n_others_per_pair}).reshape(-1);
        torch::Tensor r1_final = r1_expanded.expand({-1, n_others_per_pair}).reshape(-1);
        torch::Tensor v2_final = v2_expanded.expand({-1, n_others_per_pair}).reshape(-1);
        torch::Tensor r2_final = r2_expanded.expand({-1, n_others_per_pair}).reshape(-1);

        // 4. Call crossing_in_other with all flattened inputs
        torch::Tensor check_results = crossing_in_other(
            coords,
            v1_final, r1_final,
            v2_final, r2_final,
            v3_flat,
            rbegin3_flat, rend3_flat
            );

        // 5. Reshape results back
        torch::Tensor reshaped_check_results = check_results.reshape({nblobs, npairs, num_crossing_types, n_others_per_pair});

        // 6. Combine results: A crossing is "inside" if it's inside ALL other strips.
        torch::Tensor insides = torch::all(reshaped_check_results, /*dim=*/-1);

        return insides;
    }
    

    std::tuple<torch::Tensor, torch::Tensor>
    blob_bounds(const Coordinates& coords, const torch::Tensor& crossings,
                long nviews, const torch::Tensor& insides) {
        long nblobs = crossings.size(0);
        // long num_crossing_types = crossings.size(2); // Should be 4

        // 1. Flatten crossings to get v1, r1, v2, r2 for pitch_location
        torch::Tensor v1, r1, v2, r2;
        std::tie(v1, r1, v2, r2) = flatten_crossings(crossings, nviews);

        // 2. Determine v3 (the index of the new view)
        torch::Tensor v3 = torch::full_like(v1, nviews, torch::kLong);

        // 3. Calculate pitches for all crossings
        torch::Tensor pitches = coords.pitch_location(v1, r1, v2, r2, v3);

        // 4. Reshape pitches to (nblobs, npairs * 4) to align with flattened 'insides'
        torch::Tensor blob_pitches = pitches.reshape({nblobs, -1});

        // 5. Flatten 'insides' to (nblobs, npairs * 4) to match blob_pitches
        torch::Tensor insides_flat = insides.reshape({nblobs, -1});

        // 6. Apply the 'insides' mask to blob_pitches for min/max calculation
        torch::Tensor masked_pitches_for_min = torch::where(
            insides_flat,
            blob_pitches,
            torch::full_like(blob_pitches, std::numeric_limits<double>::infinity())
            );

        torch::Tensor masked_pitches_for_max = torch::where(
            insides_flat,
            blob_pitches,
            torch::full_like(blob_pitches, -std::numeric_limits<double>::infinity())
            );

        // 7. Calculate pmin and pmax for each blob across the (npairs * 4) dimension
        torch::Tensor pmin = std::get<0>(torch::min(masked_pitches_for_min, /*dim=*/1));
        torch::Tensor pmax = std::get<0>(torch::max(masked_pitches_for_max, /*dim=*/1));

        // 8. Handle cases where a blob has NO 'inside' crossings
        // torch::Tensor no_valid_crossings_mask = !torch::any(insides_flat, /*dim=*/1);
        torch::Tensor no_valid_crossings_mask = torch::logical_not(torch::any(insides_flat, /*dim=*/1));

        pmin = torch::where(no_valid_crossings_mask, torch::tensor(0.0, pmin.dtype()), pmin);
        pmax = torch::where(no_valid_crossings_mask, torch::tensor(0.0, pmax.dtype()), pmax);

        // 9. Convert pitches (pmin, pmax) to ray indices using coords.pitch_index
        long next_layer = nviews;

        torch::Tensor lo = coords.pitch_index(pmin, torch::tensor(next_layer, torch::kLong));
        torch::Tensor hi = coords.pitch_index(pmax, torch::tensor(next_layer, torch::kLong)) + 1;

        return std::make_tuple(lo, hi);
    }

    std::tuple<torch::Tensor, torch::Tensor>
    bounds_clamp(torch::Tensor lo, torch::Tensor hi, long nmeasures) {
        lo.index_put_({lo < 0}, 0);
        lo.index_put_({lo > nmeasures}, nmeasures);
        hi.index_put_({hi < 0}, 0);
        hi.index_put_({hi > nmeasures}, nmeasures);
        return std::make_tuple(lo, hi);
    }

    torch::Tensor trivial_blobs() {
        return torch::tensor({{{0, 1}, {0, 1}}}, torch::kLong);
    }

// --- Tiling class implementation ---

    Tiling::Tiling(const torch::Tensor& blobs_in, const torch::Tensor& crossings_in, const torch::Tensor& insides_in)
        : blobs(blobs_in), crossings(crossings_in), insides(insides_in) {}

    long Tiling::nviews() const {
        return blobs.size(1);
    }

    long Tiling::nblobs() const {
        return blobs.size(0);
    }

    Tiling Tiling::select(const torch::Tensor& ind) const {
        return Tiling(blobs.index({ind}),
                      crossings.index({ind}),
                      insides.index({ind}));
    }

    std::string Tiling::as_string(int depth) const {
        std::stringstream ss;
        ss << "blobs shape: " << blobs.sizes();
        std::string lines = ss.str();
        if (depth > 1) {
            lines += "\n\tblobs=\n" + torch::str(blobs);
        }
        if (depth > 2) {
            lines += "\n\tcrossings=\n" + torch::str(crossings);
        }
        return lines;
    }

    Tiling trivial() {
        torch::Tensor blobs = trivial_blobs();
        torch::Tensor crossings = blob_crossings(blobs);
        // Non-trivial case, call blob_insides() but the trivial case can be made
        // without requiring a coords object.
        torch::Tensor insides = torch::ones({1, 1, 4}, torch::kBool);

        return Tiling(blobs, crossings, insides);
    }

    torch::Tensor
    apply_activity(const Coordinates& coords, const torch::Tensor& blobs_in,
                   const torch::Tensor& activity) {
        torch::Tensor lo, hi;
        std::tie(lo, hi) = retrieve_new_blob_bounds(coords, blobs_in, activity);
        torch::Tensor new_blobs = expand_blobs_with_activity(blobs_in, lo, hi, activity);

        return new_blobs;
    }

    std::tuple<torch::Tensor, torch::Tensor>
    retrieve_new_blob_bounds(const Coordinates& coords, const torch::Tensor& blobs_in,
                        const torch::Tensor& activity) {
        long nviews = blobs_in.size(1);

        torch::Tensor crossings = blob_crossings(blobs_in);
        torch::Tensor insides = blob_insides(coords, blobs_in, crossings);
        torch::Tensor lo, hi;
        std::tie(lo, hi) = blob_bounds(coords, crossings, nviews, insides);

        auto results = bounds_clamp(lo, hi, activity.size(0));

        return results;
    }


} // namespace WireCell::SPNG::RayGrid

