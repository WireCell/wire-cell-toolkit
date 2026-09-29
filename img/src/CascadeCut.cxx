// The refinement operator of CascadeDeghosting (wcfm doc 14): BlobCutting's bisection, duplicated (fork rule:
// img/src/BlobCutting.cxx is a production component and stays untouched) with one addition, the depth of every
// node counted from the uncut blob.  BlobCutting's max_depth bounds the bisections of one parent->child chain
// within ONE call; cutting level by level (32, 16, 8, 4) would restart that count at every level and cut deeper
// than one direct cut to 4 wires, which is the tier the models were trained on (doc 09 `_sub4f`: BlobCutting
// length 4, max_depth 10).  Carrying the depth across levels makes the cascade's cells exactly the direct cut's
// (doc 12 T0: the widths are nested cuts of one bisection tree; checked in doctest_cascade_deghosting).
//
// split_blob_once and the recursion are copied from BlobCutting.cxx (lines 56-173 at toolkit af93ac68).

#include "WireCellImg/CascadeGraph.h"

#include "WireCellUtil/RayClustering.h"  // surrounding()
#include "WireCellUtil/RayTiling.h"

using namespace WireCell;

namespace {

    bool is_wire_layer(const RayGrid::Strip& strip) { return strip.layer >= 2 && strip.layer <= 4; }

    int strip_width(const RayGrid::Strip& strip) { return strip.bounds.second - strip.bounds.first; }

    bool is_blob_covered_by(const RayGrid::Blob& small_blob, const RayGrid::Blob& big_blob)
    {
        RayGrid::blobs_t big_blobs = {big_blob};
        RayGrid::blobs_t small_blobs = {small_blob};
        auto big_ref = big_blobs.begin();
        auto small_ref = small_blobs.begin();
        return RayGrid::surrounding(big_ref, small_ref);
    }

    // BlobCutting.cxx split_blob_once, verbatim
    std::vector<RayGrid::Blob> split_blob_once(const RayGrid::Coordinates& coords, const RayGrid::Blob& original_blob,
                                               double nudge, int min_length)
    {
        std::vector<RayGrid::Blob> result;

        const auto& strips = original_blob.strips();
        if (strips.size() < 3) {
            result.push_back(original_blob);
            return result;
        }

        int longest_strip_index = -1;
        int max_length = 0;
        for (size_t i = 0; i < strips.size(); ++i) {
            const auto& strip = strips[i];
            if (!is_wire_layer(strip)) continue;
            const int length = strip_width(strip);
            if (length > max_length) {
                max_length = length;
                longest_strip_index = (int) i;
            }
        }
        if (longest_strip_index == -1 || max_length <= min_length) {
            result.push_back(original_blob);
            return result;
        }

        const auto& longest_strip = strips[longest_strip_index];
        const int mid_point = (longest_strip.bounds.first + longest_strip.bounds.second) / 2;
        const RayGrid::Strip first_half_strip{longest_strip.layer, {longest_strip.bounds.first, mid_point}};
        const RayGrid::Strip second_half_strip{longest_strip.layer, {mid_point, longest_strip.bounds.second}};

        RayGrid::Blob first_blob;
        RayGrid::Blob second_blob;
        for (int i = 0; i < (int) strips.size(); ++i) {
            if (i != longest_strip_index) {
                first_blob.add(coords, strips[i], nudge);
                second_blob.add(coords, strips[i], nudge);
            }
            else {
                first_blob.add(coords, first_half_strip, nudge);
                second_blob.add(coords, second_half_strip, nudge);
            }
        }

        RayGrid::blobs_t new_blobs;
        if (!first_blob.corners().empty()) new_blobs.push_back(first_blob);
        if (!second_blob.corners().empty()) new_blobs.push_back(second_blob);

        RayGrid::drop_invalid(new_blobs);
        RayGrid::prune(coords, new_blobs, nudge);
        RayGrid::drop_invalid(new_blobs);

        for (const auto& processed_blob : new_blobs) {
            if (processed_blob.valid() && is_blob_covered_by(processed_blob, original_blob)) {
                result.push_back(processed_blob);
            }
        }

        if (result.empty()) {
            result.push_back(original_blob);
        }
        return result;
    }

    // BlobCutting.cxx split_blob_recursively, with the node depth carried along
    void split_recursively(const RayGrid::Coordinates& coords, const RayGrid::Blob& blob, int depth,
                           int length_threshold, const Img::Cascade::CutParams& par,
                           std::vector<std::pair<RayGrid::Blob, int>>& out)
    {
        if (!Img::Cascade::needs_cutting(blob, length_threshold) || par.max_depth - depth <= 0) {
            out.emplace_back(blob, depth);
            return;
        }
        auto split_result = split_blob_once(coords, blob, par.nudge, par.min_length);
        if (split_result.size() <= 1) {
            out.emplace_back(blob, depth);
            return;
        }
        for (const auto& sub_blob : split_result) {
            split_recursively(coords, sub_blob, depth + 1, length_threshold, par, out);
        }
    }
}  // namespace

bool Img::Cascade::needs_cutting(const RayGrid::Blob& blob, int length_threshold)
{
    for (const auto& strip : blob.strips()) {
        if (is_wire_layer(strip) && strip_width(strip) > length_threshold) {
            return true;
        }
    }
    return false;
}

std::vector<std::pair<RayGrid::Blob, int>> Img::Cascade::cut_shape(const RayGrid::Coordinates& coords,
                                                                   const RayGrid::Blob& blob, int depth, int width,
                                                                   const CutParams& par)
{
    std::vector<std::pair<RayGrid::Blob, int>> out;
    split_recursively(coords, blob, depth, width, par, out);
    return out;
}
