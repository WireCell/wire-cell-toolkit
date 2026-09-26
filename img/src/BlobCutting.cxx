// BlobCutting: split wide blobs into sub-blobs.  See the header for the
// algorithm and for what was changed relative to Xuyang Ning's original
// (XN_apply-pointcloud, toolkit commit 19a7f8f5).

#include "WireCellImg/BlobCutting.h"
#include "WireCellAux/SimpleBlob.h"
#include "WireCellUtil/RayTiling.h"
#include "WireCellUtil/RayClustering.h"  // surrounding()
#include "WireCellUtil/NamedFactory.h"

#include <vector>

WIRECELL_FACTORY(BlobCutting, WireCell::Img::BlobCutting,
                 WireCell::INamed,
                 WireCell::IFunctionNode<WireCell::IBlobSet, WireCell::IBlobSet>,
                 WireCell::IConfigurable)

using namespace WireCell;

Img::BlobCutting::BlobCutting()
    : Aux::Logger("BlobCutting", "img")
{
}

Img::BlobCutting::~BlobCutting()
{
}

WireCell::Configuration Img::BlobCutting::default_configuration() const
{
    Configuration cfg;
    cfg["length_threshold"] = m_length_threshold;  // wires; a wider U/V/W strip triggers a cut
    cfg["nudge"] = m_nudge;                        // RayGrid::Blob::add / prune tolerance
    cfg["max_depth"] = m_max_depth;                // bisections per parent->child chain
    cfg["min_length"] = m_min_length;              // never cut a strip this narrow (or narrower)
    cfg["ident_base"] = m_ident_base;              // first sub-blob ident of each frame
    return cfg;
}

void Img::BlobCutting::configure(const WireCell::Configuration& cfg)
{
    m_length_threshold = get(cfg, "length_threshold", m_length_threshold);
    m_nudge = get(cfg, "nudge", m_nudge);
    m_max_depth = get(cfg, "max_depth", m_max_depth);
    m_min_length = get(cfg, "min_length", m_min_length);
    m_ident_base = get(cfg, "ident_base", m_ident_base);
    m_next_ident = m_ident_base;
    m_have_frame = false;
    log->debug("length_threshold={} nudge={} max_depth={} min_length={} ident_base={}",
               m_length_threshold, m_nudge, m_max_depth, m_min_length, m_ident_base);
}

namespace {

    // Wire-plane strips are RayGrid layers 2, 3, 4 (0 and 1 bound the active area).
    bool is_wire_layer(const RayGrid::Strip& strip)
    {
        return strip.layer >= 2 && strip.layer <= 4;
    }

    int strip_width(const RayGrid::Strip& strip)
    {
        return strip.bounds.second - strip.bounds.first;
    }

    bool blob_needs_cutting(const RayGrid::Blob& blob, int length_threshold)
    {
        for (const auto& strip : blob.strips()) {
            if (is_wire_layer(strip) && strip_width(strip) > length_threshold) {
                return true;
            }
        }
        return false;
    }

    // RayGrid::surrounding(a, b) is true when either blob contains the other.
    // After prune() a half never exceeds its parent, so this is a belt-and-braces
    // check kept from the original; it is expected to always pass.
    bool is_blob_covered_by(const RayGrid::Blob& small_blob, const RayGrid::Blob& big_blob)
    {
        RayGrid::blobs_t big_blobs = {big_blob};
        RayGrid::blobs_t small_blobs = {small_blob};
        auto big_ref = big_blobs.begin();
        auto small_ref = small_blobs.begin();
        return RayGrid::surrounding(big_ref, small_ref);
    }

    // Split a blob in two at the mid-point of its widest wire-plane strip.
    // Returns {original} when the blob cannot be split.
    std::vector<RayGrid::Blob> split_blob_once(const RayGrid::Coordinates& coords,
                                               const RayGrid::Blob& original_blob,
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
        // half-open [lo, mid) and [mid, hi): disjoint, union = parent strip
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

        // Same post-processing as tiling: drop degenerate halves, prune the
        // untouched strips to the new corner set, drop again.
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

    std::vector<RayGrid::Blob> split_blob_recursively(const RayGrid::Coordinates& coords,
                                                      const RayGrid::Blob& blob,
                                                      int length_threshold, double nudge,
                                                      int min_length, int max_depth)
    {
        std::vector<RayGrid::Blob> result;
        if (!blob_needs_cutting(blob, length_threshold) || max_depth <= 0) {
            result.push_back(blob);
            return result;
        }
        auto split_result = split_blob_once(coords, blob, nudge, min_length);
        if (split_result.size() <= 1) {
            result.push_back(blob);
            return result;
        }
        for (const auto& sub_blob : split_result) {
            auto sub_results = split_blob_recursively(coords, sub_blob, length_threshold, nudge,
                                                      min_length, max_depth - 1);
            result.insert(result.end(), sub_results.begin(), sub_results.end());
        }
        return result;
    }

}  // namespace

bool Img::BlobCutting::operator()(const input_pointer& in, output_pointer& out)
{
    if (!in) {
        out = nullptr;
        m_next_ident = m_ident_base;
        m_have_frame = false;
        log->debug("EOS at call={}", m_count);
        ++m_count;
        return true;
    }

    // Restart the sub-blob ident sequence at every frame boundary (as
    // GridTiling does for its own idents) so idents are per-event and
    // identical between runs.
    const auto slice = in->slice();
    const int frame_ident = (slice && slice->frame()) ? slice->frame()->ident() : -1;
    if (!m_have_frame || frame_ident != m_last_frame_ident) {
        m_next_ident = m_ident_base;
        m_last_frame_ident = frame_ident;
        m_have_frame = true;
    }

    auto sbs = std::make_shared<Aux::SimpleBlobSet>(in->ident(), slice);
    const auto& input_blobs = in->blobs();
    size_t ncut = 0, nsub = 0;

    for (const auto& iblob : input_blobs) {
        const auto& shape = iblob->shape();
        if (!blob_needs_cutting(shape, m_length_threshold)) {
            sbs->m_blobs.push_back(iblob);
            continue;
        }
        auto split_shapes = split_blob_recursively(iblob->face()->raygrid(), shape,
                                                   m_length_threshold, m_nudge,
                                                   m_min_length, m_max_depth);
        SPDLOG_LOGGER_TRACE(log, "blob ident={} slice={} cut into {} sub-blobs, idents from {}",
                            iblob->ident(), in->ident(), split_shapes.size(), m_next_ident);
        ++ncut;
        nsub += split_shapes.size();
        const float value = iblob->value() / split_shapes.size();  // parent charge shared equally
        for (const auto& new_shape : split_shapes) {
            sbs->m_blobs.push_back(std::make_shared<Aux::SimpleBlob>(
                m_next_ident++, value, iblob->uncertainty(), new_shape,
                iblob->slice(), iblob->face()));
        }
    }

    if (ncut) {
        SPDLOG_LOGGER_DEBUG(log, "call={} blobset={} blobs in={} cut={} -> sub={} out={}",
                            m_count, in->ident(), input_blobs.size(), ncut, nsub, sbs->m_blobs.size());
    }
    ++m_count;
    out = sbs;
    return true;
}
