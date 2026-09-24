/** BlobCutting (wcfm doc 03): a blob with a wire-plane strip wider than
    `length_threshold` is bisected recursively into sub-blobs whose strips lie
    inside the parent's and are at most `length_threshold` wide; narrow blobs
    pass through as the same IBlob; sub-blob idents start at `ident_base` and
    restart at every frame boundary and at EOS; the parent's value is shared
    equally.  Revert-proven: with the per-frame ident reset removed, the
    second frame's idents continue from the first and CHECK_EQ(ids1, ids2)
    fails; with `length_threshold` ignored, the width check fails.
*/
#include "WireCellImg/GridTiling.h"
#include "WireCellImg/BlobCutting.h"
#include "WireCellAux/SimpleSlice.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleBlob.h"
#include "WireCellAux/Testing.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IWirePlane.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <algorithm>
#include <set>

using namespace WireCell;

// A slice whose activity is `nchan` contiguous channels in each plane of the face.
static ISlice::pointer make_slice(IAnodeFace::pointer face, int frame_ident, int slice_ident, size_t nchan)
{
    auto frame = std::make_shared<Aux::SimpleFrame>(frame_ident, 0.0);
    ISlice::map_t activity;
    for (const auto& plane : face->planes()) {
        const auto& chans = plane->channels();
        const size_t mid = chans.size() / 2;
        for (size_t i = mid; i < std::min(mid + nchan, chans.size()); ++i) {
            activity[chans[i]] = ISlice::value_t(1000.0, 10.0);
        }
    }
    return std::make_shared<Aux::SimpleSlice>(frame, slice_ident, 0.0,
                                              0.5 * units::microsecond, activity);
}

static IBlobSet::pointer tile(Img::GridTiling& gt, ISlice::pointer slice)
{
    IBlobSet::pointer out;
    REQUIRE(gt(slice, out));
    REQUIRE(out != nullptr);
    return out;
}

// Re-issue the tiled blobs with a known value so charge sharing is testable
// (GridTiling emits value 0).
static IBlobSet::pointer with_value(IBlobSet::pointer in, float value)
{
    IBlob::vector blobs;
    for (const auto& b : in->blobs()) {
        blobs.push_back(std::make_shared<Aux::SimpleBlob>(b->ident(), value, b->uncertainty(),
                                                          b->shape(), b->slice(), b->face()));
    }
    return std::make_shared<Aux::SimpleBlobSet>(in->ident(), in->slice(), blobs);
}

static int max_wire_width(const RayGrid::Blob& blob)
{
    int w = 0;
    for (const auto& s : blob.strips()) {
        if (s.layer < 2) continue;
        w = std::max(w, s.bounds.second - s.bounds.first);
    }
    return w;
}

// true if every wire-plane strip of `sub` lies inside the same-layer strip of `parent`
static bool inside(const RayGrid::Blob& sub, const RayGrid::Blob& parent)
{
    for (const auto& s : sub.strips()) {
        if (s.layer < 2) continue;
        bool found = false;
        for (const auto& p : parent.strips()) {
            if (p.layer != s.layer) continue;
            found = (s.bounds.first >= p.bounds.first && s.bounds.second <= p.bounds.second);
        }
        if (!found) return false;
    }
    return true;
}

static std::vector<int> idents_of(IBlobSet::pointer bs)
{
    std::vector<int> ids;
    for (const auto& b : bs->blobs()) ids.push_back(b->ident());
    return ids;
}

TEST_CASE("blobcutting splits wide blobs and passes narrow ones through")
{
    auto anodes = Testing::anodes("uboone");
    REQUIRE(anodes.size() > 0);
    auto face = anodes[0]->face(0);
    REQUIRE(face != nullptr);

    Img::GridTiling gt;
    auto gcfg = gt.default_configuration();
    gcfg["anode"] = "AnodePlane:0";
    gcfg["face"] = 0;
    gt.configure(gcfg);

    const int threshold = 20;
    const int ident_base = 1 << 20;
    Img::BlobCutting bc;
    auto ccfg = bc.default_configuration();
    CHECK_EQ(ccfg["length_threshold"].asInt(), 20);
    CHECK_EQ(ccfg["max_depth"].asInt(), 10);
    CHECK_EQ(ccfg["min_length"].asInt(), 2);
    CHECK_EQ(ccfg["ident_base"].asInt(), ident_base);
    ccfg["length_threshold"] = threshold;
    bc.configure(ccfg);

    // --- wide activity: 60 channels per plane -> at least one blob wider than the threshold
    auto wide = with_value(tile(gt, make_slice(face, 100, 0, 60)), 120.0f);
    REQUIRE(wide->blobs().size() > 0);
    int widest = 0;
    float value_in = 0;
    for (const auto& b : wide->blobs()) {
        widest = std::max(widest, max_wire_width(b->shape()));
        value_in += b->value();
    }
    REQUIRE(widest > threshold);

    IBlobSet::pointer cut;
    REQUIRE(bc(wide, cut));
    REQUIRE(cut != nullptr);
    CHECK_EQ(cut->ident(), wide->ident());
    CHECK(cut->slice() == wide->slice());
    CHECK(cut->blobs().size() > wide->blobs().size());

    std::set<int> seen;
    float value_out = 0;
    size_t nsub = 0;
    for (const auto& b : cut->blobs()) {
        CHECK(b->shape().valid());
        CHECK(max_wire_width(b->shape()) <= threshold);
        CHECK(seen.insert(b->ident()).second);  // unique idents
        value_out += b->value();
        // every output blob is inside some input blob (a sub-blob) or IS an input blob
        bool parent_found = false;
        for (const auto& p : wide->blobs()) {
            if (b == p) { parent_found = true; break; }
            if (inside(b->shape(), p->shape())) { parent_found = true; break; }
        }
        CHECK(parent_found);
        if (b->ident() >= ident_base) {
            ++nsub;
            CHECK(b->slice() == wide->slice());
            CHECK(b->face() == face);
        }
    }
    CHECK(nsub > 1);
    CHECK(value_out == doctest::Approx(value_in).epsilon(1e-4));  // charge shared, not lost

    // --- narrow activity: 12 channels -> nothing to cut, same IBlob pointers out
    auto narrow = tile(gt, make_slice(face, 100, 1, 12));
    REQUIRE(narrow->blobs().size() > 0);
    IBlobSet::pointer same;
    REQUIRE(bc(narrow, same));
    REQUIRE(same != nullptr);
    REQUIRE_EQ(same->blobs().size(), narrow->blobs().size());
    for (size_t i = 0; i < narrow->blobs().size(); ++i) {
        CHECK(same->blobs()[i] == narrow->blobs()[i]);
    }

    // --- a huge threshold leaves the wide set untouched too
    Img::BlobCutting bc_off;
    auto ocfg = bc_off.default_configuration();
    ocfg["length_threshold"] = 1000;
    bc_off.configure(ocfg);
    IBlobSet::pointer untouched;
    REQUIRE(bc_off(wide, untouched));
    REQUIRE_EQ(untouched->blobs().size(), wide->blobs().size());
    for (size_t i = 0; i < wide->blobs().size(); ++i) {
        CHECK(untouched->blobs()[i] == wide->blobs()[i]);
    }
}

TEST_CASE("blobcutting restarts sub-blob idents at each frame and at EOS")
{
    auto anodes = Testing::anodes("uboone");
    REQUIRE(anodes.size() > 0);
    auto face = anodes[0]->face(0);

    Img::GridTiling gt;
    auto gcfg = gt.default_configuration();
    gcfg["anode"] = "AnodePlane:0";
    gcfg["face"] = 0;
    gt.configure(gcfg);

    Img::BlobCutting bc;
    auto ccfg = bc.default_configuration();
    ccfg["length_threshold"] = 20;
    ccfg["ident_base"] = 5000;
    bc.configure(ccfg);

    auto s1 = make_slice(face, /*frame*/ 100, 0, 60);
    auto s2 = make_slice(face, /*frame*/ 200, 0, 60);   // a DIFFERENT event, no EOS between

    IBlobSet::pointer o1, o1b, o2, o3, eos;
    REQUIRE(bc(tile(gt, s1), o1));
    REQUIRE(bc(tile(gt, s1), o1b));   // same frame: idents continue
    REQUIRE(bc(tile(gt, s2), o2));    // new frame: idents restart
    const auto ids1 = idents_of(o1), ids1b = idents_of(o1b), ids2 = idents_of(o2);
    REQUIRE(ids1.size() > 1);
    CHECK_EQ(ids1.front(), 5000);
    CHECK_EQ(ids1, ids2);
    CHECK(ids1b != ids1);
    CHECK(ids1b.front() > ids1.back());

    // EOS: nullptr in, nullptr out, counter reset
    CHECK(bc(nullptr, eos));
    CHECK(eos == nullptr);
    REQUIRE(bc(tile(gt, s1), o3));
    CHECK_EQ(idents_of(o3), ids1);
}
