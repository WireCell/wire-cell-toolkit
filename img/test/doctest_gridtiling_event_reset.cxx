/** doc 81 -- GridTiling must restart its blob-ident sequence at every frame
    (event) boundary, not only at EOS.

    A process that streams MANY events through one graph (doc 76 round 2 group
    mode) gets no EOS between them, so the `m_blobs_seen` counter used to keep
    counting and event N's blobs were identified from an offset instead of from
    0.  That is not cosmetic: blob idents key the unordered_map<int,...> /
    unordered_set<int> containers in InSliceDeghosting, ProjectionDeghosting,
    BlobGrouping and LocalGeomClustering, whose ITERATION ORDER depends on the
    key values, so an offset flips order-dependent deghosting decisions and
    changes blob COUNTS and CHARGES on some events (doc 81 round 0 measured 126
    of 672 differing npz members over a 16-event nueCC48 group).

    The test feeds the SAME activity twice through ONE GridTiling instance
    under two different frame idents and requires the two blob-ident sets to be
    identical.  Revert-proven: with the reset in GridTiling.cxx removed, the
    second frame's idents are offset by the first frame's blob count and
    CHECK_EQ(id1, id2) fails.
*/
#include "WireCellImg/GridTiling.h"
#include "WireCellAux/SimpleSlice.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/Testing.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IWirePlane.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell;

// A slice whose activity is a short contiguous run of channels in each plane
// of face 0 -- enough for the three planes to cross and make blobs, small
// enough that the test stays fast.
static ISlice::pointer make_slice(IAnodeFace::pointer face, int frame_ident, int slice_ident)
{
    auto frame = std::make_shared<Aux::SimpleFrame>(frame_ident, 0.0);
    ISlice::map_t activity;
    for (const auto& plane : face->planes()) {
        const auto& chans = plane->channels();
        const size_t mid = chans.size() / 2;
        for (size_t i = mid; i < std::min(mid + 12, chans.size()); ++i) {
            activity[chans[i]] = ISlice::value_t(1000.0, 10.0);
        }
    }
    return std::make_shared<Aux::SimpleSlice>(frame, slice_ident, 0.0,
                                              0.5 * units::microsecond, activity);
}

static std::vector<int> tile(Img::GridTiling& gt, ISlice::pointer slice)
{
    IBlobSet::pointer out;
    REQUIRE(gt(slice, out));
    REQUIRE(out != nullptr);
    std::vector<int> idents;
    for (const auto& iblob : out->blobs()) idents.push_back(iblob->ident());
    return idents;
}

TEST_CASE("doc81 gridtiling restarts blob idents at each frame boundary")
{
    auto anodes = Testing::anodes("uboone");
    REQUIRE(anodes.size() > 0);
    auto anode = anodes[0];
    auto face = anode->face(0);
    REQUIRE(face != nullptr);

    // Register the anode under a known type:name so GridTiling can find it.
    Img::GridTiling gt;
    auto cfg = gt.default_configuration();
    cfg["anode"] = "AnodePlane:0";
    cfg["face"] = 0;
    gt.configure(cfg);

    auto s1 = make_slice(face, /*frame*/ 100, /*slice*/ 0);
    auto s2 = make_slice(face, /*frame*/ 200, /*slice*/ 0);   // a DIFFERENT event

    const auto id1 = tile(gt, s1);
    const auto id2 = tile(gt, s2);

    // Non-vacuous: the activity really does tile into blobs.
    REQUIRE(id1.size() > 0);
    CHECK_EQ(id1.size(), id2.size());

    // The whole point: identical input under a new frame ident gives identical
    // blob idents.  Without the per-frame reset, id2 == id1 + id1.size().
    CHECK_EQ(id1, id2);
    CHECK_EQ(id1.front(), 0);
    CHECK_EQ(id2.front(), 0);

    // EOS still resets, as before.
    IBlobSet::pointer out;
    CHECK(gt(nullptr, out));
    const auto id3 = tile(gt, s1);
    CHECK_EQ(id3, id1);
}
