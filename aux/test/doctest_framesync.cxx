// FrameSync is a hydra.  It must consume the input it buffers so that a
// caller honoring the IHydraNode contract (leftover input persists, eg the
// Pgraph Hydra wrapper) does not hand it the same frames again.

#include "WireCellAux/FrameSync.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellUtil/doctest.h"

#include <boost/any.hpp>

using namespace WireCell;

TEST_CASE("aux framesync consumes buffered input")
{
    Aux::FrameSync fsync(2);
    IHydraNodeBase& fs = fsync;  // the "any" interface a graph engine calls

    IHydraNodeBase::any_queue_vector inqs(2), outqs;
    IFrame::pointer f1 = std::make_shared<Aux::SimpleFrame>(1);
    IFrame::pointer f2 = std::make_shared<Aux::SimpleFrame>(2);
    inqs[0].push_back(f1);
    inqs[1].push_back(f2);

    REQUIRE(fs(inqs, outqs));
    CHECK(inqs[0].empty());
    CHECK(inqs[1].empty());
    REQUIRE(outqs.size() == 1);
    REQUIRE(outqs[0].size() == 2);
    CHECK(boost::any_cast<IFrame::pointer>(outqs[0][0])->ident() == 1);
    CHECK(boost::any_cast<IFrame::pointer>(outqs[0][1])->ident() == 2);

    // Call again with whatever input the caller still holds: nothing new.
    IHydraNodeBase::any_queue_vector outqs2;
    REQUIRE(fs(inqs, outqs2));
    REQUIRE(outqs2.size() == 1);
    CHECK(outqs2[0].empty());
}
