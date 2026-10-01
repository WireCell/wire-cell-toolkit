// wcp-porting-img wcfm/docs/25: TrackFitting electron_lifetime_allow_unmatched
// lets the electron-lifetime correction (icarus/docs/04) apply to clusters with
// no matched flash -- light-less simulation, where cluster_t0 is the true t0.
// Off (0) by default, so every detector that sets electron_lifetime keeps the
// flash requirement.  The correction itself is the pure
// electron_lifetime_factor tested in doctest_icarus_pr_knobs.cxx; the
// unmatched path reads the same anode and drift speed (smoke-tested on FD-HD
// event 1311, wcfm/docs/25).

#include "WireCellUtil/doctest.h"

#include "WireCellClus/TrackFitting.h"
#include "WireCellUtil/Units.h"

using namespace WireCell;
using namespace WireCell::Clus;

TEST_CASE("wcfm/25 electron_lifetime_allow_unmatched: off by default, round-trips")
{
    TrackFitting tf;
    CHECK(tf.get_parameter("electron_lifetime_allow_unmatched") == 0.0);
    tf.set_parameter("electron_lifetime_allow_unmatched", 1.0);
    CHECK(tf.get_parameter("electron_lifetime_allow_unmatched") == 1.0);
    tf.set_parameter("electron_lifetime_allow_unmatched", 0.0);
    CHECK(tf.get_parameter("electron_lifetime_allow_unmatched") == 0.0);
    // the lifetime itself stays off until set
    CHECK(tf.get_parameter("electron_lifetime") == 0.0);
}
