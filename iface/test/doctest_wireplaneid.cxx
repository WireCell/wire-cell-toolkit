#include "WireCellIface/WirePlaneId.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell;

TEST_CASE("iface wireplaneid") {

    WirePlaneId wpid(WirePlaneLayer_t::kAllLayers, 0, 0);
    CHECK(wpid.layer() == WirePlaneLayer_t::kAllLayers);
    CHECK(wpid == true);

    auto wpid_u = wpid.to_u();
    CHECK(wpid_u.index() == 0);
    CHECK(wpid_u.valid());
    CHECK(wpid_u == true);
}
