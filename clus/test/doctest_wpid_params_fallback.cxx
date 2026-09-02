// doc pdvd/25 M3: the path-skeleton builder's wire-plane parameter lookup
// must not abort when an interpolated point lands in a volume the cluster has
// no blobs in.  PDVD tracks routinely cross CRP / face boundaries, and
// IDetectorVolumes::contained_by() then names a wpid that is not a key of the
// cluster's own wpid_params; the former std::map::at() threw std::out_of_range
// and the whole PR job aborted (wire-cell rc=134 on PDVD run 039252 evt 0).
//
// The first case documents the failure the fix removes; the others pin the
// helper's contract: exact key, then the same (apa, face) under kAllLayers
// (how blob wpids are stored), then the caller's fallback.  On any run that
// never hit a missing key the result is the exact-key value, i.e. unchanged.

#include "WireCellUtil/doctest.h"
#include "WireCellClus/DynamicPointCloud.h"

#include <map>
#include <stdexcept>
#include <tuple>

using namespace WireCell;
using namespace WireCell::Clus;
using params_t = std::map<WirePlaneId, std::tuple<Facade::geo_point_t, double, double, double>>;

static params_t make_params()
{
    params_t p;
    p[WirePlaneId(kAllLayers, 0, 0)] = std::make_tuple(Facade::geo_point_t(1, 0, 0), 0.10, 0.20, 0.30);
    p[WirePlaneId(kAllLayers, 1, 0)] = std::make_tuple(Facade::geo_point_t(1, 0, 0), 0.11, 0.21, 0.31);
    p[WirePlaneId(kUlayer, 0, 3)] = std::make_tuple(Facade::geo_point_t(-1, 0, 0), 0.13, 0.23, 0.33);
    return p;
}

TEST_CASE("wpid params: the old exact-key lookup throws on a foreign volume")
{
    auto p = make_params();
    const WirePlaneId foreign(kUlayer, 1, 5);   // apa 5: no blobs of this cluster
    CHECK_THROWS_AS((void) p.at(foreign), std::out_of_range);
}

TEST_CASE("wpid params: resolve_wpid_params never throws and prefers the exact key")
{
    auto p = make_params();
    const WirePlaneId fb(kAllLayers, 0, 0);

    SUBCASE("exact key") {
        const auto& t = Facade::resolve_wpid_params(p, WirePlaneId(kUlayer, 0, 3), fb);
        CHECK(std::get<1>(t) == doctest::Approx(0.13));
    }
    SUBCASE("same (apa, face) under kAllLayers when contained_by answers another layer") {
        const auto& t = Facade::resolve_wpid_params(p, WirePlaneId(kWlayer, 1, 0), fb);
        CHECK(std::get<1>(t) == doctest::Approx(0.11));
    }
    SUBCASE("foreign volume falls back to the caller's key") {
        const auto& t = Facade::resolve_wpid_params(p, WirePlaneId(kUlayer, 1, 5), fb);
        CHECK(std::get<1>(t) == doctest::Approx(0.10));
    }
    SUBCASE("unknown apa/face falls back too") {
        const auto& t = Facade::resolve_wpid_params(p, WirePlaneId(kUnknownLayer, -1, -1), fb);
        CHECK(std::get<1>(t) == doctest::Approx(0.10));
    }
}
