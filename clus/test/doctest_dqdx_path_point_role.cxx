// doc pdvd/25 M3: TrackFitting::dQ_dx_fit's dx rule at an (apa, face) run
// boundary.  The former branch order treated a boundary at the very LAST
// trajectory point as a "first point" and read fine_tracking_path.at(i+1),
// which threw std::out_of_range (PDVD run 039252 evt 0: a fitted track whose
// last point sits in another CRP face; SBND's two one-sided anodes rarely
// place a boundary there).  The helper pins the rule: 0 first / 1 last /
// 2 middle / 3 isolated, identical to the old order wherever i+1 exists.

#include "WireCellUtil/doctest.h"
#include "WireCellClus/TrackFitting.h"

#include <utility>
#include <vector>

using namespace WireCell;
using namespace WireCell::Clus;
using paf_t = std::vector<std::pair<int, int>>;

TEST_CASE("dqdx path point role: the crash case -- a run boundary at the last point is 'last'")
{
    const paf_t paf{{0, 0}, {0, 0}, {1, 0}};
    CHECK(TrackFitting::dqdx_path_point_role(2, 3, paf) == 1);   // was 0 -> at(3) on a size-3 path
}

TEST_CASE("dqdx path point role: unchanged wherever i+1 exists")
{
    const paf_t same{{0, 0}, {0, 0}, {0, 0}, {0, 0}};
    CHECK(TrackFitting::dqdx_path_point_role(0, 4, same) == 0);
    CHECK(TrackFitting::dqdx_path_point_role(1, 4, same) == 2);
    CHECK(TrackFitting::dqdx_path_point_role(2, 4, same) == 2);
    CHECK(TrackFitting::dqdx_path_point_role(3, 4, same) == 1);

    const paf_t split{{0, 0}, {0, 0}, {1, 1}, {1, 1}};
    CHECK(TrackFitting::dqdx_path_point_role(1, 4, split) == 1);   // run ends before the boundary
    CHECK(TrackFitting::dqdx_path_point_role(2, 4, split) == 0);   // run starts at the boundary, i+1 exists
    CHECK(TrackFitting::dqdx_path_point_role(3, 4, split) == 1);

    const paf_t lone{{0, 0}, {2, 0}, {0, 0}};
    CHECK(TrackFitting::dqdx_path_point_role(1, 3, lone) == 0);    // single-point run in the middle: first, reads i+1
    CHECK(TrackFitting::dqdx_path_point_role(0, 1, lone) == 0);    // n == 1: the caller's n > 1 guard applies
}
