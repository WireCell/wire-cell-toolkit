// doc pdvd/120: the three selection rules of the beam-particle PR stage
// (which bundle, which cluster, which end) on plain rows -- no grouping,
// fitter or detector.  Written before BeamParticleFunctions.cxx existed
// (toolkit/CLAUDE.md: the test fails first).

#include "WireCellClus/BeamParticleFunctions.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell;
using namespace WireCell::Clus::PR;
using WireCell::units::cm;
using WireCell::units::us;

TEST_CASE("beam_particle_pick_bundle: window and brightest-flash rule")
{
    std::vector<BeamBundleRow> mains = {
        {10, 3, 100.0 * us, 200 * cm},   // in window, gid 3
        {11, 5, 100.5 * us, 50 * cm},    // in window, gid 5
        {12, 7, 200.0 * us, 900 * cm},   // out of window
        {13, -1, 100.2 * us, 30 * cm},   // in window but unmatched: ignored
    };
    std::vector<BeamFlashInfo> flashes = {{3, 1000.0, true}, {5, 5000.0, true}};

    SUBCASE("empty input")
    {
        auto p = beam_particle_pick_bundle({}, flashes, 99 * us, 101 * us);
        CHECK(p.gid == -1);
        CHECK(p.n_in_window == 0);
        CHECK(p.why == "none");
    }
    SUBCASE("gate off (lo >= hi) never picks")
    {
        auto p = beam_particle_pick_bundle(mains, flashes, 0, 0);
        CHECK(p.gid == -1);
        CHECK(p.why == "window-off");
    }
    SUBCASE("brighter flash beats longer main")
    {
        auto p = beam_particle_pick_bundle(mains, flashes, 99 * us, 101 * us);
        CHECK(p.gid == 5);
        CHECK(p.n_in_window == 2);
        CHECK(p.n_gids == 2);
        CHECK(p.why == "brightest");
    }
    SUBCASE("half-open window: t0 == hi is excluded")
    {
        auto p = beam_particle_pick_bundle(mains, flashes, 99 * us, 100.5 * us);
        CHECK(p.gid == 3);
        CHECK(p.n_in_window == 1);
    }
    SUBCASE("out-of-window mains are ignored")
    {
        auto p = beam_particle_pick_bundle(mains, flashes, 150 * us, 250 * us);
        CHECK(p.gid == 7);
        CHECK(p.n_in_window == 1);
        CHECK(p.why == "longest");   // gid 7 has no flash row
    }
    SUBCASE("no valid flash row -> longest main")
    {
        auto p = beam_particle_pick_bundle(mains, {}, 99 * us, 101 * us);
        CHECK(p.gid == 3);
        CHECK(p.why == "longest");
    }
    SUBCASE("equal pe -> longest main")
    {
        std::vector<BeamFlashInfo> eq = {{3, 1000.0, true}, {5, 1000.0, true}};
        auto p = beam_particle_pick_bundle(mains, eq, 99 * us, 101 * us);
        CHECK(p.gid == 3);
        CHECK(p.why == "longest");
    }
    SUBCASE("equal pe and length -> smallest gid")
    {
        std::vector<BeamBundleRow> same = {{20, 9, 100.0 * us, 10 * cm}, {21, 4, 100.0 * us, 10 * cm}};
        std::vector<BeamFlashInfo> eq = {{9, 7.0, true}, {4, 7.0, true}};
        auto p = beam_particle_pick_bundle(same, eq, 99 * us, 101 * us);
        CHECK(p.gid == 4);
    }
    SUBCASE("one bundle, two mains: the bundle's longest main is what counts")
    {
        std::vector<BeamBundleRow> two = {{30, 2, 100.0 * us, 10 * cm}, {31, 2, 100.0 * us, 80 * cm},
                                          {32, 6, 100.0 * us, 40 * cm}};
        auto p = beam_particle_pick_bundle(two, {}, 99 * us, 101 * us);
        CHECK(p.gid == 2);
        CHECK(p.n_in_window == 3);
        CHECK(p.n_gids == 2);
    }
}

TEST_CASE("beam_particle_pick_main: nearest to the nominal entry")
{
    std::vector<BeamMainCand> cands = {
        {4000151, 8 * cm, 167 * cm},
        {4000115, 50 * cm, 125 * cm},
        {4000082, 400 * cm, 98 * cm},
    };
    CHECK(beam_particle_pick_main({}, 0) == -1);
    CHECK(beam_particle_pick_main(cands, 0) == 0);
    SUBCASE("min_length filter")
    {
        CHECK(beam_particle_pick_main(cands, 150 * cm) == 0);
        CHECK(beam_particle_pick_main(cands, 170 * cm) == -1);
        std::vector<BeamMainCand> c2 = {{1, 1 * cm, 5 * cm}, {2, 30 * cm, 100 * cm}};
        CHECK(beam_particle_pick_main(c2, 10 * cm) == 1);
    }
    SUBCASE("tie on distance -> longer, then smaller id")
    {
        std::vector<BeamMainCand> t = {{9, 5 * cm, 10 * cm}, {3, 5 * cm, 20 * cm}, {1, 5 * cm, 20 * cm}};
        CHECK(beam_particle_pick_main(t, 0) == 2);
    }
}

TEST_CASE("beam_particle_choose_entry: the end nearer the nominal, direction on a tie")
{
    const Point nominal(110 * cm, 159 * cm, 0.6 * cm);
    const Vector beam_dir(-0.095, -0.704, 0.704);
    const Point a(118 * cm, 172 * cm, 2 * cm);     // the face end
    const Point b(53 * cm, 101 * cm, 138 * cm);    // deep inside

    SUBCASE("nearer end wins, either argument order")
    {
        auto c = beam_particle_choose_entry(a, b, nominal, beam_dir, 5 * cm);
        CHECK(c.entry == a);
        CHECK(c.exit == b);
        CHECK(c.dist == doctest::Approx((a - nominal).magnitude()));
        CHECK(c.cos_beam > 0.9);
        CHECK_FALSE(c.tie_by_dir);
        auto d = beam_particle_choose_entry(b, a, nominal, beam_dir, 5 * cm);
        CHECK(d.entry == a);
        CHECK(d.cos_beam == doctest::Approx(c.cos_beam));
    }
    SUBCASE("within tie_tol the direction decides")
    {
        // both ends equidistant from the nominal; only p->q runs along the beam
        const Point p(110 * cm, 169 * cm, 0.6 * cm);
        const Point q(110 * cm, 149 * cm, 0.6 * cm);
        auto c = beam_particle_choose_entry(q, p, nominal, beam_dir, 5 * cm);
        CHECK(c.tie_by_dir);
        CHECK(c.entry == p);   // p -> q is toward -y = along beam_dir
        CHECK(c.cos_beam > 0);
        auto d = beam_particle_choose_entry(q, p, nominal, beam_dir, 0);
        CHECK_FALSE(d.tie_by_dir);   // tolerance 0: pure distance, q listed first keeps q on an exact tie
        CHECK(d.entry == q);
        CHECK(d.cos_beam < 0);
    }
    SUBCASE("degenerate a == b")
    {
        auto c = beam_particle_choose_entry(a, a, nominal, beam_dir, 5 * cm);
        CHECK(c.entry == a);
        CHECK(c.cos_beam == doctest::Approx(0.0));
        CHECK(c.dist == doctest::Approx((a - nominal).magnitude()));
    }
}
