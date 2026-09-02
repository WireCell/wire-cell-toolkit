// Tests for inset_scope_fv(), the separation-scoped y/z fiducial inset behind
// ClusteringSeparate's `fv_inset_yz` knob (sbnd_xin/docs/97).
//
// WHY THE KNOB EXISTS.  JudgeSeparateDec_2 decides whether a cluster reaches a
// detector surface by comparing its convex-hull extremes against the scope's
// FV bounds -- for the y and z low/high faces, a bare `p.y() < FV_ymin` style
// test with no margin.  On SBND the configured FV is inset from the physical
// wall by only 0.65-0.85 cm, so that test effectively asks for a point ON the
// wall: a cosmic that stops a few centimetres short contributes no surface
// contact, and a two-track overcluster never shows the two independent
// surfaces the judge needs.  PDHD and PDVD avoid this by insetting
// FV_y*/FV_z* ~15 cm in the DetectorVolumes metadata
// (clus/docs/clustering-separate-fv.md), but that block is shared -- via
// select_scope_fv -- with clustering_neutrino and the containment taggers.
// inset_scope_fv applies the same inset to one pass and nothing else.
//
// Covers, in order:
//  * OFF is bit-identical: inset <= 0 returns every field unchanged, so the
//    knob's "byte-identical when off" guarantee is a property of the code and
//    not only of the A/B gates;
//  * the inset moves exactly the four y/z bounds, by exactly the amount asked,
//    and leaves x, every *_margin and the direction vectors alone;
//  * the mechanism, on the real SBND numbers and the real measured endpoints
//    of run 105 subrun 23 event 21 (docs/96 sec 5.2): with production FV that
//    cluster's two wall exits are invisible, with a 15 cm inset both are seen;
//  * the CAUSAL NEGATIVE CONTROL: feeding the same points back through a
//    zero-inset FV removes the effect entirely.  Corrupt exactly the quantity
//    the knob supplies and the new behavior must disappear.

#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Units.h"

#include "WireCellClus/ClusteringFuncs.h"

using namespace WireCell;
using namespace WireCell::Clus::Facade;

namespace {

    // The SBND per-APA fiducial block as compiled today (mm), read out of a
    // production Q/L config: dv-apa0-0 / a0f0pA.  The active volume is
    // y in [-1999.65, 1999.65], z in [0, 5010], so this FV is inset 6.53 mm in
    // y and 8.5 mm in z -- the "already inset" claim docs/96 sec 6.1 measures
    // to be too small to matter.
    ScopeFV sbnd_fv()
    {
        ScopeFV fv;
        fv.xmin = -2010.5;      fv.xmax = -25.0;
        fv.ymin = -1993.12;     fv.ymax = 1993.12;
        fv.zmin = 8.5;          fv.zmax = 5001.5;
        fv.xmin_margin = 20;    fv.xmax_margin = 20;
        fv.ymin_margin = 25;    fv.ymax_margin = 25;
        fv.zmin_margin = 30;    fv.zmax_margin = 30;
        return fv;
    }

    // The y/z surface-contact tests JudgeSeparateDec_2 actually performs on a
    // hull extreme (clustering_separate.cxx: `hy_points.at(0).y() > det_FV_ymax`
    // and the ymin/zmin/zmax twins).  No margin term on these four faces.
    int yz_contacts(const ScopeFV& fv, const geo_point_t& p)
    {
        int n = 0;
        if (p.y() > fv.ymax) ++n;
        if (p.y() < fv.ymin) ++n;
        if (p.z() > fv.zmax) ++n;
        if (p.z() < fv.zmin) ++n;
        return n;
    }

}  // namespace

TEST_CASE("clus sep fv inset off is bit identical")
{
    const ScopeFV fv = sbnd_fv();

    for (double off : {0.0, -1.0, -15 * units::cm}) {
        const ScopeFV out = inset_scope_fv(fv, off);
        CHECK(out.xmin == fv.xmin);
        CHECK(out.xmax == fv.xmax);
        CHECK(out.ymin == fv.ymin);
        CHECK(out.ymax == fv.ymax);
        CHECK(out.zmin == fv.zmin);
        CHECK(out.zmax == fv.zmax);
        CHECK(out.xmin_margin == fv.xmin_margin);
        CHECK(out.xmax_margin == fv.xmax_margin);
        CHECK(out.ymin_margin == fv.ymin_margin);
        CHECK(out.ymax_margin == fv.ymax_margin);
        CHECK(out.zmin_margin == fv.zmin_margin);
        CHECK(out.zmax_margin == fv.zmax_margin);
        CHECK(out.vertical_dir == fv.vertical_dir);
        CHECK(out.beam_dir == fv.beam_dir);
    }
}

TEST_CASE("clus sep fv inset moves exactly the four y z bounds")
{
    const ScopeFV fv = sbnd_fv();
    const double d = 15 * units::cm;
    const ScopeFV out = inset_scope_fv(fv, d);

    CHECK(out.ymin == doctest::Approx(fv.ymin + d));
    CHECK(out.ymax == doctest::Approx(fv.ymax - d));
    CHECK(out.zmin == doctest::Approx(fv.zmin + d));
    CHECK(out.zmax == doctest::Approx(fv.zmax - d));

    // x is a drift coordinate: an inset there would change the out-of-time
    // apparent-x test, which is a different mechanism.  It must not move.
    CHECK(out.xmin == fv.xmin);
    CHECK(out.xmax == fv.xmax);
    // The margins are added outward to reach the physical boundary; the inset
    // is about the inner bound only.
    CHECK(out.ymin_margin == fv.ymin_margin);
    CHECK(out.zmax_margin == fv.zmax_margin);
}

TEST_CASE("clus sep fv inset exposes the wall exits Dec_2 cannot see")
{
    const ScopeFV prod = sbnd_fv();

    // The two real wall exits of the 439 cm overcluster in run 105 subrun 23
    // event 21 (docs/96 sec 5.2), in mm: one 3.0 cm short of the y-min wall,
    // one 1.4 cm short of the z-max wall.  Both are genuine cosmic endpoints.
    const geo_point_t exit_y(-1000.0, -1969.3, 2500.0);
    const geo_point_t exit_z(-1000.0, 0.0, 4996.5);

    SUBCASE("production FV: neither exit is a surface contact")
    {
        CHECK(yz_contacts(prod, exit_y) == 0);
        CHECK(yz_contacts(prod, exit_z) == 0);
    }

    SUBCASE("15 cm inset: both exits count, so the two-surface test is reachable")
    {
        const ScopeFV inset = inset_scope_fv(prod, 15 * units::cm);
        CHECK(yz_contacts(inset, exit_y) == 1);
        CHECK(yz_contacts(inset, exit_z) == 1);
        CHECK(yz_contacts(inset, exit_y) + yz_contacts(inset, exit_z) == 2);
    }

    SUBCASE("causal negative control: remove the inset, lose the effect")
    {
        // Same points, same code path, inset set to the value that makes the
        // knob a no-op.  If this still reported contacts, the SUBCASE above
        // would be measuring something other than the inset.
        const ScopeFV zero = inset_scope_fv(prod, 0.0);
        CHECK(yz_contacts(zero, exit_y) == 0);
        CHECK(yz_contacts(zero, exit_z) == 0);
    }

    SUBCASE("a point already outside is unaffected by the inset")
    {
        const geo_point_t beyond(-1000.0, -2100.0, 2500.0);   // past the y-min wall
        CHECK(yz_contacts(prod, beyond) == 1);
        CHECK(yz_contacts(inset_scope_fv(prod, 15 * units::cm), beyond) == 1);
    }
}
