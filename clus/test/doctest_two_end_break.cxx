// doc sbnd_xin/docs/pr/48 sec 6 -- segment_two_end_break_scan: the two-end
// residual-range back-to-back junction scan.  Synthetic fit trajectories with
// the production-like 0.6 cm spacing and synthetic-but-shaped stopping
// templates (LinterpFunction components + a ParticleDataSet, wired through
// the factory exactly as production does).
//
// The physics contract under test: a single particle has at most one Bragg
// end; dQ/dx rising at BOTH ends of one segment implies a junction at the
// interior dip, located by the argmin of the joint two-arm stopping score.

#include "WireCellClus/PRSegment.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/ParticleDataSet.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <cmath>

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::PR;

namespace {

const double STEP = 0.6 * units::cm;

// dQ/dx (e/cm) vs residual range (cm) -- muon-like: Bragg rise toward r=0
// over a ~50e3 plateau.  The scan divides by units::cm itself.
double muon_like(double r_cm)
{
    return 50e3 + 130e3 * std::exp(-r_cm / 2.5);
}
// proton-like: steeper and hotter.
double proton_like(double r_cm)
{
    return 60e3 + 240e3 * std::exp(-r_cm / 1.5);
}

Clus::ParticleDataSet::pointer make_pdata()
{
    static Clus::ParticleDataSet::pointer pdata;
    if (pdata) return pdata;

    PluginManager& pm = PluginManager::instance();
    pm.add("WireCellAux");
    pm.add("WireCellClus");

    auto make_fn = [](const std::string& name, double (*shape)(double)) {
        auto icfg = Factory::lookup_tn<IConfigurable>("LinterpFunction:" + name);
        Configuration cfg = icfg->default_configuration();
        cfg["start"] = 0.0;
        cfg["step"] = 0.25;   // r in cm
        Configuration values(Json::arrayValue);
        for (int i = 0; i < 400; i++) values.append(shape(0.25 * i));
        cfg["values"] = values;
        icfg->configure(cfg);
    };
    make_fn("teb_muon", &muon_like);
    make_fn("teb_proton", &proton_like);
    // electron template: flat MIP-ish (never the best stopper).
    auto icfg = Factory::lookup_tn<IConfigurable>("LinterpFunction:teb_electron");
    Configuration cfg = icfg->default_configuration();
    cfg["start"] = 0.0;
    cfg["step"] = 50.0;
    Configuration values(Json::arrayValue);
    values.append(50e3);
    values.append(50e3);
    cfg["values"] = values;
    icfg->configure(cfg);

    auto pcfg_holder = Factory::lookup_tn<IConfigurable>("ParticleDataSet:teb_pdata");
    Configuration pcfg = pcfg_holder->default_configuration();
    pcfg["dedx_functions"]["muon"] = "LinterpFunction:teb_muon";
    pcfg["dedx_functions"]["proton"] = "LinterpFunction:teb_proton";
    pcfg["dedx_functions"]["electron"] = "LinterpFunction:teb_electron";
    pcfg_holder->configure(pcfg);
    pdata = std::dynamic_pointer_cast<Clus::ParticleDataSet>(pcfg_holder);
    return pdata;
}

// Straight segment along +z of total length L, with dQ/dx given per point by
// `profile(arc)` in e/cm.  Owned by a Graph so SegmentPtr semantics hold.
struct TestTrack {
    Graph graph;
    SegmentPtr seg;
};

std::shared_ptr<TestTrack> make_track(double L, double (*profile)(double, double), double par)
{
    auto tt = std::make_shared<TestTrack>();
    auto v1 = make_vertex(tt->graph);
    v1->wcpt().point = Point(0, 0, 0);
    auto v2 = make_vertex(tt->graph);
    v2->wcpt().point = Point(0, 0, L);
    tt->seg = make_segment(tt->graph, v1, v2);
    const int n = static_cast<int>(L / STEP) + 1;
    for (int i = 0; i < n; i++) {
        Fit f;
        const double arc = i * STEP;
        f.point = Point(0, 0, arc);
        f.index = i;
        f.dx = STEP;
        f.dQ = profile(arc, par) / units::cm * STEP;   // e/cm -> internal dQ over dx
        tt->seg->fits().push_back(f);
    }
    return tt;
}

// Total synthetic track length, set by each test before building its track.
double g_teb_total_L = 0;

// Two-Bragg profile: junction at arc = par.  Arm A stops at the segment
// START (residual range = arc), arm B stops at the END.  The junction
// itself carries a genuine local dQ/dx dip spanning ~2 fit samples -- the
// signature every measured event shows (57903: 0.81 MIP right at truth vs
// a ~2 MIP plateau; 56211: 2.37; 51513: 1.03) and the marker route R1
// localizes on.
double two_bragg(double arc, double junction)
{
    // 30e3 e/cm = 0.70x the 43e3 MIP median: a physical junction dip, above
    // the 0.6x dip_floor that vetoes instrumental (dead-region) dips.
    if (std::abs(arc - junction) < 0.9 * units::cm) return 30e3;
    if (arc <= junction) return muon_like(arc / units::cm);
    return proton_like((g_teb_total_L - arc) / units::cm);
}

// Single-Bragg (ordinary stopping muon entering at arc 0, stopping at the END).
double one_bragg(double arc, double /*unused*/)
{
    return muon_like((g_teb_total_L - arc) / units::cm);
}

// Flat MIP.
double flat_mip(double /*arc*/, double /*unused*/)
{
    return 50e3;
}

} // namespace

TEST_CASE("clus pr48 two-end break: locates a mid-track junction and fires R1")
{
    auto pdata = make_pdata();
    REQUIRE(pdata);
    const double L = 60 * units::cm;
    const double junction = 20 * units::cm;
    g_teb_total_L = L;
    auto tt = make_track(L, &two_bragg, junction);

    TwoEndBreakOptions opt;   // doc pr/48 defaults
    auto res = segment_two_end_break_scan(tt->seg, pdata, opt);

    REQUIRE(res.break_idx >= 0);
    const double located = tt->seg->fits()[res.break_idx].point.z();
    CHECK(std::abs(located - junction) < 2.0 * units::cm);
    CHECK(res.ratio_lo > 1.3);
    CHECK(res.ratio_hi > 1.3);
    CHECK(res.flagA);
    CHECK(res.flagB);
    CHECK(res.route1);
    CHECK(res.found);
    // Straight track: route R2's turn angle must NOT be what fired.
    CHECK(res.turn_deg < 5.0);
}

TEST_CASE("clus pr48 two-end break: junction near the far end still located")
{
    auto pdata = make_pdata();
    const double L = 25 * units::cm;
    const double junction = 21.5 * units::cm;   // far arm 3.5 cm -- the 56211 shape
    g_teb_total_L = L;
    auto tt = make_track(L, &two_bragg, junction);

    TwoEndBreakOptions opt;
    auto res = segment_two_end_break_scan(tt->seg, pdata, opt);

    REQUIRE(res.break_idx >= 0);
    const double located = tt->seg->fits()[res.break_idx].point.z();
    CHECK(std::abs(located - junction) < 2.0 * units::cm);
    CHECK(res.found);
}

TEST_CASE("clus pr48 two-end break: kinked weak-rise junction fires R2 at the bend")
{
    // The 57485 shape: a genuine wide-baseline bend at the junction, one
    // end's rise too weak for route R1's 1.3 floor but above R2's 1.15.
    auto pdata = make_pdata();
    const double L = 60 * units::cm;
    const double junction = 25 * units::cm;
    const double th = 45.0 * M_PI / 180.0;

    auto tt = std::make_shared<TestTrack>();
    auto v1 = make_vertex(tt->graph);
    v1->wcpt().point = Point(0, 0, 0);
    auto v2 = make_vertex(tt->graph);
    v2->wcpt().point = Point(0, (L - junction) * std::sin(th), junction + (L - junction) * std::cos(th));
    tt->seg = make_segment(tt->graph, v1, v2);
    const int n = static_cast<int>(L / STEP) + 1;
    for (int i = 0; i < n; i++) {
        Fit f;
        const double arc = i * STEP;
        if (arc <= junction) {
            f.point = Point(0, 0, arc);
        } else {
            f.point = Point(0, (arc - junction) * std::sin(th), junction + (arc - junction) * std::cos(th));
        }
        f.index = i;
        f.dx = STEP;
        // Arm A: strong muon Bragg toward the start.  Arm B: weak soft rise
        // toward the end (~1.25x interior -- fails R1's 1.3, passes R2's 1.15).
        double val;
        if (arc <= junction) val = muon_like(arc / units::cm);
        else val = 50e3 + 20e3 * std::exp(-(L - arc) / units::cm / 2.0);
        f.dQ = val / units::cm * STEP;
        tt->seg->fits().push_back(f);
    }

    TwoEndBreakOptions opt;
    auto res = segment_two_end_break_scan(tt->seg, pdata, opt);

    MESSAGE("r2case found=" << res.found << " r1=" << res.route1 << " r2=" << res.route2
            << " idx_dip=" << res.idx_dip << " idx_turn=" << res.idx_turn
            << " turn=" << res.turn_deg << " ratio=(" << res.ratio_lo << "," << res.ratio_hi
            << ") sA=" << res.sA << " sB=" << res.sB
            << " flagA=" << res.flagA << " flagB=" << res.flagB);

    CHECK(!res.route1);
    CHECK(res.route2);
    CHECK(res.found);
    REQUIRE(res.break_idx >= 0);
    // Located at the bend (the turn maximum), not at some plateau dip.
    const double arc_located = res.break_idx * STEP;
    CHECK(std::abs(arc_located - junction) < 3.0 * units::cm);
    CHECK(res.turn_deg > 40.0);
}

TEST_CASE("clus pr48 two-end break: single-Bragg stopping muon never fires")
{
    auto pdata = make_pdata();
    const double L = 60 * units::cm;
    g_teb_total_L = L;
    auto tt = make_track(L, &one_bragg, 0);

    TwoEndBreakOptions opt;
    auto res = segment_two_end_break_scan(tt->seg, pdata, opt);
    CHECK(!res.found);
}

TEST_CASE("clus pr48 two-end break: flat MIP track never fires")
{
    auto pdata = make_pdata();
    const double L = 60 * units::cm;
    g_teb_total_L = L;
    auto tt = make_track(L, &flat_mip, 0);

    TwoEndBreakOptions opt;
    auto res = segment_two_end_break_scan(tt->seg, pdata, opt);
    CHECK(!res.found);
    // The cheap pre-gate must have declined before any scan ran.
    CHECK(res.break_idx == -1);
}

TEST_CASE("clus pr48 two-end break: too-short segment and null inputs decline")
{
    auto pdata = make_pdata();
    const double L = 8 * units::cm;   // below min_len 10 cm
    g_teb_total_L = L;
    auto tt = make_track(L, &two_bragg, 4 * units::cm);

    TwoEndBreakOptions opt;
    CHECK(!segment_two_end_break_scan(tt->seg, pdata, opt).found);
    CHECK(!segment_two_end_break_scan(nullptr, pdata, opt).found);
    CHECK(!segment_two_end_break_scan(tt->seg, nullptr, opt).found);
}
