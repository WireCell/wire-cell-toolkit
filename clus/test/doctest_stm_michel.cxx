// doc pdvd/48 -- the graph-only logic behind CheckSTM_Michel on synthetic
// PR graphs: chain walk, residual-range profile, Bragg contrast and the
// Michel / continuation / delta / hadron arm classification.  No detector,
// fitter or particle dataset: the muon table is injected as a function.

#include "WireCellClus/PRSegment.h"
#include "WireCellClus/PRVertex.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/StmMichelFunctions.h"

#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <cmath>

using namespace WireCell;
using namespace WireCell::Clus::PR;

namespace {

const double STEP = 0.6 * units::cm;
const double MIP_MED = 43000 / units::cm;   // internal units, as PatternAlgorithms::m_mip_dqdx_median

VertexPtr make_vtx(Graph& g, double x_cm, double y_cm, double z_cm)
{
    auto v = make_vertex(g);
    v->wcpt().point = Point(x_cm * units::cm, y_cm * units::cm, z_cm * units::cm);
    return v;
}

// Straight segment v1 -> v2, flat dQ/dx = ratio x MIP_MED, fits every 0.6 cm.
// reverse_storage stores the fits from v2 to v1 (the profile must not care).
// bragg_last_cm > 0 multiplies dQ by bragg_factor over the last bragg_last_cm
// before v2.
SegmentPtr make_track(Graph& g, VertexPtr v1, VertexPtr v2, double ratio,
                      bool reverse_storage = false, double bragg_last_cm = 0, double bragg_factor = 1)
{
    auto seg = make_segment(g, v1, v2);
    const auto d = v2->wcpt().point - v1->wcpt().point;
    const double L = d.magnitude();
    const int n = static_cast<int>(L / STEP) + 1;
    std::vector<Fit> fits;
    for (int i = 0; i < n; i++) {
        Fit f;
        const double frac = (n > 1) ? double(i) / (n - 1) : 0.0;
        f.point = v1->wcpt().point + d * frac;
        f.index = i;
        f.dx = STEP;
        double q = ratio * MIP_MED * STEP;
        const double rr = L * (1 - frac);
        if (bragg_last_cm > 0 && rr <= bragg_last_cm * units::cm) q *= bragg_factor;
        f.dQ = q;
        fits.push_back(f);
    }
    if (reverse_storage) std::reverse(fits.begin(), fits.end());
    for (size_t i = 0; i < fits.size(); ++i) fits[i].index = static_cast<int>(i);
    seg->fits() = fits;
    return seg;
}

StmMichelArmThresholds thresholds()
{
    StmMichelArmThresholds th;
    th.mip_dqdx_median = MIP_MED;
    return th;
}

}  // namespace

TEST_CASE("stm_michel chain: shortest path takes the straight route over a delta detour")
{
    Graph g;
    auto entry = make_vtx(g, 0, 0, 0);
    auto mid   = make_vtx(g, 0, 0, 50);
    auto stop  = make_vtx(g, 0, 0, 100);
    auto side  = make_vtx(g, 0, 20, 60);
    auto s1 = make_track(g, entry, mid, 1.0);
    auto s2 = make_track(g, mid, stop, 1.0);
    auto d1 = make_track(g, mid, side, 1.0);     // a detour
    auto d2 = make_track(g, side, stop, 1.0);    // that also reaches the stop, longer
    auto chain = stm_michel_shortest_chain(g, entry, stop);
    REQUIRE(chain.size() == 2);
    CHECK(chain[0] == s1);
    CHECK(chain[1] == s2);
    auto vtxs = stm_michel_chain_vertices(g, chain, entry);
    REQUIRE(vtxs.size() == 3);
    CHECK(vtxs[0] == entry); CHECK(vtxs[1] == mid); CHECK(vtxs[2] == stop);

    // Insertion order of the detour must not matter.
    Graph g2;
    auto e2 = make_vtx(g2, 0, 0, 0);
    auto m2 = make_vtx(g2, 0, 0, 50);
    auto t2 = make_vtx(g2, 0, 0, 100);
    auto x2 = make_vtx(g2, 0, 20, 60);
    auto dd1 = make_track(g2, m2, x2, 1.0);
    auto dd2 = make_track(g2, x2, t2, 1.0);
    auto ss2 = make_track(g2, m2, t2, 1.0);
    auto ss1 = make_track(g2, e2, m2, 1.0);
    auto chain2 = stm_michel_shortest_chain(g2, e2, t2);
    REQUIRE(chain2.size() == 2);
    CHECK(chain2[0] == ss1);
    CHECK(chain2[1] == ss2);
}

TEST_CASE("stm_michel chain: unreachable stop returns empty")
{
    Graph g;
    auto a = make_vtx(g, 0, 0, 0);
    auto b = make_vtx(g, 0, 0, 30);
    auto c = make_vtx(g, 100, 0, 0);
    auto d = make_vtx(g, 100, 0, 30);
    make_track(g, a, b, 1.0);
    make_track(g, c, d, 1.0);
    CHECK(stm_michel_shortest_chain(g, a, d).empty());
    CHECK(stm_michel_shortest_chain(g, a, a).empty());
}

TEST_CASE("stm_michel farthest vertex: longest route out of the entry, predicate-restricted, ties on index")
{
    // entry -- 30 cm -- j -- 40 cm -- far ; j -- 10 cm -- delta ; a detached 100 cm piece
    Graph g;
    auto entry = make_vtx(g, 0, 0, 0);
    auto j     = make_vtx(g, 0, 0, 30);
    auto far   = make_vtx(g, 0, 0, 70);
    auto dl    = make_vtx(g, 0, 10, 30);
    auto x1    = make_vtx(g, 500, 0, 0);
    auto x2    = make_vtx(g, 500, 0, 100);
    auto s1 = make_track(g, entry, j, 1.0);
    auto s2 = make_track(g, j, far, 1.0);
    make_track(g, j, dl, 1.0);
    make_track(g, x1, x2, 1.0);
    auto v = stm_michel_farthest_vertex(g, entry);
    REQUIRE(v);
    CHECK(v == far);                       // 70 cm beats 40 cm; the detached piece is unreachable
    auto chain = stm_michel_shortest_chain(g, entry, v);
    REQUIRE(chain.size() == 2);
    CHECK(chain[0] == s1); CHECK(chain[1] == s2);
    // the predicate can veto the far vertex: the delta end is then the answer
    auto v2 = stm_michel_farthest_vertex(g, entry, [far](const VertexPtr& c) { return c != far; });
    CHECK(v2 == dl);
    // nothing reachable but the entry itself
    CHECK(!stm_michel_farthest_vertex(g, x1, [x2](const VertexPtr& c) { return c != x2; }));
    CHECK(!stm_michel_farthest_vertex(g, nullptr));
}

TEST_CASE("stm_michel live profile: dead points removed, geometry kept")
{
    Graph g;
    auto entry = make_vtx(g, 0, 0, 0);
    auto stop  = make_vtx(g, 0, 0, 30);
    auto s1 = make_track(g, entry, stop, 1.0);
    // starve every point with rr in [10, 20] cm to 5 % of MIP: a 10 cm stretch of
    // cells the fit could not read (the PDHD APA-edge / dead-channel pattern)
    for (auto& f : s1->fits()) {
        const double rr = (stop->wcpt().point - f.point).magnitude();
        if (rr >= 10 * units::cm && rr <= 20 * units::cm) f.dQ *= 0.05;
    }
    std::vector<SegmentPtr> chain{s1};
    auto prof = stm_michel_profile(g, chain, entry);
    REQUIRE(!prof.empty());
    int n_dead = -1;
    auto live = stm_michel_profile_live(prof, 0.15 * 43000.0, n_dead);
    CHECK(n_dead > 10);
    CHECK(live.L.size() + n_dead == prof.L.size());
    CHECK(live.total_length == doctest::Approx(prof.total_length));
    for (double q : live.dQdx) CHECK(q >= 0.15 * 43000.0);
    for (double rr : live.rr) CHECK((rr < 10 * units::cm || rr > 20 * units::cm));
    // the medians: with the dead stretch inside the plateau window [10, 20]
    // (short-track halving of [20, 40]) the RAW contrast is a fake; the live one is 1
    auto flat = [](double) { return 43000.0; };
    auto raw = stm_michel_bragg_contrast(prof, flat, 0.5 * units::cm, 3 * units::cm, 20 * units::cm, 40 * units::cm);
    auto lv  = stm_michel_bragg_contrast(live, flat, 0.5 * units::cm, 3 * units::cm, 20 * units::cm, 40 * units::cm);
    REQUIRE(raw.valid);
    CHECK(raw.contrast > 5.0);        // plateau median ~0
    CHECK(!lv.valid);                 // no live plateau points at all -> cannot judge, not "no Bragg"
    // frac = 0 keeps everything (the doc pdvd/48 behaviour)
    int nd0 = -1;
    auto same = stm_michel_profile_live(prof, 0.0, nd0);
    CHECK(nd0 == 0);
    CHECK(same.L.size() == prof.L.size());
}

TEST_CASE("stm_michel profile: L monotone, rr = 0 at the stop, storage order irrelevant")
{
    Graph g;
    auto entry = make_vtx(g, 0, 0, 0);
    auto mid   = make_vtx(g, 0, 0, 30);
    auto stop  = make_vtx(g, 0, 0, 60);
    auto s1 = make_track(g, entry, mid, 1.0, /*reverse_storage=*/true);
    auto s2 = make_track(g, mid, stop, 1.0, /*reverse_storage=*/false);
    std::vector<SegmentPtr> chain{s1, s2};
    auto prof = stm_michel_profile(g, chain, entry);
    REQUIRE(!prof.empty());
    for (size_t i = 1; i < prof.L.size(); ++i) CHECK(prof.L[i] >= prof.L[i - 1]);
    CHECK(prof.rr.back() == doctest::Approx(0.0));
    CHECK(prof.total_length == doctest::Approx(60 * units::cm).epsilon(0.02));
    // the first point sits at the entry, the last at the stop
    CHECK((prof.pts.front() - entry->wcpt().point).magnitude() < 0.1 * units::cm);
    CHECK((prof.pts.back() - stop->wcpt().point).magnitude() < 0.1 * units::cm);
    // dQ/dx is in e/cm: 1.0 x 43000
    CHECK(prof.dQdx.front() == doctest::Approx(43000.0).epsilon(1e-6));
}

TEST_CASE("stm_michel bragg contrast: ramped tail vs flat plateau; short track flagged")
{
    Graph g;
    auto entry = make_vtx(g, 0, 0, 0);
    auto stop  = make_vtx(g, 0, 0, 60);
    auto s = make_track(g, entry, stop, 1.0, false, /*bragg_last_cm=*/3.0, /*factor=*/3.0);
    auto prof = stm_michel_profile(g, {s}, entry);
    // an injected "muon table": 2.5x over the last 3 cm, else flat
    auto table = [](double rr_cm) { return rr_cm <= 3.0 ? 2.5 * 43000.0 : 43000.0; };
    auto b = stm_michel_bragg_contrast(prof, table, 0.5 * units::cm, 3 * units::cm, 20 * units::cm, 40 * units::cm);
    REQUIRE(b.valid);
    CHECK_FALSE(b.short_track);
    CHECK(b.contrast == doctest::Approx(3.0).epsilon(0.01));
    CHECK(b.expected == doctest::Approx(2.5).epsilon(0.01));

    // a flat track has contrast ~1 and fails a 60 % expected-rise bar
    auto flat = make_track(g, entry, stop, 1.0);
    auto pf = stm_michel_profile(g, {flat}, entry);
    auto bf = stm_michel_bragg_contrast(pf, table, 0.5 * units::cm, 3 * units::cm, 20 * units::cm, 40 * units::cm);
    REQUIRE(bf.valid);
    CHECK(bf.contrast == doctest::Approx(1.0).epsilon(0.01));
    CHECK(bf.contrast < 0.6 * bf.expected);

    // 25 cm track: the plateau window halves and short_track is set
    Graph g2;
    auto e2 = make_vtx(g2, 0, 0, 0);
    auto t2 = make_vtx(g2, 0, 0, 25);
    auto s2 = make_track(g2, e2, t2, 1.0, false, 3.0, 3.0);
    auto p2 = stm_michel_profile(g2, {s2}, e2);
    auto b2 = stm_michel_bragg_contrast(p2, table, 0.5 * units::cm, 3 * units::cm, 20 * units::cm, 40 * units::cm);
    CHECK(b2.short_track);
    CHECK(b2.valid);
}

TEST_CASE("stm_michel stop arm: short wide-angle low-MIP arm is Michel, collinear MIP arm is continuation")
{
    Graph g;
    auto entry = make_vtx(g, 0, 0, 0);
    auto stop  = make_vtx(g, 0, 0, 60);
    auto muon = make_track(g, entry, stop, 1.0, false, 3.0, 3.0);
    // Michel: 12 cm at 60 deg off the muon axis, 0.9 MIP
    auto vm = make_vtx(g, 0, 12 * std::sin(60 * M_PI / 180), 60 + 12 * std::cos(60 * M_PI / 180));
    auto michel = make_track(g, stop, vm, 0.9);
    // continuation: 20 cm straight on, 1.0 MIP
    auto vc = make_vtx(g, 0, 0, 80);
    auto cont = make_track(g, stop, vc, 1.0);
    auto th = thresholds();
    auto am = stm_michel_classify_stop_arm(g, muon, michel, stop, th);
    CHECK(am.kind == StmMichelArm::kMichel);
    CHECK(am.kink_deg == doctest::Approx(60.0).epsilon(0.05));
    CHECK(am.terminal);
    auto ac = stm_michel_classify_stop_arm(g, muon, cont, stop, th);
    CHECK(ac.kind == StmMichelArm::kContinuation);
    CHECK(ac.kink_deg < 5.0);
    // a 6 cm collinear MIP piece is a continuation, never a Michel (doc pdvd/42 sec 4.4)
    Graph g3;
    auto e3 = make_vtx(g3, 0, 0, 0);
    auto s3 = make_vtx(g3, 0, 0, 60);
    auto mu3 = make_track(g3, e3, s3, 1.0, false, 3.0, 3.0);
    auto c3 = make_vtx(g3, 0, 0, 66);
    auto piece = make_track(g3, s3, c3, 0.9);
    auto ap = stm_michel_classify_stop_arm(g3, mu3, piece, s3, th);
    CHECK(ap.kind == StmMichelArm::kContinuation);
    CHECK(ap.kind != StmMichelArm::kMichel);
    // doc pdhd/03: a collinear 20 cm MIP arm stays a continuation even when the
    // track/shower separation flagged it shower-like
    cont->set_flags(SegmentFlags::kShowerTrajectory);
    auto acs = stm_michel_classify_stop_arm(g, muon, cont, stop, th);
    CHECK(acs.kind == StmMichelArm::kContinuation);
    CHECK(acs.shower_like);
    cont->unset_flags(SegmentFlags::kShowerTrajectory);
    // ... while a shower-flagged arm at 0.4 MIP (outside the continuation band) is a Michel
    auto vlo = make_vtx(g, 0, 0, 78);
    auto lo = make_track(g, stop, vlo, 0.4);
    lo->set_flags(SegmentFlags::kShowerTrajectory);
    auto alo = stm_michel_classify_stop_arm(g, muon, lo, stop, th);
    CHECK(alo.kind == StmMichelArm::kMichel);
    // doc pdhd/03 sec 6.8: a collinear 5 cm 1.6-MIP shower-flagged stub is a Michel under the
    // doc pdvd/48 rule (flag alone) and Other once michel_shower_min_kink_deg = 15
    auto vst = make_vtx(g, 0, 0, 65);
    auto stubseg = make_track(g, stop, vst, 1.6);
    stubseg->set_flags(SegmentFlags::kShowerTrajectory);
    auto ast0 = stm_michel_classify_stop_arm(g, muon, stubseg, stop, th);
    CHECK(ast0.kind == StmMichelArm::kMichel);
    auto th15 = th; th15.michel_shower_min_kink_deg = 15;
    auto ast1 = stm_michel_classify_stop_arm(g, muon, stubseg, stop, th15);
    CHECK(ast1.kind == StmMichelArm::kOther);
    CHECK(stm_michel_classify_stop_arm(g, muon, michel, stop, th15).kind == StmMichelArm::kMichel);   // 60 deg still a Michel
    // a 60 cm straight track at 60 deg is neither (too long for a Michel, not collinear)
    auto vo = make_vtx(g, 0, 60 * std::sin(60 * M_PI / 180), 60 + 60 * std::cos(60 * M_PI / 180));
    auto other = make_track(g, stop, vo, 1.0);
    auto ao = stm_michel_classify_stop_arm(g, muon, other, stop, th);
    CHECK(ao.kind == StmMichelArm::kOther);
}

TEST_CASE("stm_michel chain arm: 5 cm terminal arm is delta, 15 cm 1.6x MIP branching arm is hadron")
{
    Graph g;
    auto entry = make_vtx(g, 0, 0, 0);
    auto mid   = make_vtx(g, 0, 0, 40);
    auto stop  = make_vtx(g, 0, 0, 80);
    auto s1 = make_track(g, entry, mid, 1.0);
    auto s2 = make_track(g, mid, stop, 1.0);
    auto vd = make_vtx(g, 0, 5, 40);
    auto delta = make_track(g, mid, vd, 1.2);
    auto vh = make_vtx(g, 0, -15, 40);
    auto vh2 = make_vtx(g, 0, -20, 45);
    auto hadron = make_track(g, mid, vh, 1.6);
    make_track(g, vh, vh2, 1.5);   // it branches on
    auto th = thresholds();
    auto ad = stm_michel_classify_chain_arm(g, s1, delta, mid, th);
    CHECK(ad.kind == StmMichelArm::kDelta);
    CHECK(ad.terminal);
    auto ah = stm_michel_classify_chain_arm(g, s1, hadron, mid, th);
    CHECK(ah.kind == StmMichelArm::kHadron);
    CHECK_FALSE(ah.terminal);
    (void)s2;
}
