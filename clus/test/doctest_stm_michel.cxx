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
#include <limits>

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

TEST_CASE("stm_michel michel gate (doc pdvd/73 P2): all P2 fields off = the doc pdvd/48 expression")
{
    for (double smk : {-1.0, 15.0}) {
        auto th = thresholds();
        th.michel_shower_min_kink_deg = smk;
        int n = 0, mismatch = 0;
        for (double len : {0.0, 5.0, 10.0, 20.0, 25.0, 30.0})
        for (double far : {0.0, 5.0, 15.0, 30.0, 70.0})
        for (double mip : {0.1, 0.2, 0.3, 0.31, 1.0, 1.99, 2.0, 2.5})
        for (double kink : {-1.0, 10.0, 14.9, 15.0, 29.9, 30.0, 59.9, 60.0, 90.0})
        for (int sh : {0, 1}) {
            const double L = len * units::cm, F = far * units::cm;
            const bool kink_ok = kink >= 0;
            const bool shower_admits = sh && (th.michel_shower_min_kink_deg < 0 || !kink_ok ||
                                              kink >= th.michel_shower_min_kink_deg);
            const bool legacy = L + F <= th.michel_max_len && mip > th.michel_mip_lo && mip < th.michel_mip_hi &&
                                (shower_admits || (kink_ok && kink >= th.michel_min_kink_deg));
            // kink_w is ignored while (c) is off
            if (stm_michel_michel_gate(L, F, mip, kink, 45.0, sh, th) != legacy) ++mismatch;
            ++n;
        }
        CHECK(n == 6 * 5 * 8 * 9 * 2);
        CHECK(mismatch == 0);
    }
}

TEST_CASE("stm_michel michel gate (doc pdvd/73 P2): the three operating points and their boundaries")
{
    auto th = thresholds();
    th.michel_shower_min_kink_deg = 15;   // PDVD production
    const double cm = units::cm;
    // (a) 039349_11/19's arm 19004: 8.7 cm, 0.22 MIP, 81 deg -- under the 0.3 floor
    CHECK_FALSE(stm_michel_michel_gate(8.7 * cm, 0, 0.22, 81.1, -1, false, th));
    auto ta = th; ta.michel_mip_lo_turned = 0.15;
    CHECK(stm_michel_michel_gate(8.7 * cm, 0, 0.22, 81.1, -1, false, ta));
    CHECK(stm_michel_michel_gate(8.7 * cm, 0, 0.22, 60.0, -1, false, ta));         // the turn is inclusive
    CHECK_FALSE(stm_michel_michel_gate(8.7 * cm, 0, 0.22, 59.9, -1, false, ta));   // not turned: the 0.3 floor stands
    CHECK_FALSE(stm_michel_michel_gate(8.7 * cm, 0, 0.22, -1, -1, true, ta));      // an unmeasurable turn never lowers it
    CHECK_FALSE(stm_michel_michel_gate(8.7 * cm, 0, 0.15, 81.1, -1, false, ta));   // the lowered floor is strict
    // (b) a shower-flagged 8.7 cm arm carrying a 54 cm subtree
    CHECK_FALSE(stm_michel_michel_gate(8.7 * cm, 54 * cm, 0.5, 81.1, -1, true, th));
    auto tb = th; tb.michel_far_len_shower_max = 60 * cm;
    CHECK(stm_michel_michel_gate(8.7 * cm, 54 * cm, 0.5, 81.1, -1, true, tb));
    CHECK(stm_michel_michel_gate(8.7 * cm, 60 * cm, 0.5, 81.1, -1, true, tb));         // inclusive
    CHECK_FALSE(stm_michel_michel_gate(8.7 * cm, 60.1 * cm, 0.5, 81.1, -1, true, tb));
    CHECK_FALSE(stm_michel_michel_gate(26 * cm, 0, 0.5, 81.1, -1, true, tb));          // the arm itself stays <= michel_max_len
    CHECK_FALSE(stm_michel_michel_gate(8.7 * cm, 54 * cm, 0.5, 81.1, -1, false, tb));  // not shower-flagged: len + far_len
    // (c) 039253_3/61's arm 61007: 7.8 cm, 0.62 MIP, 18 deg over the classifier window
    CHECK_FALSE(stm_michel_michel_gate(7.8 * cm, 7.1 * cm, 0.62, 18.0, 40.0, false, th));   // off: kink_w ignored
    auto tc = th; tc.michel_kink_window = 5 * cm;
    CHECK(stm_michel_michel_gate(7.8 * cm, 7.1 * cm, 0.62, 18.0, 40.0, false, tc));
    CHECK(stm_michel_michel_gate(7.8 * cm, 7.1 * cm, 0.62, 18.0, 30.0, false, tc));          // inclusive
    CHECK_FALSE(stm_michel_michel_gate(7.8 * cm, 7.1 * cm, 0.62, 18.0, 29.9, false, tc));
    CHECK_FALSE(stm_michel_michel_gate(7.8 * cm, 7.1 * cm, 0.62, 18.0, -1, false, tc));
    // none of them rescues a hot arm
    auto tall = th;
    tall.michel_mip_lo_turned = 0.15; tall.michel_far_len_shower_max = 60 * cm; tall.michel_kink_window = 5 * cm;
    CHECK_FALSE(stm_michel_michel_gate(8.7 * cm, 0, 2.0, 81.1, 81.1, true, tall));
}

TEST_CASE("stm_michel stop arm (doc pdvd/73 P2): a turned 0.2 MIP arm; a shower arm's subtree, fenced at the stop")
{
    Graph g;
    auto entry = make_vtx(g, 0, 0, 0);
    auto stop  = make_vtx(g, 0, 0, 60);
    auto muon = make_track(g, entry, stop, 1.0, false, 3.0, 3.0);
    auto th = thresholds();
    th.michel_shower_min_kink_deg = 15;
    // (a) 8 cm at 90 deg, 0.2 MIP, terminal
    auto va = make_vtx(g, 0, 8, 60);
    auto arm_a = make_track(g, stop, va, 0.2);
    CHECK(stm_michel_classify_stop_arm(g, muon, arm_a, stop, th).kind == StmMichelArm::kOther);
    auto ta = th; ta.michel_mip_lo_turned = 0.15;
    auto aa = stm_michel_classify_stop_arm(g, muon, arm_a, stop, ta);
    CHECK(aa.kind == StmMichelArm::kMichel);
    CHECK(aa.kink_deg == doctest::Approx(90.0).epsilon(0.02));
    CHECK(aa.kink_w_deg == -1);   // (c) off: not measured
    auto tc = th; tc.michel_kink_window = 5 * units::cm;
    CHECK(stm_michel_classify_stop_arm(g, muon, arm_a, stop, tc).kink_w_deg == doctest::Approx(90.0).epsilon(0.02));
    // (b) a shower-flagged 8 cm arm at 90 deg the other way; its far vertex
    // carries a 40 cm branch and a two-segment loop back into the stop
    auto vb = make_vtx(g, 0, -8, 60);
    auto arm_b = make_track(g, stop, vb, 0.6);
    arm_b->set_flags(SegmentFlags::kShowerTrajectory);
    auto vb2 = make_vtx(g, 0, -48, 60);
    make_track(g, vb, vb2, 0.6);
    auto vl = make_vtx(g, 0, -4, 64);
    make_track(g, vb, vl, 0.6);
    make_track(g, vl, stop, 0.6);
    const double hop = std::sqrt(32.0) * units::cm;   // vb -> vl and vl -> stop
    // fenced: the branch and the first hop of the loop, never the muon behind the stop
    CHECK(stm_michel_far_subtree_len(g, vb, arm_b, stop, 1e9) == doctest::Approx(40 * units::cm + hop));
    // the unfenced walk goes through the stop into the muon chain
    CHECK(segment_far_subtree_track_length(g, vb, arm_b, 1e9) > 40 * units::cm + hop + 50 * units::cm);
    auto b0 = stm_michel_classify_stop_arm(g, muon, arm_b, stop, th);
    CHECK_FALSE(b0.terminal);
    CHECK(b0.kind == StmMichelArm::kOther);        // len + far_len > 25 cm
    auto tb = th; tb.michel_far_len_shower_max = 60 * units::cm;
    auto b1 = stm_michel_classify_stop_arm(g, muon, arm_b, stop, tb);
    CHECK(b1.kind == StmMichelArm::kMichel);
    CHECK(b1.far_len == doctest::Approx(40 * units::cm + hop));
    auto tb40 = th; tb40.michel_far_len_shower_max = 40 * units::cm;
    CHECK(stm_michel_classify_stop_arm(g, muon, arm_b, stop, tb40).kind == StmMichelArm::kOther);
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

// doc pdhd/15: charge -> energy for a Michel piece the fitter never reached.
TEST_CASE("stm_michel charge to energy: the KineChargeOptions arithmetic, guarded")
{
    // The one line every charge-based energy in this repo ends on
    // (NeutrinoEnergyReco.cxx:509): Q / recom / fudge * W / 1e6 * MeV.
    const double W = 23.6;            // eV per ion pair
    const double Q = 3.85e5;          // PDVD 039252_15 cluster 309, the piece the
                                      // doc-14 code added to the shower and never counted
    const double track = stm_michel_charge_to_energy(Q, 0.7, 0.95, W) / units::MeV;
    CHECK(track == doctest::Approx(Q / 0.7 / 0.95 * W / 1e6));
    CHECK(track == doctest::Approx(13.67).epsilon(0.01));
    // that piece's own segment_cal_kine_dQdx on the production arm was 13.00 MeV,
    // so the TRACK factor pair puts an unfitted piece on the fitted scale
    CHECK(std::abs(track - 13.00) / 13.00 < 0.10);

    // the SHOWER pair is the one that overshoots (doc pdhd/15 sec 6: median
    // ratio 1.98 PDVD / 1.67 PDHD against the chain's own dQ/dx)
    const double shower = stm_michel_charge_to_energy(Q, 0.5, 0.8, W) / units::MeV;
    CHECK(shower / track == doctest::Approx(0.7 * 0.95 / (0.5 * 0.8)));
    CHECK(shower > 1.6 * track);

    // linear in the charge
    CHECK(stm_michel_charge_to_energy(2 * Q, 0.7, 0.95, W)
          == doctest::Approx(2 * stm_michel_charge_to_energy(Q, 0.7, 0.95, W)));

    // guards: nothing here may return a NaN or a negative energy, because a
    // non-finite value passes no gate and fails every one silently
    CHECK(stm_michel_charge_to_energy(0.0, 0.7, 0.95, W) == 0.0);
    CHECK(stm_michel_charge_to_energy(-1.0, 0.7, 0.95, W) == 0.0);
    CHECK(stm_michel_charge_to_energy(std::nan(""), 0.7, 0.95, W) == 0.0);
    CHECK(stm_michel_charge_to_energy(Q, 0.0, 0.95, W) == 0.0);
    CHECK(stm_michel_charge_to_energy(Q, 0.7, 0.0, W) == 0.0);
    CHECK(stm_michel_charge_to_energy(Q, 0.7, 0.95, 0.0) == 0.0);
}

namespace {
    // Local mirrors of the two gen models, transcribed from
    // gen/src/PracticalRecombinationModels.cxx:37-42 and
    // gen/src/RecombinationModels.cxx:152-160.  Local because a clus test must
    // not depend on the gen plugin (the doctest_cal_kine_dqdx_zero_dx pattern).
    struct BoxFwd : public IRecombinationModel {
        double m_efield, m_a{0.93}, m_b{0.212}, m_rho{1.38}, m_wi{23.6e-6};
        explicit BoxFwd(double E) : m_efield(E) {}
        virtual ~BoxFwd() {}
        double operator()(double dE, double dX) override
        {
            const double tmp = (dE / units::MeV * units::cm / dX) * m_b / (m_efield * m_rho);
            return std::log(m_a + tmp) / tmp * dE / m_wi;
        }
        double dE(double, double) override { return 0.0; }   // unused here
    };
    struct PowerFwd : public IRecombinationModel {
        double m_a{0.93}, m_k, m_p{1.0}, m_c, m_pivot{2.1}, m_wi{23.6e-6};
        PowerFwd(double k, double C) : m_k(k), m_c(C) {}
        virtual ~PowerFwd() {}
        double operator()(double dE, double dX) override
        {
            const double dedx = dE / units::MeV * units::cm / dX;
            if (dedx <= 0) return 0.0;
            const double u = m_k * std::pow(dedx / m_pivot, m_p);
            return m_c * (std::log(m_a + u) / u) * dedx / m_wi * (dX / units::cm);
        }
        double dE(double, double) override { return 0.0; }
    };
}

// doc pdhd/17 sec 9: the unfitted-charge conversion read out of the bound
// recombination model instead of the hard-coded 0.7 x 0.95 pair.
TEST_CASE("stm_michel_charge_to_energy_model: MeV per electron comes from the bound model")
{
    const double Q = 1e6;                       // electrons
    const double MIP = 2.1;                     // MeV/cm, the assumed dE/dx
    const double flat = stm_michel_charge_to_energy(Q, 0.7, 0.95, 23.6) / units::MeV;

    // ABSOLUTE numbers, not just ratios: the practical-unit models carry a
    // units::cm/units::MeV factor of 10 that a wrong conversion silently eats
    // (doc 88 / e6fb7ef3), and a ratio test would not see it.
    IRecombinationModel::pointer pdvd_box = std::make_shared<BoxFwd>(0.45);
    IRecombinationModel::pointer pdhd_box = std::make_shared<BoxFwd>(0.4959);
    // pr.jsonnet {pdvd,pdhd}_stm_recomb: p = 1, k = beta'*pivot, C measured.
    IRecombinationModel::pointer pdvd_pow = std::make_shared<PowerFwd>(0.7169082125603865, 0.7941);
    IRecombinationModel::pointer pdhd_pow = std::make_shared<PowerFwd>(0.6505519170239442, 0.8120);

    const double e_pdvd_box = stm_michel_charge_to_energy_model(Q, pdvd_box, MIP) / units::MeV;
    const double e_pdhd_box = stm_michel_charge_to_energy_model(Q, pdhd_box, MIP) / units::MeV;
    const double e_pdvd_pow = stm_michel_charge_to_energy_model(Q, pdvd_pow, MIP) / units::MeV;
    const double e_pdhd_pow = stm_michel_charge_to_energy_model(Q, pdhd_pow, MIP) / units::MeV;

    // 1e6 electrons at MIP, in MeV = 1e6 * Wi / (C * R(2.1))
    CHECK(flat == doctest::Approx(35.48872).epsilon(1e-6));
    CHECK(e_pdvd_box == doctest::Approx(33.91269).epsilon(1e-6));
    CHECK(e_pdhd_box == doctest::Approx(33.53843).epsilon(1e-6));
    CHECK(e_pdvd_pow == doctest::Approx(42.70582).epsilon(1e-6));
    CHECK(e_pdhd_pow == doctest::Approx(41.30349).epsilon(1e-6));

    // The point of the knob, stated as the four ratios doc pdhd/17 sec 6 quotes:
    // the flat pair was within 5 % of the UNCALIBRATED inverse and is 16-20 %
    // from the calibrated one.  Both models are seen here, so "it follows
    // whatever is bound" is tested rather than asserted.
    CHECK(e_pdvd_box / flat == doctest::Approx(0.9556).epsilon(1e-3));
    CHECK(e_pdhd_box / flat == doctest::Approx(0.9450).epsilon(1e-3));
    CHECK(e_pdvd_pow / flat == doctest::Approx(1.2034).epsilon(1e-3));
    CHECK(e_pdhd_pow / flat == doctest::Approx(1.1638).epsilon(1e-3));

    // dx cancels: the function fixes dx = 1 cm, and any other choice agrees.
    // (Checked through the model directly, since the function takes no dx.)
    const double dx2 = 7.3 * units::cm;
    const double dE2 = MIP * units::MeV / units::cm * dx2;
    CHECK(Q * dE2 / (*pdhd_pow)(dE2, dx2) / units::MeV == doctest::Approx(e_pdhd_pow).epsilon(1e-9));

    // MIP-EQUIVALENT, and the direction of the bias is not the obvious one:
    // quenching RISES with dE/dx, so a denser deposit needs MORE MeV per
    // electron and assuming MIP UNDER-estimates it.
    const double e_dense = stm_michel_charge_to_energy_model(Q, pdhd_pow, 5.0) / units::MeV;
    CHECK(e_dense > e_pdhd_pow);
    CHECK(e_dense / e_pdhd_pow == doctest::Approx(1.201).epsilon(1e-2));

    // linear in the charge
    CHECK(stm_michel_charge_to_energy_model(2 * Q, pdhd_pow, MIP)
          == doctest::Approx(2 * stm_michel_charge_to_energy_model(Q, pdhd_pow, MIP)));

    // guards: a non-finite or negative energy passes no gate and fails every
    // one silently, so every bad input must give exactly 0.
    IRecombinationModel::pointer null_model;
    CHECK(stm_michel_charge_to_energy_model(0.0, pdhd_pow, MIP) == 0.0);
    CHECK(stm_michel_charge_to_energy_model(-1.0, pdhd_pow, MIP) == 0.0);
    CHECK(stm_michel_charge_to_energy_model(std::nan(""), pdhd_pow, MIP) == 0.0);
    CHECK(stm_michel_charge_to_energy_model(Q, null_model, MIP) == 0.0);
    CHECK(stm_michel_charge_to_energy_model(Q, pdhd_pow, 0.0) == 0.0);
    CHECK(stm_michel_charge_to_energy_model(Q, pdhd_pow, -2.1) == 0.0);
    // Below the Modified Box's A < 1 zero crossing the forward charge goes
    // NEGATIVE; that must give 0, not a negative energy.  The crossing is where
    // A + u = 1, i.e. u = 0.07: at p = 1 that is dE/dx = 0.07*pivot/k, so
    // 0.205 MeV/cm on PDVD and 0.226 on PDHD.  (RecombinationModels.cxx:156's
    // "~0.75 MeV/cm" is the SBND fit's crossing at p = 1.362179, NOT these.)
    CHECK((*pdhd_pow)(0.30 * units::MeV / units::cm * units::cm, 1.0 * units::cm) > 0.0);
    CHECK((*pdhd_pow)(0.15 * units::MeV / units::cm * units::cm, 1.0 * units::cm) < 0.0);
    CHECK(stm_michel_charge_to_energy_model(Q, pdhd_pow, 0.15) == 0.0);
    CHECK(stm_michel_charge_to_energy_model(Q, pdvd_pow, 0.15) == 0.0);
}

// doc pdvd/51: the capture-gamma predicates.  Both are trivial arithmetic and
// that is exactly why they are pinned: the ring's inner edge has to abut the
// Michel's admission test with NO overlap and NO gap, and the whole reason the
// gamma is a separate object class is that a Michel radius wide enough to reach
// it would also feed the Michel piece assembly.
TEST_CASE("stm_michel stop gamma: the ring abuts the Michel radius, and the caps bite")
{
    using WireCell::Clus::PR::stm_michel_stop_gamma_ring;
    const double inner = 15.0, outer = 35.0, maxlen = 10.0;

    // 039252_0 cluster 77's gamma: 6 points, 2.17 cm, 26.5 cm from the stop.
    CHECK(stm_michel_stop_gamma_ring(26.5, 2.17, inner, outer, maxlen));

    // The inner edge is EXCLUSIVE and the Michel's own test is `<=`, so a
    // cluster at exactly the Michel radius belongs to the Michel and to nothing
    // else -- the two partitions must not overlap...
    CHECK_FALSE(stm_michel_stop_gamma_ring(15.0, 2.0, inner, outer, maxlen));
    // ...and must not leave a gap either: anything past it is the gamma's.
    CHECK(stm_michel_stop_gamma_ring(15.0001, 2.0, inner, outer, maxlen));

    // The outer edge is inclusive.
    CHECK(stm_michel_stop_gamma_ring(35.0, 2.0, inner, outer, maxlen));
    CHECK_FALSE(stm_michel_stop_gamma_ring(35.0001, 2.0, inner, outer, maxlen));

    // Compactness.  039252_0 cluster 77's other five same-bundle neighbours sit
    // at 99.3 .. 168.9 cm; the one thing that keeps a NEAR long object out is
    // this cap.
    CHECK(stm_michel_stop_gamma_ring(20.0, 10.0, inner, outer, maxlen));
    CHECK_FALSE(stm_michel_stop_gamma_ring(20.0, 10.0001, inner, outer, maxlen));
    CHECK_FALSE(stm_michel_stop_gamma_ring(99.3, 2.68, inner, outer, maxlen));

    // A NaN admits nothing.  It passes no gate and fails every one silently
    // (doc pdhd/15 sec 7 lost an entire object's energy to one).
    const double nan = std::numeric_limits<double>::quiet_NaN();
    CHECK_FALSE(stm_michel_stop_gamma_ring(nan, 2.0, inner, outer, maxlen));
    CHECK_FALSE(stm_michel_stop_gamma_ring(20.0, nan, inner, outer, maxlen));

    // The feature-off configuration cannot admit anything: with the outer edge
    // at or below the inner one the ring is empty for every distance.
    for (double d : {0.0, 1.0, 15.0, 20.0, 35.0, 400.0})
        CHECK_FALSE(stm_michel_stop_gamma_ring(d, 1.0, inner, /*outer*/ inner, maxlen));
}

TEST_CASE("stm_michel stop gamma: the energy window is a SEPARATE, post-fit stage")
{
    using WireCell::Clus::PR::stm_michel_stop_gamma_energy;
    const double lo = 0.2, hi = 20.0;

    // The measured population: p10 0.37, p50 0.97, p90 4.37 MeV (d16vnu).
    CHECK(stm_michel_stop_gamma_energy(0.80, lo, hi));   // 039252_0 cluster 77
    CHECK(stm_michel_stop_gamma_energy(0.97, lo, hi));

    // Inclusive at both ends.
    CHECK(stm_michel_stop_gamma_energy(lo, lo, hi));
    CHECK(stm_michel_stop_gamma_energy(hi, lo, hi));
    CHECK_FALSE(stm_michel_stop_gamma_energy(0.199, lo, hi));
    CHECK_FALSE(stm_michel_stop_gamma_energy(20.001, lo, hi));

    // Zero and negative are not "small": they are a failed measurement.
    CHECK_FALSE(stm_michel_stop_gamma_energy(0.0, lo, hi));
    CHECK_FALSE(stm_michel_stop_gamma_energy(-1.0, lo, hi));
    CHECK_FALSE(stm_michel_stop_gamma_energy(std::numeric_limits<double>::quiet_NaN(), lo, hi));
}

TEST_CASE("stm_michel survey: the admission radius is a maximum over live stages")
{
    using WireCell::Clus::PR::stm_michel_admit_radius;
    const double michel = 15.0, gamma = 35.0, survey = 60.0;

    // doc pdvd/53.  The shipped ProtoDUNE setting: all three stages on.
    CHECK(stm_michel_admit_radius(michel, gamma, true, survey, true) == doctest::Approx(60.0));

    // Each stage off drops exactly its own contribution, and nothing else.
    CHECK(stm_michel_admit_radius(michel, gamma, true,  survey, false) == doctest::Approx(35.0));
    CHECK(stm_michel_admit_radius(michel, gamma, false, survey, true)  == doctest::Approx(60.0));
    CHECK(stm_michel_admit_radius(michel, gamma, false, survey, false) == doctest::Approx(15.0));

    // THE GATE THIS ROUND RESTS ON: survey off reproduces the doc pdvd/51
    // admission exactly, for every survey radius -- including one larger than
    // the gamma's, which is the only case that could widen anything.
    for (double r : {0.0, 10.0, 35.0, 60.0, 1e4})
        CHECK(stm_michel_admit_radius(michel, gamma, true, r, false) == doctest::Approx(35.0));

    // The survey never NARROWS: a radius inside the gamma ring leaves the gamma
    // stage's reach untouched, so turning the survey on cannot remove a
    // companion the doc pdvd/51 chain admitted.
    for (double r : {0.0, 1.0, 15.0, 34.9})
        CHECK(stm_michel_admit_radius(michel, gamma, true, r, true) == doctest::Approx(35.0));

    // A non-finite radius from a live stage contributes nothing rather than
    // poisoning the maximum -- a NaN max() is implementation-defined and would
    // silently admit everything or nothing (feedback: a NaN fails every gate).
    const double nan = std::numeric_limits<double>::quiet_NaN();
    CHECK(stm_michel_admit_radius(michel, nan, true, survey, true) == doctest::Approx(60.0));
    CHECK(stm_michel_admit_radius(michel, gamma, true, nan, true) == doctest::Approx(35.0));
    CHECK(stm_michel_admit_radius(nan, gamma, true, survey, true) == doctest::Approx(60.0));
}

TEST_CASE("stm_michel survey: the ring and the survey partition the companions")
{
    using WireCell::Clus::PR::stm_michel_stop_gamma_ring;
    using WireCell::Clus::PR::stm_michel_admit_radius;
    const double michel = 15.0, gamma = 35.0, survey = 60.0, maxlen = 10.0;
    const double admit = stm_michel_admit_radius(michel, gamma, true, survey, true);

    // The five same-bundle blobs of 039253_14 cluster 49, at their measured
    // distances.  Every one is now ADMITTED (so it is fitted and gets a role-6
    // row), and exactly the two inside the ring are offered to the gamma stage.
    struct Case { double d; bool admitted; bool offered_to_gamma; };
    const Case cases[] = {
        {14.51, true,  false},   // cluster 185 -- inside the Michel radius
        {33.98, true,  true },   // cluster 184 -- in the ring
        {35.67, true,  false},   // cluster 186 -- past the ring, survey only
        {48.89, true,  false},   // cluster 183 -- survey only
        {54.44, true,  false},   // cluster 187 -- survey only
        {33.09, true,  true },   // 039253_13 cluster 431, the owner's gamma
        {99.30, false, false},   // 039252_0's next neighbour: beyond the survey too
    };
    for (const auto& c : cases) {
        CHECK((c.d <= admit) == c.admitted);
        CHECK(stm_michel_stop_gamma_ring(c.d, 2.0, michel, gamma, maxlen) == c.offered_to_gamma);
    }
}


// doc pdvd/57 -- the STOP RETREAT.  Three chain segments: a plateau body
// (seg1), a plateau segment whose last 10 cm is a Bragg-boosted rise (seg2),
// and a trailing segment (seg3) that stands in for the Michel/other object
// the fit's Steiner path reached beyond the muon's real stop.  seg3's charge
// and reverse_storage vary case by case; seg1/seg2 are shared.
namespace {
StmMichelRetreatThresholds retreat_thresholds(int max_drop = 2, double max_drop_len_cm = 25.0,
                                              double peak_window_cm = 15.0)
{
    StmMichelRetreatThresholds th;
    th.max_drop = max_drop;
    th.collapse_frac = 0.5;
    th.peak_frac = 1.4;
    th.peak_window = peak_window_cm * units::cm;
    th.max_drop_len = max_drop_len_cm * units::cm;
    th.plateau_lo = 20 * units::cm;
    th.plateau_hi = 40 * units::cm;
    th.min_dqdx_live = 0;
    th.min_tail_pts = 2;
    return th;
}
}

TEST_CASE("stm_michel stop retreat: fires on a collapsed tail behind a Bragg rise")
{
    Graph g;
    auto entry = make_vtx(g, 0, 0, 0);
    auto mid   = make_vtx(g, 0, 0, 50);
    auto near  = make_vtx(g, 0, 0, 80);
    auto stop  = make_vtx(g, 0, 0, 88);
    auto s1 = make_track(g, entry, mid, 1.0);
    auto s2 = make_track(g, mid, near, 1.0, false, /*bragg_last_cm=*/10, /*bragg_factor=*/3.0);
    auto s3 = make_track(g, near, stop, 0.15);   // the collapsed tail
    std::vector<SegmentPtr> chain{s1, s2, s3};
    auto prof = stm_michel_profile(g, chain, entry);
    REQUIRE(!prof.empty());
    auto rr = stm_michel_stop_retreat(prof, static_cast<int>(chain.size()), retreat_thresholds());
    CHECK(rr.n_drop == 1);
    CHECK(rr.drop_len == doctest::Approx(8 * units::cm).epsilon(0.05));
    CHECK(rr.plateau > 0);
}

TEST_CASE("stm_michel stop retreat: does not fire on a 1.0 MIP continuation tail")
{
    Graph g;
    auto entry = make_vtx(g, 0, 0, 0);
    auto mid   = make_vtx(g, 0, 0, 50);
    auto near  = make_vtx(g, 0, 0, 80);
    auto stop  = make_vtx(g, 0, 0, 88);
    auto s1 = make_track(g, entry, mid, 1.0);
    auto s2 = make_track(g, mid, near, 1.0, false, 10, 3.0);
    auto s3 = make_track(g, near, stop, 1.0);   // the muon carrying on, not collapsing
    std::vector<SegmentPtr> chain{s1, s2, s3};
    auto prof = stm_michel_profile(g, chain, entry);
    auto rr = stm_michel_stop_retreat(prof, static_cast<int>(chain.size()), retreat_thresholds());
    CHECK(rr.n_drop == 0);
}

TEST_CASE("stm_michel stop retreat: does not fire when nothing peaks before the collapse")
{
    Graph g;
    auto entry = make_vtx(g, 0, 0, 0);
    auto mid   = make_vtx(g, 0, 0, 50);
    auto near  = make_vtx(g, 0, 0, 80);
    auto stop  = make_vtx(g, 0, 0, 88);
    auto s1 = make_track(g, entry, mid, 1.0);
    auto s2 = make_track(g, mid, near, 1.0);   // flat -- no Bragg rise to retreat TO
    auto s3 = make_track(g, near, stop, 0.15);
    std::vector<SegmentPtr> chain{s1, s2, s3};
    auto prof = stm_michel_profile(g, chain, entry);
    auto rr = stm_michel_stop_retreat(prof, static_cast<int>(chain.size()), retreat_thresholds());
    CHECK(rr.n_drop == 0);
}

TEST_CASE("stm_michel stop retreat: max_drop = 0 is off")
{
    Graph g;
    auto entry = make_vtx(g, 0, 0, 0);
    auto mid   = make_vtx(g, 0, 0, 50);
    auto near  = make_vtx(g, 0, 0, 80);
    auto stop  = make_vtx(g, 0, 0, 88);
    auto s1 = make_track(g, entry, mid, 1.0);
    auto s2 = make_track(g, mid, near, 1.0, false, 10, 3.0);
    auto s3 = make_track(g, near, stop, 0.15);
    std::vector<SegmentPtr> chain{s1, s2, s3};
    auto prof = stm_michel_profile(g, chain, entry);
    auto rr = stm_michel_stop_retreat(prof, static_cast<int>(chain.size()), retreat_thresholds(/*max_drop=*/0));
    CHECK(rr.n_drop == 0);
}

TEST_CASE("stm_michel stop retreat: a collapsed tail longer than max_drop_len is refused")
{
    Graph g;
    auto entry = make_vtx(g, 0, 0, 0);
    auto mid   = make_vtx(g, 0, 0, 50);
    auto near  = make_vtx(g, 0, 0, 80);
    auto stop  = make_vtx(g, 0, 0, 88);   // an 8 cm tail
    auto s1 = make_track(g, entry, mid, 1.0);
    auto s2 = make_track(g, mid, near, 1.0, false, 10, 3.0);
    auto s3 = make_track(g, near, stop, 0.15);
    std::vector<SegmentPtr> chain{s1, s2, s3};
    auto prof = stm_michel_profile(g, chain, entry);
    auto rr = stm_michel_stop_retreat(prof, static_cast<int>(chain.size()), retreat_thresholds(2, /*max_drop_len_cm=*/2.0));
    CHECK(rr.n_drop == 0);   // 8 cm tail > 2 cm cap
}

TEST_CASE("stm_michel stop retreat: fires the same way when the last segment is stored reversed")
{
    Graph g;
    auto entry = make_vtx(g, 0, 0, 0);
    auto mid   = make_vtx(g, 0, 0, 50);
    auto near  = make_vtx(g, 0, 0, 80);
    auto stop  = make_vtx(g, 0, 0, 88);
    auto s1 = make_track(g, entry, mid, 1.0);
    auto s2 = make_track(g, mid, near, 1.0, false, 10, 3.0);
    auto s3 = make_track(g, near, stop, 0.15, /*reverse_storage=*/true);
    std::vector<SegmentPtr> chain{s1, s2, s3};
    auto prof = stm_michel_profile(g, chain, entry);
    auto rr = stm_michel_stop_retreat(prof, static_cast<int>(chain.size()), retreat_thresholds());
    CHECK(rr.n_drop == 1);
    CHECK(rr.drop_len == doctest::Approx(8 * units::cm).epsilon(0.05));
}

// doc pdvd/74 (P3): how the retreat reads the dropped tail.  Built from
// arrays -- the overshoot's charge changes row by row, which make_track
// cannot state: a 1.0 plateau (seg 0, L 0-44), a body (seg 1, L 46-80) whose
// last 10 cm rise to 3.0, and a 3 cm overshoot (seg 2) whose first row is the
// vertex again, at the same position -- as stm_michel_profile writes it
// (039252_16/32 carries 92720 on both copies).  The live floor is 0.15, so
// rows below it are what the doc 57 reading calls dead.
namespace {
StmMichelProfile overshoot_profile(double q_vtx, double q_a, double q_mid, double q_e)
{
    std::vector<double> L, q; std::vector<Point> pts; std::vector<int> seg;
    auto add = [&](double l_cm, double c, int s) {
        L.push_back(l_cm * units::cm); q.push_back(c); pts.push_back(Point(0, 0, l_cm) * units::cm); seg.push_back(s);
    };
    for (double l = 0; l <= 44; l += 2) add(l, 1.0, 0);
    for (double l = 46; l <= 70; l += 2) add(l, 1.0, 1);
    add(72, 1.6, 1); add(74, 2.0, 1); add(76, 2.4, 1); add(78, 2.8, 1); add(80, 3.0, 1);
    add(80, q_vtx, 2);                                  // the vertex again: seg 2's first row
    add(80.6, q_a, 2); add(81.2, q_mid, 2); add(81.8, q_mid, 2); add(82.4, q_mid, 2); add(83.0, q_e, 2);
    StmMichelProfile p;
    p.L = L; p.dQdx = q; p.pts = pts; p.seg_idx = seg;
    p.total_length = L.back();
    p.rr.resize(L.size());
    for (size_t i = 0; i < L.size(); ++i) p.rr[i] = p.total_length - L[i];
    return p;
}

StmMichelRetreatThresholds p3_thresholds(bool strict, bool sublive)
{
    auto th = retreat_thresholds();
    th.min_dqdx_live = 0.15;
    th.tail_strict = strict;
    th.tail_sublive = sublive;
    return th;
}

// doc pdvd/57's reading, copied literally, for the one drop these shapes
// allow (seg 2 alone; seg 1 as well would be 37 cm > max_drop_len).
int legacy_one_drop(const StmMichelProfile& p, const StmMichelRetreatThresholds& th)
{
    std::vector<double> pl;
    for (size_t i = 0; i < p.rr.size(); ++i)
        if (p.dQdx[i] >= th.min_dqdx_live && p.rr[i] >= th.plateau_lo && p.rr[i] <= th.plateau_hi) pl.push_back(p.dQdx[i]);
    if (pl.size() < 3) return 0;
    const double plateau = stm_michel_median(pl);
    double bL = -1;
    for (size_t i = 0; i < p.L.size(); ++i) if (p.seg_idx[i] >= 2) { bL = p.L[i]; break; }
    std::vector<double> tail;
    for (size_t i = 0; i < p.L.size(); ++i)
        if (p.L[i] >= bL && p.dQdx[i] >= th.min_dqdx_live) tail.push_back(p.dQdx[i]);
    if (static_cast<int>(tail.size()) < th.min_tail_pts) return 0;
    if (!(stm_michel_median(tail) < th.collapse_frac * plateau)) return 0;
    std::vector<double> kq, kr;
    for (size_t i = 0; i < p.L.size(); ++i)
        if (p.L[i] <= bL && p.dQdx[i] >= th.min_dqdx_live) { kq.push_back(p.dQdx[i]); kr.push_back(bL - p.L[i]); }
    double peak = -1; int nw = 0;
    for (size_t i = 0; i < kq.size(); ++i) {
        const size_t lo = i ? i - 1 : 0, hi = std::min(kq.size() - 1, i + 1);
        std::vector<double> w(kq.begin() + lo, kq.begin() + hi + 1);
        std::sort(w.begin(), w.end());
        if (kr[i] > th.peak_window) continue;
        ++nw;
        peak = std::max(peak, w[w.size() / 2]);
    }
    return (nw >= 3 && peak >= th.peak_frac * plateau) ? 1 : 0;
}
}

TEST_CASE("stm_michel stop retreat P3: both readings off = the doc 57 reading, on a grid of overshoot shapes")
{
    int n = 0, fires = 0;
    for (double q_vtx : {3.0, 0.2})
        for (double q_a : {0.1, 0.4, 0.6, 1.0, 3.0})
            for (double q_mid : {0.0, 0.1, 0.3, 1.0})
                for (double q_e : {0.05, 0.2, 1.0}) {
                    auto p = overshoot_profile(q_vtx, q_a, q_mid, q_e);
                    auto th = p3_thresholds(false, false);
                    const int got = stm_michel_stop_retreat(p, 3, th).n_drop;
                    CHECK(got == legacy_one_drop(p, th));
                    ++n;
                    fires += got;
                }
    CHECK(n == 120);
    CHECK(fires > 0);   // the grid exercises both answers
    CHECK(fires < n);
}

TEST_CASE("stm_michel stop retreat P3: the vertex row sets a short tail's median -- tail_strict reads past it")
{
    // 039252_16/32's shape: the doc 57 tail is {3, 3, 0.6, 0.2} (median 1.8),
    // the strict one {0.6, 0.2} (median 0.4 < 0.5 x plateau).
    auto p = overshoot_profile(3.0, 0.6, 0.1, 0.2);
    CHECK(stm_michel_stop_retreat(p, 3, p3_thresholds(false, false)).n_drop == 0);
    auto rr = stm_michel_stop_retreat(p, 3, p3_thresholds(true, false));
    CHECK(rr.n_drop == 1);
    CHECK(rr.drop_len == doctest::Approx(3.0 * units::cm));
}

TEST_CASE("stm_michel stop retreat P3: a collapse below the live floor is read only with tail_sublive")
{
    // 039253_3/61's shape: nothing live past the vertex.
    auto p = overshoot_profile(3.0, 0.1, 0.1, 0.1);
    CHECK(stm_michel_stop_retreat(p, 3, p3_thresholds(false, false)).n_drop == 0);   // {3, 3}
    CHECK(stm_michel_stop_retreat(p, 3, p3_thresholds(true, false)).n_drop == 0);    // {} -- too few rows to judge
    CHECK(stm_michel_stop_retreat(p, 3, p3_thresholds(true, true)).n_drop == 1);     // {0.1 x 5}
    CHECK(stm_michel_stop_retreat(p, 3, p3_thresholds(false, true)).n_drop == 1);    // {3, 3, 0.1 x 5}: median 0.1
}

TEST_CASE("stm_michel stop retreat P3: tail_sublive takes a dead stretch for a collapse (the known hazard)")
{
    // dQ/dx 0 past a live Bragg rise is what a run of dead channels reads.
    // Nothing at this site knows the channel map, so the sub-live reading
    // retreats off it exactly as off a collapsed Michel -- pinned so the
    // hazard is stated behaviour, not a surprise.
    auto p = overshoot_profile(3.0, 0.0, 0.0, 0.0);
    CHECK(stm_michel_stop_retreat(p, 3, p3_thresholds(false, false)).n_drop == 0);
    CHECK(stm_michel_stop_retreat(p, 3, p3_thresholds(true, true)).n_drop == 1);
}

TEST_CASE("stm_michel stop retreat P3: a live continuation past the vertex is refused under every reading")
{
    auto p = overshoot_profile(3.0, 1.0, 1.0, 1.0);
    for (bool s : {false, true})
        for (bool u : {false, true}) CHECK(stm_michel_stop_retreat(p, 3, p3_thresholds(s, u)).n_drop == 0);
}

// ---------------------------------------------------------------------------
// doc pdvd/58 (T1c): stm_michel_row_kink_deg and stm_michel_stop_split.
// These two build the StmMichelProfile directly rather than through a graph:
// row_kink_deg only ever reads L/pts, and stop_split needs points bent at an
// exact, hand-picked row inside a single segment -- easier to state exactly
// with parallel arrays than to coax out of make_track's straight geometry.
namespace {

StmMichelProfile mkprof(const std::vector<double>& L_cm, const std::vector<double>& dQdx,
                        const std::vector<Point>& pts, const std::vector<int>& seg_idx)
{
    REQUIRE(L_cm.size() == dQdx.size());
    REQUIRE(L_cm.size() == pts.size());
    REQUIRE(L_cm.size() == seg_idx.size());
    StmMichelProfile p;
    p.dQdx = dQdx;
    p.pts = pts;
    p.seg_idx = seg_idx;
    p.L.resize(L_cm.size());
    for (size_t i = 0; i < L_cm.size(); ++i) p.L[i] = L_cm[i] * units::cm;
    p.total_length = p.L.empty() ? 0 : p.L.back();
    p.rr.resize(p.L.size());
    for (size_t i = 0; i < p.L.size(); ++i) p.rr[i] = p.total_length - p.L[i];
    return p;
}

StmMichelSplitThresholds split_thresholds(int max_split = 1, double kink_min_deg = 15.0,
                                          double min_drop_cm = 3.0, double peak_window_cm = 15.0,
                                          double dir_window_cm = 5.0, double max_drop_len_cm = 25.0)
{
    StmMichelSplitThresholds th;
    th.max_split = max_split;
    th.kink_min_deg = kink_min_deg;
    th.min_drop = min_drop_cm * units::cm;
    th.collapse_frac = 0.5;
    th.peak_frac = 1.4;
    th.peak_window = peak_window_cm * units::cm;
    th.dir_window = dir_window_cm * units::cm;
    th.max_drop_len = max_drop_len_cm * units::cm;
    th.plateau_lo = 20 * units::cm;
    th.plateau_hi = 40 * units::cm;
    th.min_dqdx_live = 0;
    th.min_tail_pts = 3;
    return th;
}

// The single-bend profile every basic stop_split test starts from: a flat
// plateau (seg 0, 100/cm, L 0-35), a warmup peak inside the LAST chain
// segment (seg 1, 200/cm, L 40-50, still on the z axis), a bend exactly at
// L=50 (still 200/cm there -- the bend is a position change, not yet the
// collapse), then the collapsed tail (15/cm, L 50-60) that turns 90 degrees
// off the z axis at that same row.  Candidate row = index 10 (L=50): the
// only one where seg_idx==1 (last segment), rr in range, tail beyond it
// collapsed, and a Bragg peak survives within dir_window/peak_window of it.
StmMichelProfile bent_profile(bool bend)
{
    std::vector<double> L; std::vector<double> dQdx; std::vector<Point> pts; std::vector<int> seg;
    auto add = [&](double l_cm, double q, double x, double y, double z, int s) {
        L.push_back(l_cm); dQdx.push_back(q); pts.push_back(Point(x, y, z) * units::cm); seg.push_back(s);
    };
    for (double l = 0; l <= 35; l += 5) add(l, 100, 0, 0, l, 0);        // seg 0: the plateau reference
    add(40, 200, 0, 0, 40, 1);                                        // seg 1: warmup peak, still +z
    add(45, 200, 0, 0, 45, 1);
    add(50, 200, 0, 0, 50, 1);                                        // the bend row
    if (bend) {
        add(55, 15, 0, 5, 50, 1);                                    // leg turns +y here
        add(60, 15, 0, 10, 50, 1);
    }
    else {
        add(55, 15, 0, 0, 55, 1);                                    // no bend: straight on down +z
        add(60, 15, 0, 0, 60, 1);
    }
    return mkprof(L, dQdx, pts, seg);
}

}  // namespace

TEST_CASE("stm_michel_row_kink_deg: 0 deg on a straight run, unmeasurable past either end")
{
    std::vector<double> L; std::vector<Point> pts;
    for (double l = 0; l <= 30; l += 5) { L.push_back(l); pts.push_back(Point(0, 0, l) * units::cm); }
    std::vector<double> zeros_d(L.size(), 0); std::vector<int> zeros_i(L.size(), 0);
    auto prof = mkprof(L, zeros_d, pts, zeros_i);
    CHECK(stm_michel_row_kink_deg(prof, 3, 5 * units::cm) == doctest::Approx(0).epsilon(1e-6));
    CHECK(stm_michel_row_kink_deg(prof, 0, 5 * units::cm) == -1);   // no room behind
    CHECK(stm_michel_row_kink_deg(prof, 6, 5 * units::cm) == -1);   // no room ahead
    CHECK(stm_michel_row_kink_deg(prof, 3, 100 * units::cm) == -1); // window past both ends
}

TEST_CASE("stm_michel_row_kink_deg: 90 deg at a right-angle bend")
{
    std::vector<double> L{0, 5, 10, 15, 20};
    std::vector<Point> pts{Point(0, 0, 0) * units::cm, Point(0, 0, 5) * units::cm, Point(0, 0, 10) * units::cm,
                           Point(0, 5, 10) * units::cm, Point(0, 10, 10) * units::cm};
    std::vector<double> zeros_d(L.size(), 0); std::vector<int> zeros_i(L.size(), 0);
    auto prof = mkprof(L, zeros_d, pts, zeros_i);
    CHECK(stm_michel_row_kink_deg(prof, 2, 5 * units::cm) == doctest::Approx(90).epsilon(1e-6));
}

TEST_CASE("stm_michel stop split: fires on a collapsed tail behind a Bragg rise, WITH a kink")
{
    auto prof = bent_profile(/*bend=*/true);
    auto sp = stm_michel_stop_split(prof, /*n_chain_segs=*/2, split_thresholds());
    REQUIRE(sp.ok);
    CHECK(sp.index == 10);
    CHECK(sp.cut_rr == doctest::Approx(10 * units::cm).epsilon(0.05));
    CHECK(sp.drop_len == doctest::Approx(10 * units::cm).epsilon(0.05));
    CHECK(sp.kink_deg == doctest::Approx(90).epsilon(0.5));
    CHECK(sp.plateau > 0);
}

TEST_CASE("stm_michel stop split: the SAME collapsed tail does NOT fire with no kink -- the discriminator")
{
    auto prof = bent_profile(/*bend=*/false);
    auto sp = stm_michel_stop_split(prof, /*n_chain_segs=*/2, split_thresholds());
    CHECK_FALSE(sp.ok);
}

TEST_CASE("stm_michel stop split: max_split = 0 is off")
{
    auto prof = bent_profile(true);
    auto sp = stm_michel_stop_split(prof, 2, split_thresholds(/*max_split=*/0));
    CHECK_FALSE(sp.ok);
}

TEST_CASE("stm_michel stop split: a qualifying row outside the LAST chain segment is refused")
{
    auto prof = bent_profile(true);   // every point carries seg_idx 0 or 1
    auto sp = stm_michel_stop_split(prof, /*n_chain_segs=*/3, split_thresholds());   // last_seg = 2
    CHECK_FALSE(sp.ok);
}

TEST_CASE("stm_michel stop split: a drop shorter than min_drop is refused")
{
    auto prof = bent_profile(true);   // the only qualifying row has rr = 10 cm
    auto sp = stm_michel_stop_split(prof, 2, split_thresholds(1, 15.0, /*min_drop_cm=*/11.0));
    CHECK_FALSE(sp.ok);
}

TEST_CASE("stm_michel stop split: picks the LARGER of two qualifying kinks")
{
    // Two genuine bends inside the last chain segment, both landing on a
    // collapsed tail with a Bragg rise still in reach: a shallow 45 deg turn
    // at L=44, then a sharper 90 deg turn at L=48.  dir_window/peak_window
    // are loosened from the shipped defaults purely so both bends fit inside
    // one short synthetic segment; the shipped operating point is what doc
    // 58's tests above pin, not this one.
    const double s = 1.0 / std::sqrt(2.0);
    std::vector<double> L; std::vector<double> dQdx; std::vector<Point> pts; std::vector<int> seg;
    auto add = [&](double l_cm, double q, double x, double y, double z, int sidx) {
        L.push_back(l_cm); dQdx.push_back(q); pts.push_back(Point(x, y, z) * units::cm); seg.push_back(sidx);
    };
    for (double l = 0; l <= 30; l += 5) add(l, 100, 0, 0, l, 0);   // seg 0: plateau reference
    add(34, 200, 0, 0, 34, 1); add(36, 200, 0, 0, 36, 1);          // seg 1: warmup peak
    add(38, 200, 0, 0, 38, 1); add(40, 200, 0, 0, 40, 1);
    add(42, 15, 0, 0, 42, 1);
    add(44, 15, 0, 0, 44, 1);                                     // row A: bend to (s,0,s) starts here
    add(46, 15, 2 * s, 0, 44 + 2 * s, 1);                         // +z -> unit_A, 45 deg at row A
    add(48, 15, 4 * s, 0, 44 + 4 * s, 1);                         // row B: bend to (s,0,-s) starts here
    add(50, 15, 4 * s + 2 * s, 0, 44 + 4 * s - 2 * s, 1);         // unit_A -> unit_B, 90 deg at row B
    add(52, 15, 4 * s + 4 * s, 0, 44 + 4 * s - 4 * s, 1);
    auto prof = mkprof(L, dQdx, pts, seg);
    auto th = split_thresholds(/*max_split=*/1, /*kink_min_deg=*/30.0, /*min_drop_cm=*/3.0,
                               /*peak_window_cm=*/25.0, /*dir_window_cm=*/2.0, /*max_drop_len_cm=*/30.0);
    auto sp = stm_michel_stop_split(prof, /*n_chain_segs=*/2, th);
    REQUIRE(sp.ok);
    CHECK(sp.index == 14);   // L=48, the 90 deg row -- NOT index 12 (L=44, 45 deg), even though it is found first
    CHECK(sp.kink_deg == doctest::Approx(90).epsilon(0.5));
}

// ---- doc pdvd/70 (P1): topology-first stop evidence -----------------------
// The rule the verdict applies and the census predicted offline are the same
// expression: clear R_NO_BRAGG | R_SHAPE_FLAT (| R_PROFILE_SPARSE when asked)
// on a Michel that exists, attached or bridged, at or above both minima.

TEST_CASE("stm_michel topology clear: a good Michel clears the two shape bits and nothing else")
{
    const unsigned shape = R_NO_BRAGG | R_SHAPE_FLAT;
    CHECK(stm_michel_topology_clear(shape, 1, 1, 24.0, 6.6, 10.0, 3.0, false) == shape);
    CHECK(stm_michel_topology_clear(R_SHAPE_FLAT, 1, 2, 12.9, 4.8, 10.0, 3.0, false) == R_SHAPE_FLAT);
    CHECK(stm_michel_topology_clear(R_NO_BRAGG, 1, 1, 10.0, 3.0, 10.0, 3.0, false) == R_NO_BRAGG);   // both minima inclusive
    // any other bit survives, so the caller's is_stm still rejects
    for (unsigned b : {unsigned(R_NO_CHAIN), unsigned(R_STOP_UNMATCHED), unsigned(R_NOT_MUON_PID),
                       unsigned(R_CONTINUATION), unsigned(R_STOP_NEAR_BOUNDARY), unsigned(R_VERTEX_HADRON),
                       unsigned(R_SHORT), unsigned(R_PROFILE_SPARSE), unsigned(R_PLATEAU_OFF_MIP),
                       unsigned(R_STOP_INTO_DEAD), unsigned(R_CLUSTER_NOT_TRACK), unsigned(R_PROFILE_GEOMETRY)}) {
        const unsigned clr = stm_michel_topology_clear(b | shape, 1, 1, 30.0, 8.0, 10.0, 3.0, false);
        CHECK(clr == shape);
        CHECK(((b | shape) & ~clr) == b);
    }
    CHECK(stm_michel_topology_clear(0u, 1, 1, 30.0, 8.0, 10.0, 3.0, true) == 0u);   // nothing to clear
}

TEST_CASE("stm_michel topology clear: R_PROFILE_SPARSE only with clears_sparse")
{
    const unsigned bits = R_SHAPE_FLAT | R_PROFILE_SPARSE;
    CHECK(stm_michel_topology_clear(bits, 1, 1, 25.4, 6.5, 10.0, 3.0, false) == unsigned(R_SHAPE_FLAT));
    CHECK(stm_michel_topology_clear(bits, 1, 1, 25.4, 6.5, 10.0, 3.0, true) == bits);
    CHECK(stm_michel_topology_clear(R_PROFILE_SPARSE, 1, 2, 27.1, 22.1, 10.0, 3.0, true) == unsigned(R_PROFILE_SPARSE));
}

TEST_CASE("stm_michel topology clear: no Michel, a charge-only object, a small one or a NaN clears nothing")
{
    const unsigned shape = R_NO_BRAGG | R_SHAPE_FLAT;
    CHECK(stm_michel_topology_clear(shape, 0, 1, 30.0, 8.0, 10.0, 3.0, true) == 0u);   // michel_found 0 (e.g. a T2c/T3c veto)
    CHECK(stm_michel_topology_clear(shape, 1, 0, 30.0, 8.0, 10.0, 3.0, true) == 0u);
    CHECK(stm_michel_topology_clear(shape, 1, 3, 30.0, 8.0, 10.0, 3.0, true) == 0u);   // charge only: no topology
    CHECK(stm_michel_topology_clear(shape, 1, 1, 9.966, 4.2, 10.0, 3.0, true) == 0u);  // 039349_77/52's KE
    CHECK(stm_michel_topology_clear(shape, 1, 1, 21.6, 2.3, 10.0, 3.0, true) == 0u);   // short
    CHECK(stm_michel_topology_clear(shape, 1, 1, std::numeric_limits<double>::quiet_NaN(), 8.0, 10.0, 3.0, true) == 0u);
    CHECK(stm_michel_topology_clear(shape, 1, 1, 30.0, std::numeric_limits<double>::quiet_NaN(), 10.0, 3.0, true) == 0u);
}

// ---- doc pdvd/71 (P4): the Michel's isolated gamma blobs --------------------
// The owner's three criteria: along the Michel electron, a dot near the stop,
// and an energy that over-clustering cannot inflate.  Arguments: d_stop, len,
// cos, d_mich, d_body, ke; then radius 35, max_len 10, cos_min 0.5, max_ke 20.

TEST_CASE("stm_michel gamma gate: each gate, its boundary, and the order they fire in")
{
    CHECK(stm_michel_gamma_gate(26.5, 2.2, 0.90, 20.0, 40.0, 3.1, 35.0, 10.0, 0.5, 20.0) == 0);
    CHECK(stm_michel_gamma_gate(35.0, 10.0, 0.5, 20.0, 20.0001, 20.0, 35.0, 10.0, 0.5, 20.0) == 0);   // every boundary inclusive
    CHECK(stm_michel_gamma_gate(35.0001, 2.0, 0.9, 20.0, 40.0, 3.0, 35.0, 10.0, 0.5, 20.0) == 1);     // radius
    CHECK(stm_michel_gamma_gate(20.0, 10.0001, 0.9, 20.0, 40.0, 3.0, 35.0, 10.0, 0.5, 20.0) == 2);    // length
    CHECK(stm_michel_gamma_gate(20.0, 2.0, 0.4999, 20.0, 40.0, 3.0, 35.0, 10.0, 0.5, 20.0) == 3);     // cone
    CHECK(stm_michel_gamma_gate(20.0, 2.0, -0.9, 20.0, 40.0, 3.0, 35.0, 10.0, 0.5, 20.0) == 3);       // behind the stop
    CHECK(stm_michel_gamma_gate(20.0, 2.0, 0.9, 12.0, 12.0, 3.0, 35.0, 10.0, 0.5, 20.0) == 4);        // a tie goes to the body
    CHECK(stm_michel_gamma_gate(20.0, 2.0, 0.9, 12.0, 11.0, 3.0, 35.0, 10.0, 0.5, 20.0) == 4);
    CHECK(stm_michel_gamma_gate(20.0, 2.0, 0.9, 12.0, 40.0, 20.0001, 35.0, 10.0, 0.5, 20.0) == 5);    // an over-clustered lump
    // the first failing gate is the one reported
    CHECK(stm_michel_gamma_gate(50.0, 30.0, -1.0, 12.0, 1.0, 99.0, 35.0, 10.0, 0.5, 20.0) == 1);
    CHECK(stm_michel_gamma_gate(20.0, 30.0, -1.0, 12.0, 1.0, 99.0, 35.0, 10.0, 0.5, 20.0) == 2);
    // the knobs move the boundaries: a 60 cm radius admits the 50 cm blob
    CHECK(stm_michel_gamma_gate(50.0, 2.0, 0.9, 30.0, 60.0, 3.0, 60.0, 10.0, 0.5, 20.0) == 0);
}

TEST_CASE("stm_michel gamma gate: a non-finite input is rejected before any other gate")
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    CHECK(stm_michel_gamma_gate(nan, 2.0, 0.9, 20.0, 40.0, 3.0, 35.0, 10.0, 0.5, 20.0) == 6);
    CHECK(stm_michel_gamma_gate(20.0, 2.0, nan, 20.0, 40.0, 3.0, 35.0, 10.0, 0.5, 20.0) == 6);
    CHECK(stm_michel_gamma_gate(20.0, 2.0, 0.9, 20.0, 40.0, nan, 35.0, 10.0, 0.5, 20.0) == 6);
    CHECK(stm_michel_gamma_gate(20.0, 2.0, 0.9, 20.0, inf, 3.0, 35.0, 10.0, 0.5, 20.0) == 6);
}

TEST_CASE("stm_michel gamma take: the running total never passes the cap; a blob that does not fit is skipped")
{
    using V = std::vector<int>;
    CHECK(stm_michel_gamma_take(30.0, {2.0, 5.0, 1.0}, 60.0) == V{1, 1, 1});
    CHECK(stm_michel_gamma_take(50.0, {5.0, 12.0, 3.0}, 60.0) == V{1, 0, 1});   // 55, 67 skipped, 58
    CHECK(stm_michel_gamma_take(50.0, {10.0}, 60.0) == V{1});                   // exactly at the cap is taken
    CHECK(stm_michel_gamma_take(76.0, {0.5, 1.0}, 60.0) == V{0, 0});            // a core already over the cap gets nothing
    CHECK(stm_michel_gamma_take(30.0, {}, 60.0).empty());
    const double nan = std::numeric_limits<double>::quiet_NaN();
    CHECK(stm_michel_gamma_take(nan, {1.0}, 60.0) == V{0});
    CHECK(stm_michel_gamma_take(30.0, {nan, -1.0, 2.0}, 60.0) == V{0, 0, 1});
}

// ---------------------------------------------------------------------------
// doc pdvd/82 (doc 78 action item 2): the peak-relative collapsed-tail
// admission.  Every profile here has a HOT tail -- 0.9 x the plateau, which is
// nothing like a collapse and is exactly what the 14 items of doc 78 sec 2.3
// carry: the fit rides through the Michel, so the charge after the Bragg peak
// falls to 0.5-1.9 x plateau, never under it.  0.9 x plateau is 0.45 x the peak.
// ---------------------------------------------------------------------------
namespace {

// The split's bent_profile above with the tail raised from 15 (0.15 x plateau,
// a real collapse) to 90 (0.9 x plateau, 0.45 x peak).  Fork by duplication:
// bent_profile is what the doc 58 cases pin and must not move.
StmMichelProfile hot_tail_split_profile(bool bend)
{
    std::vector<double> L; std::vector<double> dQdx; std::vector<Point> pts; std::vector<int> seg;
    auto add = [&](double l_cm, double q, double x, double y, double z, int s) {
        L.push_back(l_cm); dQdx.push_back(q); pts.push_back(Point(x, y, z) * units::cm); seg.push_back(s);
    };
    for (double l = 0; l <= 35; l += 5) add(l, 100, 0, 0, l, 0);        // seg 0: the plateau reference
    add(40, 200, 0, 0, 40, 1);                                        // seg 1: the Bragg peak, still +z
    add(45, 200, 0, 0, 45, 1);
    add(50, 200, 0, 0, 50, 1);                                        // the candidate row
    if (bend) {
        add(55, 90, 0, 5, 50, 1);                                     // the Michel, turning +y
        add(60, 90, 0, 10, 50, 1);
    }
    else {
        add(55, 90, 0, 0, 55, 1);                                     // straight on: a muon still going
        add(60, 90, 0, 0, 60, 1);
    }
    return mkprof(L, dQdx, pts, seg);
}

// The same shape as three CHAIN SEGMENTS, so stm_michel_stop_retreat can drop
// the last one: the bend sits exactly on the segment boundary (L = 55), which
// is the row the retreat measures its kink at.
StmMichelProfile hot_tail_retreat_profile(bool bend)
{
    std::vector<double> L; std::vector<double> dQdx; std::vector<Point> pts; std::vector<int> seg;
    auto add = [&](double l_cm, double q, double x, double y, double z, int s) {
        L.push_back(l_cm); dQdx.push_back(q); pts.push_back(Point(x, y, z) * units::cm); seg.push_back(s);
    };
    for (double l = 0; l <= 35; l += 5) add(l, 100, 0, 0, l, 0);       // seg 0: plateau
    add(40, 200, 0, 0, 40, 1); add(45, 200, 0, 0, 45, 1);             // seg 1: the Bragg peak
    add(50, 200, 0, 0, 50, 1);
    if (bend) {                                                       // seg 2: the hot tail
        add(55, 90, 0, 0, 55, 2);                                     // the boundary row, still +z
        add(60, 90, 0, 5, 55, 2);                                     // turns +y past it -> 90 deg AT L=55
        add(65, 90, 0, 10, 55, 2);
    }
    else {
        add(55, 90, 0, 0, 55, 2);
        add(60, 90, 0, 0, 60, 2);
        add(65, 90, 0, 0, 65, 2);
    }
    return mkprof(L, dQdx, pts, seg);
}

StmMichelRetreatThresholds doc82_retreat_th(double tail_peak_frac, double kink_min = 25.0)
{
    auto th = retreat_thresholds();
    th.tail_peak_frac = tail_peak_frac;
    th.tail_peak_kink_min = kink_min;
    th.dir_window = 5 * units::cm;
    return th;
}

StmMichelSplitThresholds doc82_split_th(double tail_peak_frac, double kink_min = 25.0)
{
    auto th = split_thresholds();
    th.tail_peak_frac = tail_peak_frac;
    th.tail_peak_kink_min = kink_min;
    return th;
}

}  // namespace

TEST_CASE("stm_michel doc82 retreat: a hot tail is refused with the knob off and taken with it on")
{
    auto prof = hot_tail_retreat_profile(/*bend=*/true);
    // off: the doc 57 reading alone -- 0.9 x plateau is not a collapse
    auto off = stm_michel_stop_retreat(prof, 3, doc82_retreat_th(0.0));
    CHECK(off.n_drop == 0);
    CHECK_FALSE(off.by_tail_peak);
    // on: 0.9 x plateau IS 0.45 x the surviving peak, and the boundary row turns 90 deg
    auto on = stm_michel_stop_retreat(prof, 3, doc82_retreat_th(0.5));
    REQUIRE(on.n_drop == 1);
    CHECK(on.by_tail_peak);
    CHECK(on.drop_len == doctest::Approx(10 * units::cm).epsilon(0.05));
    CHECK(on.plateau == doctest::Approx(100).epsilon(0.01));
    CHECK(on.last_peak == doctest::Approx(200).epsilon(0.01));
    CHECK(on.last_tail_med == doctest::Approx(90).epsilon(0.01));
    CHECK(on.last_kink_deg == doctest::Approx(90).epsilon(0.5));
}

TEST_CASE("stm_michel doc82 retreat: the bend is what carries the looser tail -- straight is refused")
{
    auto prof = hot_tail_retreat_profile(/*bend=*/false);
    auto on = stm_michel_stop_retreat(prof, 3, doc82_retreat_th(0.5));
    CHECK(on.n_drop == 0);
    // and the same profile is refused by a kink bar above the bend it does have
    auto bent = hot_tail_retreat_profile(true);
    CHECK(stm_michel_stop_retreat(bent, 3, doc82_retreat_th(0.5, /*kink_min=*/95.0)).n_drop == 0);
}

TEST_CASE("stm_michel doc82 retreat: a tail above tail_peak_frac x peak is still refused")
{
    auto prof = hot_tail_retreat_profile(true);   // tail / peak = 0.45
    CHECK(stm_michel_stop_retreat(prof, 3, doc82_retreat_th(0.4)).n_drop == 0);
    CHECK(stm_michel_stop_retreat(prof, 3, doc82_retreat_th(0.5)).n_drop == 1);
}

TEST_CASE("stm_michel doc82 retreat: a genuine collapse still fires, and NOT by the peak rule")
{
    Graph g;
    auto entry = make_vtx(g, 0, 0, 0);
    auto mid   = make_vtx(g, 0, 0, 50);
    auto near  = make_vtx(g, 0, 0, 80);
    auto stop  = make_vtx(g, 0, 0, 88);
    auto s1 = make_track(g, entry, mid, 1.0);
    auto s2 = make_track(g, mid, near, 1.0, false, 10, 3.0);
    auto s3 = make_track(g, near, stop, 0.15);   // the doc 57 collapse
    std::vector<SegmentPtr> chain{s1, s2, s3};
    auto prof = stm_michel_profile(g, chain, entry);
    auto legacy = stm_michel_stop_retreat(prof, 3, doc82_retreat_th(0.0));
    auto with82 = stm_michel_stop_retreat(prof, 3, doc82_retreat_th(0.5));
    CHECK(with82.n_drop == legacy.n_drop);            // the knob only ADDS admissions
    CHECK(with82.drop_len == legacy.drop_len);
    CHECK_FALSE(with82.by_tail_peak);                 // this one is the plateau rule's, straight chain and all
}

TEST_CASE("stm_michel doc82 split: a hot tail is refused with the knob off and taken with it on")
{
    auto prof = hot_tail_split_profile(/*bend=*/true);
    auto off = stm_michel_stop_split(prof, 2, doc82_split_th(0.0));
    CHECK_FALSE(off.ok);
    auto on = stm_michel_stop_split(prof, 2, doc82_split_th(0.5));
    REQUIRE(on.ok);
    CHECK(on.by_tail_peak);
    CHECK(on.index == 10);
    CHECK(on.cut_rr == doctest::Approx(10 * units::cm).epsilon(0.05));
    CHECK(on.kink_deg == doctest::Approx(90).epsilon(0.5));
    CHECK(on.tail_med == doctest::Approx(90).epsilon(0.01));
    CHECK(on.peak == doctest::Approx(200).epsilon(0.01));
}

TEST_CASE("stm_michel doc82 split: no bend, or a bend under tail_peak_kink_min, is refused")
{
    CHECK_FALSE(stm_michel_stop_split(hot_tail_split_profile(false), 2, doc82_split_th(0.5)).ok);
    CHECK_FALSE(stm_michel_stop_split(hot_tail_split_profile(true), 2, doc82_split_th(0.5, 95.0)).ok);
}

TEST_CASE("stm_michel doc82 split: the doc 58 collapse cases are untouched by the knob")
{
    // bent_profile's 0.15 x plateau tail: the plateau rule admits it, so the
    // bend bar the peak rule would impose never applies -- tail_peak_kink_min
    // is set absurdly high here to prove exactly that.
    auto legacy = stm_michel_stop_split(bent_profile(true), 2, split_thresholds());
    auto with82 = stm_michel_stop_split(bent_profile(true), 2, doc82_split_th(0.5, 179.0));
    REQUIRE(legacy.ok);
    REQUIRE(with82.ok);
    CHECK(with82.index == legacy.index);
    CHECK(with82.cut_rr == legacy.cut_rr);
    CHECK_FALSE(with82.by_tail_peak);
    CHECK_FALSE(stm_michel_stop_split(bent_profile(false), 2, doc82_split_th(0.5, 179.0)).ok);
}

// doc pdvd/83 (doc 78 action item 4): a Michel leaving the chain BEFORE the
// stop.  The shape of all four owner-scan items: entry -> pen (the muon) ->
// stop (a short stub), with the Michel hanging off `pen`, turned back.
namespace {
struct NearChain {
    std::vector<SegmentPtr> chain;
    std::vector<VertexPtr> vtxs;
    VertexPtr pen;
};
// muon along +z to z = z_pen, then a stub of `stub_cm` to the stop
NearChain near_chain(Graph& g, double z_pen, double stub_cm)
{
    NearChain c;
    auto entry = make_vtx(g, 0, 0, 0);
    c.pen = make_vtx(g, 0, 0, z_pen);
    auto stop = make_vtx(g, 0, 0, z_pen + stub_cm);
    c.chain = {make_track(g, entry, c.pen, 1.0), make_track(g, c.pen, stop, 1.0, false, stub_cm, 2.5)};
    c.vtxs = {entry, c.pen, stop};
    return c;
}
// an arm `len_cm` long leaving v at `kink_deg` from +z, in the y-z plane
SegmentPtr arm_at(Graph& g, VertexPtr v, double len_cm, double kink_deg, double ratio)
{
    const auto p = v->wcpt().point;
    const double a = kink_deg * M_PI / 180;
    auto far = make_vtx(g, p.x() / units::cm, p.y() / units::cm + len_cm * std::sin(a),
                        p.z() / units::cm + len_cm * std::cos(a));
    return make_track(g, v, far, ratio);
}
const auto no_skip = [](const SegmentPtr&) { return false; };
}  // namespace

TEST_CASE("stm_michel near-stop arm (doc pdvd/83): a turned-back arm off the penultimate vertex is the Michel")
{
    Graph g;
    auto c = near_chain(g, 60, 2.5);
    auto arm = arm_at(g, c.pen, 14.4, 146, 1.1);    // 039349_64/24 s24003's numbers
    auto th = thresholds();
    auto r = stm_michel_near_stop_arms(g, c.chain, c.vtxs, 5 * units::cm, th, no_skip);
    REQUIRE(r.vtx_index == 1);
    CHECK(r.dist == doctest::Approx(2.5 * units::cm).epsilon(0.02));
    CHECK(r.n_examined == 1);
    REQUIRE(r.michel.size() == 1);
    CHECK(r.michel[0].seg == arm);
    CHECK(r.michel[0].kind == StmMichelArm::kMichel);
    CHECK(r.michel[0].kink_deg == doctest::Approx(146).epsilon(0.05));   // against the INCOMING segment, not the stub
    // off, or a one-segment chain: nothing
    CHECK(stm_michel_near_stop_arms(g, c.chain, c.vtxs, 0, th, no_skip).vtx_index == -1);
    CHECK(stm_michel_near_stop_arms(g, {c.chain[0]}, {c.vtxs[0], c.vtxs[1]}, 5 * units::cm, th, no_skip).vtx_index == -1);
    // the skip predicate (the caller's chain set + claimed ids) is honoured
    auto rs = stm_michel_near_stop_arms(g, c.chain, c.vtxs, 5 * units::cm, th,
                                        [&](const SegmentPtr& s) { return s == arm; });
    CHECK(rs.vtx_index == -1);
    CHECK(rs.n_examined == 0);
}

TEST_CASE("stm_michel near-stop arm (doc pdvd/83): the distance gate")
{
    Graph g;
    auto c = near_chain(g, 51, 9.0);
    arm_at(g, c.pen, 14.4, 146, 1.1);
    auto th = thresholds();
    auto r5 = stm_michel_near_stop_arms(g, c.chain, c.vtxs, 5 * units::cm, th, no_skip);
    CHECK(r5.vtx_index == -1);
    CHECK(r5.n_examined == 0);                       // the vertex is never reached
    auto r10 = stm_michel_near_stop_arms(g, c.chain, c.vtxs, 10 * units::cm, th, no_skip);
    CHECK(r10.vtx_index == 1);
    CHECK(r10.dist == doctest::Approx(9.0 * units::cm).epsilon(0.02));
}

TEST_CASE("stm_michel near-stop arm (doc pdvd/83): deltas, hadrons and a continuation are not taken")
{
    Graph g;
    auto c = near_chain(g, 60, 2.5);
    arm_at(g, c.pen, 5, 120, 1.2);                   // <= 8 cm terminal: the interior loop's delta
    auto h = arm_at(g, c.pen, 15, 100, 1.6);         // > 8 cm at 1.6 MIP: the interior loop's hadron
    make_track(g, find_other_vertex(g, h, c.pen), make_vtx(g, 0, 20, 50), 1.5);
    auto th = thresholds();
    auto r = stm_michel_near_stop_arms(g, c.chain, c.vtxs, 5 * units::cm, th, no_skip);
    CHECK(r.vtx_index == -1);
    CHECK(r.n_examined == 0);
    // a collinear MIP arm is offered the gate and read as a continuation -- never taken
    Graph g2;
    auto c2 = near_chain(g2, 60, 2.5);
    arm_at(g2, c2.pen, 10, 5, 1.0);
    auto r2 = stm_michel_near_stop_arms(g2, c2.chain, c2.vtxs, 5 * units::cm, th, no_skip);
    CHECK(r2.vtx_index == -1);
    CHECK(r2.n_examined == 1);
    // below the Michel charge floor (0.3 MIP): offered, refused -- 039349_69/56 s56006 reads 0.29
    Graph g3;
    auto c3 = near_chain(g3, 60, 6.5);
    arm_at(g3, c3.pen, 13.6, 90, 0.25);
    auto r3 = stm_michel_near_stop_arms(g3, c3.chain, c3.vtxs, 7 * units::cm, th, no_skip);
    CHECK(r3.vtx_index == -1);
    CHECK(r3.n_examined == 1);
}

TEST_CASE("stm_michel near-stop arm (doc pdvd/83): the nearest vertex wins, longest arm first")
{
    Graph g;
    auto entry = make_vtx(g, 0, 0, 0);
    auto v1 = make_vtx(g, 0, 0, 40);
    auto v2 = make_vtx(g, 0, 0, 56);
    auto stop = make_vtx(g, 0, 0, 58.5);
    std::vector<SegmentPtr> chain{make_track(g, entry, v1, 1.0), make_track(g, v1, v2, 1.0), make_track(g, v2, stop, 1.0)};
    std::vector<VertexPtr> vtxs{entry, v1, v2, stop};
    arm_at(g, v1, 12, 140, 1.0);                     // 18.5 cm from the stop
    auto a_short = arm_at(g, v2, 9, 130, 0.9);       // 2.5 cm from the stop
    auto a_long = arm_at(g, v2, 13, 150, 1.0);
    auto th = thresholds();
    auto r = stm_michel_near_stop_arms(g, chain, vtxs, 20 * units::cm, th, no_skip);
    REQUIRE(r.vtx_index == 2);
    REQUIRE(r.michel.size() == 2);
    CHECK(r.michel[0].seg == a_long);
    CHECK(r.michel[1].seg == a_short);
    CHECK(r.n_examined == 2);                        // v1 is never reached: v2 already won
}

// ---- doc pdvd/84 (doc 78 action item 3): the moved-stop veto's exemptions ----

TEST_CASE("stm_michel moved-stop spare: both off vetoes, the kink wins first, then the reach")
{
    using S = StmMichelMovedStopSpare;
    CHECK(stm_michel_moved_stop_spare(132.6, 11.9, -1, -1) == S::kVeto);   // both tests off
    CHECK(stm_michel_moved_stop_spare(132.6, 11.9, 60, -1) == S::kKink);   // doc 72 alone
    CHECK(stm_michel_moved_stop_spare(132.6, 11.9, 60, 6.5) == S::kKink);  // both pass: the kink is counted, as before
    CHECK(stm_michel_moved_stop_spare(59.58, 9.3, 60, 6.5) == S::kReach);  // under the kink, over the reach
    CHECK(stm_michel_moved_stop_spare(59.58, 9.3, 60, -1) == S::kVeto);    // reach off: doc 72's answer
    CHECK(stm_michel_moved_stop_spare(59.58, 9.3, -1, 6.5) == S::kReach);  // kink off, reach alone
    CHECK(stm_michel_moved_stop_spare(58.68, 5.4, 60, 6.5) == S::kVeto);   // neither
}

TEST_CASE("stm_michel moved-stop spare: boundaries inclusive, unmeasured and non-finite never spare")
{
    using S = StmMichelMovedStopSpare;
    CHECK(stm_michel_moved_stop_spare(60.0, 0.0, 60, 6.5) == S::kKink);    // kink == min
    CHECK(stm_michel_moved_stop_spare(10.0, 6.5, 60, 6.5) == S::kReach);   // reach == min
    CHECK(stm_michel_moved_stop_spare(-1.0, 0.0, 0, -1) == S::kVeto);      // kink -1 (unmeasurable), even at min 0
    CHECK(stm_michel_moved_stop_spare(-1.0, 7.0, 0, 6.5) == S::kReach);    // ... but the reach still reads
    const double nan = std::numeric_limits<double>::quiet_NaN();
    CHECK(stm_michel_moved_stop_spare(nan, nan, 60, 6.5) == S::kVeto);
    CHECK(stm_michel_moved_stop_spare(nan, 9.3, 60, 6.5) == S::kReach);
    CHECK(stm_michel_moved_stop_spare(132.6, nan, 60, 6.5) == S::kKink);
}

TEST_CASE("stm_michel moved-stop spare: the record's 12 veto instances at 60 deg / 6.5 cm spare only the owner's Michels")
{
    // d84_t2c_census.py sec 2: every distinct T2c instance on 28 PDVD arms,
    // (kink deg, len + far cm, owner-confirmed Michel?).
    struct I { double kink, reach; bool michel; };
    const I rec[] = {
        {59.58, 9.3, true},  {59.58, 9.3, true},  {89.34, 7.5, true},  {89.34, 7.5, true},    // 039252_2/79
        {132.64, 11.9, true},                                                                 // 039349_48/21
        {17.19, 5.1, false}, {43.47, 4.9, false},                                             // 039252_4/55, 039349_20/41
        {48.21, 5.8, false}, {48.21, 5.8, false}, {58.68, 5.4, false},                        // 039349_61/62
        {33.93, 5.6, false}, {23.89, 5.4, false},                                             // 039349_72/54, 039349_75/73
    };
    for (const auto& i : rec) {
        const bool spared = stm_michel_moved_stop_spare(i.kink, i.reach, 60, 6.5) != StmMichelMovedStopSpare::kVeto;
        CHECK(spared == i.michel);
    }
    // The kink alone (doc 72's 60) leaves 039252_2/79 at 59.58 vetoed.
    CHECK(stm_michel_moved_stop_spare(59.58, 9.3, 60, -1) == StmMichelMovedStopSpare::kVeto);
}
