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
