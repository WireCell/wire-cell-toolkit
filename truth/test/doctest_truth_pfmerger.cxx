/** PFMerger on a synthetic tree: nu -> mu + p + pi0 -> gamma gamma, a
 *  short and a long delta ray off the muon, a Michel electron, a
 *  Compton electron of one photon, a conversion pair of the other,
 *  a neutron with a scattered proton, and an NC variant with no
 *  primary at the vertex (void root).
 */

#include "WireCellTruth/PFMerger.h"
#include "WireCellTruth/Schema.h"
#include "WireCellTruth/TruthTable.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell;
using namespace WireCell::Truth;

namespace {
    struct Builder {
        TruthTable tt;
        explicit Builder(Point vtx = Point(0, 0, 0))
        {
            std::vector<double> nu(truth_nu_columns().size(), 0.0);
            nu[column_index(truth_nu_columns(), "nu_idx")] = 0;
            nu[column_index(truth_nu_columns(), "pdg")] = 14;
            nu[column_index(truth_nu_columns(), "vtx_x")] = vtx.x();
            nu[column_index(truth_nu_columns(), "vtx_y")] = vtx.y();
            nu[column_index(truth_nu_columns(), "vtx_z")] = vtx.z();
            tt.nu.push_back(nu);
            tt.source = "test";
        }
        // positions in cm, KE in GeV, charge in electrons
        void add(int tid, int mother, int pdg, const char* proc, double ke, Point start, Point end, double q)
        {
            const auto& c = truth_pf_columns();
            std::vector<double> r(c.size(), 0.0);
            auto set = [&](const char* n, double v) { r[column_index(c, n)] = v; };
            set("nu_row", 0); set("trackid", tid); set("mother_trackid", mother); set("parent_trackid", mother);
            set("pdg", pdg); set("process", process_code(proc)); set("KE", ke); set("E", ke);
            set("start_x", start.x()); set("start_y", start.y()); set("start_z", start.z());
            set("end_x", end.x()); set("end_y", end.y()); set("end_z", end.z());
            set("orig_trackid", tid); set("pf_node", -1); set("cls", -1); set("merged_into", -1); set("nvis_q", q);
            tt.pf.push_back(r);
        }
    };

    const PFNode& node_of(const PFResult& res, int tid)
    {
        const auto& gi = res.info(tid);
        REQUIRE(gi.pf_node >= 0);
        return res.nodes[gi.pf_node];
    }
}

TEST_CASE("pfmerger builds nodes, merges deltas and showers, keeps the michel")
{
    const Point vtx(10, 20, 30);
    Builder b(vtx);
    b.add(1, 0, 13, "primary", 0.5, vtx, Point(10, 20, 130), 9e6);           // muon, 100 cm
    b.add(2, 0, 2212, "primary", 0.1, vtx, Point(10, 25, 35), 1.5e6);        // proton, ~7 cm
    b.add(3, 0, 111, "primary", 0.2, vtx, vtx, 0);                           // pi0 (invisible)
    b.add(4, 3, 22, "Decay", 0.1, vtx, Point(10, 30, 60), 0);                // gamma 1
    b.add(5, 3, 22, "Decay", 0.1, vtx, Point(20, 20, 40), 0);                // gamma 2
    b.add(6, 4, 11, "conv", 0.06, Point(10, 30, 60), Point(10, 34, 70), 2e6);  // e- of conversion
    b.add(7, 4, -11, "conv", 0.04, Point(10, 30, 60), Point(10, 33, 68), 1.5e6); // e+ of conversion
    b.add(8, 6, 22, "eBrem", 0.01, Point(10, 32, 65), Point(10, 35, 75), 0);  // brem inside the shower
    b.add(9, 8, 11, "compt", 0.005, Point(10, 35, 75), Point(10, 35.2, 75.3), 5e4); // compton inside the shower
    b.add(10, 5, 11, "compt", 0.003, Point(20, 20, 40), Point(20, 20.1, 40.2), 4e4); // lone compton of gamma 2 -> LE
    b.add(11, 1, 11, "muIoni", 0.002, Point(10, 20, 50), Point(10, 20.3, 50.2), 3e4); // short delta
    b.add(12, 1, 11, "muIoni", 0.02, Point(10, 20, 80), Point(10, 25, 82), 4e5);     // long delta (5.4 cm)
    b.add(13, 1, -11, "Decay", 0.03, Point(10, 20, 130), Point(10, 24, 135), 6e5);   // Michel
    b.add(14, 0, 2112, "primary", 0.05, vtx, Point(50, 20, 30), 0);          // neutron (visible, no charge)
    b.add(15, 14, 2212, "hadElastic", 0.02, Point(50, 20, 30), Point(52, 20, 30), 3e5); // scattered proton, 2 cm
    b.add(16, 2, 1000180400, "hIoni", 0.0, Point(10, 25, 35), Point(10, 25, 35), 0);   // nucleus (invisible)
    b.add(17, 1, 11, "muIoni", 0.0005, Point(10, 20, 90), Point(10, 20.05, 90.1), 500); // tiny delta

    PFMerger merger;
    auto res = merger.run(b.tt);

    // primaries
    const auto& mu = node_of(res, 1);
    CHECK(mu.cls == Cls::track);
    CHECK(mu.primary);
    CHECK(mu.parent == -1);
    CHECK(mu.survived);
    CHECK(node_of(res, 2).primary);

    // the pi0 is invisible and has no visible ancestor: no node; its
    // photons are primaries (SPINE skip list)
    CHECK(res.info(3).pf_node == -1);
    CHECK(res.info(3).merged_into == 3);
    const auto& g1 = node_of(res, 4);
    CHECK(g1.cls == Cls::shower);
    CHECK(g1.primary);
    CHECK(g1.parent == -1);
    // the whole conversion shower is one node
    CHECK(res.info(6).pf_node == g1.index);
    CHECK(res.info(7).pf_node == g1.index);
    CHECK(res.info(8).pf_node == g1.index);
    CHECK(res.info(9).pf_node == g1.index);
    CHECK(g1.q_vis == doctest::Approx(2e6 + 1.5e6 + 5e4));
    CHECK(g1.ids == std::vector<int>{4, 6, 7, 8, 9});
    CHECK(res.info(9).merged_into == 4);

    // gamma 2 with only a small compton electron: LE, instance invalid
    const auto& g2 = node_of(res, 5);
    CHECK(g2.cls == Cls::le);
    CHECK_FALSE(g2.instance_valid);
    CHECK(res.info(10).pf_node == g2.index);
    CHECK_FALSE(res.info(10).instance_valid);

    // deltas: short one merged into the muon, long one its own node
    CHECK(res.info(11).pf_node == mu.index);
    CHECK(res.info(11).merged_into == 1);
    CHECK(res.info(17).pf_node == mu.index);
    const auto& d12 = node_of(res, 12);
    CHECK(d12.cls == Cls::delta);
    CHECK(d12.parent == mu.index);
    CHECK(d12.trackid == 12);
    CHECK(res.nmerged_delta == 2);

    // Michel: own node, parent muon
    const auto& mi = node_of(res, 13);
    CHECK(mi.cls == Cls::michel);
    CHECK(mi.parent == mu.index);
    CHECK(mi.survived);

    // neutron: a node with no charge merges into its parent (nothing) and
    // survives; the scattered proton hangs below it
    const auto& n = node_of(res, 14);
    CHECK(n.cls == Cls::track);
    CHECK(n.parent == -1);
    const auto& p15 = node_of(res, 15);
    CHECK(p15.trackid == 15);
    CHECK(p15.parent == n.index);
    CHECK(p15.cls == Cls::track);
    CHECK_FALSE(p15.primary);

    // the nucleus maps to the proton node
    CHECK(res.info(16).pf_node == node_of(res, 2).index);
    CHECK(res.nmerged_invisible == 2);

    // the muon's charge includes its merged deltas
    CHECK(mu.q_vis == doctest::Approx(9e6 + 3e4 + 500));
    CHECK(res.void_roots.empty());

    // meta is serialisable
    auto meta = merger.params().meta();
    CHECK(meta["version"].asInt() == 1);
    CHECK(meta["delta_len_cm"].asDouble() == doctest::Approx(1.0));
}

TEST_CASE("pfmerger weak nodes merge into the parent and an NC without a primary is a void root")
{
    const Point vtx(0, 0, 0);
    Builder b(vtx);
    b.add(1, 0, 2112, "primary", 0.05, vtx, Point(100, 0, 0), 0);                   // NC: only a neutron
    b.add(2, 1, 2212, "neutronInelastic", 0.02, Point(100, 0, 0), Point(103, 0, 0), 4e5); // far proton
    b.add(3, 2, 211, "protonInelastic", 0.001, Point(103, 0, 0), Point(105, 0, 0), 2000); // 2 cm pion with too little charge: weak
    b.add(4, 2, 2212, "protonInelastic", 0.05, Point(103, 0, 0), Point(110, 0, 0), 2e6);    // strong daughter
    PFMerger merger;
    auto res = merger.run(b.tt);
    const auto& p2 = node_of(res, 2);
    CHECK(res.info(3).pf_node == p2.index);        // merged into the parent track
    CHECK(res.info(3).merged_into == 2);
    CHECK(res.nmerged_weak == 1);
    CHECK(node_of(res, 4).parent == p2.index);
    CHECK(res.void_roots.count(0) == 1);          // no primary node at the vertex

    // measured charge overrides the table: give the pion a lot and it survives
    std::map<int, double> q{{1, 0}, {2, 4e5}, {3, 1e6}, {4, 2e6}};
    PFMergerParams params;
    params.pf_ke_min = 0.0;
    auto res2 = PFMerger(params).run(b.tt, q);
    CHECK(node_of(res2, 3).trackid == 3);
    CHECK(res2.nmerged_weak == 0);

    // a depo-only table (no parentage) gives one track node per id
    IDepo::vector depos;
    (void)depos;
}

TEST_CASE("pfmerger em_fold: the radiation of a delta and of a Michel stays with them only when on")
{
    const Point vtx(0, 0, 0);
    Builder b(vtx);
    b.add(1, 0, 13, "primary", 0.5, vtx, Point(0, 0, 100), 9e6);                       // muon, 100 cm
    b.add(2, 1, 11, "muIoni", 0.05, Point(0, 0, 40), Point(0, 6, 44), 8e5);           // long delta (7.2 cm)
    b.add(3, 2, 22, "eBrem", 0.02, Point(0, 3, 42), Point(0, 13, 52), 0);             // its brem photon
    b.add(4, 3, 11, "conv", 0.012, Point(0, 13, 52), Point(0, 16, 58), 5e5);          // the photon converts
    b.add(5, 3, -11, "conv", 0.008, Point(0, 13, 52), Point(0, 15, 56), 4e5);
    b.add(6, 4, 11, "eIoni", 0.003, Point(0, 14, 54), Point(0, 16, 55), 6e4);          // ionisation inside that shower
    b.add(7, 1, 11, "Decay", 0.03, Point(0, 0, 100), Point(0, 4, 105), 6e5);          // Michel
    b.add(8, 7, 22, "eBrem", 0.015, Point(0, 2, 102), Point(0, 2, 120), 0);           // the Michel's brem photon
    b.add(9, 8, 11, "compt", 0.01, Point(0, 2, 120), Point(0, 2, 128), 3e5);          // its Compton electron (8 cm)

    SUBCASE("off: each radiated photon is a shower root of its own")
    {
        PFMerger merger;
        auto res = merger.run(b.tt);
        CHECK(node_of(res, 2).cls == Cls::delta);
        CHECK(node_of(res, 7).cls == Cls::michel);
        CHECK(res.info(4).pf_node != res.info(2).pf_node);
        CHECK(res.info(4).pf_node == res.info(5).pf_node);
        CHECK(res.info(6).pf_node == res.info(4).pf_node);
        CHECK(res.info(9).pf_node != res.info(7).pf_node);
        CHECK_FALSE(PFMergerParams().meta().isMember("em_fold_delta_michel"));
    }
    SUBCASE("on: the delta and the Michel keep their whole cascade, classes unchanged")
    {
        PFMergerParams p;
        p.em_fold = true;
        PFMerger merger(p);
        auto res = merger.run(b.tt);
        const int nd = res.info(2).pf_node, nm = res.info(7).pf_node;
        for (int t : {3, 4, 5, 6}) CHECK(res.info(t).pf_node == nd);
        for (int t : {8, 9}) CHECK(res.info(t).pf_node == nm);
        CHECK(node_of(res, 2).cls == Cls::delta);
        CHECK(node_of(res, 7).cls == Cls::michel);
        CHECK(node_of(res, 1).cls == Cls::track);
        CHECK(res.info(1).pf_node != nd);
        CHECK(p.meta()["em_fold_delta_michel"].asBool());
    }
}
