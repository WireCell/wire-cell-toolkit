/** BlobDepoLabel must reproduce BlobDepoFill bit for bit (the fill
 *  value of every blob) while DepoBlobAssoc keeps, per blob, the
 *  charge of every G4 track that contributed, and BlobDepoLabelTensor
 *  emits that table with the schema of Truth::Schema.
 *
 *  Setup as img/test/doctest_blobdepofill_faces.cxx: one slice of the
 *  dune10kt-1x2x6 anode tiled with GridTiling (loaded as a plugin, no
 *  link to WireCellImg), two drifted depos with G4 priors (track 7,
 *  pdg 13, energy 2.5 and track 8, pdg 2212, energy 4.0) at the
 *  centroids of two blobs.
 */

#include "WireCellTruth/DepoBlobAssoc.h"
#include "WireCellTruth/BlobDepoLabel.h"
#include "WireCellTruth/BlobDepoLabelTensor.h"
#include "WireCellTruth/Schema.h"

#include "WireCellAux/SimpleCluster.h"
#include "WireCellAux/SimpleDepo.h"
#include "WireCellAux/SimpleDepoSet.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleSlice.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IBlobDepoFill.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/ITiling.h"
#include "WireCellIface/IWirePlane.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Testing.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/GraphTools.h"
#include "WireCellUtil/doctest.h"

#include <cmath>
#include <map>

using namespace WireCell;

namespace {

    IAnodePlane::pointer make_fdhd_anode(int ident)
    {
        Testing::load_plugins({"WireCellGen", "WireCellSigProc", "WireCellImg"});
        {
            auto icfg = Factory::lookup<IConfigurable>("WireSchemaFile");
            auto cfg = icfg->default_configuration();
            cfg["filename"] = "dune10kt-1x2x6-wires-larsoft-v1.json.bz2";
            icfg->configure(cfg);
        }
        const std::string tn = String::format("AnodePlane:%d", ident);
        auto icfg = Factory::lookup_tn<IConfigurable>(tn);
        auto cfg = icfg->default_configuration();
        cfg["ident"] = ident;
        cfg["wire_schema"] = "WireSchemaFile";
        const double apa_plane = 39.5355 * units::mm, res_plane = 130.0155 * units::mm, cpa = 3629.1625 * units::mm;
        cfg["faces"][0]["anode"] = apa_plane;
        cfg["faces"][0]["response"] = res_plane;
        cfg["faces"][0]["cathode"] = cpa;
        cfg["faces"][1]["anode"] = -apa_plane;
        cfg["faces"][1]["response"] = -res_plane;
        cfg["faces"][1]["cathode"] = -cpa;
        icfg->configure(cfg);
        return Factory::find_tn<IAnodePlane>(tn);
    }

    ISlice::pointer make_slice(IAnodeFace::pointer face, double start, double span)
    {
        auto frame = std::make_shared<Aux::SimpleFrame>(100, 0.0);
        ISlice::map_t activity;
        // two separated activity groups per plane so the slice tiles
        // into several blobs
        for (const auto& plane : face->planes()) {
            const auto& chans = plane->channels();
            const size_t mid = chans.size() / 2;
            for (size_t i = mid; i < std::min(mid + 12, chans.size()); ++i) {
                activity[chans[i]] = ISlice::value_t(1000.0, 10.0);
            }
            for (size_t i = mid + 40; i < std::min(mid + 52, chans.size()); ++i) {
                activity[chans[i]] = ISlice::value_t(1000.0, 10.0);
            }
        }
        return std::make_shared<Aux::SimpleSlice>(frame, 0, start, span, activity);
    }

    Point blob_centroid(IAnodeFace::pointer face, const IBlob::pointer& blob)
    {
        const auto& coords = face->raygrid();
        Point cen(0, 0, 0);
        int ncorner = 0;
        for (const auto& [a, b] : blob->shape().corners()) {
            cen += coords.ray_crossing(a, b);
            ++ncorner;
        }
        REQUIRE(ncorner >= 3);
        return cen / ncorner;
    }

    struct Setup {
        IAnodePlane::pointer anode;
        IAnodeFace::pointer face;
        ICluster::pointer cluster;
        IDepoSet::pointer depos;
        int ident_a{-1}, ident_b{-1};   // blob idents the two depos sit in
        double qa{5000.0}, qb{3000.0};
    };

    Setup make_setup()
    {
        Setup s;
        s.anode = make_fdhd_anode(8);
        REQUIRE(s.anode);
        s.face = s.anode->face(0);
        REQUIRE(s.face);
        const double span = 2.0 * units::us;
        auto slice = make_slice(s.face, 0.0, span);

        auto tiling = Factory::lookup_tn<ITiling>("GridTiling:t19");
        {
            auto icfg = Factory::lookup_tn<IConfigurable>("GridTiling:t19");
            auto gcfg = icfg->default_configuration();
            gcfg["anode"] = String::format("AnodePlane:%d", s.anode->ident());
            gcfg["face"] = 0;
            icfg->configure(gcfg);
        }
        IBlobSet::pointer bs;
        REQUIRE((*tiling)(slice, bs));
        REQUIRE(bs);
        const auto blobs = bs->blobs();
        REQUIRE(blobs.size() >= 2);

        const double xresp = s.face->planes()[2]->pimpos()->origin()[0];
        const double behind = 3e-6 * units::mm;
        auto make_depo = [&](const IBlob::pointer& blob, double q, int id, int pdg, double energy) {
            const auto cen = blob_centroid(s.face, blob);
            const Point pos(xresp - s.face->dirx() * behind, cen.y(), cen.z());
            // the undrifted G4 depo (the prior) carries id, pdg and energy;
            // the drifted one carries the id only, as Drifter does.
            const Point pos0(pos.x() + 100 * units::cm, cen.y(), cen.z());
            auto prior = std::make_shared<Aux::SimpleDepo>(-1.0 * units::ms, pos0, q, nullptr, 0.0, 0.0, id, pdg, energy);
            return std::make_shared<Aux::SimpleDepo>(0.5 * span, pos, q, prior, 0.1 * units::mm, 1.0 * units::mm, id, 0, 1.0);
        };
        const size_t ib = blobs.size() - 1;
        s.ident_a = blobs[0]->ident();
        s.ident_b = blobs[ib]->ident();
        IDepo::vector dv{make_depo(blobs[0], s.qa, 7, 13, 2.5), make_depo(blobs[ib], s.qb, 8, 2212, 4.0)};
        s.depos = std::make_shared<Aux::SimpleDepoSet>(0, dv);

        cluster_graph_t g;
        for (const auto& b : blobs) {
            boost::add_vertex(cluster_node_t(IBlob::pointer(b)), g);
        }
        s.cluster = std::make_shared<Aux::SimpleCluster>(g, 100);
        return s;
    }

    Configuration fill_cfg(const Configuration& base)
    {
        Configuration cfg = base;
        cfg["speed"] = 1.6 * units::mm / units::us;
        cfg["time_offset"] = 0.0;
        return cfg;
    }

    std::map<int, float> blob_values(const ICluster::pointer& cl)
    {
        std::map<int, float> out;
        for (const auto& v : GraphTools::mir(boost::vertices(cl->graph()))) {
            const auto& node = cl->graph()[v];
            if (node.code() != 'b') continue;
            auto b = std::get<IBlob::pointer>(node.ptr);
            out[b->ident()] = b->value();
        }
        return out;
    }
}

TEST_CASE("blobdepolabel reproduces blobdepofill bit for bit and keeps per-track charge")
{
    auto s = make_setup();

    // reference: the production component, loaded as a plugin
    ICluster::pointer ref;
    {
        auto fill = Factory::lookup_tn<IBlobDepoFill>("BlobDepoFill:t19ref");
        auto icfg = Factory::lookup_tn<IConfigurable>("BlobDepoFill:t19ref");
        icfg->configure(fill_cfg(icfg->default_configuration()));
        REQUIRE((*fill)(std::make_tuple(s.cluster, s.depos), ref));
        REQUIRE(ref);
    }
    const auto refvals = blob_values(ref);
    REQUIRE(refvals.size() >= 2);
    double refsum = 0;
    for (const auto& [id, v] : refvals) refsum += v;
    CHECK(refsum > 0.5 * (s.qa + s.qb));

    // the new component, same interface
    ICluster::pointer out;
    Truth::BlobDepoLabel label;
    label.configure(fill_cfg(label.default_configuration()));
    REQUIRE(label(std::make_tuple(s.cluster, s.depos), out));
    REQUIRE(out);
    const auto vals = blob_values(out);
    REQUIRE(vals.size() == refvals.size());
    for (const auto& [id, v] : refvals) {
        REQUIRE(vals.count(id) == 1);
        CHECK(vals.at(id) == v);     // bitwise, not approximate
    }

    // the association itself
    Truth::DepoBlobAssocParams p;
    p.speed = 1.6 * units::mm / units::us;
    auto res = Truth::depo_blob_assoc(s.cluster, *s.depos->depos(), p);
    CHECK(res.nblobs == refvals.size());
    CHECK(res.ndepos == 2);
    std::map<int, double> per_track;
    size_t nfilled = 0;
    const auto& gr = s.cluster->graph();
    for (const auto& [vd, e] : res.blobs) {
        double sum = 0;
        for (const auto& [tid, qe] : e.by_track) {
            sum += qe.first;
            per_track[tid] += qe.first;
            CHECK(qe.second == doctest::Approx(qe.first / (tid == 7 ? s.qa : s.qb) * (tid == 7 ? 2.5 : 4.0)));
        }
        CHECK(sum == doctest::Approx(e.q_total).epsilon(1e-9));
        auto b = std::get<IBlob::pointer>(gr[vd].ptr);
        CHECK((float)e.q_total == refvals.at(b->ident()));
        if (!e.by_track.empty()) ++nfilled;
        if (b->ident() == s.ident_a) {
            REQUIRE(e.by_track.count(7) == 1);
            CHECK(e.by_track.at(7).first > 0.5 * s.qa);
        }
        if (b->ident() == s.ident_b) {
            REQUIRE(e.by_track.count(8) == 1);
            CHECK(e.by_track.at(8).first > 0.5 * s.qb);
        }
    }
    CHECK(nfilled >= 2);
    CHECK(per_track.size() == 2);
    CHECK(per_track[7] <= s.qa * (1 + 1e-9));
    CHECK(per_track[8] <= s.qb * (1 + 1e-9));
    CHECK(per_track[7] + per_track[8] == doctest::Approx(refsum).epsilon(1e-6));

    // the tensor output
    Truth::BlobDepoLabelTensor bt;
    {
        auto cfg = fill_cfg(bt.default_configuration());
        cfg["datapath"] = "clusters/%d";
        bt.configure(cfg);
    }
    ITensorSet::pointer ts;
    REQUIRE(bt(std::make_tuple(s.cluster, s.depos), ts));
    REQUIRE(ts);
    CHECK(ts->ident() == 100);
    ITensor::pointer table, summary;
    size_t ncluster = 0;
    for (const auto& t : *ts->tensors()) {
        const auto& md = t->metadata();
        const std::string dt = md["datatype"].asString();
        if (dt == "blobtruth") table = t;
        else if (dt == "blobtruth_sum") summary = t;
        else ++ncluster;
    }
    REQUIRE(table);
    REQUIRE(summary);
    CHECK(ncluster >= 1);
    CHECK(table->metadata()["datapath"].asString() == "clusters/100/blobtruth");
    REQUIRE(table->shape().size() == 2);
    CHECK(table->shape()[1] == 6);
    CHECK(table->metadata()["columns"].size() == 6);
    CHECK(table->metadata()["columns"][3].asString() == "trackid");
    REQUIRE(summary->shape().size() == 2);
    CHECK(summary->shape()[0] == refvals.size());
    CHECK(summary->shape()[1] == 7);
    CHECK(summary->metadata()["ndepos"].asUInt64() == 2);
    CHECK(summary->metadata()["q_depos"].asDouble() == doctest::Approx(s.qa + s.qb));
    CHECK(summary->metadata()["q_depos_window"].asDouble() == doctest::Approx(s.qa + s.qb));
    const double* srow = reinterpret_cast<const double*>(summary->data());
    size_t n_a = 0, n_ghost = 0;
    for (size_t r = 0; r < summary->shape()[0]; ++r) {
        const int bident = (int)srow[7 * r];
        CHECK((int)srow[7 * r + 1] == s.face->ident());
        const double q_true = srow[7 * r + 3];
        const int tid_top = (int)srow[7 * r + 5];
        CHECK((float)q_true == refvals.at(bident));
        if (bident == s.ident_a) { ++n_a; CHECK(tid_top == 7); }
        if (tid_top == Truth::tid_ghost) { ++n_ghost; CHECK(q_true == 0.0); }
    }
    CHECK(n_a == 1);
    const double* trow = reinterpret_cast<const double*>(table->data());
    double tsum = 0;
    for (size_t r = 0; r < table->shape()[0]; ++r) {
        tsum += trow[6 * r + 4];
        CHECK((trow[6 * r + 3] == 7 || trow[6 * r + 3] == 8));
    }
    CHECK(tsum == doctest::Approx(refsum).epsilon(1e-6));

    // extra smearing changes the numbers but keeps the bookkeeping consistent
    p.extra_smear_tran = 2.0 * units::mm;
    auto res2 = Truth::depo_blob_assoc(s.cluster, *s.depos->depos(), p);
    double sum2 = 0, tot2 = 0;
    for (const auto& [vd, e] : res2.blobs) {
        tot2 += e.q_total;
        for (const auto& [tid, qe] : e.by_track) sum2 += qe.first;
    }
    CHECK(sum2 == doctest::Approx(tot2).epsilon(1e-9));
    CHECK(tot2 > 0.0);
}

TEST_CASE("blobdepolabel folds larg4's negative EM-daughter ids to the ancestor when asked")
{
    // the two depos keyed on the LArSoft "-ancestor" convention: the proton's
    // deposit is a dropped EM daughter of track 1 (-1, the ghost sentinel of
    // the per-point labels) and the muon's of track 7 (-7).
    auto s = make_setup();
    IDepo::vector dv;
    for (const auto& d : *s.depos->depos()) {
        const int id = d->id() == 7 ? -7 : -1;
        auto prior = std::make_shared<Aux::SimpleDepo>(d->prior()->time(), d->prior()->pos(), d->charge(), nullptr,
                                                       0.0, 0.0, id, d->prior()->pdg(), d->prior()->energy());
        dv.push_back(std::make_shared<Aux::SimpleDepo>(d->time(), d->pos(), d->charge(), prior,
                                                       d->extent_long(), d->extent_tran(), id, 0, 1.0));
    }
    Truth::DepoBlobAssocParams p;
    p.speed = 1.6 * units::mm / units::us;

    // off (the R1 behaviour): the keys stay negative
    auto res = Truth::depo_blob_assoc(s.cluster, dv, p);
    std::set<int> keys;
    for (const auto& [vd, e] : res.blobs) for (const auto& [tid, qe] : e.by_track) keys.insert(tid);
    CHECK(keys.count(-1) == 1);
    CHECK(keys.count(-7) == 1);
    CHECK(keys.count(1) == 0);

    // on: the charge is keyed on the ancestor, totals unchanged
    p.fold_negative_ids = true;
    auto res2 = Truth::depo_blob_assoc(s.cluster, dv, p);
    std::set<int> keys2;
    double q1 = 0.0, q7 = 0.0;
    for (const auto& [vd, e] : res2.blobs) {
        for (const auto& [tid, qe] : e.by_track) {
            keys2.insert(tid);
            if (tid == 1) q1 += qe.first;
            if (tid == 7) q7 += qe.first;
        }
    }
    CHECK(keys2.count(1) == 1);
    CHECK(keys2.count(7) == 1);
    CHECK(keys2.count(-1) == 0);
    CHECK(keys2.count(-7) == 0);
    double qm1 = 0.0, qm7 = 0.0;
    for (const auto& [vd, e] : res.blobs) {
        for (const auto& [tid, qe] : e.by_track) {
            if (tid == -1) qm1 += qe.first;
            if (tid == -7) qm7 += qe.first;
        }
    }
    CHECK(q1 == qm1);
    CHECK(q7 == qm7);
    CHECK(q1 > 0.0);
    CHECK(q7 > 0.0);

    // the config key round-trips through the component
    Truth::BlobDepoLabel label;
    auto cfg = label.default_configuration();
    CHECK(cfg["fold_negative_ids"].asBool() == false);
}
