/** TruthPointLabel on the tiled FD-HD slice of doctest_truth_blobdepolabel:
 *  a live tree built the way PointTreeBuilding::sample_live does
 *  (blob_clusters order, "3d" + "scalar" per blob node) gets a
 *  "truth" local point cloud on every blob node with the fixed key set
 *  of Schema, the blob with the muon depo is labelled tid 7 / pdg 13,
 *  the blob with the proton depo tid 8 / pdg 2212, every other blob is
 *  a ghost, the polyline fields follow the depo trajectory, and the
 *  tree survives the TensorDM round trip with the new cloud
 *  key-homogeneous.  The pgraph node TruthPointLabeler runs the same
 *  on tensors.
 */

#include "WireCellTruth/DepoBlobAssoc.h"
#include "WireCellTruth/PFMerger.h"
#include "WireCellTruth/Schema.h"
#include "WireCellTruth/TruthPointLabel.h"
#include "WireCellTruth/TruthPointLabeler.h"
#include "WireCellTruth/TruthTable.h"

#include "WireCellAux/ClusterHelpers.h"
#include "WireCellAux/SamplingHelpers.h"
#include "WireCellAux/SimpleCluster.h"
#include "WireCellAux/SimpleDepo.h"
#include "WireCellAux/SimpleDepoSet.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleSlice.h"
#include "WireCellAux/SimpleTensorSet.h"
#include "WireCellAux/TensorDMpointtree.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/ITiling.h"
#include "WireCellIface/IWirePlane.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/Testing.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/GraphTools.h"
#include "WireCellUtil/doctest.h"

#include <cmath>
#include <map>

using namespace WireCell;
using WireCell::PointCloud::Array;
using WireCell::PointCloud::Dataset;
using Points = WireCell::PointCloud::Tree::Points;

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

    const double speed = 1.6 * units::mm / units::us;
    const double tick = 0.5 * units::us;
    const double span = 2.0 * units::us;

    struct Setup {
        IAnodePlane::pointer anode;
        IAnodeFace::pointer face;
        ICluster::pointer cluster;
        IDepoSet::pointer depos;
        int ident_a{-1}, ident_b{-1};
        Point cen_a, cen_b;
        double xresp{0};
        std::map<int, Point> centroids;
        double qa{5000.0}, qb{3000.0};
    };

    Setup make_setup()
    {
        Setup s;
        s.anode = make_fdhd_anode(8);
        REQUIRE(s.anode);
        s.face = s.anode->face(0);
        REQUIRE(s.face);
        auto slice = make_slice(s.face, 0.0, span);

        auto tiling = Factory::lookup_tn<ITiling>("GridTiling:t19pl");
        {
            auto icfg = Factory::lookup_tn<IConfigurable>("GridTiling:t19pl");
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

        s.xresp = s.face->planes()[2]->pimpos()->origin()[0];
        const double behind = 3e-6 * units::mm;
        // a track of 5 depos per particle: 4 mm steps along +z through
        // the blob centroid (only the central one lands in the blob's
        // slice; the others are in neighbouring slices and carry the
        // trajectory)
        IDepo::vector dv;
        auto make_track = [&](const IBlob::pointer& blob, double q, int id, int pdg, double energy, Point& cen_out) {
            const auto cen = blob_centroid(s.face, blob);
            cen_out = cen;
            for (int k = -2; k <= 2; ++k) {
                const Point pos(s.xresp - s.face->dirx() * behind, cen.y(), cen.z() + k * 4 * units::mm);
                const Point pos0(pos.x() + 100 * units::cm, pos.y(), pos.z());
                const double t0 = -1.0 * units::ms + k * 0.1 * units::us;
                auto prior = std::make_shared<Aux::SimpleDepo>(t0, pos0, q, nullptr, 0.0, 0.0, id, pdg, energy);
                // one slice is 2 us: the k=0 depo is in slice 0, the others in slices k
                dv.push_back(std::make_shared<Aux::SimpleDepo>(0.5 * span + k * span, pos, q, prior,
                                                              0.1 * units::mm, 1.0 * units::mm, id, 0, 1.0));
            }
        };
        const size_t ib = blobs.size() - 1;
        s.ident_a = blobs[0]->ident();
        s.ident_b = blobs[ib]->ident();
        make_track(blobs[0], s.qa, 7, 13, 2.5, s.cen_a);
        make_track(blobs[ib], s.qb, 8, 2212, 4.0, s.cen_b);
        s.depos = std::make_shared<Aux::SimpleDepoSet>(0, dv);

        cluster_graph_t g;
        for (const auto& b : blobs) {
            boost::add_vertex(cluster_node_t(IBlob::pointer(b)), g);
            s.centroids[b->ident()] = blob_centroid(s.face, b);
        }
        s.cluster = std::make_shared<Aux::SimpleCluster>(g, 100);
        return s;
    }

    // The reco-frame position used by both the tree and the labels:
    // x = response plane + drift of the arrival time (t=0 at the plane).
    Point reco_position(const Setup& s, const IDepo& d)
    {
        return Point(s.xresp + s.face->dirx() * speed * d.time(), d.pos().y(), d.pos().z());
    }

    // A live tree as PointTreeBuilding::sample_live makes it: a node per
    // cluster, a node per blob with "3d" (3 points around the centroid)
    // and "scalar".
    Points::node_ptr make_tree(const Setup& s)
    {
        auto root = std::make_unique<Points::node_t>();
        const auto& gr = s.cluster->graph();
        auto clusters = Aux::blob_clusters(gr);
        for (const auto& [cid, vdescs] : clusters) {
            auto cnode = root->insert();
            for (const auto& vdesc : vdescs) {
                if (gr[vdesc].code() != 'b') continue;
                auto iblob = std::get<IBlob::pointer>(gr[vdesc].ptr);
                const Point cen = s.centroids.at(iblob->ident());
                const double x = s.xresp + s.face->dirx() * speed * (iblob->slice()->start() + 0.5 * span);
                std::vector<double> xs{x, x, x}, ys{cen.y(), cen.y() + 2 * units::mm, cen.y() - 2 * units::mm},
                    zs{cen.z(), cen.z() + 1 * units::mm, cen.z() - 1 * units::mm};
                Dataset pc3d;
                pc3d.add("x", Array(xs));
                pc3d.add("y", Array(ys));
                pc3d.add("z", Array(zs));
                Dataset scalar;
                Aux::fill_scalar_blob(scalar, *iblob, tick);
                Aux::fill_scalar_aux(scalar);
                Aux::fill_scalar_center(scalar, pc3d);
                PointCloud::Tree::named_pointclouds_t pcs;
                pcs.emplace("scalar", std::move(scalar));
                pcs.emplace("3d", std::move(pc3d));
                cnode->insert(Points(std::move(pcs)));
            }
        }
        return root;
    }

    Truth::TruthTable make_table()
    {
        Truth::TruthTable tt;
        const auto& nuc = Truth::truth_nu_columns();
        std::vector<double> nu(nuc.size(), 0.0);
        nu[Truth::column_index(nuc, "nu_idx")] = 0;
        nu[Truth::column_index(nuc, "pdg")] = 14;
        tt.nu.push_back(nu);
        const auto& pfc = Truth::truth_pf_columns();
        auto row = [&](int tid, int pdg, double ke, double q) {
            std::vector<double> r(pfc.size(), 0.0);
            r[Truth::column_index(pfc, "trackid")] = tid;
            r[Truth::column_index(pfc, "orig_trackid")] = tid * 100;
            r[Truth::column_index(pfc, "pdg")] = pdg;
            r[Truth::column_index(pfc, "process")] = Truth::process_code("primary");
            r[Truth::column_index(pfc, "KE")] = ke;
            r[Truth::column_index(pfc, "end_z")] = 50;    // 50 cm long
            r[Truth::column_index(pfc, "nvis_q")] = q;
            r[Truth::column_index(pfc, "nu_row")] = 0;
            tt.pf.push_back(r);
        };
        row(7, 13, 0.5, 25000);
        row(8, 2212, 0.1, 15000);
        tt.source = "test";
        return tt;
    }

    // every blob node of the tree, depth first
    std::vector<Points::node_t*> blob_nodes(Points::node_t& root)
    {
        std::vector<Points::node_t*> out;
        for (auto* c : root.children()) {
            for (auto* b : c->children()) out.push_back(b);
        }
        return out;
    }
}

TEST_CASE("truth point labels on a sampled live tree")
{
    auto s = make_setup();
    auto root = make_tree(s);
    const size_t nblobs = blob_nodes(*root).size();
    REQUIRE(nblobs >= 2);

    Truth::DepoBlobAssocParams ap;
    ap.speed = speed;
    ap.keep_depos = true;
    auto res = Truth::depo_blob_assoc(s.cluster, *s.depos->depos(), ap);
    auto tt = make_table();
    std::map<int, double> q_vis;
    for (const auto& [vd, e] : res.blobs) {
        for (const auto& [tid, qe] : e.by_track) q_vis[tid] += qe.first;
    }
    Truth::PFMergerParams pp;
    pp.q_min = 100;
    auto pf = Truth::PFMerger(pp).run(tt, q_vis);
    REQUIRE(pf.nodes.size() == 2);

    Truth::depo_position_f position = [&](const IDepo& d) { return reco_position(s, d); };
    auto polylines = Truth::build_polylines(*s.depos->depos(), ap, position);
    REQUIRE(polylines.size() == 2);
    CHECK(polylines.at(7)->nvertices() == 5);
    // in the reco frame the depos of neighbouring slices sit one slice
    // of drift apart in x: the trajectory is diagonal
    const double step = std::hypot(4 * units::mm, speed * span);
    CHECK(polylines.at(7)->length() == doctest::Approx(4 * step));

    Truth::PointLabelParams lp;
    lp.tick = tick;
    auto sum = Truth::label_live_tree(*root, s.cluster, *s.depos->depos(), res, tt, pf, polylines, lp, position);
    CHECK(sum.nnodes == nblobs);
    CHECK(sum.npoints == 3 * nblobs);
    CHECK(sum.nblobs_skipped == 0);
    CHECK(sum.nghost + 2 <= nblobs);   // at least the two filled blobs are not ghosts

    // every node has the truth PC with the fixed key set and the 3d size
    const auto& keys = Truth::truth_pc_keys();
    size_t n7 = 0, n8 = 0, nghost = 0;
    for (auto* node : blob_nodes(*root)) {
        const auto& lpcs = node->value.local_pcs();
        REQUIRE(lpcs.count(Truth::truth_pc_name) == 1);
        const auto& tpc = lpcs.at(Truth::truth_pc_name);
        CHECK(tpc.size_major() == lpcs.at("3d").size_major());
        auto got = tpc.keys();
        std::vector<std::string> want(keys.begin(), keys.end());
        std::sort(want.begin(), want.end());
        CHECK(got == want);
        const auto tid = tpc.get("tid")->elements<int>();
        const auto pdg = tpc.get("pdg")->elements<int>();
        const auto cls = tpc.get("cls")->elements<int>();
        const auto orig = tpc.get("orig_tid")->elements<int>();
        const auto valid = tpc.get("instance_valid")->elements<int>();
        const auto frac = tpc.get("pid_frac")->elements<double>();
        const auto dtraj = tpc.get("d_traj")->elements<double>();
        const auto dirz = tpc.get("dir_z")->elements<double>();
        const auto st = tpc.get("s_t")->elements<double>();
        const auto dedx = tpc.get("dedx_true")->elements<double>();
        const auto vvtx = tpc.get("v_heat_vtx")->elements<double>();
        const auto iid = tpc.get("iid")->elements<int>();
        const auto pfn = tpc.get("pf_node")->elements<int>();
        for (size_t i = 0; i < tid.size(); ++i) {
            if (tid[i] == Truth::tid_ghost) {
                ++nghost;
                CHECK(valid[i] == 0);
                CHECK(pdg[i] == 0);
                CHECK(pfn[i] == -1);
                CHECK(dtraj[i] == doctest::Approx(lp.d_cap));
                continue;
            }
            CHECK(valid[i] == 1);
            CHECK(frac[i] > 0.0);
            CHECK(frac[i] <= 1.0);
            CHECK(cls[i] == (int)Truth::Cls::track);
            CHECK(iid[i] == 0);
            CHECK(pfn[i] >= 0);
            CHECK(orig[i] == tid[i] * 100);
            // the trajectory runs diagonally in (x, z) through the point
            CHECK(std::abs(dirz[i]) == doctest::Approx(4 * units::mm / step));
            CHECK(dtraj[i] < 5 * units::mm);
            CHECK(st[i] >= 0.0);
            CHECK(st[i] <= 1.0);
            CHECK(dedx[i] > 0.0);
            CHECK(vvtx[i] >= 0.0);
            if (tid[i] == 7) { ++n7; CHECK(pdg[i] == 13); }
            else if (tid[i] == 8) { ++n8; CHECK(pdg[i] == 2212); }
            else { CHECK(false); }
        }
    }
    CHECK(n7 >= 3);
    CHECK(n8 >= 3);
    CHECK(nghost == 3 * sum.nghost);

    // round trip through the tensor data model keeps every array
    auto tens = Aux::TensorDM::as_tensors(*root, "pointtrees/0/live");
    auto back = Aux::TensorDM::as_pctree(tens, "pointtrees/0/live");
    REQUIRE(back);
    for (auto* node : blob_nodes(*back)) {
        const auto& lpcs = node->value.local_pcs();
        REQUIRE(lpcs.count("truth") == 1);
        CHECK(lpcs.at("truth").keys().size() == keys.size());
        CHECK(lpcs.at("truth").size_major() == 3);
    }

    // the dedx: 2.5 energy per step for the muon
    auto hit = polylines.at(7)->nearest(s.cen_a);
    CHECK(hit.dedx == doctest::Approx(2.5 / step));
}

TEST_CASE("truth point labeler node relabels a tensor set")
{
    auto s = make_setup();
    auto root = make_tree(s);
    auto tens = Aux::TensorDM::as_tensors(*root, "pointtrees/100/live");
    Configuration md;
    md["ident"] = 100;
    auto ints = std::make_shared<Aux::SimpleTensorSet>(100, md, std::make_shared<ITensor::vector>(tens));

    Truth::TruthPointLabeler node;
    auto cfg = node.default_configuration();
    cfg["speed"] = speed;
    cfg["anode"] = String::format("AnodePlane:%d", s.anode->ident());
    cfg["face"] = 0;
    cfg["drift_speed"] = speed;
    cfg["reco_time_offset"] = 0.0;
    cfg["tick"] = tick;
    cfg["pfmerger"]["q_min"] = 100;
    node.configure(cfg);

    ITensorSet::pointer out;
    REQUIRE(node(std::make_tuple(ints, s.cluster, s.depos), out));
    REQUIRE(out);
    CHECK(out->ident() == 100);
    CHECK(out->metadata().isMember("labels_meta"));
    CHECK(out->metadata()["labels_meta"]["pf_source"].asString() == "depos");
    CHECK(out->metadata()["labels_meta"]["schema_version"].asInt() == Truth::schema_version);

    std::map<std::string, ITensor::pointer> bytype;
    for (const auto& t : *out->tensors()) bytype[t->metadata()["datatype"].asString()] = t;
    REQUIRE(bytype.count("pctree"));
    REQUIRE(bytype.count("blobtruth"));
    REQUIRE(bytype.count("truth_nu"));
    REQUIRE(bytype.count("truth_pf"));
    CHECK(bytype.at("truth_pf")->metadata()["datapath"].asString() == "truth/100/pf");
    CHECK(bytype.at("truth_pf")->shape()[0] == 2);   // the two depo keys
    CHECK(bytype.at("blobtruth")->metadata()["datapath"].asString() == "clusters/100/blobtruth");

    auto back = Aux::TensorDM::as_pctree(*out->tensors(), "pointtrees/100/live");
    REQUIRE(back);
    size_t nlab = 0;
    for (auto* n : blob_nodes(*back)) {
        REQUIRE(n->value.local_pcs().count("truth") == 1);
        const auto tid = n->value.local_pcs().at("truth").get("tid")->elements<int>();
        for (int t : tid) if (t >= 0) ++nlab;
    }
    CHECK(nlab >= 6);

    // EOS passes through
    REQUIRE(node(std::make_tuple(nullptr, nullptr, nullptr), out));
    CHECK(out == nullptr);

    // A set without a point tree at inpath (a face PointTreeBuilding never
    // sampled, PDHD wall faces) passes through unchanged: same object, no
    // truth tensors added.
    auto dead = Aux::TensorDM::as_tensors(*make_tree(s), "pointtrees/101/dead");
    auto ints2 = std::make_shared<Aux::SimpleTensorSet>(101, md, std::make_shared<ITensor::vector>(dead));
    ITensorSet::pointer out2;
    REQUIRE(node(std::make_tuple(ints2, s.cluster, s.depos), out2));
    REQUIRE(out2);
    CHECK(out2.get() == ints2.get());
    CHECK(out2->tensors()->size() == dead.size());
    CHECK(!out2->metadata().isMember("labels_meta"));
}

// pr_ml doc 23 sec 7.4: the two default-off trajectory knobs of the labeler.
namespace {

    ITensorSet::pointer run_node(const Setup& s, IDepoSet::pointer depos, bool drift_volume_only, bool per_particle)
    {
        auto root = make_tree(s);
        auto tens = Aux::TensorDM::as_tensors(*root, "pointtrees/100/live");
        Configuration md;
        md["ident"] = 100;
        auto ints = std::make_shared<Aux::SimpleTensorSet>(100, md, std::make_shared<ITensor::vector>(tens));
        Truth::TruthPointLabeler node;
        auto cfg = node.default_configuration();
        cfg["speed"] = speed;
        cfg["anode"] = String::format("AnodePlane:%d", s.anode->ident());
        cfg["face"] = 0;
        cfg["drift_speed"] = speed;
        // Aux::time2drift places x from the collection wires; put the
        // depos on the tree's x (response plane + drift of the time)
        const double xorig = s.face->planes()[2]->wires().front()->center().x();
        cfg["reco_time_offset"] = s.face->dirx() * (s.xresp - xorig) / speed;
        cfg["tick"] = tick;
        cfg["pfmerger"]["q_min"] = 100;
        if (drift_volume_only) cfg["polyline_drift_volume_only"] = true;
        if (per_particle) cfg["polyline_per_particle"] = true;
        node.configure(cfg);
        ITensorSet::pointer out;
        REQUIRE(node(std::make_tuple(ints, s.cluster, depos), out));
        REQUIRE(out);
        return out;
    }

    // every truth column of the relabelled tree, concatenated over the blob nodes
    std::map<std::string, std::vector<double>> truth_columns(const ITensorSet::pointer& out)
    {
        auto back = Aux::TensorDM::as_pctree(*out->tensors(), "pointtrees/100/live");
        REQUIRE(back);
        std::map<std::string, std::vector<double>> cols;
        for (auto* n : blob_nodes(*back)) {
            const auto& tpc = n->value.local_pcs().at("truth");
            for (const auto& k : tpc.keys()) {
                const auto arr = tpc.get(k);
                auto& c = cols[k];
                if (arr->is_type<int>()) {
                    for (int v : arr->elements<int>()) c.push_back(v);
                }
                else {
                    for (double v : arr->elements<double>()) c.push_back(v);
                }
            }
        }
        return cols;
    }
}

TEST_CASE("truth point labels: per-particle polylines keep folded daughters apart")
{
    // key 5 = an ancestor along z at y = 0 and its unsaved EM daughter
    // (larg4 id -5, folded) along z at y = 20 mm, their depos interleaved
    // in time: one polyline per key zigzags between the two lines.
    IDepo::vector dv;
    for (int k = 0; k < 4; ++k) {
        const double z = k * 10 * units::mm;
        auto pa = std::make_shared<Aux::SimpleDepo>(2 * k * units::ns, Point(0, 0, z), 1.0, nullptr, 0, 0, 5, 11, 1.0);
        auto pd = std::make_shared<Aux::SimpleDepo>((2 * k + 1) * units::ns, Point(0, 20 * units::mm, z), 1.0, nullptr, 0, 0, -5, 22, 1.0);
        dv.push_back(std::make_shared<Aux::SimpleDepo>(0, pa->pos(), 1.0, pa, 0, 0, 5));
        dv.push_back(std::make_shared<Aux::SimpleDepo>(0, pd->pos(), 1.0, pd, 0, 0, -5));
    }
    Truth::DepoBlobAssocParams ap;
    ap.fold_negative_ids = true;
    Truth::depo_position_f position = [](const IDepo& d) { return d.pos(); };

    auto keyed = Truth::build_polylines(dv, ap, position);
    REQUIRE(keyed.size() == 1);
    CHECK(keyed.at(5)->nvertices() == 8);
    auto parts = Truth::build_particle_polylines(dv, ap, position);
    REQUIRE(parts.size() == 1);
    REQUIRE(parts.at(5).size() == 2);   // particles -5 and 5, in id order
    CHECK(parts.at(5)[0]->start().y() == doctest::Approx(20 * units::mm));
    CHECK(parts.at(5)[1]->start().y() == doctest::Approx(0));
    CHECK(parts.at(5)[1]->length() == doctest::Approx(30 * units::mm));

    // a point on the ancestor between two of its depos: on the particle,
    // off the zigzag's chords
    const Point p(0, 0, 15 * units::mm);
    CHECK(keyed.at(5)->nearest(p).dist > 4 * units::mm);   // nearest chord D1 -> A2 at 4.47 mm
    double dmin = 1e9;
    for (const auto& q : parts.at(5)) dmin = std::min(dmin, q->nearest(p).dist);
    CHECK(dmin == doctest::Approx(0).epsilon(1e-9));

    // Stage A index mode: depo ids are rows of sed_ids, id_map gives the
    // (folded) key and orig_ids the particle
    IDepo::vector dvi;
    const std::vector<int> rows{0, 1, 2, 3};
    for (int r : rows) {
        auto pr = std::make_shared<Aux::SimpleDepo>(r * units::ns, Point(0, 0, r * 10 * units::mm), 1.0, nullptr, 0, 0, r);
        dvi.push_back(std::make_shared<Aux::SimpleDepo>(0, pr->pos(), 1.0, pr, 0, 0, r));
    }
    ap.id_map = std::make_shared<std::vector<int>>(std::vector<int>{5, 5, -5, -5});
    const std::vector<int> orig{5, 5, 11, 12};
    auto pi = Truth::build_particle_polylines(dvi, ap, position, &orig);
    REQUIRE(pi.size() == 1);
    CHECK(pi.at(5).size() == 3);
    CHECK(Truth::build_polylines(dvi, ap, position).at(5)->nvertices() == 4);
}

TEST_CASE("truth point labeler: polyline_drift_volume_only drops the other drift volume's depos, per_particle is neutral on single-particle keys")
{
    auto s = make_setup();
    const auto base = truth_columns(run_node(s, s.depos, false, false));
    REQUIRE(base.count("d_traj"));

    // a depo of track 7 in the other drift volume (face 1 side), in time
    // between two of the track's depos and 30 mm off in z: the legacy
    // polyline places it with face 0's drift relation and detours
    IDepo::vector dv(s.depos->depos()->begin(), s.depos->depos()->end());
    {
        const Point pos(-s.xresp, s.cen_a.y(), s.cen_a.z() + 30 * units::mm);
        auto prior = std::make_shared<Aux::SimpleDepo>(-1.0 * units::ms + 0.05 * units::us, pos, s.qa, nullptr,
                                                      0.0, 0.0, 7, 13, 2.5);
        dv.push_back(std::make_shared<Aux::SimpleDepo>(span, pos, s.qa, prior, 0.1 * units::mm, 1.0 * units::mm, 7, 0, 1.0));
    }
    auto mirrored = std::make_shared<Aux::SimpleDepoSet>(0, dv);

    const auto legacy = truth_columns(run_node(s, mirrored, false, false));
    CHECK(legacy.at("tid") == base.at("tid"));          // the association never sees it
    CHECK(legacy.at("s_t") != base.at("s_t"));          // the legacy trajectory does
    for (size_t i = 0; i < base.at("tid").size(); ++i) {
        if (base.at("tid")[i] >= 0) CHECK(base.at("d_traj")[i] < 5 * units::mm);   // the frames line up
    }
    const auto fixed = truth_columns(run_node(s, mirrored, true, false));
    CHECK(fixed == base);

    // the same depo on this face's drift side but beyond the face's z
    // extent (over a neighbouring anode) is part of the trajectory and stays
    IDepo::vector dv2(s.depos->depos()->begin(), s.depos->depos()->end());
    {
        const auto sens = s.face->sensitive();
        const auto& bb = sens.bounds();
        const double zout = std::max(bb.first.z(), bb.second.z()) + 100 * units::mm;
        const Point pos(s.xresp - s.face->dirx() * 3e-6 * units::mm, s.cen_a.y(), zout);
        auto prior = std::make_shared<Aux::SimpleDepo>(-1.0 * units::ms + 0.05 * units::us, pos, s.qa, nullptr,
                                                      0.0, 0.0, 7, 13, 2.5);
        dv2.push_back(std::make_shared<Aux::SimpleDepo>(span, pos, s.qa, prior, 0.1 * units::mm, 1.0 * units::mm, 7, 0, 1.0));
    }
    auto beside = std::make_shared<Aux::SimpleDepoSet>(0, dv2);
    const auto beside_legacy = truth_columns(run_node(s, beside, false, false));
    CHECK(beside_legacy.at("s_t") != base.at("s_t"));
    CHECK(truth_columns(run_node(s, beside, true, false)) == beside_legacy);

    // keys 7 and 8 are one particle each: the per-particle path gives the
    // legacy bytes
    auto out_pp = run_node(s, s.depos, false, true);
    CHECK(truth_columns(out_pp) == base);
    const auto pl = out_pp->metadata()["labels_meta"]["point_label"];
    CHECK(pl["polyline_per_particle"].asBool());
    CHECK(!pl.isMember("polyline_drift_volume_only"));
    // off: the metadata keeps its legacy keys
    const auto pl0 = run_node(s, s.depos, false, false)->metadata()["labels_meta"]["point_label"];
    CHECK(!pl0.isMember("polyline_per_particle"));
    CHECK(!pl0.isMember("polyline_drift_volume_only"));
}
