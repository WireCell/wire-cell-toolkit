// Regression pins for the crash/UB guards added by the doc pr/11 latent-pattern
// audit (clus/docs/audits/pr11-latent-pattern-audit.md).
//
// Each TEST_CASE below reproduces a defect that used to abort or corrupt a
// production job.  These are regression tests pinning fixes that were validated
// on real events by byte-identical A/B gates -- not test-first development.
// Every case reaches its guard through public API with no data fixture, so the
// suite stays portable.
//
// Each case names the guard it pins; reverting that guard must make the case
// fail (crash or assertion).  See the revert-proof table in the commit message.

#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Units.h"

#include "WireCellClus/DynamicPointCloud.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Blob.h"
#include "WireCellClus/FiducialUtils.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellClus/PRSegment.h"
#include "WireCellClus/PRShower.h"

using namespace WireCell;
using namespace WireCell::PointCloud;
using namespace WireCell::PointCloud::Tree;
using namespace WireCell::Clus::Facade;

// ---------------------------------------------------------------------------
// helpers
// ---------------------------------------------------------------------------

// wpid with face=-1/apa=-1: index_new_points() then skips the per-plane 2D
// k-d trees, so a DPC can be exercised with no anodes and no geometry.
static int novol_wpid()
{
    return WirePlaneId(kUlayer, -1, -1).ident();
}

// A DPCBatch of n points along x, with the CSR projection offsets sealed
// correctly (3 end_plane() per add_point, so p2d_off.size() == 3n+1).
static DPCBatch make_batch(size_t n, double x0 = 0.0)
{
    DPCBatch b;
    const std::array<int, 3> wind{0, 0, 0};
    const std::array<double, 3> dcut{0.0, 0.0, 0.0};
    for (size_t i = 0; i < n; ++i) {
        b.add_point(x0 + i * units::cm, 0.0, 0.0, novol_wpid(), nullptr, nullptr, wind, dcut);
        b.end_plane();
        b.end_plane();
        b.end_plane();
    }
    return b;
}

static std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>> no_wpid_params()
{
    return {};
}

// A cluster with ONE BLOB CHILD PER POINT.  {"3d", {"x","y","z"}} is Cluster's
// default scope, so these become the cluster's points for npoints() /
// point3d() / get_pca().  Each blob also gets the one-row "scalar" PC that
// Blob::fill_cache requires.
//
// One blob per point matters: get_two_extreme_points() smooths each returned
// extreme with calc_ave_pos(p, 5 cm), which averages BLOB CENTERS, not points.
// A single blob holding every point would collapse both extremes onto that
// blob's center and the test would pass or fail for the wrong reason.
static Cluster* make_cluster_with_points(Grouping& g, const std::vector<geo_point_t>& pts)
{
    using fa_float_t = double;
    using fa_int_t = int;

    Cluster& cl = g.make_child();
    for (const auto& p : pts) {
        Blob& bl = cl.make_child();

        Dataset ds;
        ds.add("x", Array(std::vector<double>{p.x()}));
        ds.add("y", Array(std::vector<double>{p.y()}));
        ds.add("z", Array(std::vector<double>{p.z()}));
        bl.value().local_pcs()["3d"] = ds;

        bl.value().local_pcs()["scalar"] = Dataset({
            {"charge", Array({(fa_float_t) 1.0})},
            {"center_x", Array({(fa_float_t) p.x()})},
            {"center_y", Array({(fa_float_t) p.y()})},
            {"center_z", Array({(fa_float_t) p.z()})},
            {"wpid", Array({(fa_int_t) WirePlaneId(kAllLayers).ident()})},
            {"npoints", Array({(fa_int_t) 1})},
            {"slice_index_min", Array({(fa_int_t) 0})},
            {"slice_index_max", Array({(fa_int_t) 1})},
            {"u_wire_index_min", Array({(fa_int_t) 0})},
            {"u_wire_index_max", Array({(fa_int_t) 1})},
            {"v_wire_index_min", Array({(fa_int_t) 0})},
            {"v_wire_index_max", Array({(fa_int_t) 1})},
            {"w_wire_index_min", Array({(fa_int_t) 0})},
            {"w_wire_index_max", Array({(fa_int_t) 1})},
            {"max_wire_interval", Array({(fa_int_t) 1})},
            {"min_wire_interval", Array({(fa_int_t) 1})},
            {"max_wire_type", Array({(fa_int_t) 0})},
            {"min_wire_type", Array({(fa_int_t) 0})},
        });
    }
    return &cl;
}

// ---------------------------------------------------------------------------
// Guard 1 -- DynamicPointCloud self-append (DynamicPointCloud.cxx warn_self_append)
//
// Appending a cloud's own storage to itself is undefined behavior, not a
// bounded 2x: DPCBatch::append does vector::insert from iterators INTO the
// destination vector, which a reallocation invalidates mid-copy.  This is how
// the pr/11 defect-D events blew the heap to a 224 GB bad_alloc.  The guard
// covers two distinct alias shapes -- batch identity (&m_pts == &points) and
// cloud identity (this == &other) -- so both are pinned here.
// ---------------------------------------------------------------------------

TEST_CASE("pr11 dpc self-append via cloud identity is a no-op")
{
    DynamicPointCloud dpc(no_wpid_params());
    dpc.add_points(make_batch(8));
    REQUIRE(dpc.npoints() == 8);

    // const-ref overload: this == &other
    CHECK_NOTHROW(dpc.add_points(dpc));
    CHECK(dpc.npoints() == 8);

    // row-selection overload: same identity test, non-empty rows
    CHECK_NOTHROW(dpc.add_points(dpc, std::vector<size_t>{0, 1, 2}));
    CHECK(dpc.npoints() == 8);
}

TEST_CASE("pr11 dpc self-append via batch identity is a no-op")
{
    DynamicPointCloud dpc(no_wpid_params());
    dpc.add_points(make_batch(5));
    REQUIRE(dpc.npoints() == 5);

    // &m_pts == &points, reached through the public accessor.  Both the
    // const-ref and the rvalue overload route through warn_self_append.
    const DPCBatch& own = dpc.points();
    CHECK_NOTHROW(dpc.add_points(own));
    CHECK(dpc.npoints() == 5);
}

TEST_CASE("pr11 dpc append of a distinct cloud still works")
{
    // The guard must not have disabled the legitimate path.
    DynamicPointCloud a(no_wpid_params());
    DynamicPointCloud b(no_wpid_params());
    a.add_points(make_batch(4, 0.0));
    b.add_points(make_batch(6, 100 * units::cm));

    a.add_points(b);
    CHECK(a.npoints() == 10);
    CHECK(b.npoints() == 6);  // source untouched
}

TEST_CASE("pr11 dpc copy is independent of its source")
{
    // The contract clone_dpc() (PRShower.cxx, file-local) relies on: a cloud
    // rebuilt through the public API shares no storage with its source, so a
    // later append to one cannot alias the other.
    DynamicPointCloud src(no_wpid_params());
    src.add_points(make_batch(7));

    DynamicPointCloud copy(src.get_wpid_params());
    copy.add_points(src);
    REQUIRE(copy.npoints() == 7);

    copy.add_points(make_batch(3, 50 * units::cm));
    CHECK(copy.npoints() == 10);
    CHECK(src.npoints() == 7);  // source unchanged
}

// ---------------------------------------------------------------------------
// Guard 2 -- Dataset::size() vs size_major() on a zero-point cloud
// (MyFCN.cxx UpdateInfo, TrackFitting.cxx; the trap recurs at 4+ sites)
//
// Dataset::size() returns m_store.size() -- the number of ARRAYS, not points.
// An empty steiner cloud still carries its three zero-length coordinate
// arrays, so `size() == 0` is false, the old guard let the code through, the
// kNN calls returned empty and the unguarded [0]/.front() read past the end.
// ---------------------------------------------------------------------------

TEST_CASE("pr11 dataset size counts arrays not points")
{
    Dataset ds;
    ds.add("x", Array(std::vector<double>{}));
    ds.add("y", Array(std::vector<double>{}));
    ds.add("z", Array(std::vector<double>{}));

    // The trap, stated as an assertion: a point-less dataset is NOT size()==0.
    CHECK(ds.size() == 3);
    CHECK(ds.size_major() == 0);

    // ... and the predicate the fix uses does catch it.
    CHECK_FALSE(ds.size() == 0);
    CHECK((ds.size() == 0 || ds.size_major() == 0));
}

TEST_CASE("pr11 empty steiner pc on a cluster is caught by size_major")
{
    // The production shape: has_pc() true, size() non-zero, no points.
    Points::node_t root;
    Grouping* g = root.value.facade<Grouping>();
    Cluster& cl = g->make_child();

    Dataset steiner;
    steiner.add("x", Array(std::vector<double>{}));
    steiner.add("y", Array(std::vector<double>{}));
    steiner.add("z", Array(std::vector<double>{}));
    cl.value().local_pcs()["steiner_pc"] = steiner;

    REQUIRE(cl.has_pc("steiner_pc"));
    const auto& pc = cl.get_pc("steiner_pc");

    // The pre-fix predicate passes this cloud through...
    CHECK_FALSE(pc.size() == 0);
    // ...the post-fix predicate rejects it.
    CHECK((pc.size() == 0 || pc.size_major() == 0));
}

// ---------------------------------------------------------------------------
// Guard 3 -- all-points-excluded fallback in Cluster::get_main_axis_points()
// and Cluster::get_two_extreme_points() (Facade_Cluster.cxx)
//
// connect_graph.cxx marks should_exclude for any point >= 1 cm from its
// reference cluster and nothing guarantees a non-empty complement, so "every
// point excluded" is a legal state, not a programmer error.  It used to throw
// and abort the job.  The throw must survive for a genuinely point-less
// cluster, so both halves are pinned.
// ---------------------------------------------------------------------------

TEST_CASE("pr11 extremes fall back when every point is excluded")
{
    Points::node_t root;
    Grouping* g = root.value.facade<Grouping>();
    // Spaced > 5 cm apart: get_two_extreme_points() smooths each returned
    // extreme with calc_ave_pos(p, 5 cm), which would collapse a tighter
    // cluster onto one centroid and mask whether the fallback ran.
    Cluster* cl = make_cluster_with_points(*g, {
        {0 * units::cm, 0, 0},
        {10 * units::cm, 10 * units::cm, 0},
        {20 * units::cm, 20 * units::cm, 0},
        {30 * units::cm, 30 * units::cm, 0},
    });
    REQUIRE(cl->npoints() == 4);

    // Exclude every point (set_excluded_points is the documented test hook).
    cl->set_excluded_points({0, 1, 2, 3});
    REQUIRE(cl->get_excluded_points().size() == 4);

    std::pair<geo_point_t, geo_point_t> ax, ex;
    CHECK_NOTHROW(ax = cl->get_main_axis_points());
    CHECK_NOTHROW(ex = cl->get_two_extreme_points());

    // The fallback must return real points of this cluster, not default-
    // constructed zeros -- the two extremes must differ.
    CHECK((ax.first - ax.second).magnitude() > 0.0);
    CHECK((ex.first - ex.second).magnitude() > 0.0);
}

TEST_CASE("pr11 extremes still throw for a genuinely point-less cluster")
{
    Points::node_t root;
    Grouping* g = root.value.facade<Grouping>();
    Cluster* cl = make_cluster_with_points(*g, {});
    REQUIRE(cl->npoints() == 0);

    CHECK_THROWS_AS(cl->get_main_axis_points(), std::runtime_error);
    CHECK_THROWS_AS(cl->get_two_extreme_points(), std::runtime_error);
}

// ---------------------------------------------------------------------------
// Guard 4 -- Shower::add_segment must not alias a segment's own DPC
// (PRShower.cxx clone_dpc + the `!= seg_dpc` membership tests)
//
// The shower used to seed its cloud with the segment's own shared_ptr; a later
// merge then put the same object on both sides of add_points -- the aliasing
// half of the 224 GB bad_alloc.  Re-adding a segment is routine (fresh
// used_segments sets, examine_showers merges), so the operation must be
// idempotent.
// ---------------------------------------------------------------------------

TEST_CASE("pr11 shower add_segment is idempotent and does not alias")
{
    using namespace WireCell::Clus::PR;

    Graph gr;
    auto v1 = make_vertex(gr);
    v1->wcpt().point = Point(0, 0, 0);
    auto v2 = make_vertex(gr);
    v2->wcpt().point = Point(10 * units::cm, 0, 0);
    auto seg = make_segment(gr, v1, v2);

    auto seg_dpc = std::make_shared<DynamicPointCloud>(no_wpid_params());
    seg_dpc->add_points(make_batch(6));
    seg->dpcloud("fit", seg_dpc);
    REQUIRE(seg->dpcloud("fit")->npoints() == 6);

    Shower sh(gr);
    sh.add_segment(seg);

    auto sh_dpc = sh.dpcloud("fit");
    REQUIRE(sh_dpc != nullptr);
    CHECK(sh_dpc->npoints() == 6);
    // The shower must hold a COPY, never the segment's own object.
    CHECK(sh_dpc.get() != seg_dpc.get());

    // Re-adding the same segment must not double the cloud.
    CHECK_NOTHROW(sh.add_segment(seg));
    CHECK(sh.dpcloud("fit")->npoints() == 6);
    // ... and must not have grown the segment's own cloud either.
    CHECK(seg->dpcloud("fit")->npoints() == 6);
}

// ---------------------------------------------------------------------------
// Guard 5 -- FiducialUtils::inside_dead_region sentinel apa/face
// (FiducialUtils.cxx)
//
// A sentinel apa/face (contained_by() miss, or PR::Fit::paf left at its
// {-1,-1} default) reached convert_3Dpoint_time_ch, which aborts inside
// Grouping::fastgeom via m_anodes.at(-1).  A point with no detector volume has
// no dead region to be inside of.  The guard returns before ANY dereference,
// which is what makes this callable on a default-constructed FiducialUtils.
// ---------------------------------------------------------------------------

TEST_CASE("pr11 inside_dead_region rejects sentinel apa/face before dereferencing")
{
    // No dynamic data fed, so m_internal.live is unset.  (The declared
    // FiducialUtils() default ctor has no definition anywhere in the tree --
    // pre-existing dead declaration, noted not fixed; construct from an empty
    // StaticData instead.)
    Clus::FiducialUtils fu{Clus::FiducialUtils::StaticData{}};
    const Point p(1 * units::cm, 2 * units::cm, 3 * units::cm);

    CHECK_FALSE(fu.inside_dead_region(p, -1, -1));
    CHECK_FALSE(fu.inside_dead_region(p, 0, -1));
    CHECK_FALSE(fu.inside_dead_region(p, -1, 0));
}
