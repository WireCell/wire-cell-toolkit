/** doc pr/98 -- regression pins for TrackFitting::update_association's
 *  -1.0-sentinel guard (the fit_exclusion / doc pr/30 F1 path).
 *
 * REPRODUCER FIRST.  update_association arbitrates each candidate 2-D cell
 * between the segment being fit and every competitor segment via
 * segment_get_closest_2d_distances, whose underlying
 * DynamicPointCloud::get_closest_2d_point_info returns a -1.0 SENTINEL when
 * the segment has no points at the queried (plane, face, apa) -- e.g. a
 * competitor on the other drift face of a cathode-crossing cluster.  The
 * prototype's ToyPointCloud::get_closest_2d_dis returns 1e9 there ("no
 * measurement", prototype data/src/ToyPointCloud.cxx:415-437).  Unguarded:
 *
 *   - a single cross-face competitor poisons min_dis1_track to -1.0, so
 *     min_dis_track < min_dis1_track is false for EVERY cell and the whole
 *     association set beyond the 0.3 cm keep-floor is stripped (case 1);
 *   - a fitting segment with no points at the queried face gets
 *     min_dis_track = -1.0, which wins every comparison and keeps every
 *     cell unconditionally (case 2).
 *
 * The guard (TrackFitting.cxx update_association, all three plane loops)
 * maps negative distances to 1e9 for prototype parity.  Reverting the guard
 * makes cases 1 and 2 fail; case 3 pins the untouched keep rule.
 *
 * The path is reachable only with the fit_exclusion knob on (C++ default
 * false), so these pins guard a knob-on path and prove nothing about the
 * knob-off production path (which the pr/98 hash gate covers).
 *
 * Fixture notes: geometry is injected through TrackFittingTestHarness (the
 * doc pr/98 friend seam) because BuildGeometry needs live IDetectorVolumes.
 * Synthetic single-APA geometry, faces 0 and 1, all offsets 0:
 * angle_u = +90 deg so the U projection is (x, -y) and a wire index w sits
 * at transverse coordinate w * pitch; pitch = 1 cm; slope_t = 1/cm so time
 * index t sits at drift x = t cm.  All test points share x = 0, so 2-D
 * distances are purely transverse and readable by eye.
 */

#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/PRSegment.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellClus/DynamicPointCloud.h"
#include "WireCellUtil/Units.h"

#include "WireCellUtil/doctest.h"

#include <cmath>
#include <memory>
#include <algorithm>
#include <random>
#include <set>
#include <vector>

using namespace WireCell;
using namespace WireCell::Clus;
using WireCell::Clus::Facade::DPCBatch;
using WireCell::Clus::Facade::DynamicPointCloud;

namespace WireCell::Clus {

    // The doc pr/98 test seam declared as a friend in TrackFitting.h.
    struct TrackFittingTestHarness {
        static void set_geometry(TrackFitting& tf, int apa, double angle_u, double angle_v,
                                 double angle_w, double pitch, double slope_t)
        {
            for (int face = 0; face < 2; ++face) {
                WirePlaneId wpid(kAllLayers, face, apa);
                tf.wpid_offsets[wpid] = std::make_tuple(0.0, 0.0, 0.0, 0.0);
                tf.wpid_slopes[wpid] =
                    std::make_tuple(slope_t, std::make_pair(-sin(angle_u) / pitch, cos(angle_u) / pitch),
                                    std::make_pair(-sin(angle_v) / pitch, cos(angle_v) / pitch),
                                    std::make_pair(-sin(angle_w) / pitch, cos(angle_w) / pitch));
            }
        }
        // doc pdvd/45: the per-volume drift-frame shift form_map_graph would
        // compute from dirx * cluster_t0 * v_drift when excl_t0_frame > 0.
        static void set_excl_x_shift(TrackFitting& tf, int apa, int face, double shift_x)
        {
            tf.m_excl_x_shift[WirePlaneId(kAllLayers, face, apa).ident()] = shift_x;
        }
        static void clear_excl_x_shift(TrackFitting& tf) { tf.m_excl_x_shift.clear(); }
        // pr_ml doc 87b: the exact competitor index, as form_map_graph scopes it.
        static void excl_index_begin(TrackFitting& tf, const std::vector<std::shared_ptr<PR::Segment>>& segs)
        {
            tf.exclusion_index_begin(segs);
        }
        static void excl_index_end(TrackFitting& tf) { tf.exclusion_index_end(); }
        static size_t excl_index_used(const TrackFitting& tf) { return tf.m_excl_index_used; }
    };

}  // namespace WireCell::Clus

namespace {

    const double kPitch = 1.0 * units::cm;
    const double kAngleU = M_PI / 2;   // U projection = (x, -y)
    const double kAngleV = -M_PI / 2;  // V projection = (x, +y)
    const double kAngleW = 0.0;        // W projection = (x, +z)

    // Shared wpid params (both faces of apa 0) for every DynamicPointCloud:
    // get_closest_2d_point_info raises if the queried volume has no entry,
    // and the sentinel under test is the EMPTY-TREE return, not that raise.
    std::map<WirePlaneId, std::tuple<Facade::geo_point_t, double, double, double>> both_face_params()
    {
        std::map<WirePlaneId, std::tuple<Facade::geo_point_t, double, double, double>> params;
        for (int face = 0; face < 2; ++face) {
            params[WirePlaneId(kAllLayers, face, 0)] =
                std::make_tuple(Facade::geo_point_t(0, 0, 0), kAngleU, kAngleV, kAngleW);
        }
        return params;
    }

    // A "fit" cloud whose points all live on `face` of apa 0, with the
    // standard per-plane projections (the same cos*z - sin*y convention
    // DynamicPointCloud::get_closest_2d_point_info queries with).
    std::shared_ptr<DynamicPointCloud> make_face_cloud(int face, const std::vector<Point>& pts)
    {
        DPCBatch batch;
        const double angles[3] = {kAngleU, kAngleV, kAngleW};
        const WirePlaneLayer_t layers[3] = {kUlayer, kVlayer, kWlayer};
        for (const auto& p : pts) {
            batch.add_point(p.x(), p.y(), p.z(), WirePlaneId(kAllLayers, face, 0).ident(), nullptr,
                            nullptr, {0, 0, 0}, {0.0, 0.0, 0.0});
            for (int plane = 0; plane < 3; ++plane) {
                batch.add_proj(p.x(), cos(angles[plane]) * p.z() - sin(angles[plane]) * p.y(),
                               WirePlaneId(layers[plane], face, 0).ident());
                batch.end_plane();
            }
        }
        auto dpc = std::make_shared<DynamicPointCloud>(both_face_params());
        dpc->add_points(std::move(batch));
        return dpc;
    }

    struct Rig {
        TrackFitting tf;
        std::shared_ptr<PR::Graph> graph = std::make_shared<PR::Graph>();

        Rig()
        {
            // add_graph BEFORE any segment exists: sync_from_graph would
            // dereference segment->cluster()->grouping() otherwise.
            tf.add_graph(graph);
            TrackFittingTestHarness::set_geometry(tf, 0, kAngleU, kAngleV, kAngleW, kPitch,
                                                  1.0 / units::cm);
        }

        // A graph segment whose "fit" cloud holds `pts` on `face` of apa 0.
        std::shared_ptr<PR::Segment> segment_on_face(int face, const std::vector<Point>& pts)
        {
            auto v1 = PR::make_vertex(*graph);
            auto v2 = PR::make_vertex(*graph);
            auto seg = PR::make_segment(*graph, v1, v2);
            seg->dpcloud("fit", make_face_cloud(face, pts));
            return seg;
        }

        // Run update_association on one U-plane cell at (face 0, apa 0,
        // time 0, wire 2) -- i.e. drift x = 0, transverse 2 cm -- and report
        // whether the cell survived.
        bool cell_kept(std::shared_ptr<PR::Segment> seg,
                       const std::vector<std::shared_ptr<PR::Segment>>& all)
        {
            TrackFitting::PlaneData ut, vt, wt;
            ut.associated_2d_points.insert(TrackFitting::Coord2D(0, 0, 0, 2, 0, kUlayer));
            tf.update_association(seg, all, ut, vt, wt);
            return ut.associated_2d_points.size() == 1;
        }
    };

}  // namespace

// Case 1 -- competitor sentinel.  The fitted segment is 1 cm from the cell
// (beyond the 0.3 cm floor); the only competitor lives entirely on the OTHER
// face, so its face-0 U tree is empty and it answers -1.0.  Unguarded, that
// poisons min_dis1_track and strips the cell; with the guard the competitor
// contributes no measurement (1e9, prototype parity) and the cell is kept.
TEST_CASE("pr98 update_association cross-face competitor must not strip")
{
    Rig rig;
    auto fitted = rig.segment_on_face(0, {Point(0, -1.0 * units::cm, 0)});     // 1 cm from cell
    auto crossface = rig.segment_on_face(1, {Point(0, -1.0 * units::cm, 0)});  // face 1 only

    CHECK(rig.cell_kept(fitted, {fitted, crossface}));
}

// Case 2 -- own-segment sentinel.  The fitted segment has no points at the
// queried face (min_dis_track = -1.0 unguarded, which wins every comparison
// and keeps the cell); the competitor is a real 0.5 cm away.  Prototype
// behaviour: no measurement (1e9) loses to any finite competitor => drop.
TEST_CASE("pr98 update_association own-segment sentinel must not always-keep")
{
    Rig rig;
    auto fitted = rig.segment_on_face(1, {Point(0, -1.0 * units::cm, 0)});     // face 1 only
    auto competitor = rig.segment_on_face(0, {Point(0, -1.5 * units::cm, 0)});  // 0.5 cm from cell

    CHECK_FALSE(rig.cell_kept(fitted, {fitted, competitor}));
}

// Case 3 -- the keep rule itself, untouched by the guard: a strictly closer
// real competitor strips the cell; the 0.3 cm floor keeps it regardless.
TEST_CASE("pr98 update_association keep rule with real measurements")
{
    Rig rig;

    // Fitted 1 cm away, competitor 0.5 cm away => competitor wins, drop.
    auto fitted = rig.segment_on_face(0, {Point(0, -1.0 * units::cm, 0)});
    auto competitor = rig.segment_on_face(0, {Point(0, -1.5 * units::cm, 0)});
    CHECK_FALSE(rig.cell_kept(fitted, {fitted, competitor}));

    // Fitted 0.2 cm away (inside the 0.3 cm floor) => kept even though the
    // competitor at 0.1 cm is closer.
    auto hugging = rig.segment_on_face(0, {Point(0, -1.8 * units::cm, 0)});
    auto closer = rig.segment_on_face(0, {Point(0, -1.9 * units::cm, 0)});
    CHECK(rig.cell_kept(hugging, {hugging, closer}));

    // No competitor at all: min_dis1_track stays 1e9 and any finite own
    // distance keeps the cell (the NeutrinoVertexFinder single-segment
    // local-graph no-op).
    auto lone = rig.segment_on_face(0, {Point(0, -1.0 * units::cm, 0)});
    CHECK(rig.cell_kept(lone, {lone}));
}

// doc pdvd/45 -- the drift-frame shift.  update_association builds the cell's
// test point in the RAW drift frame (time index -> x through the geometric
// offset_t; here slope_t = 1/cm so time 90 is x = 90 cm) while the segment
// clouds hold t0-CORRECTED fit points.  Model a cluster whose
// dirx * t0 * v_drift = +90 cm: both clouds sit at corrected x ~ 0, the cell at
// raw x = 90 cm.  The fitted segment is transversely 1 cm from the cell; the
// competitor sits 20 cm further along x and 0.5 cm transversely from the cell.
//   mixed frame : own = hypot(90, 1) = 90.0   competitor = hypot(70, 0.5) = 70.0  -> DROPPED
//                 (the segment whose cloud reaches farthest toward the raw cell
//                 wins -- the PDVD "one winner per cluster" symptom)
//   same frame  : own = 1.0                   competitor = hypot(20, 0.5) = 20.0  -> KEPT
// A shift for another volume, or a zero shift, must leave the legacy decision.
TEST_CASE("pdvd45 update_association drift-frame shift restores proximity arbitration")
{
    Rig rig;
    auto fitted = rig.segment_on_face(0, {Point(0, -1.0 * units::cm, 0)});
    auto far_competitor = rig.segment_on_face(0, {Point(20.0 * units::cm, -1.5 * units::cm, 0)});
    auto near_competitor = rig.segment_on_face(0, {Point(0, -1.5 * units::cm, 0)});

    auto kept = [&](std::shared_ptr<PR::Segment> seg,
                    const std::vector<std::shared_ptr<PR::Segment>>& all, int time) {
        TrackFitting::PlaneData ut, vt, wt;
        ut.associated_2d_points.insert(TrackFitting::Coord2D(0, 0, time, 2, 0, kUlayer));
        rig.tf.update_association(seg, all, ut, vt, wt);
        return ut.associated_2d_points.size() == 1;
    };

    // Legacy (empty shift map): mixed frame, far competitor wins, cell dropped.
    TrackFittingTestHarness::clear_excl_x_shift(rig.tf);
    CHECK_FALSE(kept(fitted, {fitted, far_competitor}, 90));

    // Shift = +90 cm for (apa 0, face 0): same frame, cell kept.
    TrackFittingTestHarness::set_excl_x_shift(rig.tf, 0, 0, 90.0 * units::cm);
    CHECK(kept(fitted, {fitted, far_competitor}, 90));

    // A shift registered only for face 1 must not move a face-0 cell.
    TrackFittingTestHarness::clear_excl_x_shift(rig.tf);
    TrackFittingTestHarness::set_excl_x_shift(rig.tf, 0, 1, 90.0 * units::cm);
    CHECK_FALSE(kept(fitted, {fitted, far_competitor}, 90));

    // Zero shift == legacy on a same-frame cell (time 0): near competitor at
    // 0.5 cm strips it, far competitor at 20 cm does not.
    TrackFittingTestHarness::clear_excl_x_shift(rig.tf);
    TrackFittingTestHarness::set_excl_x_shift(rig.tf, 0, 0, 0.0);
    CHECK_FALSE(kept(fitted, {fitted, near_competitor}, 0));
    CHECK(kept(fitted, {fitted, far_competitor}, 0));
    TrackFittingTestHarness::clear_excl_x_shift(rig.tf);
    CHECK_FALSE(kept(fitted, {fitted, near_competitor}, 0));
    CHECK(kept(fitted, {fitted, far_competitor}, 0));
}


// ---------------------------------------------------------------------------
// pr_ml doc 87b -- the exact competitor index (TrackFitting::ExclusionIndex).
// update_association must take the same keep/drop decision on every cell with
// the index as with the full competitor scan.  The index only filters which
// competitors are asked; each candidate is still asked with the unchanged
// exclusion distance and the unchanged '<=' rule.

namespace {

    std::vector<Point> random_points(std::mt19937& rng, int n, double cx, double cy, double cz, double spread)
    {
        std::uniform_real_distribution<double> u(-1.0, 1.0);
        std::vector<Point> pts;
        for (int i = 0; i < n; ++i) {
            pts.emplace_back(cx + spread * u(rng), cy + spread * u(rng), cz + spread * u(rng));
        }
        return pts;
    }

    // U / V cells: wire w <-> transverse w cm; W cells: wire w <-> z = w cm;
    // time t <-> drift x = t cm (the Rig's synthetic geometry).
    void fill_cells(TrackFitting::PlaneData& ut, TrackFitting::PlaneData& vt, TrackFitting::PlaneData& wt)
    {
        for (int face = 0; face < 2; ++face) {
            for (int t = 0; t <= 20; t += 2) {
                for (int w = -10; w <= 10; w += 2) {
                    ut.associated_2d_points.insert(TrackFitting::Coord2D(0, face, t, w, 0, kUlayer));
                    vt.associated_2d_points.insert(TrackFitting::Coord2D(0, face, t, w, 0, kVlayer));
                }
                for (int w = 0; w <= 20; w += 2) {
                    wt.associated_2d_points.insert(TrackFitting::Coord2D(0, face, t, w, 0, kWlayer));
                }
            }
        }
    }

    struct Kept {
        std::set<TrackFitting::Coord2D> u, v, w;
        // Coord2D defines only operator<: equal sets = same size and pairwise equivalent.
        static bool same(const std::set<TrackFitting::Coord2D>& a, const std::set<TrackFitting::Coord2D>& b)
        {
            return a.size() == b.size() &&
                   std::equal(a.begin(), a.end(), b.begin(), [](const TrackFitting::Coord2D& x, const TrackFitting::Coord2D& y) {
                       return !(x < y) && !(y < x);
                   });
        }
        bool operator==(const Kept& o) const { return same(u, o.u) && same(v, o.v) && same(w, o.w); }
        size_t size() const { return u.size() + v.size() + w.size(); }
    };

    Kept run(Rig& rig, std::shared_ptr<PR::Segment> seg, const std::vector<std::shared_ptr<PR::Segment>>& all,
             bool with_index)
    {
        TrackFitting::PlaneData ut, vt, wt;
        fill_cells(ut, vt, wt);
        if (with_index) TrackFittingTestHarness::excl_index_begin(rig.tf, all);
        rig.tf.update_association(seg, all, ut, vt, wt);
        if (with_index) TrackFittingTestHarness::excl_index_end(rig.tf);
        return Kept{ut.associated_2d_points, vt.associated_2d_points, wt.associated_2d_points};
    }

}  // namespace

TEST_CASE("pr_ml 87b exclusion index: same decisions as the full scan on random segments")
{
    Rig rig;
    std::mt19937 rng(8787);
    std::uniform_real_distribution<double> ux(0.0, 20.0), uy(-10.0, 10.0), uz(0.0, 20.0), us(0.2, 3.0);
    std::uniform_int_distribution<int> un(1, 30);

    std::vector<std::shared_ptr<PR::Segment>> all;
    for (int k = 0; k < 40; ++k) {
        const int kind = k % 10;
        const auto pts = random_points(rng, un(rng), ux(rng) * units::cm, uy(rng) * units::cm,
                                       uz(rng) * units::cm, us(rng) * units::cm);
        if (kind <= 6) {
            all.push_back(rig.segment_on_face(0, pts));             // a "fit" cloud on face 0
        }
        else if (kind == 7) {
            auto seg = rig.segment_on_face(0, {});                  // empty "fit" ...
            seg->dpcloud("main", make_face_cloud(0, pts));          // ... falls back to "main"
            all.push_back(seg);
        }
        else if (kind == 8) {
            all.push_back(rig.segment_on_face(1, pts));             // face 1 only
        }
        else {
            auto v1 = PR::make_vertex(*rig.graph);
            auto v2 = PR::make_vertex(*rig.graph);
            all.push_back(PR::make_segment(*rig.graph, v1, v2));    // no cloud at all
        }
    }

    size_t kept_total = 0, cells_total = 0;
    for (const auto& seg : all) {
        const Kept full = run(rig, seg, all, false);
        const Kept indexed = run(rig, seg, all, true);
        CHECK(full == indexed);
        kept_total += full.size();
        cells_total += 3 * 2 * 11 * 11;
    }
    // The index answered scans (not only the fallbacks), and the cases are
    // not degenerate (some cells kept, some dropped).
    CHECK(TrackFittingTestHarness::excl_index_used(rig.tf) > 0);
    CHECK(kept_total > 0);
    CHECK(kept_total < cells_total);
}

TEST_CASE("pr_ml 87b exclusion index: an exact tie still drops the cell")
{
    Rig rig;
    // The U cell (time 0, wire 2) sits at transverse 2 cm.  The fitted segment
    // is at 1 cm and the competitor at 3 cm: equal distances, and the keep rule
    // drops on '<='.
    auto fitted = rig.segment_on_face(0, {Point(0, -1.0 * units::cm, 0)});
    auto twin = rig.segment_on_face(0, {Point(0, -3.0 * units::cm, 0)});
    std::vector<std::shared_ptr<PR::Segment>> all{fitted, twin};

    auto kept = [&](bool with_index) {
        TrackFitting::PlaneData ut, vt, wt;
        ut.associated_2d_points.insert(TrackFitting::Coord2D(0, 0, 0, 2, 0, kUlayer));
        if (with_index) TrackFittingTestHarness::excl_index_begin(rig.tf, all);
        rig.tf.update_association(fitted, all, ut, vt, wt);
        if (with_index) TrackFittingTestHarness::excl_index_end(rig.tf);
        return ut.associated_2d_points.size() == 1;
    };
    const size_t used0 = TrackFittingTestHarness::excl_index_used(rig.tf);
    CHECK_FALSE(kept(false));
    CHECK_FALSE(kept(true));
    CHECK(TrackFittingTestHarness::excl_index_used(rig.tf) == used0 + 1);
}

TEST_CASE("pr_ml 87b exclusion index: clouds with different wpid_params fall back to the full scan")
{
    Rig rig;
    std::mt19937 rng(87);
    std::vector<std::shared_ptr<PR::Segment>> all;
    for (int k = 0; k < 6; ++k) {
        all.push_back(rig.segment_on_face(0, random_points(rng, 10, 5.0 * units::cm, 0.0, 10.0 * units::cm,
                                                           3.0 * units::cm)));
    }
    // One competitor whose cloud carries face-0 params only: one projection no
    // longer serves every segment, so the index must step aside.
    std::map<WirePlaneId, std::tuple<Facade::geo_point_t, double, double, double>> face0;
    face0[WirePlaneId(kAllLayers, 0, 0)] = std::make_tuple(Facade::geo_point_t(0, 0, 0), kAngleU, kAngleV, kAngleW);
    auto odd = std::make_shared<DynamicPointCloud>(face0);
    all.push_back(rig.segment_on_face(0, {Point(5.0 * units::cm, 0, 10.0 * units::cm)}));
    all.back()->dpcloud("fit", odd);
    all.back()->dpcloud("main", make_face_cloud(0, {Point(5.0 * units::cm, 0, 10.0 * units::cm)}));

    const size_t used0 = TrackFittingTestHarness::excl_index_used(rig.tf);
    for (const auto& seg : all) {
        CHECK(run(rig, seg, all, false) == run(rig, seg, all, true));
    }
    // The odd segment's "fit" cloud is empty, so its chosen cloud is "main"
    // (both_face_params) and the params agree: the index stays usable.  Give it
    // a non-empty odd cloud instead and the index must not be used.
    DPCBatch b;
    b.add_point(5.0 * units::cm, 0, 10.0 * units::cm, WirePlaneId(kAllLayers, 0, 0).ident(), nullptr, nullptr, {0, 0, 0},
                {0.0, 0.0, 0.0});
    const double angles[3] = {kAngleU, kAngleV, kAngleW};
    const WirePlaneLayer_t layers[3] = {kUlayer, kVlayer, kWlayer};
    for (int plane = 0; plane < 3; ++plane) {
        b.add_proj(5.0 * units::cm, cos(angles[plane]) * 10.0 * units::cm, WirePlaneId(layers[plane], 0, 0).ident());
        b.end_plane();
    }
    odd->add_points(std::move(b));
    const size_t used1 = TrackFittingTestHarness::excl_index_used(rig.tf);
    for (const auto& seg : all) {
        CHECK(run(rig, seg, all, false) == run(rig, seg, all, true));
    }
    CHECK(TrackFittingTestHarness::excl_index_used(rig.tf) == used1);
    CHECK(used1 >= used0);
}
