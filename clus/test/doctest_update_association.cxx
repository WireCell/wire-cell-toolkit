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
