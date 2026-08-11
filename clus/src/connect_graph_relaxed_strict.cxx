// doc pr/53 sec 16: "relaxed_strict" -- a stricter fork of
// connect_graph_relaxed (fork by full duplication, production file
// untouched), selected only by ClusteringProtectBundle via its graph_name
// config key.  It closes the sub-3cm arithmetic blind spot of the relaxed
// path test (doc pr/53 sec 6): the 1cm sampling loop's last sample is p2
// itself and good by construction, so num_bad <= num_steps-1 and the
// legacy "num_bad > 2" floor is unreachable for necks under ~3cm.
// Differences from connect_graph_relaxed, all inside the kill tests
// (relaxed_strict_bad below); every other line (candidate generation, MST
// construction, <3cm restore, edge emission) is a verbatim copy:
//   S1  ratio denominator becomes the interior step count (num_steps-1);
//   S2  the "> 2" floors become ">= 2" (sum-bad floors "> 2" likewise);
//   S3  S1+S2 applied to all three path-check blocks (closest-pair,
//       dir1, dir2), since their edges are tested and MST'd independently.
// Caps (>7 single-counter, >9 sum) and the num_bad[3] >= 3 dead-W veto are
// unchanged.
//
// doc pr/53 round 7 sec 18: image_check (default false, "relaxed_strict"
// behavior byte-for-byte unchanged) adds S5, an independent OR-kill next to
// relaxed_strict_bad in all three path-check blocks: relaxed_img_bad() on
// the longest contiguous run of "ghost" steps -- 1cm samples with no 3D
// image point (this cluster's own pt_clouds, any component) within
// m_img_radius and not dead-channel-excused on any plane.  Root cause:
// relaxed_strict_bad only ever sees three INDEPENDENT 2D per-plane
// projections (Facade_Grouping.cxx has_closest_point), never intersected in
// 3D, so three planes can each see charge from a different nearby track and
// the whole straight path reads "good" with nothing physically there.  The
// predicate, its run-length floor and its edge-length cap are justified in
// WireCellClus/Graphs.h and doc pr/53 round 7 sec 18.2 (offline scan against
// the owner's flagged edges + all round-6 mover edges) -- see there before
// changing the operating point.  "relaxed_strict_img"
// (make_graph_relaxed_strict_img) is the only caller passing image_check=true.
//
// doc pr/56 round 2: two_d_check (default false, "relaxed_strict_img"
// behavior byte-for-byte unchanged) adds S6, an independent additional
// OR-kill alongside relaxed_strict_bad/relaxed_img_bad in all three
// path-check blocks: two_d_connectivity_bad() on a bounded per-plane 2D
// flood-fill in (wire_index, time_slice) space between the two candidate
// components' own real fired-pixel footprints near the closest-approach
// region (Facade::Grouping::wire_charge_row/is_wire_dead), independent of
// the straight-3D-line assumption S1-S5 all share. Root cause this closes:
// even S5's per-1cm-sample "any 3D image point nearby" test can be
// satisfied at every sample by unrelated nearby charge without the two
// components' own footprints ever actually touching at the discrete
// wire/tick level -- see doc pr/56 (design) and doc pr/56 round 2 (this
// implementation) for the motivating gap (269774 c5-c10/c10-c11, both
// individually "legitimate" under S1-S5 yet bridging a gap a direct,
// correctly-killed candidate could not). Owner's kill rule: >= 1
// non-excused plane gap (stricter than S1-S5's implicit ">= 2 views"
// convention, by explicit instruction). Dead channels are folded directly
// into the flood-fill's passability predicate (a dead cell is always
// passable); U/V prolonged/isochronous gaps are excused by reusing the
// already-computed angle1/angle2 wire-parallel test, never re-detected.
// "relaxed_strict_img_2d" (make_graph_relaxed_strict_img_2d) is the only
// caller passing two_d_check=true.
//
// doc pr/56 round 3: repairs two round-2 defects and adds one relief
// valve, all still inside this function -- (1) the BFS adjacency is now a
// true Chebyshev box (dw wires x ds slice-steps), not round 2's
// accidental 4-connected/Von-Neumann stencil; (2) the slice-adjacency
// stride is now the minimum of the seed blobs' own
// slice_index_max()-slice_index_min() (= tick_span by construction), not
// round 2's inferred gap between sorted slice_index_min values, which
// measured wrong strides on a real event; (3) s6_dis_floor abstains S6
// entirely below a tuned 3D gap, where components are close enough that a
// 2D "gap" is almost certainly a coarse-seeding/BFS artifact rather than a
// real disconnection. Per-plane stencil half-widths and the floor are
// tuned from a 117-event scan against both round-2's kill targets and a
// round-3 rescue target (an isochronous track wrongly split by a
// U-plane-only gap, attributed to prolonged induction signal from a large
// EM shower) -- see doc pr/56 round 3 before changing any operating-point
// constant in the two_d_gap_kill lambda below.
//
// doc pr/57: adds an env-gated (WCT_OC56_SCAN_DUMP) per-edge JSONL dump,
// feeding sbnd_xin/overclustering_display -- a hand-scan tool for judging
// S6 removals as real gaps vs. induction-plane signal inefficiency. Unset
// (the default) => Oc56DumpWriter::instance().enabled() is false and every
// dump-only branch below is skipped; two_d_connectivity_bad, every shipped
// constant, and the verdict path are untouched. See that doc for the
// record schema.
//
// doc pr/58: when the dump IS on, also evaluate (and dump) candidates below
// s6_dis_floor instead of abstaining outright -- owner question: is the
// floor wrongly suppressing a real, unexcused W-plane gap near a vertex
// (evt 21073). With the dump alone this is diagnostic only: the function's
// return value (and every already-shipped S6-ON arm's kill set) still
// forces `false` below the floor exactly as before; see the `below_floor`
// local and `killed_ignore_floor` dump field below. floor_w_override (new
// function parameter, default false) is the real production-behavior
// counterpart -- owner-requested, evt 21073 j=6/k=7 as the validating case:
// when on, an unexcused W-plane-only gap kills even below the floor. Its
// own flavor, "relaxed_strict_img_2d_wfloor" (make_graphs.cxx /
// Facade_Cluster.cxx), is the only caller that passes true; every existing
// flavor including "relaxed_strict_img_2d" itself is untouched.
//
// doc pr/57 round 4: adds a third dump record type, "connectivity" -- one per
// graph call, written after the emit loop, carrying the final connected-
// component label of every starting component plus the list of component
// edges this call emitted (with endpoints).  The per-candidate `killed` bit
// alone cannot answer the question the hand scan is actually asking ("should
// these two components be separated"): one pair gets up to three independent
// candidates and can also stay joined through a third component.  Dump-only,
// same WCT_OC56_SCAN_DUMP gate; nothing here reads back into the graph.

#include "WireCellClus/Graphs.h"
#include "WireCellClus/IPCTransform.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Persist.h"

#include "connect_graphs.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <map>
#include <mutex>
#include <set>
#include <string>
#include <unordered_set>

// doc pr/57 round 6: 3x3 PCA for the S6 rescue's local-direction and
// component-extent features.  Eigen is already a WCT-wide dependency
// (WireCellUtil/Array.h); no new external dependency is introduced.
#include <Eigen/Dense>

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Graphs;
using namespace WireCell::Clus::Facade;

namespace {

    // doc pr/56 round 2 S6: one grid cell in a plane's (wire_index,
    // time_slice) space.
    struct S6Cell {
        int wind;
        int slice;
    };
    inline bool operator==(const S6Cell& a, const S6Cell& b) {
        return a.wind == b.wind && a.slice == b.slice;
    }
    struct S6CellHash {
        size_t operator()(const S6Cell& c) const {
            return std::hash<long long>()(
                (static_cast<long long>(c.wind) << 32) ^ static_cast<unsigned int>(c.slice));
        }
    };

    // A cell is passable if it carries a real hit (wire_charge_row has this
    // wire_index as a key at this exact time_slice -- ctpc_* rows are
    // unfiltered, so key presence IS the live predicate, doc pr/56 round 2)
    // or is dead-channel-excused (is_wire_dead, independent x-range check).
    // Folding dead-channel-excused into "passable" is deliberate: dead
    // channels are part of what "connected" means here, not a bolt-on rule.
    bool s6_cell_passable(const Facade::Grouping* grouping, int apa, int face, int plane,
                          int wind, int slice)
    {
        const auto* row = grouping->wire_charge_row(apa, face, plane, slice);
        if (row && row->count(wind)) return true;
        return grouping->is_wire_dead(apa, face, plane, wind, slice);
    }

    // doc pr/56 round 3: bounded bidirectional-seeded BFS between seed
    // cell-sets A and B in one plane's (wire_index, time_slice) grid.
    // slice_step is the real slice-adjacency stride (see caller -- NOT
    // assumed to be 1 tick, since wire_charge_row's keys are only populated
    // every tick_span ticks). Adjacency is a Chebyshev box: |d(wind)|<=dw
    // wires, |d(slice)|<=ds slice-steps, excluding the origin -- dw=ds=1 is
    // true 8-connectivity. (Round 2 shipped an unintentional 4-connected/
    // Von-Neumann stencil, {(+-1,0),(0,+-slice_step)} only -- see file
    // header.) Returns true if A and B are connected within the window and
    // cell budget; false (a gap) otherwise, INCLUDING when the budget is
    // exhausted before a connection is found -- fails closed, an
    // inconclusive result must not silently pass. If budget_hit is
    // non-null, *budget_hit reports whether the circuit breaker actually
    // fired, so callers can log "budget-exhausted-gap" distinctly from a
    // genuine gap rather than let it silently enter tuning data as one.
    bool s6_planes_connected(
        const Facade::Grouping* grouping, int apa, int face, int plane,
        const std::vector<S6Cell>& seeds_a, const std::vector<S6Cell>& seeds_b,
        int wind_lo, int wind_hi, int slice_lo, int slice_hi,
        int slice_step, int dw, int ds, size_t cell_budget, bool* budget_hit = nullptr)
    {
        if (budget_hit) *budget_hit = false;
        if (seeds_a.empty() || seeds_b.empty()) return false;

        std::unordered_set<S6Cell, S6CellHash> target(seeds_b.begin(), seeds_b.end());
        std::unordered_set<S6Cell, S6CellHash> visited;
        std::vector<S6Cell> frontier;
        for (const auto& s : seeds_a) {
            if (target.count(s)) return true;
            if (visited.insert(s).second) frontier.push_back(s);
        }

        std::vector<std::pair<int, int>> steps;
        steps.reserve(static_cast<size_t>(2 * dw + 1) * (2 * ds + 1) - 1);
        for (int a = -dw; a <= dw; ++a) {
            for (int b = -ds; b <= ds; ++b) {
                if (a == 0 && b == 0) continue;
                steps.emplace_back(a, b * slice_step);
            }
        }

        while (!frontier.empty()) {
            if (visited.size() > cell_budget) {
                if (budget_hit) *budget_hit = true;
                return false;  // circuit breaker: fails closed
            }
            std::vector<S6Cell> next;
            for (const auto& c : frontier) {
                for (const auto& d : steps) {
                    S6Cell nb{c.wind + d.first, c.slice + d.second};
                    if (nb.wind < wind_lo || nb.wind > wind_hi) continue;
                    if (nb.slice < slice_lo || nb.slice > slice_hi) continue;
                    if (visited.count(nb)) continue;
                    if (target.count(nb)) return true;
                    if (!s6_cell_passable(grouping, apa, face, plane, nb.wind, nb.slice)) continue;
                    visited.insert(nb);
                    next.push_back(nb);
                }
            }
            frontier.swap(next);
        }
        return false;
    }

    // doc pr/62 S7 ("long-edge corridor connectivity"): the seam with S6's
    // s6_dis_cap in the two_d_gap_kill lambda below (currently 30.0cm). MUST
    // equal that literal -- S6 is defined identically false at or above this
    // value, S7 (long_gap_kill below) is defined identically false below it,
    // so no candidate is ever evaluated by neither test or both. Mirrored as
    // Graphs::long_corridor_bad's dis_floor_cm default (same value, bare cm)
    // in the public header; doctest_long_corridor.cxx pins it. If this ever
    // changes, change s6_dis_cap and both defaults together.
    constexpr double s7_dis_floor_cm = 30.0;  // cm, NOT internal units

    // doc pr/62 S7: corridor geometry for the long-edge connectivity check.
    // A capsule -- a band of half-width hw around the straight segment from
    // the candidate's own near-point lattice cell to its far-point lattice
    // cell, plus an endpoint disk of radius cap at each end -- rather than
    // S6's open rectangle. Coordinates are normalized so one wire-step and
    // one slice-step (v = (slice - s_anchor)/slice_step) are both
    // approximately physical-mm-scale, so an isotropic half-width in these
    // units is approximately a physical half-width, and it is the same
    // metric the BFS's own Chebyshev box (dw wires x ds slice-steps)
    // already uses. Defined in the plane's own lattice, anchored on the
    // candidate's own two endpoint cells, rather than by projecting an
    // arbitrary 3D point: Grouping::convert_3Dpoint_time_ch exists but is a
    // different implementation of the wire/tick map than the one that
    // produced cluster.wire_index (Facade_Cluster.cxx, ultimately
    // pimpos->closest), and the two disagree by +-1 wire at tie boundaries
    // -- not something the kill verdict should depend on.
    struct S7Corridor {
        double u1 = 0, v1 = 0, u2 = 0, v2 = 0;  // endpoints, lattice units
        double hw = 3.0;                          // band half-width, lattice units
        double cap = 8.0;                         // endpoint disk radius, lattice units
        int slice_step = 1;
        int s_anchor = 0;                          // the tick value that maps to v=0

        // Perpendicular distance from (wind,slice) to the segment, and the
        // segment parameter 'along' in [0,1] (clamped) at the closest point
        // on the segment -- 'along' is a fraction of the corridor's own
        // length, not a cm distance; census code turns it into cm via
        // edge_dis.
        void geom(int wind, int slice, double& perp, double& along) const
        {
            const double u = wind;
            const double v = (slice - s_anchor) / static_cast<double>(std::max(1, slice_step));
            const double dx = u2 - u1, dy = v2 - v1;
            const double len2 = dx * dx + dy * dy;
            double t = 0.0;
            if (len2 > 1e-9) t = ((u - u1) * dx + (v - v1) * dy) / len2;
            const double tc = std::max(0.0, std::min(1.0, t));
            const double px = u1 + tc * dx, py = v1 + tc * dy;
            perp = std::sqrt((u - px) * (u - px) + (v - py) * (v - py));
            along = tc;
        }

        bool inside(int wind, int slice) const
        {
            double perp = 0, along = 0;
            geom(wind, slice, perp, along);
            if (perp <= hw) return true;
            // Endpoint disks catch seeds/cells sitting off-axis near either
            // terminus, which the band's perpendicular test alone would clip.
            const double u = wind;
            const double v = (slice - s_anchor) / static_cast<double>(std::max(1, slice_step));
            const double d1 = std::hypot(u - u1, v - v1);
            const double d2 = std::hypot(u - u2, v - v2);
            return d1 <= cap || d2 <= cap;
        }
    };

    // doc pr/62 S7: s6_cell_passable (:159-165) forked, not modified (M10 --
    // it has a live production consumer). The exact-key wire_charge_row
    // lookup there is safe at S6's <30cm range because S6's edge_dis-scaled
    // seed radius pulls in enough blobs that both sides span several
    // residues mod slice_step. S7 uses a FIXED 2cm seed radius (the change
    // that keeps the search O(D) instead of S6's O(D^2)), which collapses
    // each side toward a single blob -- an exact-key lookup would then read
    // a residue mismatch as "no charge" independent of whether charge is
    // actually there, on every long candidate. Probing a half-open window of
    // exactly slice_step ticks centred on the cell tiles the tick axis
    // without overlap and is residue-blind; it widens acceptance by at most
    // half a slice (~1.6mm at the SBND default tick_span). Dead-channel
    // excusal is unchanged (is_wire_dead is an x-range test, tick-exactness
    // irrelevant).
    bool s7_cell_passable(const Facade::Grouping* grouping, int apa, int face, int plane,
                          int wind, int slice, int slice_step)
    {
        const int lo = slice - slice_step / 2;
        for (int t = lo; t < lo + slice_step; ++t) {
            const auto* row = grouping->wire_charge_row(apa, face, plane, t);
            if (row && row->count(wind)) return true;
        }
        return grouping->is_wire_dead(apa, face, plane, wind, slice);
    }

    // doc pr/62 S7: s6_planes_connected (:182-229) forked, not modified
    // (M10). Three differences from it, each existing to make the search
    // O(D) instead of S6's O(D^2) and to make an inconclusive result
    // harmless at long range instead of a silent kill:
    //   (1) a cell must additionally satisfy corr.inside(wind,slice) -- the
    //       search cannot roam the open rectangle, only the corridor;
    //   (2) passability is s7_cell_passable (residue-blind, see above), not
    //       s6_cell_passable;
    //   (3) the circuit breaker FAILS OPEN. s6_planes_connected's
    //       fail-closed posture is correct at S6's <30cm range, where the
    //       window is small and exhaustion is pathological; at S7's
    //       30cm-and-up range an exhausted budget on a large sparse real
    //       shower must NOT be laundered into a kill -- the caller abstains
    //       via *budget_hit instead.
    // reach_along, if non-null, is set to the largest corridor-fraction
    // (corr.geom's 'along', in [0,1]) actually reached from seeds_a --
    // census-only bookkeeping; the verdict never reads it.
    bool s7_corridor_connected(
        const Facade::Grouping* grouping, int apa, int face, int plane,
        const std::vector<S6Cell>& seeds_a, const std::vector<S6Cell>& seeds_b,
        const S7Corridor& corr,
        int wind_lo, int wind_hi, int slice_lo, int slice_hi,
        int slice_step, int dw, int ds, size_t cell_budget,
        bool* budget_hit = nullptr, double* reach_along = nullptr)
    {
        if (budget_hit) *budget_hit = false;
        if (reach_along) *reach_along = 0.0;
        if (seeds_a.empty() || seeds_b.empty()) return false;

        std::unordered_set<S6Cell, S6CellHash> target(seeds_b.begin(), seeds_b.end());
        std::unordered_set<S6Cell, S6CellHash> visited;
        std::vector<S6Cell> frontier;
        for (const auto& s : seeds_a) {
            if (reach_along) {
                double p = 0, a = 0;
                corr.geom(s.wind, s.slice, p, a);
                if (a > *reach_along) *reach_along = a;
            }
            if (target.count(s)) return true;
            if (visited.insert(s).second) frontier.push_back(s);
        }

        std::vector<std::pair<int, int>> steps;
        steps.reserve(static_cast<size_t>(2 * dw + 1) * (2 * ds + 1) - 1);
        for (int a = -dw; a <= dw; ++a) {
            for (int b = -ds; b <= ds; ++b) {
                if (a == 0 && b == 0) continue;
                steps.emplace_back(a, b * slice_step);
            }
        }

        while (!frontier.empty()) {
            if (visited.size() > cell_budget) {
                if (budget_hit) *budget_hit = true;
                return true;  // fails OPEN: caller abstains via *budget_hit
            }
            std::vector<S6Cell> next;
            for (const auto& c : frontier) {
                for (const auto& d : steps) {
                    S6Cell nb{c.wind + d.first, c.slice + d.second};
                    if (nb.wind < wind_lo || nb.wind > wind_hi) continue;
                    if (nb.slice < slice_lo || nb.slice > slice_hi) continue;
                    if (visited.count(nb)) continue;
                    if (target.count(nb)) return true;
                    if (!corr.inside(nb.wind, nb.slice)) continue;
                    if (!s7_cell_passable(grouping, apa, face, plane, nb.wind, nb.slice, slice_step)) continue;
                    visited.insert(nb);
                    if (reach_along) {
                        double p = 0, a = 0;
                        corr.geom(nb.wind, nb.slice, p, a);
                        if (a > *reach_along) *reach_along = a;
                    }
                    next.push_back(nb);
                }
            }
            frontier.swap(next);
        }
        return false;
    }

    // doc pr/57: env-gated (WCT_OC56_SCAN_DUMP) singleton JSONL writer for
    // the overclustering-separation hand-scan display. One compact JSON
    // object per line, written under a mutex -- ClusteringProtectBundle
    // processes clusters concurrently, and doc pr/56 round 3 sec 8.3
    // already learned the hard way that a dump format needing cross-line
    // correlation by a locally-scoped index (like a per-cluster component
    // number) breaks silently when lines from different clusters
    // interleave. Here the correlation key ("graph_call", from
    // next_graph_call()) is a globally unique atomic counter instead, so a
    // component record and the edge records that reference it stay
    // correctly paired no matter how writes interleave across threads.
    class Oc56DumpWriter {
       public:
        static Oc56DumpWriter& instance() {
            static Oc56DumpWriter w;
            return w;
        }
        bool enabled() const { return m_enabled; }
        long next_graph_call() { return m_graph_call.fetch_add(1); }
        void write(const Json::Value& rec) {
            if (!m_enabled) return;
            const std::string line = Persist::dumps(rec, 0) + "\n";
            std::lock_guard<std::mutex> lk(m_mu);
            m_ofs << line;
        }

       private:
        Oc56DumpWriter() {
            const char* path = std::getenv("WCT_OC56_SCAN_DUMP");
            if (!path || !*path) {
                m_enabled = false;
                return;
            }
            m_ofs.open(path, std::ios::out | std::ios::trunc);
            m_enabled = m_ofs.is_open();
        }
        bool m_enabled = false;
        std::ofstream m_ofs;
        std::mutex m_mu;
        std::atomic<long> m_graph_call{0};
    };

}  // namespace

bool Graphs::relaxed_strict_bad(int nbad, int num_steps, int cap)
{
    // S1: the sampling loop tests fractions (ii+1)/num_steps, so the last
    // sample is the far endpoint p2 -- a cluster charge point, good by
    // construction.  Only num_steps-1 samples can possibly be bad.
    const int interior = num_steps - 1;
    // S2: floor ">= 2" (legacy "> 2"), cap unchanged.
    return nbad > cap || (nbad >= 2 && nbad >= 0.75 * interior);
}

bool Graphs::relaxed_img_bad(int max_ghost_run, double dis_cm, int run_floor, double dis_cap_cm)
{
    // doc pr/53 round 7 sec 18.2: operating point chosen from an offline
    // scan against all 27 round-6 mover events' emitted closest-pair edges,
    // not guessed -- see WireCellClus/Graphs.h for the full justification.
    return max_ghost_run >= run_floor && dis_cm < dis_cap_cm;
}

bool Graphs::two_d_connectivity_bad(bool gap_u, bool gap_v, bool gap_w,
                                     bool excuse_u, bool excuse_v)
{
    // doc pr/56 round 2 S6: owner's kill rule is >= 1 non-excused gap
    // plane. W is never excused -- the prolonged/isochronous exception is
    // an induction-plane-only physical effect (see WireCellClus/Graphs.h).
    return (gap_u && !excuse_u) || (gap_v && !excuse_v) || gap_w;
}

bool Graphs::two_d_rescue_ok(const S6RescueInput& in)
{
    // doc pr/57 round 6: line-by-line port of oc56_fit.py rescue() -- the
    // rule fitted against the owner's full hand scan.  Keep the two in
    // lockstep: any constant changed here must be re-fitted there first,
    // and the pair-level validation table (doc pr/57 sec 14) re-made.
    const bool gu = in.gap[0], gv = in.gap[1], gw = in.gap[2];
    if (in.dis_cm > 5.0) return false;
    if (!(gu || gv || gw)) return false;

    const bool coll = in.ab_local_deg >= 0 && in.ab_local_deg < 20.0;
    const bool collw = in.ab_local_deg >= 0 && in.ab_local_deg < 15.0;
    const bool dead_w = in.dead_w >= 3;
    const bool sub = in.lmin_cm > 4.0 && in.npmin >= 50;

    // Co-location: the two components' local footprints OVERLAP in the
    // (connected) W view plus at least one more connected plane -- two
    // views agree these are pieces of one object at the same place, not
    // end-to-end fragments.  W-anchored on purpose: the same overlap in
    // U+V alone is common at shower junctions the owner labels good.
    if (!gw && in.has_plane[2] && in.ov[2] >= 0.6 && in.npmin >= 50) {
        for (int p = 0; p < 2; ++p) {
            if (!in.gap[p] && in.has_plane[p] && in.ov[p] >= 0.6) return true;
        }
    }

    if (gw) {
        // W-robustness (owner case c): a W gap is only forgivable when it is
        // tiny (closes at a (3,3) stencil), well covered and the break is
        // collinear -- or when a dead-W band explains it (case a).
        const bool w_tiny = in.close_mx[2] <= 3 && in.cov[2] >= 0.8;
        const bool w_ok = (dead_w && in.npmin >= 20) || (w_tiny && collw && sub);
        if (!w_ok) return false;
    }

    // Every unexcused gapped induction plane must be explained.
    bool any_unexc = false;
    for (int p = 0; p < 2; ++p) {
        const bool g = in.gap[p];
        const bool e = (p == 0) ? in.excuse_u : in.excuse_v;
        if (!g || e) continue;
        any_unexc = true;
        const double c = in.cov[p];
        if (c < 0) return false;
        const double sl = std::max(in.slope[p], 0.0);
        const double ex = std::max(in.ext_med[p], 0.0);
        const int cl = in.close_mx[p];
        bool ok = false;
        // dead-W band distorts the induction response (owner case a); the
        // coverage floor is LOW because the distortion itself makes holes
        if (dead_w && in.npmin >= 20 && c >= 0.65) ok = true;
        // direction-consistent substantial pair over a covered gap
        if (sub && coll && c >= 0.90) ok = true;
        // ...even at moderate coverage when clearly track-like (W clean)
        if (!gw && in.npmin >= 55 && collw && c >= 0.55) ok = true;
        // prolonged-signal tiers (owner case b)
        if (in.lmin_cm > 5.5 && in.npmin >= 50 && (sl >= 5.0 || ex >= 50.0) && c >= 0.82) ok = true;
        if (in.npmin >= 15 && (sl >= 8.0 || ex >= 50.0) && c >= 0.82) ok = true;
        if (in.npmin >= 5 && collw && sl >= 15.0 && c >= 0.95) ok = true;
        // big pair, tight gap, small closure, full coverage (W clean)
        if (!gw && in.npmin >= 150 && in.dis_cm < 2.0 && cl <= 4 && c >= 0.90) ok = true;
        if (!ok) return false;
    }
    // Either every unexcused induction gap was explained, or there was none
    // (W side already vouched for, induction gaps all excused).
    (void)any_unexc;
    return true;
}

bool Graphs::two_d_w_track_ok(const S6WTrackInput& in)
{
    // doc pr/64 round 4: line-by-line port of wgate_sweep.py r2w_fires()/
    // r2d_fires() at the owner-selected tightened operating point.  Keep the
    // two in lockstep: any constant changed here must be re-scored there
    // against the 899-label truth table (docs/pr/pr57r6-truth.tsv) first.
    // All thresholds live in this one block on purpose.
    constexpr double L_cut_cm    = 6.0;   // R2w: min principal-axis extent
    constexpr int    N_cut       = 50;    // R2w: min component points
    constexpr double T_cut_cm    = 2.0;   // R2w: max transverse RMS
    constexpr double A_cut_deg   = 25.0;  // R2w: max global-axis angle
    constexpr double T_tight_cm  = 1.7;   // tightening (owner choice, evt122660)
    constexpr double A_tight_deg = 6.0;   // tightening (owner choice, evt122660)
    constexpr int    Nd_cut      = 20;    // R2d: min component points
    constexpr int    dead_w_min  = 3;     // R2d: min dead W wires across the gap
    constexpr double d_dead_cm   = 3.0;   // R2d: max candidate length

    // R2d: dead-W band across a W gap.  NOT sole-vote-gated -- the fitted
    // population (evt60669) has an induction plane also voting; the dump's
    // dead-W count is the same s6_rin.dead_w feature the offline fit used.
    if (in.w_gap && in.dead_w >= dead_w_min && in.npmin >= Nd_cut && in.dis_cm < d_dead_cm) {
        return true;
    }
    // R2w: W is the sole voting plane and both components are long, thin,
    // substantial and mutually collinear by their GLOBAL principal axes.
    if (!in.w_sole_vote) return false;
    if (!(in.lmin_cm > L_cut_cm) || in.npmin < N_cut) return false;
    if (!(in.tmax_cm < T_cut_cm) || !(in.ab_global_deg < A_cut_deg)) return false;
    if (!(in.tmax_cm < T_tight_cm || in.ab_global_deg < A_tight_deg)) return false;
    return true;
}

bool Graphs::long_corridor_bad(const S7CorridorInput& in, int min_gapped_planes,
                               double dis_floor_cm, double gap_floor_cm)
{
    // doc pr/62 S7: seam with S6 -- identically false below the floor. See
    // WireCellClus/Graphs.h and s7_dis_floor_cm above.
    if (in.dis_cm < dis_floor_cm) return false;
    const bool excused[3] = {in.excuse_u, in.excuse_v, false};  // W never excused
    int nvote = 0;
    for (int p = 0; p < 3; ++p) {
        if (!in.has_plane[p]) continue;    // no evidence, no verdict
        if (in.budget_hit[p]) continue;    // inconclusive: S7 fails OPEN, so abstain
        if (!in.gap[p]) continue;
        if (excused[p]) continue;
        if (in.gap_cm[p] < gap_floor_cm) continue;   // inactive at the default 0
        ++nvote;
    }
    return nvote >= min_gapped_planes;
}

void Graphs::connect_graph_relaxed_strict(
    const Facade::Cluster& cluster,
    IDetectorVolumes::pointer dv,
    IPCTransformSet::pointer pcts,
    Weighted::Graph& graph,
    bool image_check,
    bool two_d_check,
    bool floor_w_override,
    bool two_d_rescue,
    bool long_check,
    int long_min_planes,
    bool w_track_excuse)
{
    const bool use_ctpc = true;
    const auto* grouping = cluster.grouping();
    // doc pr/53 round 7 sec 18: S5 operating point (see Graphs.h for the
    // radius justification -- matches the OFFLINE oc53_probe.Loader img_dis
    // radius used in the threshold scan). Only read when image_check.
    const double m_img_radius = 1.0 * units::cm;

    // Get all the wire plane IDs from the grouping
    const auto& wpids = grouping->wpids();

    // Key: {apa, face} pair.  Each apa/face has fixed U/V/W wire directions that
    // are the same regardless of which layer's WirePlaneId is used for the lookup.
    // Using a full WirePlaneId (which includes the layer) as the map key would cause
    // std::out_of_range when get_wireplaneid() returns a wpid whose layer differs
    // from those populated here.  Keying by {apa,face} avoids this entirely.
    using af_pair_t = std::pair<int,int>;
    std::map<af_pair_t, geo_point_t> af_U_dir;
    std::map<af_pair_t, geo_point_t> af_V_dir;
    std::map<af_pair_t, geo_point_t> af_W_dir;
    for (const auto& wpid : wpids) {
        int apa = wpid.apa();
        int face = wpid.face();
        af_pair_t af{apa, face};
        if (af_U_dir.count(af)) continue;  // already computed for this apa/face

        // Create canonical wpids for all three planes with this APA and face
        WirePlaneId wpid_u(kUlayer, face, apa);
        WirePlaneId wpid_v(kVlayer, face, apa);
        WirePlaneId wpid_w(kWlayer, face, apa);

        // Get wire directions for all planes
        Vector wire_dir_u = dv->wire_direction(wpid_u);
        Vector wire_dir_v = dv->wire_direction(wpid_v);
        Vector wire_dir_w = dv->wire_direction(wpid_w);

        // Calculate angles
        double angle_u = std::atan2(wire_dir_u.z(), wire_dir_u.y());
        double angle_v = std::atan2(wire_dir_v.z(), wire_dir_v.y());
        double angle_w = std::atan2(wire_dir_w.z(), wire_dir_w.y());

        af_U_dir[af] = geo_point_t(0, cos(angle_u), sin(angle_u));
        af_V_dir[af] = geo_point_t(0, cos(angle_v), sin(angle_v));
        af_W_dir[af] = geo_point_t(0, cos(angle_w), sin(angle_w));
    }

    // this drift direction is only used to calculate isochronous case, so this is OK ...
    const geo_vector_t drift_dir_abs(1, 0, 0);


    // Form connected components
    std::vector<int> component(num_vertices(graph));
    const size_t num = connected_components(graph, &component[0]);

    if (num <= 1) return;

    // doc pr/53 round 6 F0 (relaxed_edge_census): log-only diagnostic, default
    // OFF.  Env WCT_RELAXED_EDGE_CENSUS unset => no log lines and no behavior
    // change of any kind; the graph construction below is untouched either way.
    static const bool oc53_census = std::getenv("WCT_RELAXED_EDGE_CENSUS") != nullptr;
    static auto oc53_log = WireCell::Log::logger("clus");

    // doc pr/57: per-cluster-call dump state. graph_call is a globally
    // unique id for THIS call (one per cluster's candidate graph); the
    // component-index dedup set is local to it, so re-dumping the same
    // component's full point cloud once per edge that touches it is
    // avoided without needing any cross-call correlation. No-op (and no
    // atomic increment) when the dump is disabled.
    auto& oc56_dump = Oc56DumpWriter::instance();
    const bool oc56_dump_on = oc56_dump.enabled();
    const long oc56_graph_call = oc56_dump_on ? oc56_dump.next_graph_call() : -1;
    std::set<size_t> oc56_dumped_components;

    // Allocate exactly num point clouds (one per component)
    std::vector<std::shared_ptr<Simple3DPointCloud>> pt_clouds(num);
    std::vector<std::vector<size_t>> pt_clouds_global_indices(num);
    for (size_t c = 0; c < num; ++c) {
        pt_clouds[c] = std::make_shared<Simple3DPointCloud>();
    }

    const auto& points = cluster.points();
    for (size_t i = 0; i < component.size(); ++i) {
        size_t c = component[i];
        pt_clouds[c]->add({points[0][i], points[1][i], points[2][i]});
        pt_clouds_global_indices[c].push_back(i);
    }

    // doc pr/57 round 6: per-component (point count, principal-axis extent)
    // cache for the S6 rescue -- computed at most once per component per
    // call, only when two_d_rescue is on and a killed candidate touches the
    // component.  Keyed by component index (never iterated), so fully
    // deterministic.  npmin is capped at 20000 to match the Python fit,
    // whose component clouds were the (<=20000-point strided) dump records;
    // every rescue threshold is far below the cap so the two agree exactly.
    // doc pr/64 round 4: widened to also keep the top eigenvector and the
    // transverse RMS -- both already computed by the same eigen-decomposition
    // and previously discarded.  trms uses the (n-1) normalization to match
    // the Python fit's np.cov (oc56_autoscan.transverse_rms).
    struct S6CompStat {
        int npoints = 0;
        double extent_cm = 0;
        Eigen::Vector3d axis = Eigen::Vector3d::Zero();  ///< zero => unmeasured
        double trms_cm = -1;                             ///< -1 => unmeasured
    };
    std::map<size_t, S6CompStat> s6_comp_stats;

    // doc pr/57 round 6: component pairs whose dir1/dir2 candidate the S6
    // rescue explicitly kept.  Needed at the emit loop: the dir-MST records
    // the pair's CLOSEST-candidate tuple (process_mst_deterministically), so
    // a killed closest sets the recorded tuple to -1 and the pair silently
    // loses every dir emission -- a rescued dir bridge would evaporate
    // (evt286180 0-1 was the labelled case).  Lookup-only, never iterated.
    std::set<std::pair<size_t, size_t>> s6_rescued_dir_pairs;
    auto s6_comp_stat = [&](size_t comp) -> S6CompStat {
        auto it = s6_comp_stats.find(comp);
        if (it != s6_comp_stats.end()) return it->second;
        const auto& cloud = pt_clouds.at(comp);
        const size_t n = cloud->get_num_points();
        S6CompStat st;
        st.npoints = static_cast<int>(std::min<size_t>(n, 20000));
        if (n >= 2) {
            Eigen::Vector3d mean = Eigen::Vector3d::Zero();
            for (size_t i = 0; i < n; ++i) {
                const auto p = cloud->point(i);
                mean += Eigen::Vector3d(p.x(), p.y(), p.z());
            }
            mean /= double(n);
            Eigen::Matrix3d C = Eigen::Matrix3d::Zero();
            for (size_t i = 0; i < n; ++i) {
                const auto p = cloud->point(i);
                const Eigen::Vector3d d = Eigen::Vector3d(p.x(), p.y(), p.z()) - mean;
                C += d * d.transpose();
            }
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(C);
            const Eigen::Vector3d axis = es.eigenvectors().col(2);
            double pmin = 1e30, pmax = -1e30;
            for (size_t i = 0; i < n; ++i) {
                const auto p = cloud->point(i);
                const double pr = (Eigen::Vector3d(p.x(), p.y(), p.z()) - mean).dot(axis);
                if (pr < pmin) pmin = pr;
                if (pr > pmax) pmax = pr;
            }
            st.extent_cm = (pmax - pmin) / units::cm;
            st.axis = axis;
            if (n >= 3) {
                // eigenvalues ascending: cols 0,1 are the transverse pair.
                // C is an unnormalized outer-product sum; np.cov divides by
                // (n-1), so match that here.
                const double t2 = (std::max(es.eigenvalues()(0), 0.0) +
                                   std::max(es.eigenvalues()(1), 0.0)) /
                                  double(n - 1);
                st.trms_cm = std::sqrt(t2) / units::cm;
            }
        }
        s6_comp_stats[comp] = st;
        return st;
    };

    // doc pr/57 round 6: local principal axis of a component within 6 cm
    // (widened once to 12 cm when starved) of a break point -- the rescue's
    // direction-consistency feature.  Mirrors oc56_fit.py local_dir().
    auto s6_local_axis = [&](const std::shared_ptr<Simple3DPointCloud>& cloud,
                              const geo_point_t& q) -> std::pair<bool, Eigen::Vector3d> {
        for (const double rr : {6.0 * units::cm, 12.0 * units::cm}) {
            const auto pts = cloud->get_closest_wcpoints_radius(q, rr);
            if (pts.size() < 5) continue;
            Eigen::Vector3d mean = Eigen::Vector3d::Zero();
            for (const auto& pr : pts) {
                mean += Eigen::Vector3d(pr.second.x(), pr.second.y(), pr.second.z());
            }
            mean /= double(pts.size());
            Eigen::Matrix3d C = Eigen::Matrix3d::Zero();
            for (const auto& pr : pts) {
                const Eigen::Vector3d d =
                    Eigen::Vector3d(pr.second.x(), pr.second.y(), pr.second.z()) - mean;
                C += d * d.transpose();
            }
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(C);
            return {true, es.eigenvectors().col(2)};
        }
        return {false, Eigen::Vector3d::Zero()};
    };

    if (oc53_census) {
        oc53_log->debug("OC53CENSUS-S cluster nblobs={} npoints={} ncomp={}",
                        cluster.nchildren(), cluster.npoints(), num);
        for (size_t c = 0; c < num; ++c) {
            double lo[3] = {1e12, 1e12, 1e12}, hi[3] = {-1e12, -1e12, -1e12};
            const size_t np = pt_clouds[c]->get_num_points();
            for (size_t i = 0; i < np; ++i) {
                const auto p = pt_clouds[c]->point(i);
                const double v[3] = {p.x(), p.y(), p.z()};
                for (int d = 0; d < 3; ++d) {
                    if (v[d] < lo[d]) lo[d] = v[d];
                    if (v[d] > hi[d]) hi[d] = v[d];
                }
            }
            oc53_log->debug("OC53CENSUS-S comp c={} npts={} bbox=({:.1f},{:.1f},{:.1f})-({:.1f},{:.1f},{:.1f})cm",
                            c, np, lo[0]/units::cm, lo[1]/units::cm, lo[2]/units::cm,
                            hi[0]/units::cm, hi[1]/units::cm, hi[2]/units::cm);
        }
    }

    // Initialize distance metrics — all sentinels (-1,-1,1e9) mean "no valid connection".
    // Use direct construction to avoid a redundant zero-fill followed by an overwrite pass (C.3).
    const auto sentinel = std::make_tuple(-1, -1, 1e9);
    const std::vector<std::tuple<int,int,double>> sentinel_row(num, sentinel);
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis(num, sentinel_row);
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_mst(num, sentinel_row);
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir1(num, sentinel_row);
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir2(num, sentinel_row);
    std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir_mst(num, sentinel_row);

    // Hoist scope-transform and cluster_t0 out of all per-step CTPC loops
    const bool needs_transform = (cluster.get_default_scope().hash() != cluster.get_raw_scope().hash());
    const auto ctpc_transform = needs_transform ? pcts->pc_transform(cluster.get_scope_transform()) : nullptr;
    const double cluster_t0 = needs_transform ? cluster.get_cluster_t0() : 0.0;

    // doc pr/53 round 7 sec 18 S5: only evaluated when image_check.  Re-walks
    // the same 1cm samples as the nb/nb1 loop below (independent pass, kept
    // separate from the already-shipped, extensively-validated nb/nb1
    // accumulation so this addition cannot perturb it).  A step is "ghost"
    // when it has no 3D image point (any of THIS cluster's own closely-
    // components -- in-family precedent connect_graph_relaxed.cxx:1062-1069,
    // relaxed_pid's cross-component check) within m_img_radius and is not
    // dead-channel-excused on any plane.  Returns the longest contiguous run.
    auto max_ghost_run = [&](const geo_point_t& p1, const geo_point_t& p2, int num_steps,
                              const WirePlaneId& wpid_p1, const WirePlaneId& wpid_p2) -> int {
        int run = 0, best = 0;
        for (int ii = 0; ii != num_steps; ii++) {
            geo_point_t test_p(
                p1.x() + (p2.x() - p1.x())/num_steps*(ii + 1),
                p1.y() + (p2.y() - p1.y())/num_steps*(ii + 1),
                p1.z() + (p2.z() - p1.z())/num_steps*(ii + 1)
            );
            bool ghost;
            auto test_wpid = get_wireplaneid(test_p, wpid_p1, wpid_p2, dv);
            if (test_wpid.apa() != -1) {
                geo_point_t test_p_raw = test_p;
                if (needs_transform) {
                    test_p_raw = ctpc_transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                }
                auto scores = grouping->test_good_point(test_p_raw, test_wpid.apa(), test_wpid.face());
                double img_dis = 1e9;
                for (size_t q = 0; q != num; ++q) {
                    img_dis = std::min(img_dis, pt_clouds.at(q)->get_closest_dis(test_p));
                }
                const bool dead_excused = (scores[3] + scores[4] + scores[5]) > 0;
                ghost = (img_dis > m_img_radius) && !dead_excused;
            } else {
                // Outside all APA volumes: no dead-channel excuse is
                // available, same convention as the nb/nb1 loop's apa==-1
                // branch (unconditionally counted bad there too).
                ghost = true;
            }
            run = ghost ? run + 1 : 0;
            if (run > best) best = run;
        }
        return best;
    };

    // doc pr/56 round 2 S6 (repaired + retuned round 3): only evaluated
    // when two_d_check (default false, no-op).  Per-plane 2D wind/tick
    // connectivity between components A and B's own real fired-pixel
    // footprints near the candidate closest-approach region -- independent
    // of the straight-3D-line assumption S1-S5 all share.  See file header
    // and WireCellClus/Graphs.h (two_d_connectivity_bad) for the design.
    //
    // doc pr/56 round 3: two defects fixed, one relief valve added, all
    // inside this lambda + s6_planes_connected above --
    //   - the BFS adjacency is now a true Chebyshev box (dw wires x ds
    //     slice-steps), not round 2's accidental 4-connected stencil
    //     ({(+-1,0),(0,+-slice_step)} only);
    //   - slice_step is now the MINIMUM of the seed blobs' own
    //     slice_index_max()-slice_index_min() (= tick_span by
    //     construction, aux/src/SamplingHelpers.cxx:90-91), not round 2's
    //     "first gap between sorted distinct slice_index_min values",
    //     which measured wrong strides (12..204 instead of 4) whenever the
    //     seed region spanned non-adjacent blobs;
    //   - s6_dis_floor: below this 3D gap the two components are close
    //     enough that a 2D "gap" almost certainly reflects a
    //     coarse-seeding/BFS artifact rather than a real disconnection;
    //     S6 abstains entirely (same "no evidence, no verdict" posture as
    //     the APA/face and empty-seed abstentions below).
    // Per-plane stencil half-widths, the window's fixed excess reach and
    // the distance floor are all tuned from the round-3 117-event scan
    // (doc pr/56 round 3 sec 8) against the round-2/round-3 named target
    // edges -- see that doc before changing any of these constants.
    auto two_d_gap_kill = [&](size_t j, size_t k, const char* blk,
                               const std::shared_ptr<Simple3DPointCloud>& cloud_a,
                               const std::vector<size_t>& gidx_a,
                               const std::shared_ptr<Simple3DPointCloud>& cloud_b,
                               const std::vector<size_t>& gidx_b,
                               const geo_point_t& q1, const WirePlaneId& wpid_q1,
                               const geo_point_t& q2, const WirePlaneId& wpid_q2,
                               double edge_dis, double a1, double a2) -> bool {
        if (!two_d_check) return false;

        // Distance cap: keeps the search window small; long candidates are
        // already governed by S1-S5.  Not a physics claim, just a budget.
        const double s6_dis_cap = 30.0 * units::cm;
        if (edge_dis >= s6_dis_cap) return false;

        // doc pr/56 round 3: below this 3D gap, abstain outright.  Tuned
        // from the 117-event scan (doc pr/56 round 3 sec 8): separates the
        // 71372 j=12/k=13 rescue target (0.76cm) from both 269774 kill
        // targets (2.00cm, 2.69cm) with margin on both sides, and on the
        // full sample cuts the 0-1cm S6-evaluated-edge kill rate from 53%
        // (repaired baseline, no floor) to 0% while leaving every bucket
        // at >=2cm untouched (a floor only ever REMOVES kills below it; it
        // cannot change a verdict at or above it).
        const double s6_dis_floor = 1.0 * units::cm;
        const bool below_floor = edge_dis < s6_dis_floor;
        // doc pr/58: keep evaluating below the floor whenever the scan dump
        // is on (diagnostic only -- see the `killed` line below, which
        // still forces `false` in that case) OR floor_w_override is on
        // (a real production behavior change, see connect_graphs.h). Every
        // already-shipped S6-ON arm ("relaxed_strict_img_2d") passes
        // floor_w_override=false and is byte-for-byte unaffected either way.
        if (below_floor && !oc56_dump_on && !floor_w_override) return false;

        // Only meaningful within a single APA/face -- a candidate that
        // crosses TPCs cannot be shown connected by a 2D wire/tick check.
        // Abstain (S1-S5 already carry their own APA-boundary handling for
        // this same pair) rather than force a verdict from no evidence.
        if (wpid_q1.apa() < 0 || wpid_q1.apa() != wpid_q2.apa() || wpid_q1.face() != wpid_q2.face()) {
            return false;
        }
        const int apa = wpid_q1.apa();
        const int face = wpid_q1.face();

        constexpr double s6_pi = 3.141592653589793;
        constexpr double s6_wire_angle_tol = 12.5 / 180.0 * s6_pi;
        const double seed_radius = std::max(edge_dis, 1.0 * units::cm) + 2.0 * units::cm;

        auto seeds_a_wc = cloud_a->get_closest_wcpoints_radius(q1, seed_radius);
        auto seeds_b_wc = cloud_b->get_closest_wcpoints_radius(q2, seed_radius);
        if (seeds_a_wc.empty() || seeds_b_wc.empty()) return false;  // no data to seed with; abstain

        // doc pr/56 round 3: slice-adjacency stride, from each seed blob's
        // OWN span, not round 2's inferred gap across the whole seed set
        // (see comment above -- the round-2 approach measured 11 of 97
        // edges with a wrong stride on evt 71372).
        auto blob_stride = [&](const std::vector<std::pair<size_t, geo_point_t>>& wc,
                                const std::vector<size_t>& gidx, int& best) {
            for (const auto& p : wc) {
                const auto* blob = cluster.blob_with_point(gidx.at(p.first));
                if (!blob) continue;
                const int span = blob->slice_index_max() - blob->slice_index_min();
                if (span > 0 && (best <= 0 || span < best)) best = span;
            }
        };
        int slice_step = -1;
        blob_stride(seeds_a_wc, gidx_a, slice_step);
        blob_stride(seeds_b_wc, gidx_b, slice_step);
        if (slice_step <= 0) slice_step = 4;  // SBND MaskSlices tick_span default -- fallback only

        // doc pr/56 round 3: per-plane stencil half-widths.  (1,1) is true
        // 8-connectivity (the repaired round-2 baseline).  U/V carry the
        // isochronous/prolonged-shower relaxation (owner report: a real
        // track broken by a U-plane-only gap, attributed to prolonged
        // induction signal from a large EM shower, not prolonged path
        // *topology* -- so the existing angle1/angle2 excusal below cannot
        // fire for it).  W stays at the repaired baseline: it is never
        // excused, and no case has yet justified relaxing it further.
        // doc pr/56 round 3 sec 8: adjacency widening (dw,ds > 1) was
        // tried and REJECTED -- (2,2) on U/V rescues the 71372 j=12/k=13
        // target (0.76cm) but ALSO rescues 269774 j=10/k=11 (2.69cm, a
        // MUST-STAY-KILLED control): both are U-plane-only unexcused gaps,
        // and pure adjacency widening cannot tell them apart by distance,
        // only by wire/tick proximity, which both happen to share. This is
        // exactly the owner-approved fallback case ("if widening U
        // adjacency alone can't separate them, use a 3D-gap floor
        // instead") -- so adjacency stays at the repaired (1,1) baseline
        // on every plane, and s6_dis_floor below does the separating
        // work, since it discriminates by the one axis (3D distance) on
        // which the two controls actually differ.
        const int s6_dw_uv = 1;
        const int s6_ds_uv = 1;
        const int s6_dw_w = 1;
        const int s6_ds_w = 1;
        // Fixed excess window reach beyond any operating-point stencil, so
        // the window is IDENTICAL for every (dw,ds) up to this reach --
        // connectivity is then monotone in (dw,ds), and the round-3 scan
        // matrix's entry at the shipped operating point equals the
        // verdict actually used, by construction.
        const int s6_reach = 4;
        const int margin_w = 3 + s6_reach;    // wires
        const int margin_s = 2 + s6_reach;    // slice-steps
        const size_t cell_budget = 20000;     // larger window+stencil than round 2's 4000

        bool gap[3] = {false, false, false};
        bool budget_hit[3] = {false, false, false};
        std::string matrix[3];  // doc pr/56 round 3 census + doc pr/57 dump, see below
        Json::Value oc56_dump_planes(Json::arrayValue);  // doc pr/57, dump only
        // doc pr/57 round 6: per-plane seed/window snapshot kept for the S6
        // rescue's feature measurements.  Filled only when two_d_rescue is
        // on; a plane skipped for missing seeds stays has=false and never
        // votes in the rescue either.
        struct S6PlaneData {
            bool has = false;
            std::vector<S6Cell> sa, sb;
            int wlo = 0, whi = 0, slo = 0, shi = 0;
        };
        S6PlaneData s6_pdata[3];
        for (int plane = 0; plane != 3; ++plane) {
            std::vector<S6Cell> seeds_a2d, seeds_b2d;
            int wlo = 1 << 30, whi = -(1 << 30), slo = 1 << 30, shi = -(1 << 30);
            auto collect = [&](const std::vector<std::pair<size_t, geo_point_t>>& wc,
                                const std::vector<size_t>& gidx, std::vector<S6Cell>& out) {
                for (const auto& p : wc) {
                    const size_t gi = gidx.at(p.first);
                    const auto* blob = cluster.blob_with_point(gi);
                    if (!blob) continue;
                    const int wind = cluster.wire_index(gi, plane);
                    const int slice = blob->slice_index_min();
                    out.push_back({wind, slice});
                    wlo = std::min(wlo, wind); whi = std::max(whi, wind);
                    slo = std::min(slo, slice); shi = std::max(shi, slice);
                }
            };
            collect(seeds_a_wc, gidx_a, seeds_a2d);
            collect(seeds_b_wc, gidx_b, seeds_b2d);
            if (seeds_a2d.empty() || seeds_b2d.empty()) continue;  // no data this plane; don't vote gap

            const int win_wlo = wlo - margin_w, win_whi = whi + margin_w;
            const int win_slo = slo - margin_s * slice_step, win_shi = shi + margin_s * slice_step;
            const int dw = (plane == 2) ? s6_dw_w : s6_dw_uv;
            const int ds = (plane == 2) ? s6_ds_w : s6_ds_uv;
            const bool connected = s6_planes_connected(
                grouping, apa, face, plane, seeds_a2d, seeds_b2d,
                win_wlo, win_whi, win_slo, win_shi,
                slice_step, dw, ds, cell_budget, &budget_hit[plane]);
            gap[plane] = !connected;

            if (two_d_rescue) {
                s6_pdata[plane] = {true, seeds_a2d, seeds_b2d,
                                   win_wlo, win_whi, win_slo, win_shi};
            }

            // doc pr/56 round 3 / doc pr/57: a full dw={1..4} x ds={1..4}
            // connectivity matrix at the SAME fixed window (s6_reach
            // above), so any candidate per-plane operating point is
            // evaluable offline from one run without rerunning
            // wire-cell -- either for the census log or for the scan
            // dump's near-miss badges.  Diagnostic only -- never affects
            // gap[]/killed.  matrix[i*4+j] is dw=i+1, ds=j+1 ('1'
            // connected, '0' gap).  Built into ONE log line / dump record
            // per edge (not one per plane): round-3's first attempt at
            // per-plane lines correlated by (blk,j,k) across sequential
            // log entries broke silently whenever ClusteringProtectBundle
            // processes clusters concurrently and their log lines
            // interleave -- a single self-contained record per edge has no
            // such correlation to get wrong.
            if (oc53_census || oc56_dump_on) {
                matrix[plane].reserve(16);
                for (int mdw = 1; mdw <= 4; ++mdw) {
                    for (int mds = 1; mds <= 4; ++mds) {
                        const bool c = s6_planes_connected(
                            grouping, apa, face, plane, seeds_a2d, seeds_b2d,
                            win_wlo, win_whi, win_slo, win_shi,
                            slice_step, mdw, mds, cell_budget, nullptr);
                        matrix[plane].push_back(c ? '1' : '0');
                    }
                }
            }

            if (oc56_dump_on) {
                // doc pr/57: enumerate every INTEGER slice in the window,
                // not just multiples of slice_step from win_slo. The real
                // BFS above steps by slice_step multiples starting from
                // EACH SEED'S OWN slice value, so seeds on different
                // residue classes mod slice_step query different,
                // interleaved lattices -- the set of cells the BFS can
                // touch is a union across them, not one win_slo-anchored
                // lattice. Stepping by 1 sidesteps that: it is a superset
                // of every cell any BFS invocation (any seed, any dw/ds)
                // could touch inside this window, so the offline replay in
                // oc56_dump_check.py can query it exactly as the BFS would
                // and is guaranteed not to see a spurious "no data" gap
                // that is really just an unwalked residue in the dump.
                // Capped to bound pathological (large slice_step) windows;
                // the fallback stride actually used is recorded per-plane
                // (native_step) so the checker knows when a native-1 dump
                // was not attempted rather than silently trusting it.
                // "fired" and "dead" are reported separately so the viewer
                // can shade them differently; dead is collapsed to
                // per-wire slice RANGES, not per-cell, or the payload
                // explodes.
                const long win_wires = static_cast<long>(win_whi - win_wlo) + 1;
                const long win_slices_native = static_cast<long>(win_shi - win_slo) + 1;
                const long oc56_dump_cap = 300000;
                const bool oc56_native = (win_wires * win_slices_native) <= oc56_dump_cap;
                const int dump_step = oc56_native ? 1 : slice_step;

                Json::Value fired(Json::arrayValue);
                Json::Value dead(Json::arrayValue);
                for (int wind = win_wlo; wind <= win_whi; ++wind) {
                    int run_lo = -1;
                    for (int slice = win_slo; slice <= win_shi; slice += dump_step) {
                        const auto* row = grouping->wire_charge_row(apa, face, plane, slice);
                        if (row && row->count(wind)) {
                            Json::Value c(Json::arrayValue);
                            c.append(wind);
                            c.append(slice);
                            fired.append(c);
                        }
                        const bool is_dead = grouping->is_wire_dead(apa, face, plane, wind, slice);
                        if (is_dead && run_lo < 0) run_lo = slice;
                        if (!is_dead && run_lo >= 0) {
                            Json::Value r(Json::arrayValue);
                            r.append(wind);
                            r.append(run_lo);
                            r.append(slice - dump_step);
                            dead.append(r);
                            run_lo = -1;
                        }
                    }
                    if (run_lo >= 0) {
                        Json::Value r(Json::arrayValue);
                        r.append(wind);
                        r.append(run_lo);
                        r.append(win_shi);
                        dead.append(r);
                        run_lo = -1;
                    }
                }
                Json::Value seeds_a_j(Json::arrayValue), seeds_b_j(Json::arrayValue);
                for (const auto& c : seeds_a2d) {
                    Json::Value v(Json::arrayValue);
                    v.append(c.wind);
                    v.append(c.slice);
                    seeds_a_j.append(v);
                }
                for (const auto& c : seeds_b2d) {
                    Json::Value v(Json::arrayValue);
                    v.append(c.wind);
                    v.append(c.slice);
                    seeds_b_j.append(v);
                }
                Json::Value win(Json::arrayValue);
                win.append(win_wlo);
                win.append(win_whi);
                win.append(win_slo);
                win.append(win_shi);
                Json::Value prec;
                prec["plane"] = plane;
                prec["win"] = win;
                prec["native_step"] = dump_step;
                prec["fired"] = fired;
                prec["dead"] = dead;
                prec["seeds_a"] = seeds_a_j;
                prec["seeds_b"] = seeds_b_j;
                oc56_dump_planes.append(prec);
            }
        }

        const bool excuse_u = a1 < s6_wire_angle_tol;
        const bool excuse_v = a2 < s6_wire_angle_tol;
        // doc pr/58: the hypothetical verdict from gap/excuse alone, BEFORE
        // the below_floor override -- dumped as `killed_ignore_floor` for
        // the scan display regardless of floor_w_override.
        const bool killed_if_evaluated = two_d_connectivity_bad(gap[0], gap[1], gap[2], excuse_u, excuse_v);
        // Below the floor: abstain (false) UNLESS floor_w_override is on
        // AND the gap is on W specifically -- W is never excused by
        // two_d_connectivity_bad, so gap[2] alone is exactly "unexcused W
        // gap". A U/V-only gap still can't kill below the floor even with
        // the knob on: the floor exists to guard against noisy near-floor
        // induction-plane false positives, and this knob only trusts W.
        // (doc pr/57 round 6: non-const because two_d_rescue below may
        // un-kill; with two_d_rescue false the value never changes.)
        bool killed = below_floor ? (floor_w_override && gap[2]) : killed_if_evaluated;
        const bool s6_killed_pre = killed;

        // doc pr/57 round 6: the S6 rescue.  Only ever runs on candidates
        // S6 (or the floor override) actually killed, and only within the
        // display population (<= 5 cm), so its cost and its reach are both
        // bounded.  Features mirror oc56_fit.py edge_features(); the
        // decision itself is the pure Graphs::two_d_rescue_ok().
        bool s6_rescued = false;
        Graphs::S6RescueInput s6_rin;
        bool s6_rin_valid = false;
        if (two_d_rescue && killed && edge_dis <= 5.0 * units::cm) {
            s6_rin.dis_cm = edge_dis / units::cm;
            for (int p = 0; p < 3; ++p) s6_rin.gap[p] = gap[p];
            s6_rin.excuse_u = excuse_u;
            s6_rin.excuse_v = excuse_v;
            const auto sta = s6_comp_stat(j);
            const auto stb = s6_comp_stat(k);
            s6_rin.npmin = std::min(sta.npoints, stb.npoints);
            s6_rin.lmin_cm = std::min(sta.extent_cm, stb.extent_cm);
            const auto axa = s6_local_axis(cloud_a, q1);
            const auto axb = s6_local_axis(cloud_b, q2);
            if (axa.first && axb.first) {
                const double c = std::min(1.0, std::abs(axa.second.dot(axb.second)));
                s6_rin.ab_local_deg = std::acos(c) * 180.0 / s6_pi;
            }
            for (int plane = 0; plane != 3; ++plane) {
                const auto& pd = s6_pdata[plane];
                s6_rin.has_plane[plane] = pd.has;
                if (!pd.has) continue;
                // smallest diagonal stencil that closes the plane; the
                // fixed-reach window makes connectivity monotone in
                // (dw,ds), so min over connected of max(dw,ds) is exactly
                // the smallest d with (d,d) connected -- 3 extra bounded
                // BFS runs at most, killed candidates only.
                int close_mx = 1;
                if (gap[plane]) {
                    close_mx = 5;
                    for (int d = 2; d <= 4; ++d) {
                        if (s6_planes_connected(grouping, apa, face, plane,
                                                pd.sa, pd.sb, pd.wlo, pd.whi,
                                                pd.slo, pd.shi, slice_step,
                                                d, d, cell_budget, nullptr)) {
                            close_mx = d;
                            break;
                        }
                    }
                }
                s6_rin.close_mx[plane] = close_mx;

                // seed spans, overlap and (wire,time)-view slope
                int alo = 1 << 30, ahi = -(1 << 30), blo2 = 1 << 30, bhi2 = -(1 << 30);
                int smin = 1 << 30, smax = -(1 << 30);
                for (const auto& cell : pd.sa) {
                    alo = std::min(alo, cell.wind); ahi = std::max(ahi, cell.wind);
                    smin = std::min(smin, cell.slice); smax = std::max(smax, cell.slice);
                }
                for (const auto& cell : pd.sb) {
                    blo2 = std::min(blo2, cell.wind); bhi2 = std::max(bhi2, cell.wind);
                    smin = std::min(smin, cell.slice); smax = std::max(smax, cell.slice);
                }
                const int ovw = std::min(ahi, bhi2) - std::max(alo, blo2) + 1;
                const int wminspan = std::min(ahi - alo, bhi2 - blo2) + 1;
                s6_rin.ov[plane] = std::max(0.0, double(ovw) / double(wminspan));
                if (pd.sa.size() + pd.sb.size() >= 3) {
                    const int dwspan = std::max(ahi, bhi2) - std::min(alo, blo2);
                    s6_rin.slope[plane] = (dwspan < 1)
                                              ? 999.0
                                              : double(smax - smin) / double(dwspan);
                }

                // coverage / per-wire fired extent / dead wires across the
                // combined seed span.  Same stride decision as the dump
                // (full-window size), so the Python fit and this replay see
                // identical lattices.
                const int lo = std::min(alo, blo2), hi = std::max(ahi, bhi2);
                const long win_wires2 = static_cast<long>(pd.whi - pd.wlo) + 1;
                const long win_slices2 = static_cast<long>(pd.shi - pd.slo) + 1;
                const int dstep = (win_wires2 * win_slices2 <= 300000) ? 1 : slice_step;
                int n_sig = 0, n_dead = 0;
                std::vector<int> exts;
                for (int wind = lo; wind <= hi; ++wind) {
                    int fmin = 1 << 30, fmax = -(1 << 30);
                    bool dead_any = false;
                    for (int slice = pd.slo; slice <= pd.shi; slice += dstep) {
                        const auto* row = grouping->wire_charge_row(apa, face, plane, slice);
                        if (row && row->count(wind)) {
                            fmin = std::min(fmin, slice);
                            fmax = std::max(fmax, slice);
                        }
                        if (!dead_any &&
                            grouping->is_wire_dead(apa, face, plane, wind, slice)) {
                            dead_any = true;
                        }
                    }
                    if (fmax >= fmin) {
                        ++n_sig;
                        exts.push_back(fmax - fmin);
                    }
                    if (dead_any) ++n_dead;
                }
                s6_rin.cov[plane] = double(n_sig) / double(hi - lo + 1);
                if (!exts.empty()) {
                    std::sort(exts.begin(), exts.end());
                    const size_t m = exts.size();
                    s6_rin.ext_med[plane] =
                        (m % 2) ? double(exts[m / 2])
                                : 0.5 * (exts[m / 2 - 1] + exts[m / 2]);
                }
                if (plane == 2) s6_rin.dead_w = n_dead;
            }
            s6_rin_valid = true;
            s6_rescued = Graphs::two_d_rescue_ok(s6_rin);
            if (s6_rescued) {
                killed = false;
                if (std::string(blk) != "closest") {
                    s6_rescued_dir_pairs.insert({std::min(j, k), std::max(j, k)});
                }
            }
        }

        // doc pr/64 round 4: the W-plane long-track exception.  Runs on ANY
        // candidate still killed after the rescue (no <= 5 cm cap -- the
        // exposed population is every killed S6 candidate whose sole voting
        // plane is W, e.g. 276836's 5.97 cm pair).  Features come from the
        // same memoized s6_comp_stat cache; the decision is the pure
        // Graphs::two_d_w_track_ok().  Revive-only, exactly like the rescue.
        bool s6_wtrack_revived = false;
        Graphs::S6WTrackInput s6_win;
        bool s6_win_valid = false;
        if (w_track_excuse && killed) {
            const bool vote_u = gap[0] && !excuse_u;
            const bool vote_v = gap[1] && !excuse_v;
            s6_win.w_gap = gap[2];
            s6_win.w_sole_vote = gap[2] && !vote_u && !vote_v;
            s6_win.dis_cm = edge_dis / units::cm;
            const auto sta = s6_comp_stat(j);
            const auto stb = s6_comp_stat(k);
            s6_win.npmin = std::min(sta.npoints, stb.npoints);
            s6_win.lmin_cm = std::min(sta.extent_cm, stb.extent_cm);
            if (sta.trms_cm >= 0 && stb.trms_cm >= 0) {
                s6_win.tmax_cm = std::max(sta.trms_cm, stb.trms_cm);
            }
            if (sta.npoints >= 2 && stb.npoints >= 2) {
                const double c = std::min(1.0, std::abs(sta.axis.dot(stb.axis)));
                s6_win.ab_global_deg = std::acos(c) * 180.0 / s6_pi;
            }
            // dead-W is only measured for the <= 5 cm rescue population,
            // which contains all of R2d's dis < 3 cm band.
            s6_win.dead_w = s6_rin_valid ? s6_rin.dead_w : -1;
            s6_win_valid = true;
            if (Graphs::two_d_w_track_ok(s6_win)) {
                s6_wtrack_revived = true;
                killed = false;
                if (std::string(blk) != "closest") {
                    // same emission repair as the rescue above (evt286180).
                    s6_rescued_dir_pairs.insert({std::min(j, k), std::max(j, k)});
                }
            }
        }

        if (oc53_census) {
            oc53_log->debug(
                "OC56CENSUS-2D edge blk={} j={} k={} dis={:.2f}cm apa={} face={} "
                "p1=({:.2f},{:.2f},{:.2f}) p2=({:.2f},{:.2f},{:.2f}) slice_step={} "
                "gap_u={} gap_v={} gap_w={} excuse_u={} excuse_v={} "
                "budget_u={} budget_v={} budget_w={} matrix_u={} matrix_v={} matrix_w={} killed={}",
                blk, j, k, edge_dis / units::cm, apa, face,
                q1.x() / units::cm, q1.y() / units::cm, q1.z() / units::cm,
                q2.x() / units::cm, q2.y() / units::cm, q2.z() / units::cm,
                slice_step, gap[0], gap[1], gap[2], excuse_u, excuse_v,
                budget_hit[0], budget_hit[1], budget_hit[2],
                matrix[0], matrix[1], matrix[2], killed);
        }

        if (oc56_dump_on) {
            // doc pr/57: dump each component's full point cloud once per
            // cluster-graph call (oc56_graph_call / oc56_dumped_components,
            // see above), keyed by (graph_call, comp) so repeated
            // references from other edges in the same call don't re-dump
            // it -- j/k are cluster-local and get reused across many edges.
            auto dump_component = [&](size_t comp, const std::shared_ptr<Simple3DPointCloud>& cloud) {
                if (!oc56_dumped_components.insert(comp).second) return;  // already dumped this call
                Json::Value rec;
                rec["type"] = "component";
                rec["graph_call"] = static_cast<int>(oc56_graph_call);
                rec["comp"] = static_cast<int>(comp);
                rec["cluster_id"] = cluster.get_cluster_id();
                const size_t n = cloud->get_num_points();
                const size_t cap = 20000;
                const size_t stride = (n > cap) ? (n + cap - 1) / cap : 1;
                Json::Value pts(Json::arrayValue);
                for (size_t i = 0; i < n; i += stride) {
                    const auto p = cloud->point(i);
                    Json::Value pt(Json::arrayValue);
                    pt.append(p.x() / units::cm);
                    pt.append(p.y() / units::cm);
                    pt.append(p.z() / units::cm);
                    pts.append(pt);
                }
                rec["points"] = pts;
                rec["npts"] = static_cast<int>(n);
                rec["truncated"] = (stride > 1);
                oc56_dump.write(rec);
            };
            dump_component(j, cloud_a);
            dump_component(k, cloud_b);

            Json::Value rec;
            rec["type"] = "edge";
            rec["graph_call"] = static_cast<int>(oc56_graph_call);
            rec["j"] = static_cast<int>(j);
            rec["k"] = static_cast<int>(k);
            rec["blk"] = blk;
            rec["dis"] = edge_dis / units::cm;
            rec["apa"] = apa;
            rec["face"] = face;
            Json::Value p1j(Json::arrayValue), p2j(Json::arrayValue);
            p1j.append(q1.x() / units::cm);
            p1j.append(q1.y() / units::cm);
            p1j.append(q1.z() / units::cm);
            p2j.append(q2.x() / units::cm);
            p2j.append(q2.y() / units::cm);
            p2j.append(q2.z() / units::cm);
            rec["p1"] = p1j;
            rec["p2"] = p2j;
            rec["slice_step"] = slice_step;
            Json::Value gapj(Json::arrayValue), excusej(Json::arrayValue),
                budgetj(Json::arrayValue), matrixj(Json::arrayValue);
            gapj.append(gap[0]); gapj.append(gap[1]); gapj.append(gap[2]);
            excusej.append(excuse_u); excusej.append(excuse_v);
            budgetj.append(budget_hit[0]); budgetj.append(budget_hit[1]); budgetj.append(budget_hit[2]);
            matrixj.append(matrix[0]); matrixj.append(matrix[1]); matrixj.append(matrix[2]);
            rec["gap"] = gapj;
            rec["excuse"] = excusej;
            rec["budget"] = budgetj;
            rec["matrix"] = matrixj;
            rec["killed"] = killed;
            // doc pr/58: `killed` above already reflects floor_w_override
            // when this dump was produced by the "relaxed_strict_img_2d_wfloor"
            // flavor -- these three fields let the scan display show the
            // floor-ignoring hypothetical either way, and record which mode
            // produced this record.
            rec["below_floor"] = below_floor;
            rec["killed_ignore_floor"] = killed_if_evaluated;
            rec["floor_w_override"] = floor_w_override;
            // doc pr/57 round 6: rescue provenance.  `killed` above is the
            // POST-rescue verdict actually returned; these three fields let
            // the offline checker verify Python-rule == C++-rule on every
            // candidate, feature by feature.
            if (two_d_rescue) {
                rec["two_d_rescue"] = true;
                rec["killed_pre_rescue"] = s6_killed_pre;
                rec["rescued"] = s6_rescued;
                if (s6_rin_valid) {
                    Json::Value v2;
                    v2["ab"] = s6_rin.ab_local_deg;
                    v2["npmin"] = s6_rin.npmin;
                    v2["lmin"] = s6_rin.lmin_cm;
                    v2["dead_w"] = s6_rin.dead_w;
                    Json::Value cl(Json::arrayValue), cv(Json::arrayValue),
                        slp(Json::arrayValue), exm(Json::arrayValue),
                        ovj(Json::arrayValue), hp(Json::arrayValue);
                    for (int p = 0; p < 3; ++p) {
                        cl.append(s6_rin.close_mx[p]);
                        cv.append(s6_rin.cov[p]);
                        slp.append(s6_rin.slope[p]);
                        exm.append(s6_rin.ext_med[p]);
                        ovj.append(s6_rin.ov[p]);
                        hp.append(s6_rin.has_plane[p]);
                    }
                    v2["close"] = cl;
                    v2["cov"] = cv;
                    v2["slope"] = slp;
                    v2["ext_med"] = exm;
                    v2["ov"] = ovj;
                    v2["has_plane"] = hp;
                    rec["v2"] = v2;
                }
            }
            // doc pr/64 round 4: W-track exception provenance -- absent
            // entirely unless the flavor enables it, so every existing
            // flavor's dump records are byte-identical.
            if (w_track_excuse) {
                rec["w_track"] = true;
                rec["w_track_revived"] = s6_wtrack_revived;
                if (s6_win_valid) {
                    Json::Value v3;
                    v3["w_gap"] = s6_win.w_gap;
                    v3["w_sole"] = s6_win.w_sole_vote;
                    v3["npmin"] = s6_win.npmin;
                    v3["lmin"] = s6_win.lmin_cm;
                    v3["tmax"] = s6_win.tmax_cm;
                    v3["ab_global"] = s6_win.ab_global_deg;
                    v3["dead_w"] = s6_win.dead_w;
                    rec["v3"] = v3;
                }
            }
            rec["planes"] = oc56_dump_planes;
            oc56_dump.write(rec);
        }
        return killed;
    };

    // doc pr/62 S7 ("long-edge corridor connectivity"): only evaluated when
    // long_check (default false, no-op). Independent additional OR-kill,
    // symmetric in structure to two_d_gap_kill above but covering exactly
    // the band that lambda declines to judge -- candidates at or above
    // s7_dis_floor_cm (== S6's s6_dis_cap, 30cm). See file header and
    // WireCellClus/Graphs.h (Graphs::long_corridor_bad) for the design and
    // root cause. gi1/gi2 are the GLOBAL point indices of q1/q2 themselves
    // (q1/q2 are cloud members, not arbitrary points), which is what lets
    // the corridor be anchored on their own native lattice cells
    // (cluster.wire_index + blob_with_point(...)->slice_index_min()) with
    // no 3D->lattice projection in the verdict path.
    auto long_gap_kill = [&](size_t j, size_t k, const char* blk,
                             const std::shared_ptr<Simple3DPointCloud>& cloud_a,
                             const std::vector<size_t>& gidx_a,
                             const std::shared_ptr<Simple3DPointCloud>& cloud_b,
                             const std::vector<size_t>& gidx_b,
                             size_t gi1, size_t gi2,
                             const geo_point_t& q1, const WirePlaneId& wpid_q1,
                             const geo_point_t& q2, const WirePlaneId& wpid_q2,
                             double edge_dis, double a1, double a2) -> bool {
        if (!long_check) return false;

        // No double jeopardy with S6: identically false below the seam.
        if (edge_dis < s7_dis_floor_cm * units::cm) return false;

        // Only meaningful within a single APA/face, same reasoning as S6.
        if (wpid_q1.apa() < 0 || wpid_q1.apa() != wpid_q2.apa() || wpid_q1.face() != wpid_q2.face()) {
            return false;
        }
        const int apa = wpid_q1.apa();
        const int face = wpid_q1.face();

        const auto* blob1 = cluster.blob_with_point(gi1);
        const auto* blob2 = cluster.blob_with_point(gi2);
        if (!blob1 || !blob2) return false;  // no blob data; abstain

        constexpr double s7_pi = 3.141592653589793;
        constexpr double s7_wire_angle_tol = 12.5 / 180.0 * s7_pi;
        const bool excuse_u = a1 < s7_wire_angle_tol;
        const bool excuse_v = a2 < s7_wire_angle_tol;

        // Fixed 2cm seed radius (NOT edge_dis-scaled like S6's) -- the
        // change that keeps the search O(D). q1/q2 are themselves cloud
        // members so the radius query is never empty.
        const double s7_seed_radius = 2.0 * units::cm;
        auto seeds_a_wc = cloud_a->get_closest_wcpoints_radius(q1, s7_seed_radius);
        auto seeds_b_wc = cloud_b->get_closest_wcpoints_radius(q2, s7_seed_radius);
        if (seeds_a_wc.empty() || seeds_b_wc.empty()) return false;

        auto blob_stride = [&](const std::vector<std::pair<size_t, geo_point_t>>& wc,
                                const std::vector<size_t>& gidx, int& best) {
            for (const auto& p : wc) {
                const auto* blob = cluster.blob_with_point(gidx.at(p.first));
                if (!blob) continue;
                const int span = blob->slice_index_max() - blob->slice_index_min();
                if (span > 0 && (best <= 0 || span < best)) best = span;
            }
        };
        int slice_step = -1;
        blob_stride(seeds_a_wc, gidx_a, slice_step);
        blob_stride(seeds_b_wc, gidx_b, slice_step);
        if (slice_step <= 0) slice_step = 4;  // SBND MaskSlices tick_span default -- fallback only

        // Corridor half-width/endpoint-disk radius: a starting operating
        // point, not yet fitted -- doc pr/62, unlike S6, has no hand-scan
        // labels in this distance band. cap must fully cover the 2cm seed
        // radius's footprint in lattice units; hw is a conservative interior
        // value versus S6's own fixed excess reach (7 wires / 6 slice-steps,
        // :741-748 above). See doc pr/62 for the reasoning and the census
        // fields that would inform refitting it.
        const double s7_hw = 3.0;
        const double s7_cap = 8.0;
        const int s7_margin = static_cast<int>(std::ceil(s7_cap)) + 2;
        const size_t s7_cell_budget = 20000;  // backstop only; the O(D) design keeps this slack
        constexpr int s7_dw = 1, s7_ds = 1;   // S6's repaired (1,1) baseline

        bool has_plane[3] = {false, false, false};
        bool gap[3] = {false, false, false};
        bool budget_hit[3] = {false, false, false};
        double gap_cm[3] = {0, 0, 0};
        int u1arr[3] = {0, 0, 0}, u2arr[3] = {0, 0, 0};
        const int s1_anchor = blob1->slice_index_min();
        const int s2_anchor = blob2->slice_index_min();

        for (int plane = 0; plane != 3; ++plane) {
            const int u1 = cluster.wire_index(gi1, plane);
            const int u2 = cluster.wire_index(gi2, plane);
            u1arr[plane] = u1; u2arr[plane] = u2;

            std::vector<S6Cell> seeds_a2d, seeds_b2d;
            auto collect = [&](const std::vector<std::pair<size_t, geo_point_t>>& wc,
                                const std::vector<size_t>& gidx, std::vector<S6Cell>& out) {
                for (const auto& p : wc) {
                    const size_t gi = gidx.at(p.first);
                    const auto* blob = cluster.blob_with_point(gi);
                    if (!blob) continue;
                    out.push_back({cluster.wire_index(gi, plane), blob->slice_index_min()});
                }
            };
            collect(seeds_a_wc, gidx_a, seeds_a2d);
            collect(seeds_b_wc, gidx_b, seeds_b2d);
            if (seeds_a2d.empty() || seeds_b2d.empty()) continue;  // no data this plane; don't vote

            S7Corridor corr;
            corr.u1 = u1; corr.v1 = 0.0;
            corr.u2 = u2; corr.v2 = double(s2_anchor - s1_anchor) / double(slice_step);
            corr.hw = s7_hw; corr.cap = s7_cap;
            corr.slice_step = slice_step; corr.s_anchor = s1_anchor;

            const int win_wlo = std::min(u1, u2) - s7_margin;
            const int win_whi = std::max(u1, u2) + s7_margin;
            const int win_slo = std::min(s1_anchor, s2_anchor) - s7_margin * slice_step;
            const int win_shi = std::max(s1_anchor, s2_anchor) + s7_margin * slice_step;

            has_plane[plane] = true;
            double reach = 0.0;
            const bool connected = s7_corridor_connected(
                grouping, apa, face, plane, seeds_a2d, seeds_b2d, corr,
                win_wlo, win_whi, win_slo, win_shi,
                slice_step, s7_dw, s7_ds, s7_cell_budget,
                &budget_hit[plane], &reach);
            gap[plane] = !connected;
            gap_cm[plane] = gap[plane] ? (1.0 - reach) * (edge_dis / units::cm) : 0.0;
        }

        Graphs::S7CorridorInput s7in;
        s7in.dis_cm = edge_dis / units::cm;
        for (int p = 0; p < 3; ++p) {
            s7in.has_plane[p] = has_plane[p];
            s7in.gap[p] = gap[p];
            s7in.budget_hit[p] = budget_hit[p];
            s7in.gap_cm[p] = gap_cm[p];
        }
        s7in.excuse_u = excuse_u;
        s7in.excuse_v = excuse_v;
        const bool killed = Graphs::long_corridor_bad(s7in, long_min_planes);

        if (oc53_census) {
            oc53_log->debug(
                "OC62CENSUS-S7 edge blk={} j={} k={} dis={:.2f}cm apa={} face={} "
                "cell1=[{},{},{}] cell2=[{},{},{}] s1={} s2={} slice_step={} hw={:.1f} cap={:.1f} "
                "has=[{},{},{}] gap=[{},{},{}] budget=[{},{},{}] gapcm=[{:.2f},{:.2f},{:.2f}] "
                "excuse_u={} excuse_v={} min_gapped={} killed={}",
                blk, j, k, edge_dis / units::cm, apa, face,
                u1arr[0], u1arr[1], u1arr[2], u2arr[0], u2arr[1], u2arr[2],
                s1_anchor, s2_anchor, slice_step, s7_hw, s7_cap,
                has_plane[0], has_plane[1], has_plane[2],
                gap[0], gap[1], gap[2],
                budget_hit[0], budget_hit[1], budget_hit[2],
                gap_cm[0], gap_cm[1], gap_cm[2],
                excuse_u, excuse_v, long_min_planes, killed);
        }
        return killed;
    };

    // Calculate distances between components
    for (size_t j = 0; j != num; j++) {
        for (size_t k = j + 1; k != num; k++) {
            // Get closest points between components
            index_index_dis[j][k] = pt_clouds.at(j)->get_closest_points(*pt_clouds.at(k));

            // C.4: skip the expensive Hough probes when clouds are too far apart to benefit.
            // get_closest_point_along_vec is called with max_dis=80 cm, so a closest-pair
            // distance already >= 80 cm guarantees both directional probes return nothing.
            const bool close_enough = std::get<2>(index_index_dis[j][k]) < 80 * units::cm;

            // Skip small clouds
            if (close_enough &&
                ((num < 100 && pt_clouds.at(j)->get_num_points() > 100 && pt_clouds.at(k)->get_num_points() > 100 &&
                  (pt_clouds.at(j)->get_num_points() + pt_clouds.at(k)->get_num_points()) > 400) ||
                 (pt_clouds.at(j)->get_num_points() > 500 && pt_clouds.at(k)->get_num_points() > 500))) {

                // Get closest points and calculate directions
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis[j][k]));
                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis[j][k]));

                geo_vector_t dir1 = cluster.vhough_transform(p1, 30 * units::cm, Cluster::HoughParamSpace::theta_phi, pt_clouds.at(j), 
                                                pt_clouds_global_indices.at(j));
                geo_vector_t dir2 = cluster.vhough_transform(p2, 30 * units::cm, Cluster::HoughParamSpace::theta_phi, pt_clouds.at(k),
                                                pt_clouds_global_indices.at(k)); 
                dir1 = dir1 * -1;
                dir2 = dir2 * -1;

                std::pair<int, double> result1 = pt_clouds.at(k)->get_closest_point_along_vec(p1, dir1, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm);

                if (result1.first >= 0) {
                    index_index_dis_dir1[j][k] = std::make_tuple(std::get<0>(index_index_dis[j][k]), 
                                                                result1.first, result1.second);
                }

                std::pair<int, double> result2 = pt_clouds.at(j)->get_closest_point_along_vec(p2, dir2, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm); 

                if (result2.first >= 0) {
                    index_index_dis_dir2[j][k] = std::make_tuple(result2.first,
                                                                std::get<1>(index_index_dis[j][k]), 
                                                                result2.second);
                }
            }
            // Now check the path 

            {
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis[j][k]));
                auto wpid_p1 = cluster.wire_plane_id(pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis[j][k])));
                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis[j][k]));
                auto wpid_p2 = cluster.wire_plane_id(pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis[j][k])));

                double dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
                double step_dis = 1.0 * units::cm;
                int num_steps = dis/step_dis + 1;

                

                // Track different types of "bad" points
                int num_bad[4] = {0,0,0,0};   // more than one of three are bad
                int num_bad1[4] = {0,0,0,0};  // at least one of three are bad
                int num_bad2[3] = {0,0,0};    // number of dead channels

                // Check points along path
                for (int ii = 0; ii != num_steps; ii++) {
                    geo_point_t test_p(
                        p1.x() + (p2.x() - p1.x())/num_steps*(ii + 1),
                        p1.y() + (p2.y() - p1.y())/num_steps*(ii + 1),
                        p1.z() + (p2.z() - p1.z())/num_steps*(ii + 1)
                    );

                    // Test point quality using grouping parameters
                    std::vector<int> scores;
                    if (use_ctpc) {
                        auto test_wpid = get_wireplaneid(test_p, wpid_p1, wpid_p2, dv);
                        if (test_wpid.apa()!=-1){
                            geo_point_t test_p_raw = test_p;
                            if (needs_transform) {
                                test_p_raw = ctpc_transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                            }
                            scores = grouping->test_good_point(test_p_raw, test_wpid.apa(), test_wpid.face());

                            // Check overall quality
                            if (scores[0] + scores[3] + scores[1] + scores[4] + (scores[2]+scores[5])*2 < 3) {
                                num_bad[0]++;
                            }
                            if (scores[0]+scores[3]==0) num_bad[1]++;
                            if (scores[1]+scores[4]==0) num_bad[2]++;
                            if (scores[2]+scores[5]==0) num_bad[3]++;

                            if (scores[3]!=0) num_bad2[0]++;
                            if (scores[4]!=0) num_bad2[1]++;
                            if (scores[5]!=0) num_bad2[2]++;

                            if (scores[0] + scores[3] + scores[1] + scores[4] + (scores[2]+scores[5]) < 3) {
                                num_bad1[0]++;
                            }
                            if (scores[0]+scores[3]==0) num_bad1[1]++;
                            if (scores[1]+scores[4]==0) num_bad1[2]++;
                            if (scores[2]+scores[5]==0) num_bad1[3]++;
                        } else {
                            // Step is outside all APA volumes (between APAs).
                            // Count as fully bad: no signal from any plane can validate this gap.
                            num_bad[0]++;  num_bad[1]++;  num_bad[2]++;  num_bad[3]++;
                            num_bad1[0]++; num_bad1[1]++; num_bad1[2]++; num_bad1[3]++;
                        }
                    }
                }

                auto test_wpid = get_wireplaneid(p1, wpid_p1, p2, wpid_p2, dv);

                // Calculate angles between directions
                geo_vector_t tempV1(0, p2.y() - p1.y(), p2.z() - p1.z());
                geo_vector_t tempV5;

                double angle1 = tempV1.angle(af_U_dir.at({test_wpid.apa(), test_wpid.face()})); 
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2)) * sin(angle1),
                        0);
                angle1 = tempV5.angle(drift_dir_abs);

                double angle2 = tempV1.angle(af_V_dir.at({test_wpid.apa(), test_wpid.face()}));
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2)) * sin(angle2),
                        0);
                angle2 = tempV5.angle(drift_dir_abs);

                double angle1p = tempV1.angle(af_W_dir.at({test_wpid.apa(), test_wpid.face()}));
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2)) * sin(angle1p),
                        0); 
                angle1p = tempV5.angle(drift_dir_abs);

                tempV5.set(p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z());
                double angle3 = tempV5.angle(drift_dir_abs);

                bool flag_strong_check = true;

                // Define constants for readability
                constexpr double pi = 3.141592653589793;
                constexpr double perp_angle_tol = 10.0/180.0*pi;
                constexpr double wire_angle_tol = 12.5/180.0*pi;
                constexpr double perp_angle = pi/2.0;
                constexpr double invalid_dist = 1e9;

                if (fabs(angle3 - perp_angle) < perp_angle_tol) {
                    geo_vector_t tempV2 = cluster.vhough_transform(p1, 15*units::cm);
                    geo_vector_t tempV3 = cluster.vhough_transform(p2, 15*units::cm);
                    
                    if (fabs(tempV2.angle(drift_dir_abs) - perp_angle) < perp_angle_tol &&
                        fabs(tempV3.angle(drift_dir_abs) - perp_angle) < perp_angle_tol) {
                        flag_strong_check = false;
                    }
                }
                else if (angle1 < wire_angle_tol || angle2 < wire_angle_tol || angle1p < wire_angle_tol) {
                    flag_strong_check = false;
                }

                // Helper function to invalidate distance
                auto invalidate_distance = [&]() {
                    index_index_dis[j][k] = std::make_tuple(-1, -1, invalid_dist);
                };

                // Kill tests: relaxed_strict_bad = S1+S2 (see file header);
                // branch structure identical to connect_graph_relaxed.
                if (flag_strong_check) {
                    if (relaxed_strict_bad(num_bad1[0], num_steps)) {
                        invalidate_distance();
                    }
                }
                else {
                    bool parallel_angles = (angle1 < wire_angle_tol && angle2 < wire_angle_tol) ||
                                        (angle1p < wire_angle_tol && angle1 < wire_angle_tol) ||
                                        (angle1p < wire_angle_tol && angle2 < wire_angle_tol);

                    if (parallel_angles) {
                        if (relaxed_strict_bad(num_bad[0], num_steps)) {
                            invalidate_distance();
                        }
                    }
                    else if (angle1 < wire_angle_tol) {
                        int sum_bad = num_bad[2] + num_bad[3];
                        if (relaxed_strict_bad(sum_bad, num_steps, 9) || num_bad[3] >= 3) {
                            invalidate_distance();
                        }
                    }
                    else if (angle2 < wire_angle_tol) {
                        int sum_bad = num_bad[1] + num_bad[3];
                        if (relaxed_strict_bad(sum_bad, num_steps, 9) || num_bad[3] >= 3) {
                            invalidate_distance();
                        }
                    }
                    else if (angle1p < wire_angle_tol) {
                        int sum_bad = num_bad[2] + num_bad[1];
                        if (relaxed_strict_bad(sum_bad, num_steps, 9)) {
                            invalidate_distance();
                        }
                    }
                    else if (relaxed_strict_bad(num_bad[0], num_steps)) {
                        invalidate_distance();
                    }
                }

                // doc pr/53 round 7 sec 18 S5: independent OR-kill, evaluated
                // regardless of branch (unlike relaxed_strict_bad above, the
                // 3D-image test does not depend on wire-direction angles).
                // No-op, byte-identical to round-6 "relaxed_strict" when
                // image_check is false (default).
                int ghost_run = 0;
                if (image_check) {
                    ghost_run = max_ghost_run(p1, p2, num_steps, wpid_p1, wpid_p2);
                    if (relaxed_img_bad(ghost_run, dis / units::cm)) {
                        invalidate_distance();
                    }
                    if (oc53_census) {
                        oc53_log->debug(
                            "OC53CENSUS-IMG closest j={} k={} dis={:.2f}cm nsteps={} ghost_run={} killed={}",
                            j, k, dis / units::cm, num_steps, ghost_run,
                            std::get<0>(index_index_dis[j][k]) < 0);
                    }
                }

                // doc pr/56 round 2 S6: independent additional OR-kill (owner's
                // rule: kill on >= 1 non-excused per-plane 2D gap).  No-op,
                // byte-identical to image_check-only behavior when
                // two_d_check is false (default).  Only evaluated on
                // candidates S1-S5 above have not already killed.
                if (two_d_check && std::get<0>(index_index_dis[j][k]) >= 0) {
                    if (two_d_gap_kill(j, k, "closest", pt_clouds.at(j), pt_clouds_global_indices.at(j),
                                        pt_clouds.at(k), pt_clouds_global_indices.at(k),
                                        p1, wpid_p1, p2, wpid_p2, dis, angle1, angle2)) {
                        invalidate_distance();
                    }
                }

                // doc pr/62 S7: independent additional OR-kill on the band
                // S6 declines (edge_dis >= s7_dis_floor_cm). No-op,
                // byte-identical to two_d_rescue behavior when long_check is
                // false (default). Only evaluated on candidates S1-S6 above
                // have not already killed.
                if (long_check && std::get<0>(index_index_dis[j][k]) >= 0) {
                    const size_t gi1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis[j][k]));
                    const size_t gi2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis[j][k]));
                    if (long_gap_kill(j, k, "closest", pt_clouds.at(j), pt_clouds_global_indices.at(j),
                                      pt_clouds.at(k), pt_clouds_global_indices.at(k), gi1, gi2,
                                      p1, wpid_p1, p2, wpid_p2, dis, angle1, angle2)) {
                        invalidate_distance();
                    }
                }

                if (oc53_census) {
                    const bool killed = std::get<0>(index_index_dis[j][k]) < 0;
                    oc53_log->debug(
                        "OC53CENSUS-S closest j={} k={} p1=({:.2f},{:.2f},{:.2f}) p2=({:.2f},{:.2f},{:.2f}) "
                        "dis={:.2f}cm nsteps={} strong={} nb=[{},{},{},{}] nb1=[{},{},{},{}] killed={}",
                        j, k, p1.x()/units::cm, p1.y()/units::cm, p1.z()/units::cm,
                        p2.x()/units::cm, p2.y()/units::cm, p2.z()/units::cm,
                        dis/units::cm, num_steps, flag_strong_check,
                        num_bad[0], num_bad[1], num_bad[2], num_bad[3],
                        num_bad1[0], num_bad1[1], num_bad1[2], num_bad1[3], killed);
                }
            }

            // Now check path again ...
            if (std::get<0>(index_index_dis_dir1[j][k]) >= 0) {
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis_dir1[j][k])); 
                auto wpid_p1 = cluster.wire_plane_id(pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir1[j][k])));
                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis_dir1[j][k])); 
                auto wpid_p2 = cluster.wire_plane_id(pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir1[j][k])));

                double dis = sqrt(pow(p1.x() - p2.x(), 2) + 
                                pow(p1.y() - p2.y(), 2) + 
                                pow(p1.z() - p2.z(), 2));
                double step_dis = 1.0 * units::cm;
                int num_steps = dis/step_dis + 1;
                int num_bad = 0;
                int num_bad1 = 0;

                // Check intermediate points along path
                for (int ii = 0; ii != num_steps; ii++) {
                    geo_point_t test_p(
                        p1.x() + (p2.x() - p1.x())/num_steps*(ii + 1),
                        p1.y() + (p2.y() - p1.y())/num_steps*(ii + 1),
                        p1.z() + (p2.z() - p1.z())/num_steps*(ii + 1)
                    );

                    if (use_ctpc) {
                        auto test_wpid = get_wireplaneid(test_p, wpid_p1, wpid_p2, dv);
                        if (test_wpid.apa()!=-1){
                            geo_point_t test_p_raw = test_p;
                            if (needs_transform) {
                                test_p_raw = ctpc_transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                            }
                            const bool good_point = grouping->is_good_point(test_p_raw, test_wpid.apa(), test_wpid.face());
                            if (!good_point) {
                                num_bad++;
                            }
                            if (!grouping->is_good_point(test_p_raw, test_wpid.apa(), test_wpid.face(), 0.6*units::cm, 1, 0)) {
                                num_bad1++;
                            }
                        } else {
                            // Step is outside all APA volumes — count as bad.
                            num_bad++;
                            num_bad1++;
                        }
                    }
                }

                auto test_wpid = get_wireplaneid(p1, wpid_p1, p2, wpid_p2, dv);

                // Calculate angles
                geo_vector_t tempV1(0, p2.y() - p1.y(), p2.z() - p1.z());
                geo_vector_t tempV5;

                double angle1 = tempV1.angle(af_U_dir.at({test_wpid.apa(), test_wpid.face()}));
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2))*sin(angle1),
                        0);
                angle1 = tempV5.angle(drift_dir_abs);
                
                double angle2 = tempV1.angle(af_V_dir.at({test_wpid.apa(), test_wpid.face()}));
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2))*sin(angle2),
                        0);
                angle2 = tempV5.angle(drift_dir_abs);
                
                tempV5.set(p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z());
                double angle3 = tempV5.angle(drift_dir_abs);
                
                double angle1p = tempV1.angle(af_W_dir.at({test_wpid.apa(), test_wpid.face()}));
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2))*sin(angle1p),
                        0);
                angle1p = tempV5.angle(drift_dir_abs);

                const double pi = 3.141592653589793;
                if (fabs(angle3 - pi/2) < 10.0/180.0*pi || 
                    angle1 < 12.5/180.0*pi ||
                    angle2 < 12.5/180.0*pi || 
                    angle1p < 7.5/180.0*pi) {
                    // Parallel or prolonged case
                    if (relaxed_strict_bad(num_bad, num_steps)) {
                        index_index_dis_dir1[j][k] = std::make_tuple(-1, -1, 1e9);
                    }
                }
                else {
                    if (relaxed_strict_bad(num_bad1, num_steps)) {
                        index_index_dis_dir1[j][k] = std::make_tuple(-1, -1, 1e9);
                    }
                }

                // doc pr/53 round 7 sec 18 S5 (see the closest-pair block for
                // the full comment). No-op when image_check is false.
                int dir1_ghost_run = 0;
                if (image_check) {
                    dir1_ghost_run = max_ghost_run(p1, p2, num_steps, wpid_p1, wpid_p2);
                    if (relaxed_img_bad(dir1_ghost_run, dis / units::cm)) {
                        index_index_dis_dir1[j][k] = std::make_tuple(-1, -1, 1e9);
                    }
                    if (oc53_census) {
                        oc53_log->debug("OC53CENSUS-IMG dir1 j={} k={} dis={:.2f}cm nsteps={} ghost_run={} killed={}",
                                        j, k, dis / units::cm, num_steps, dir1_ghost_run,
                                        std::get<0>(index_index_dis_dir1[j][k]) < 0);
                    }
                }

                // doc pr/56 round 2 S6 (see the closest-pair block for the
                // full comment). No-op when two_d_check is false (default).
                if (two_d_check && std::get<0>(index_index_dis_dir1[j][k]) >= 0) {
                    if (two_d_gap_kill(j, k, "dir1", pt_clouds.at(j), pt_clouds_global_indices.at(j),
                                        pt_clouds.at(k), pt_clouds_global_indices.at(k),
                                        p1, wpid_p1, p2, wpid_p2, dis, angle1, angle2)) {
                        index_index_dis_dir1[j][k] = std::make_tuple(-1, -1, 1e9);
                    }
                }

                // doc pr/62 S7 (see the closest-pair block for the full
                // comment). No-op when long_check is false (default).
                if (long_check && std::get<0>(index_index_dis_dir1[j][k]) >= 0) {
                    const size_t gi1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir1[j][k]));
                    const size_t gi2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir1[j][k]));
                    if (long_gap_kill(j, k, "dir1", pt_clouds.at(j), pt_clouds_global_indices.at(j),
                                      pt_clouds.at(k), pt_clouds_global_indices.at(k), gi1, gi2,
                                      p1, wpid_p1, p2, wpid_p2, dis, angle1, angle2)) {
                        index_index_dis_dir1[j][k] = std::make_tuple(-1, -1, 1e9);
                    }
                }

                if (oc53_census) {
                    const bool killed = std::get<0>(index_index_dis_dir1[j][k]) < 0;
                    oc53_log->debug("OC53CENSUS-S dir1 j={} k={} dis={:.2f}cm nsteps={} nb={} nb1={} killed={}",
                                    j, k, dis/units::cm, num_steps, num_bad, num_bad1, killed);
                }
            }

            //Now check path again ... 
            // Now check the path...
            if (std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis_dir2[j][k]));
                auto wpid_p1 = cluster.wire_plane_id(pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir2[j][k])));
                geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis_dir2[j][k]));
                auto wpid_p2 = cluster.wire_plane_id(pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir2[j][k])));

                double dis = sqrt(pow(p1.x() - p2.x(), 2) + 
                                pow(p1.y() - p2.y(), 2) + 
                                pow(p1.z() - p2.z(), 2));
                double step_dis = 1.0 * units::cm;
                int num_steps = dis/step_dis + 1;
                int num_bad = 0;
                int num_bad1 = 0;

                // Check points along path
                for (int ii = 0; ii != num_steps; ii++) {
                    geo_point_t test_p(
                        p1.x() + (p2.x() - p1.x())/num_steps*(ii + 1),
                        p1.y() + (p2.y() - p1.y())/num_steps*(ii + 1),
                        p1.z() + (p2.z() - p1.z())/num_steps*(ii + 1)
                    );

                    if (use_ctpc) {
                        auto test_wpid = get_wireplaneid(test_p, wpid_p1, wpid_p2, dv);
                        if (test_wpid.apa()!=-1){
                            geo_point_t test_p_raw = test_p;
                            if (needs_transform) {
                                test_p_raw = ctpc_transform->backward(test_p, cluster_t0, test_wpid.face(), test_wpid.apa());
                            }
                            const bool good_point = grouping->is_good_point(test_p_raw, test_wpid.apa(), test_wpid.face());
                            if (!good_point) {
                                num_bad++;
                            }
                            if (!grouping->is_good_point(test_p_raw, test_wpid.apa(), test_wpid.face(), 0.6*units::cm, 1, 0)) {
                                num_bad1++;
                            }
                        } else {
                            // Step is outside all APA volumes — count as bad.
                            num_bad++;
                            num_bad1++;
                        }
                    }
                }

                auto test_wpid = get_wireplaneid(p1, wpid_p1, p2, wpid_p2, dv);

                // Calculate angles between directions
                geo_vector_t tempV1(0, p2.y() - p1.y(), p2.z() - p1.z());
                geo_vector_t tempV5;

                double angle1 = tempV1.angle(af_U_dir.at({test_wpid.apa(), test_wpid.face()}));
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2))*sin(angle1),
                        0);
                angle1 = tempV5.angle(drift_dir_abs);

                double angle2 = tempV1.angle(af_V_dir.at({test_wpid.apa(), test_wpid.face()}));
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2))*sin(angle2),
                        0);
                angle2 = tempV5.angle(drift_dir_abs);

                tempV5.set(p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z());
                double angle3 = tempV5.angle(drift_dir_abs);

                double angle1p = tempV1.angle(af_W_dir.at({test_wpid.apa(), test_wpid.face()}));
                tempV5.set(fabs(p2.x() - p1.x()),
                        sqrt(pow(p2.y() - p1.y(), 2) + pow(p2.z() - p1.z(), 2))*sin(angle1p),
                        0);
                angle1p = tempV5.angle(drift_dir_abs);

                const double pi = 3.141592653589793;
                bool is_parallel = fabs(angle3 - pi/2) < 10.0/180.0*pi || 
                                angle1 < 12.5/180.0*pi ||
                                angle2 < 12.5/180.0*pi || 
                                angle1p < 7.5/180.0*pi;

                if (is_parallel) {
                    // Parallel or prolonged case
                    if (relaxed_strict_bad(num_bad, num_steps)) {
                        index_index_dis_dir2[j][k] = std::make_tuple(-1, -1, 1e9);
                    }
                }
                else {
                    if (relaxed_strict_bad(num_bad1, num_steps)) {
                        index_index_dis_dir2[j][k] = std::make_tuple(-1, -1, 1e9);
                    }
                }

                // doc pr/53 round 7 sec 18 S5 (see the closest-pair block for
                // the full comment). No-op when image_check is false.
                int dir2_ghost_run = 0;
                if (image_check) {
                    dir2_ghost_run = max_ghost_run(p1, p2, num_steps, wpid_p1, wpid_p2);
                    if (relaxed_img_bad(dir2_ghost_run, dis / units::cm)) {
                        index_index_dis_dir2[j][k] = std::make_tuple(-1, -1, 1e9);
                    }
                    if (oc53_census) {
                        oc53_log->debug("OC53CENSUS-IMG dir2 j={} k={} dis={:.2f}cm nsteps={} ghost_run={} killed={}",
                                        j, k, dis / units::cm, num_steps, dir2_ghost_run,
                                        std::get<0>(index_index_dis_dir2[j][k]) < 0);
                    }
                }

                // doc pr/56 round 2 S6 (see the closest-pair block for the
                // full comment). No-op when two_d_check is false (default).
                if (two_d_check && std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                    if (two_d_gap_kill(j, k, "dir2", pt_clouds.at(j), pt_clouds_global_indices.at(j),
                                        pt_clouds.at(k), pt_clouds_global_indices.at(k),
                                        p1, wpid_p1, p2, wpid_p2, dis, angle1, angle2)) {
                        index_index_dis_dir2[j][k] = std::make_tuple(-1, -1, 1e9);
                    }
                }

                // doc pr/62 S7 (see the closest-pair block for the full
                // comment). No-op when long_check is false (default).
                if (long_check && std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                    const size_t gi1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir2[j][k]));
                    const size_t gi2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir2[j][k]));
                    if (long_gap_kill(j, k, "dir2", pt_clouds.at(j), pt_clouds_global_indices.at(j),
                                      pt_clouds.at(k), pt_clouds_global_indices.at(k), gi1, gi2,
                                      p1, wpid_p1, p2, wpid_p2, dis, angle1, angle2)) {
                        index_index_dis_dir2[j][k] = std::make_tuple(-1, -1, 1e9);
                    }
                }

                if (oc53_census) {
                    const bool killed = std::get<0>(index_index_dis_dir2[j][k]) < 0;
                    oc53_log->debug("OC53CENSUS-S dir2 j={} k={} dis={:.2f}cm nsteps={} nb={} nb1={} killed={}",
                                    j, k, dis/units::cm, num_steps, num_bad, num_bad1, killed);
                }
            }
        }
    }

    // deal with MST of first type
    {
        Weighted::Graph temp_graph(num);
        for (size_t j = 0; j != num; j++) {
            for (size_t k = j + 1; k != num; k++) {
                int index1 = j;
                int index2 = k;
                if (std::get<0>(index_index_dis[j][k]) >= 0) {
                    if (!boost::edge(index1, index2, temp_graph).second) {
                        add_edge(index1, index2, std::get<2>(index_index_dis[j][k]), temp_graph);
                    }
                }
            }
        }

        // Process MST
        process_mst_deterministically(temp_graph, index_index_dis, index_index_dis_mst);
    }

    // MST of the direction ...
    {
        Weighted::Graph temp_graph(num);

        for (size_t j = 0; j != num; j++) {
            for (size_t k = j + 1; k != num; k++) {
                int index1 = j;
                int index2 = k;
                if (std::get<0>(index_index_dis_dir1[j][k]) >= 0 || std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                    if (!boost::edge(index1, index2, temp_graph).second) {
                        // Add edge with minimum distance from both directions
                    add_edge(
                        index1, index2,
                        std::min(std::get<2>(index_index_dis_dir1[j][k]), std::get<2>(index_index_dis_dir2[j][k])),
                        temp_graph);
                    }
                }
            }
        }

        process_mst_deterministically(temp_graph, index_index_dis, index_index_dis_dir_mst);

    }

    std::vector<std::pair<size_t, size_t>> oc53_pairs;  // emitted component pairs (census only)

    // doc pr/57 round 4: dump-only record of the edges this call actually
    // emits between components.  The hand-scan display shows one CANDIDATE at
    // a time and its `killed` bit, but the owner's judgement is about the two
    // components -- and one pair gets up to three independent candidates
    // (closest/dir1/dir2), while the pair can also stay joined through a third
    // component.  So "was this pair really separated" is only knowable here,
    // after every emit decision is made, which is why it is a separate record
    // rather than a field on the (much earlier) per-candidate edge records.
    // No-op when the dump is off; nothing below reads or perturbs the graph.
    Json::Value oc56_emitted(Json::arrayValue);
    auto oc56_note_emit = [&](size_t j, size_t k, const char* src, int gind1, int gind2,
                              double dis, bool dup) {
        if (!oc56_dump_on) return;
        Json::Value e;
        e["j"] = static_cast<int>(j);
        e["k"] = static_cast<int>(k);
        e["src"] = src;
        e["dis"] = dis / units::cm;
        // Endpoints from cluster.points() -- the same frame pt_clouds (and so
        // the dumped component point clouds) are built from, so the display
        // can draw this edge directly against them.
        Json::Value p1j(Json::arrayValue), p2j(Json::arrayValue);
        for (int d = 0; d != 3; ++d) p1j.append(points[d][gind1] / units::cm);
        for (int d = 0; d != 3; ++d) p2j.append(points[d][gind2] / units::cm);
        e["p1"] = p1j;
        e["p2"] = p2j;
        // dup: this exact vertex pair already carried an edge, so add_edge was
        // skipped.  The components are connected either way -- what changes is
        // only whether THIS emit created the link.
        e["dup"] = dup;
        oc56_emitted.append(e);
    };

    for (size_t j = 0; j != num; j++) {
        for (size_t k = j + 1; k != num; k++) {
            if (std::get<2>(index_index_dis[j][k]) < 3 * units::cm) {
                index_index_dis_mst[j][k] = index_index_dis[j][k];
            }

            // establish the path ...
            if (std::get<0>(index_index_dis_mst[j][k]) >= 0) {
                const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_mst[j][k]));
                const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_mst[j][k]));
                const float dis = std::get<2>(index_index_dis_mst[j][k]);
                const bool dup = boost::edge(gind1, gind2, graph).second;
                if (!dup) {
                    /*auto edge =*/ add_edge(gind1, gind2, dis, graph);
                }
                oc56_note_emit(j, k, "mst", gind1, gind2, dis, dup);
            }

            if (std::get<0>(index_index_dis_dir_mst[j][k]) >= 0) {
                if (std::get<0>(index_index_dis_dir1[j][k]) >= 0) {
                    const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir1[j][k]));
                    const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir1[j][k]));
                    float dis;
                    if (std::get<2>(index_index_dis_dir1[j][k]) > 5 * units::cm) {
                        dis = std::get<2>(index_index_dis_dir1[j][k]) * 1.1;
                    }
                    else {
                        dis = std::get<2>(index_index_dis_dir1[j][k]);
                    }
                    const bool dup = boost::edge(gind1, gind2, graph).second;
                    if(!dup) {
                        /*auto edge =*/ add_edge(gind1, gind2, dis, graph);
                    }
                    oc56_note_emit(j, k, "dir1", gind1, gind2, dis, dup);
                }
                if (std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                    const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir2[j][k]));
                    const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir2[j][k]));
                    float dis;
                    if (std::get<2>(index_index_dis_dir2[j][k]) > 5 * units::cm) {
                        dis = std::get<2>(index_index_dis_dir2[j][k]) * 1.1;
                    }
                    else {
                        dis = std::get<2>(index_index_dis_dir2[j][k]);
                    }
                    // }
                    const bool dup = boost::edge(gind1, gind2, graph).second;
                    if (!dup) {
                        /*auto edge =*/ add_edge(gind1, gind2, dis, graph);
                    }
                    oc56_note_emit(j, k, "dir2", gind1, gind2, dis, dup);
                }
            }
            else if ((two_d_rescue || w_track_excuse) && s6_rescued_dir_pairs.count({j, k})) {
                // doc pr/57 round 6: this pair's dir candidate was
                // explicitly RESCUED by two_d_rescue_ok, but the dir-MST
                // gate above is closed -- its recorded tuple is the
                // closest candidate's, which a killed closest sets to -1
                // (process_mst_deterministically), so every dir emission
                // for the pair silently evaporates.  A rescued candidate
                // is a validated artifact bridge; emit it anyway.  Only
                // pairs the rescue touched can enter here, so the reach
                // of this branch is exactly the fitted rescue population.
                if (std::get<0>(index_index_dis_dir1[j][k]) >= 0) {
                    const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir1[j][k]));
                    const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir1[j][k]));
                    float dis = std::get<2>(index_index_dis_dir1[j][k]);
                    if (dis > 5 * units::cm) dis *= 1.1;
                    const bool dup = boost::edge(gind1, gind2, graph).second;
                    if (!dup) {
                        add_edge(gind1, gind2, dis, graph);
                    }
                    oc56_note_emit(j, k, "dir1", gind1, gind2, dis, dup);
                }
                if (std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                    const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir2[j][k]));
                    const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir2[j][k]));
                    float dis = std::get<2>(index_index_dis_dir2[j][k]);
                    if (dis > 5 * units::cm) dis *= 1.1;
                    const bool dup = boost::edge(gind1, gind2, graph).second;
                    if (!dup) {
                        add_edge(gind1, gind2, dis, graph);
                    }
                    oc56_note_emit(j, k, "dir2", gind1, gind2, dis, dup);
                }
            }

            if (oc53_census) {
                const bool e_mst = std::get<0>(index_index_dis_mst[j][k]) >= 0;
                const bool e_dir = std::get<0>(index_index_dis_dir_mst[j][k]) >= 0 &&
                                   (std::get<0>(index_index_dis_dir1[j][k]) >= 0 ||
                                    std::get<0>(index_index_dis_dir2[j][k]) >= 0);
                if (e_mst || e_dir) oc53_pairs.emplace_back(j, k);
            }

        }  // k
    }  // j

    // doc pr/57 round 4: one connectivity record per graph call, written after
    // every emit decision above -- the graph is final here
    // (make_graph_relaxed_strict_img_2d_wfloor returns it directly).  `final`
    // is recomputed from the graph itself and `edges` is accumulated
    // independently at the three emit sites, so the display's reachability
    // answer can be cross-checked against a second, independent computation
    // (oc56_dump_check.py).  Keyed by the same globally-unique graph_call as
    // the component and edge records, so the join survives interleaved writes
    // from concurrently-processed clusters.
    if (oc56_dump_on) {
        std::vector<int> final_component(num_vertices(graph));
        const size_t final_num = connected_components(graph, &final_component[0]);
        Json::Value rec;
        rec["type"] = "connectivity";
        rec["graph_call"] = static_cast<int>(oc56_graph_call);
        rec["cluster_id"] = cluster.get_cluster_id();
        rec["ncomp"] = static_cast<int>(num);
        rec["nfinal"] = static_cast<int>(final_num);
        Json::Value finalj(Json::arrayValue);
        for (size_t c = 0; c != num; ++c) {
            // Every point of a starting component keeps a common final label,
            // so its first global index stands for the whole component.
            finalj.append(final_component.at(pt_clouds_global_indices.at(c).at(0)));
        }
        rec["final"] = finalj;
        rec["edges"] = oc56_emitted;
        oc56_dump.write(rec);
    }

    if (oc53_census) {
        oc53_log->debug("OC53CENSUS-S summary ncomp={} nedges={}", num, oc53_pairs.size());
        // Bridge status of each emitted component-pair edge: DFS reachability
        // over the remaining emitted pairs.  A bridge=true edge is the sole
        // link between its two components' sides.
        for (size_t e = 0; e < oc53_pairs.size(); ++e) {
            const size_t ja = oc53_pairs[e].first;
            const size_t ka = oc53_pairs[e].second;
            std::vector<std::vector<size_t>> adj(num);
            for (size_t f = 0; f < oc53_pairs.size(); ++f) {
                if (f == e) continue;
                adj[oc53_pairs[f].first].push_back(oc53_pairs[f].second);
                adj[oc53_pairs[f].second].push_back(oc53_pairs[f].first);
            }
            std::vector<char> seen(num, 0);
            std::vector<size_t> stack{ja};
            seen[ja] = 1;
            while (!stack.empty()) {
                const size_t v = stack.back();
                stack.pop_back();
                for (const size_t w : adj[v]) {
                    if (!seen[w]) { seen[w] = 1; stack.push_back(w); }
                }
            }
            const bool e_mst = std::get<0>(index_index_dis_mst[ja][ka]) >= 0;
            const bool e_dir1 = std::get<0>(index_index_dis_dir1[ja][ka]) >= 0;
            const bool e_dir2 = std::get<0>(index_index_dis_dir2[ja][ka]) >= 0;
            oc53_log->debug("OC53CENSUS-S edge j={} k={} dis={:.2f}cm mst={} dir1={} dir2={} bridge={}",
                            ja, ka, std::get<2>(index_index_dis[ja][ka])/units::cm,
                            e_mst, e_dir1, e_dir2, !seen[ka]);
        }
    }
}
