#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/ClusteringFuncsMixins.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"

class ClusteringSeparate;
WIRECELL_FACTORY(ClusteringSeparate, ClusteringSeparate,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)


using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Graphs;
using namespace WireCell::Clus::Facade;
using namespace WireCell::PointCloud::Tree;


// See declaration in ClusteringFuncs.h.  Chooses the fiducial volume matching the
// scope `dv` was configured for: per-APA stages get that drift volume's FV (union of
// the APA's configured faces); multi-APA (all-APA) stages get the "overall" cryostat
// FV, reproducing the legacy dv->metadata(WirePlaneId(0)) reads bit-for-bit.  The
// scope is read from dv->wpident_faces() (the configured drift volumes), NOT from the
// live grouping, so all-APA stays "overall" even when only one TPC has activity.
// Env-gated gate tracing (set WCT_SEP_DEBUG=1): prints per-cluster separation
// gate quantities.  Kept permanently because separation triggers depend on
// detector FV tuning and re-diagnosing them needs these numbers.
static bool sep_debug()
{
    static const bool v = (getenv("WCT_SEP_DEBUG") != nullptr);
    return v;
}

ScopeFV WireCell::Clus::Facade::select_scope_fv(IDetectorVolumes::pointer dv, bool common_face_x)
{
    const Configuration overall = dv->metadata(WirePlaneId(0));

    // Read a double field from a per-face block, falling back to "overall" when the
    // block does not define it (per-face blocks may define only x-bounds).
    auto field = [&](const Configuration& blk, const char* key) -> double {
        if (blk.isMember(key) && !blk[key].isNull()) return blk[key].asDouble();
        return overall[key].asDouble();
    };
    auto read_dir = [&](const char* key, geo_point_t def) -> geo_point_t {
        const Json::Value j = overall[key];
        if (!j.isNull() && j.isArray() && j.size() >= 3)
            return geo_point_t(j[0].asDouble(), j[1].asDouble(), j[2].asDouble());
        return def;
    };

    ScopeFV fv;
    // vertical_dir / beam_dir are detector-global physical constants.
    fv.vertical_dir = read_dir("vertical_dir", geo_point_t(0, 1, 0));
    fv.beam_dir = read_dir("beam_dir", geo_point_t(0, 0, 1));

    // Configured drift volumes (face-level wpids) this dv was built for.
    std::vector<WirePlaneId> faces;
    std::set<int> apas;
    for (const auto& kv : dv->wpident_faces()) {
        const WirePlaneId wpid(kv.first);
        faces.push_back(wpid);
        apas.insert(wpid.apa());
    }

    // Multi-APA or empty -> cryostat envelope (legacy behavior, bit-identical).
    if (faces.empty() || apas.size() > 1) {
        fv.xmin = field(overall, "FV_xmin");  fv.xmax = field(overall, "FV_xmax");
        fv.ymin = field(overall, "FV_ymin");  fv.ymax = field(overall, "FV_ymax");
        fv.zmin = field(overall, "FV_zmin");  fv.zmax = field(overall, "FV_zmax");
        fv.xmin_margin = field(overall, "FV_xmin_margin");  fv.xmax_margin = field(overall, "FV_xmax_margin");
        fv.ymin_margin = field(overall, "FV_ymin_margin");  fv.ymax_margin = field(overall, "FV_ymax_margin");
        fv.zmin_margin = field(overall, "FV_zmin_margin");  fv.zmax_margin = field(overall, "FV_zmax_margin");

        // Drift-side group: several APAs imaging one common drift side carry
        // identical per-face FV_x metadata -- adopt that x-range so out-of-time
        // apparent-x (past the cathode / behind the anode) is testable.  Mixed
        // faces (e.g. an all-APA scope spanning both drift sides) never agree
        // and keep the overall x.
        if (common_face_x && !faces.empty()) {
            bool common = true;
            double xmin = 0, xmax = 0, xmin_m = 0, xmax_m = 0;
            for (size_t i = 0; i != faces.size(); i++) {
                const Configuration blk = dv->metadata(faces[i]);
                const double bxmin = field(blk, "FV_xmin"), bxmax = field(blk, "FV_xmax");
                if (i == 0) {
                    xmin = bxmin;  xmax = bxmax;
                    xmin_m = field(blk, "FV_xmin_margin");  xmax_m = field(blk, "FV_xmax_margin");
                }
                else if (bxmin != xmin || bxmax != xmax) {
                    common = false;
                    break;
                }
            }
            if (common) {
                fv.xmin = xmin;  fv.xmax = xmax;
                fv.xmin_margin = xmin_m;  fv.xmax_margin = xmax_m;
            }
        }
        return fv;
    }

    // Single APA -> union (outermost envelope) over the configured per-(APA,face)
    // blocks.  For a single-face APA (e.g. SBND) this is just that one block.  Each
    // outward margin is carried from the face contributing that extreme.
    bool first = true;
    for (const auto& wpid : faces) {
        const Configuration blk = dv->metadata(wpid);
        const double xmin = field(blk, "FV_xmin"), xmax = field(blk, "FV_xmax");
        const double ymin = field(blk, "FV_ymin"), ymax = field(blk, "FV_ymax");
        const double zmin = field(blk, "FV_zmin"), zmax = field(blk, "FV_zmax");
        if (first) {
            fv.xmin = xmin;  fv.xmax = xmax;  fv.ymin = ymin;  fv.ymax = ymax;  fv.zmin = zmin;  fv.zmax = zmax;
            fv.xmin_margin = field(blk, "FV_xmin_margin");  fv.xmax_margin = field(blk, "FV_xmax_margin");
            fv.ymin_margin = field(blk, "FV_ymin_margin");  fv.ymax_margin = field(blk, "FV_ymax_margin");
            fv.zmin_margin = field(blk, "FV_zmin_margin");  fv.zmax_margin = field(blk, "FV_zmax_margin");
            first = false;
            continue;
        }
        if (xmin < fv.xmin) { fv.xmin = xmin; fv.xmin_margin = field(blk, "FV_xmin_margin"); }
        if (xmax > fv.xmax) { fv.xmax = xmax; fv.xmax_margin = field(blk, "FV_xmax_margin"); }
        if (ymin < fv.ymin) { fv.ymin = ymin; fv.ymin_margin = field(blk, "FV_ymin_margin"); }
        if (ymax > fv.ymax) { fv.ymax = ymax; fv.ymax_margin = field(blk, "FV_ymax_margin"); }
        if (zmin < fv.zmin) { fv.zmin = zmin; fv.zmin_margin = field(blk, "FV_zmin_margin"); }
        if (zmax > fv.zmax) { fv.zmax = zmax; fv.zmax_margin = field(blk, "FV_zmax_margin"); }
    }
    return fv;
}


// SBND-only two-track boundary tag (default OFF; gated by sbnd_boundary_tag).
//
// A single track touches any given planar detector face at most once, so >=2
// well-separated, charge-dense endpoints on ONE face imply >=2 tracks.  SBND's
// fiducial "outside" thresholds sit ~3.5-4 cm inside the wires (1 cm inset +
// 2.5-3 cm margin), so genuine cosmic endpoints that land right at the upstream
// (z-min) wall are captured as hull extremes but never flagged "outside" by
// JudgeSeparateDec_2 -- this tag uses a near-wall band to see them.  Per the
// agreed SBND scope it fires only when the upstream face carries two such tips
// AND a 3rd dense exit sits on a different face (the full evt-139220
// "two-on-upstream + one-elsewhere" topology), which is the lowest-false-split
// form.  On a trigger it seeds the detected tips into independent_points so
// Separate_1's endpoint selection picks the real prongs.  boundary_points is the
// convex hull already populated by JudgeSeparateDec_2.
static bool JudgeSeparateDec_SBND_boundary(const Cluster* cluster,
                                           const std::vector<geo_point_t>& boundary_points,
                                           const WireCell::Clus::Facade::ScopeFV& fv,
                                           const double cluster_length,
                                           std::vector<geo_point_t>& independent_points)
{
    // SBND tunables, hard-coded to match the existing threshold style in this file
    // (pinned against data evt 139220; see clustering-separate doc).
    const double LEN_MIN    = 250 * units::cm;  // long-cosmic regime
    const double Z_NEAR     = 10  * units::cm;  // upstream near-wall band: z < zmin + Z_NEAR
    const double OTHER_NEAR = 5   * units::cm;  // near-wall band for the other faces
    const double GROUP_DIS  = 25  * units::cm;  // merge upstream tips within this distance
    const double SEP_MIN    = 40  * units::cm;  // the two upstream tips must be this far apart
    const double GAP_MIN    = 10  * units::cm;  // tip-midpoint must be this far from charge
                                                // (two-track divergence ~14 cm on evt 139220;
                                                //  a single grazing track keeps it ~0-3 cm)
    const int    DENS_MIN   = 75;               // charge-dense endpoint (nnearby within 15 cm)

    if (cluster_length < LEN_MIN) return false;

    // 1. Dense hull points within the upstream (z-min) near-wall band.
    std::vector<geo_point_t> up;
    for (const auto& p : boundary_points) {
        if (p.z() < fv.zmin + Z_NEAR) {
            geo_point_t tp(p.x(), p.y(), p.z());
            if (cluster->nnearby(tp, 15 * units::cm) > DENS_MIN) up.push_back(p);
        }
    }
    if (up.size() < 2) return false;

    // 2. Greedy-group the upstream points (representative = first point of each group).
    std::vector<geo_point_t> reps;
    for (const auto& p : up) {
        bool merged = false;
        for (const auto& r : reps)
            if ((p - r).magnitude() < GROUP_DIS) { merged = true; break; }
        if (!merged) reps.push_back(p);
    }
    if (reps.size() < 2) return false;

    // the two most-separated upstream tips
    double best = -1; geo_point_t a, b;
    for (size_t i = 0; i + 1 < reps.size(); i++)
        for (size_t j = i + 1; j < reps.size(); j++) {
            double d = (reps[i] - reps[j]).magnitude();
            if (d > best) { best = d; a = reps[i]; b = reps[j]; }
        }
    if (best < SEP_MIN) return false;

    // 3. Gap guard: the midpoint between the two tips must lie in a charge gap; a
    //    single track grazing the wall would keep the midpoint sitting on charge.
    geo_point_t mid((a.x() + b.x()) / 2., (a.y() + b.y()) / 2., (a.z() + b.z()) / 2.);
    if (cluster->get_closest_dis(mid) < GAP_MIN) return false;

    // 4. Require >=1 dense exit on a face OTHER than upstream-zmin.
    geo_point_t other; bool have_other = false;
    for (const auto& p : boundary_points) {
        if (p.z() < fv.zmin + Z_NEAR) continue;  // skip the upstream band
        bool near = p.y() < fv.ymin + OTHER_NEAR || p.y() > fv.ymax - OTHER_NEAR ||
                    p.x() < fv.xmin + OTHER_NEAR || p.x() > fv.xmax - OTHER_NEAR ||
                    p.z() > fv.zmax - OTHER_NEAR;
        if (!near) continue;
        geo_point_t tp(p.x(), p.y(), p.z());
        if (cluster->nnearby(tp, 15 * units::cm) > DENS_MIN) { other = p; have_other = true; break; }
    }
    if (!have_other) return false;

    // Seed the detected prongs (15 cm dedup) so Separate_1 sees the real tips.
    auto seed = [&](const geo_point_t& q) {
        for (const auto& ip : independent_points)
            if ((ip - q).magnitude() < 15 * units::cm) return;
        independent_points.push_back(q);
    };
    seed(a); seed(b); seed(other);
    return true;
}


// ----------------------------------------------------------------------------
// Knob-gated refinements applied to the "family" of clusters produced by
// separating one original cluster (default OFF; gated by collinear_recover /
// band_recarve).  Tunables are hard-coded in the style of the rest of this
// file, pinned against PDVD run 39324 evt 0 group4567 -- see
// clus/docs/clustering-separate-refine.md.
// ----------------------------------------------------------------------------

// Step A (collinear_recover): recover stranded collinear track tips.
// Separate_1 carves its "path" cluster by 2D proximity to a graph shortest
// path, so a track tip sitting beyond a real imaging gap larger than
// Separate_2's 5 cm relink never joins the path and is left in the leftover
// cluster.  For each long, thin, non-isochronous family member T, move blobs
// from sibling members that continue T's PCA axis beyond its endpoints (small
// perpendicular offset, collinear local direction, contiguous axial run) into
// T.  Dead family entries are nulled in place, never erased.
// With interior_reclaim (collinear_interior knob) additionally claim sibling
// blobs sitting ON T's axis *inside* its span: at a track crossing the final
// Separate_2 connectivity relink can absorb an interior segment of T into the
// other track's cluster (the crossing region touches at <5 cm).
static void recover_collinear_tips(Grouping& live_grouping, const Tree::Scope& scope,
                                   std::vector<Cluster*>& family, bool interior_reclaim = false)
{
    const double TRACK_LEN_MIN   = 50 * units::cm;  // T must be a long ...
    const double TRACK_RATIO_MAX = 0.15;            // ... thin (pca eval1/eval0) track
    const double BAND_EXCL_ANG   = 15.;             // skip T near-isochronous (deg from perp-to-drift):
                                                    // wide bands are mutually near-collinear and would
                                                    // steal each other's far blobs
    const double PERP_MAX        = 4 * units::cm;   // blob center to T-axis distance
    const double GAP_MAX         = 15 * units::cm;  // max axial gap, endpoint->blob and blob->blob
    const double DIR_ANG_MAX     = 15.;             // blob local dir vs T axis (deg, sign-free)

    const geo_point_t drift_dir(1, 0, 0);

    for (size_t ti = 0; ti != family.size(); ti++) {
        Cluster* track = family.at(ti);
        if (!track) continue;
        if (track->get_length() < TRACK_LEN_MIN) continue;
        const auto& pca = track->get_pca();
        if (pca.axis.size() < 2 || pca.values.at(0) <= 0) continue;
        if (pca.values.at(1) / pca.values.at(0) >= TRACK_RATIO_MAX) continue;
        geo_point_t axis(pca.axis.at(0).x(), pca.axis.at(0).y(), pca.axis.at(0).z());
        if (fabs(axis.angle(drift_dir) - 3.1415926 / 2.) / 3.1415926 * 180. < BAND_EXCL_ANG) continue;

        const geo_point_t center = pca.center;
        const auto ends = track->get_main_axis_points();
        const double t_end1 = (ends.first - center).dot(axis);
        const double t_end2 = (ends.second - center).dot(axis);
        const double t_lo = std::min(t_end1, t_end2);
        const double t_hi = std::max(t_end1, t_end2);

        // Collect candidate blobs beyond either endpoint -- queries only,
        // against the frozen family; mutations happen below.
        struct Cand { size_t oi; int bi; double t; };
        std::vector<Cand> cands;
        for (size_t oi = 0; oi != family.size(); oi++) {
            if (oi == ti) continue;
            Cluster* donor = family.at(oi);
            if (!donor) continue;
            const auto& blobs = donor->children();
            for (int bi = 0; bi != (int) blobs.size(); bi++) {
                const geo_point_t bc = blobs.at(bi)->center_pos();
                const geo_point_t rel = bc - center;
                const double t = rel.dot(axis);
                if (t >= t_lo && t <= t_hi) continue;  // only beyond the endpoints
                const geo_point_t perp = rel - axis * t;
                if (perp.magnitude() >= PERP_MAX) continue;
                geo_point_t ldir = donor->vhough_transform(bc, 15 * units::cm);
                double ang = axis.angle(ldir) / 3.1415926 * 180.;
                if (ang > 90.) ang = 180. - ang;
                if (ang >= DIR_ANG_MAX) continue;
                cands.push_back({oi, bi, t});
            }
        }
        if (cands.empty() && !interior_reclaim) continue;

        // Contiguity: walk outward from each endpoint, accepting candidates
        // while the axial gap stays within GAP_MAX.  Sort by axial coordinate
        // with (family index, blob index) tie-breaks for determinism.
        std::map<size_t, std::vector<int>> moves;  // family idx -> blob indices
        auto walk = [&](bool high_end) {
            std::vector<const Cand*> side;
            for (const auto& c : cands)
                if (high_end ? c.t > t_hi : c.t < t_lo) side.push_back(&c);
            std::sort(side.begin(), side.end(), [&](const Cand* a, const Cand* b) {
                const double ta = high_end ? a->t : -a->t;
                const double tb = high_end ? b->t : -b->t;
                if (ta != tb) return ta < tb;
                if (a->oi != b->oi) return a->oi < b->oi;
                return a->bi < b->bi;
            });
            double prev = high_end ? t_hi : -t_lo;
            for (const Cand* c : side) {
                const double t = high_end ? c->t : -c->t;
                if (t - prev > GAP_MAX) break;
                prev = t;
                moves[c->oi].push_back(c->bi);
            }
        };
        walk(true);
        walk(false);

        // Interior reclaim: absorb whole SHORT sibling fragments lying along
        // T's axis inside its span.  Separate_1's carve can shed small
        // mid-track fragments of T (a crossing track's relink region pulls
        // them away) that later proximity merges then attach to the WRONG
        // track.  Fragment-level tests only -- per-blob local directions on a
        // 20 cm fragment are noise.  Tips are handled by the walk above.
        if (interior_reclaim) {
            const double FRAG_LEN_MAX     = 50 * units::cm;  // only short fragments move (claimers
                                                             // are >= TRACK_LEN_MIN, so disjoint)
            const double FRAG_PERP_MAX    = 8 * units::cm;   // every blob center near T's axis
            const double FRAG_DIR_ANG_MAX = 30.;             // fragment main axis vs T axis; the PDHD
                                                             // 27409 evt 40900 chunk reads ~19 deg
            const double FRAG_DIR_LEN_MIN = 6 * units::cm;   // below this the fragment pca direction
                                                             // is noise; geometry gates alone decide
            for (size_t oi = 0; oi != family.size(); oi++) {
                if (oi == ti) continue;
                if (moves.count(oi)) continue;  // tip walk already claims from this donor
                Cluster* donor = family.at(oi);
                if (!donor) continue;
                const double dlen = donor->get_length();
                if (dlen >= FRAG_LEN_MAX) continue;
                const auto& blobs = donor->children();
                if (blobs.empty()) continue;
                double max_perp = 0;
                bool in_span = true;
                for (const auto* b : blobs) {
                    const geo_point_t rel = geo_point_t(b->center_pos()) - center;
                    const double t = rel.dot(axis);
                    if (t <= t_lo || t >= t_hi) { in_span = false; break; }
                    const geo_point_t perp = rel - axis * t;
                    max_perp = std::max(max_perp, perp.magnitude());
                }
                double dang = -1;
                if (in_span && max_perp < FRAG_PERP_MAX && dlen >= FRAG_DIR_LEN_MIN) {
                    const auto& dpca = donor->get_pca();
                    if (!dpca.axis.empty() && dpca.values.at(0) > 0) {
                        geo_point_t daxis(dpca.axis.at(0).x(), dpca.axis.at(0).y(), dpca.axis.at(0).z());
                        dang = axis.angle(daxis) / 3.1415926 * 180.;
                        if (dang > 90.) dang = 180. - dang;
                    }
                }
                if (sep_debug())
                    std::cout << "SEPDBG intclaim T=" << ti << " donor=" << oi
                              << " len=" << dlen / units::cm << " nblob=" << blobs.size()
                              << " in_span=" << in_span << " max_perp=" << max_perp / units::cm
                              << " ang=" << dang << std::endl;
                if (!in_span || max_perp >= FRAG_PERP_MAX) continue;
                if (dang >= FRAG_DIR_ANG_MAX) continue;  // dang<0 (too short / degenerate pca) passes
                for (int bi = 0; bi != (int) blobs.size(); bi++) moves[oi].push_back(bi);
            }
        }

        if (moves.empty()) continue;

        const auto tscope = track->get_default_scope();
        const auto ttrans = track->get_scope_transform(tscope);
        for (auto& [oi, bidxs] : moves) {
            Cluster* donor = family.at(oi);
            if ((int) bidxs.size() == donor->nchildren()) {
                // whole donor moves; do not leave a childless husk behind
                track->take_children(*donor, true);
                family.at(oi) = nullptr;
                live_grouping.destroy_child(donor);
                assert(donor == nullptr);
            }
            else {
                std::vector<int> b2groupid(donor->nchildren(), -1);  // -1 = stay
                for (int bi : bidxs) b2groupid.at(bi) = 0;
                auto pieces = live_grouping.separate(donor, b2groupid, false);
                Cluster* piece = pieces.at(0);
                track->take_children(*piece, true);
                live_grouping.destroy_child(piece);
                assert(piece == nullptr);
            }
        }
        // take_children does not invalidate the facade cache; the
        // set_default_scope round-trip performs the full clear (see
        // Facade_Cluster.cxx clear_cache notes).
        track->set_default_scope(tscope);
        track->set_scope_filter(tscope, true);
        track->set_scope_transform(tscope, ttrans);

        size_t nmoved = 0;
        for (const auto& [oi, bidxs] : moves) nmoved += bidxs.size();
        std::cout << "Separate collinear_recover: track len " << track->get_length() / units::cm
                  << " cm claimed " << nmoved << " blobs from " << moves.size() << " donors" << std::endl;
    }
}

// Step A1b (collinear_member_merge): rejoin a single straight track that the
// carve cut into two (or more) long thin pieces.  collinear_recover cannot do
// this: it claims BLOBS along a track's axis, but a >=50 cm sibling is never
// an interior-reclaim donor and the 4 cm tip-walk gate loses a slightly-bent
// cosmic within a couple of metres (PDHD 27409 evt 40924: two pieces of one
// straight cosmic, axes 6.5 deg apart, touching at 0.3 cm, stayed two
// clusters).  Merge a touching pair of long thin members when their union is
// still one thin straight track.  Measured discriminators (blob-center proxy
// numbers from the display point clouds; re-pinned with SEPDBG):
//  - the 40924 pair (must merge): angle 6.5 deg, centroid-to-other-line
//    offsets 11.7/6.5 cm, union perp rms 9.1 cm, union eval ratio 0.006;
//  - 27409 evt 40904's forking pair (must stay): 20 deg / 42.6, 28.9 / 24.0;
//  - 27409 evt 40908's separated parallel tracks (tightest control, 5.1 deg):
//    offsets 26.9/28.2 cm, union rms 16.0 cm, ratio 0.014.
static void merge_collinear_members(Grouping& live_grouping, const Tree::Scope& scope,
                                    std::vector<Cluster*>& family,
                                    const double len_min_short = -1)
{
    const double LEN_MIN        = 100 * units::cm;  // longer piece (anchor)
    // Shorter piece's floor: default (<0) keeps the symmetric 100 cm gate;
    // the grouping-wide stitch lowers it (anchor + bridge, cathode_connect's
    // min_length_short pattern) -- PDVD 39252 evt 298637's 104 cm (display)
    // piece reads just under 100 cm by get_length().
    const double SHORT_MIN      = len_min_short > 0 ? len_min_short : LEN_MIN;
    const double THIN_RATIO_MAX = 0.05;             // pca eval1/eval0 per piece
    const double TOUCH_MAX      = 5 * units::cm;    // pieces of one cut track touch
    const double ANG_MAX        = 12.;              // main-axis angle (deg); a genuine
                                                    // one-track pair reads 10.01 (PDVD
                                                    // 39252 evt 298637) while the keep
                                                    // controls sit at 20 (fork) and 5.1
                                                    // (parallel pair, killed by rms)
    const double OFF_MAX        = 30 * units::cm;   // each centroid off the other's line;
                                                    // a long (321 cm) slightly-curved
                                                    // track puts its centroid 25.5 cm off
                                                    // its continuation's line -- union rms
                                                    // is the real discriminator
    const double UNION_RMS_MAX  = 7 * units::cm;    // union perp rms about its own 3D line:
                                                    // genuine one-cosmic rejoins read <=5.7
                                                    // across the PDHD+PDVD suite; the one
                                                    // fork-adjacent (kinked) pair that must
                                                    // not rejoin reads 8.1 (27409 evt 40904,
                                                    // whose downstream connect1 then fused
                                                    // the two full fork prongs)
    const double UNION_RATIO    = 0.010;            // union pca eval1/eval0

    bool merged_any = true;
    while (merged_any) {  // a track cut into 3 needs a second pass
        merged_any = false;
        for (size_t fi = 0; fi + 1 < family.size() && !merged_any; fi++) {
            Cluster* a = family.at(fi);
            if (!a) continue;
            if (a->get_length() < SHORT_MIN) continue;
            const auto& pca_a = a->get_pca();
            if (pca_a.axis.size() < 1 || pca_a.values.size() < 2 || pca_a.values.at(0) <= 0) continue;
            if (pca_a.values.at(1) / pca_a.values.at(0) >= THIN_RATIO_MAX) continue;
            const geo_point_t ca = pca_a.center;
            geo_point_t da(pca_a.axis.at(0).x(), pca_a.axis.at(0).y(), pca_a.axis.at(0).z());

            for (size_t fj = fi + 1; fj < family.size(); fj++) {
                Cluster* b = family.at(fj);
                if (!b) continue;
                if (b->get_length() < SHORT_MIN) continue;
                // at least one piece must be a long anchor
                if (a->get_length() < LEN_MIN && b->get_length() < LEN_MIN) continue;
                const auto& pca_b = b->get_pca();
                if (pca_b.axis.size() < 1 || pca_b.values.size() < 2 || pca_b.values.at(0) <= 0) continue;
                if (pca_b.values.at(1) / pca_b.values.at(0) >= THIN_RATIO_MAX) continue;
                const geo_point_t cb = pca_b.center;
                geo_point_t db(pca_b.axis.at(0).x(), pca_b.axis.at(0).y(), pca_b.axis.at(0).z());

                const double ang =
                    fabs(da.angle(db)) / 3.1415926 * 180.;
                const double ang_sf = std::min(ang, 180. - ang);  // sign-free
                // each centroid's perpendicular distance to the OTHER's line
                geo_point_t vab(cb.x() - ca.x(), cb.y() - ca.y(), cb.z() - ca.z());
                const double ta = vab.dot(da);
                geo_point_t ra(vab.x() - ta * da.x(), vab.y() - ta * da.y(), vab.z() - ta * da.z());
                geo_point_t vba(ca.x() - cb.x(), ca.y() - cb.y(), ca.z() - cb.z());
                const double tb = vba.dot(db);
                geo_point_t rb(vba.x() - tb * db.x(), vba.y() - tb * db.y(), vba.z() - tb * db.z());

                bool pass = (ang_sf < ANG_MAX) && (ra.magnitude() < OFF_MAX) && (rb.magnitude() < OFF_MAX);

                double touch = -1;
                if (pass) {
                    auto dis = a->get_closest_points(*b);
                    touch = std::get<2>(dis);
                    pass = touch < TOUCH_MAX;
                }

                double union_rms = -1, union_ratio = -1;
                if (pass) {
                    // union 3D pca of the pair's npoints-weighted blob centers:
                    // top eigenpair by power iteration + one deflation step
                    double sw = 0, sm[3] = {0, 0, 0};
                    for (Cluster* m : {a, b})
                        for (const Blob* blob : m->children()) {
                            const geo_point_t bc = blob->center_pos();
                            const double w = blob->npoints();
                            sw += w;
                            sm[0] += w * bc.x();
                            sm[1] += w * bc.y();
                            sm[2] += w * bc.z();
                        }
                    if (sw <= 0) { pass = false; }
                    else {
                        for (int i = 0; i != 3; i++) sm[i] /= sw;
                        double cov[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
                        for (Cluster* m : {a, b})
                            for (const Blob* blob : m->children()) {
                                const geo_point_t bc = blob->center_pos();
                                const double w = blob->npoints();
                                const double d[3] = {bc.x() - sm[0], bc.y() - sm[1], bc.z() - sm[2]};
                                for (int i = 0; i != 3; i++)
                                    for (int j = 0; j != 3; j++) cov[i][j] += w * d[i] * d[j];
                            }
                        for (int i = 0; i != 3; i++)
                            for (int j = 0; j != 3; j++) cov[i][j] /= sw;
                        auto power = [](const double m[3][3], const double s[3], double out[3]) {
                            double v[3] = {s[0], s[1], s[2]};
                            double n0 = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
                            if (n0 <= 0) { v[0] = 1; v[1] = 0; v[2] = 0; n0 = 1; }
                            for (int i = 0; i != 3; i++) v[i] /= n0;
                            double lam = 0;
                            for (int it = 0; it != 30; it++) {
                                double w[3] = {0, 0, 0};
                                for (int i = 0; i != 3; i++)
                                    for (int j = 0; j != 3; j++) w[i] += m[i][j] * v[j];
                                lam = v[0] * w[0] + v[1] * w[1] + v[2] * w[2];
                                const double n = std::sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
                                if (n <= 0) break;
                                for (int i = 0; i != 3; i++) v[i] = w[i] / n;
                            }
                            for (int i = 0; i != 3; i++) out[i] = v[i];
                            return lam;
                        };
                        const double seed[3] = {da.x(), da.y(), da.z()};
                        double e0[3];
                        const double lam0 = power(cov, seed, e0);
                        if (lam0 <= 0) { pass = false; }
                        else {
                            // perp rms about the principal line
                            double r2 = 0;
                            for (Cluster* m : {a, b})
                                for (const Blob* blob : m->children()) {
                                    const geo_point_t bc = blob->center_pos();
                                    const double w = blob->npoints();
                                    const double d[3] = {bc.x() - sm[0], bc.y() - sm[1], bc.z() - sm[2]};
                                    const double t = d[0] * e0[0] + d[1] * e0[1] + d[2] * e0[2];
                                    r2 += w * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2] - t * t);
                                }
                            union_rms = std::sqrt(r2 / sw);
                            // deflate and find the 2nd eigenvalue
                            double cov1[3][3];
                            for (int i = 0; i != 3; i++)
                                for (int j = 0; j != 3; j++) cov1[i][j] = cov[i][j] - lam0 * e0[i] * e0[j];
                            const double seed1[3] = {e0[1], -e0[0], e0[2]};  // anything off e0
                            double e1[3];
                            const double lam1 = power(cov1, seed1, e1);
                            union_ratio = lam1 > 0 ? lam1 / lam0 : 0;
                            pass = (union_rms < UNION_RMS_MAX) && (union_ratio < UNION_RATIO);
                        }
                    }
                }

                if (sep_debug())
                    std::cout << "SEPDBG colmerge m" << fi << "/m" << fj
                              << " lens=" << a->get_length() / units::cm << "/" << b->get_length() / units::cm
                              << " ca=(" << ca.x() / units::cm << "," << ca.y() / units::cm << "," << ca.z() / units::cm << ")"
                              << " cb=(" << cb.x() / units::cm << "," << cb.y() / units::cm << "," << cb.z() / units::cm << ")"
                              << " ang=" << ang_sf
                              << " offs=" << ra.magnitude() / units::cm << "/" << rb.magnitude() / units::cm
                              << " touch=" << (touch < 0 ? -1. : touch / units::cm)
                              << " union_rms=" << (union_rms < 0 ? -1. : union_rms / units::cm)
                              << " ratio=" << union_ratio
                              << (pass ? " MERGE" : " keep") << std::endl;
                if (!pass) continue;

                const auto tscope = a->get_default_scope();
                const auto ttrans = a->get_scope_transform(tscope);
                a->take_children(*b, true);
                family.at(fj) = nullptr;
                live_grouping.destroy_child(b);
                assert(b == nullptr);
                a->set_default_scope(tscope);
                a->set_scope_filter(tscope, true);
                a->set_scope_transform(tscope, ttrans);

                std::cout << "Separate collinear_member_merge: rejoined pieces"
                          << " len " << a->get_length() / units::cm << " cm, angle "
                          << ang_sf << " deg, union rms " << union_rms / units::cm
                          << " cm" << std::endl;
                merged_any = true;
                break;
            }
        }
    }
}

// Step A3 (band_merge_back): re-assemble a single isochronous band that the
// carve hatched into interleaved pieces.  With separation firing on band-
// topology clusters (FV insets), Separate_1's path carve can cut ONE wide band
// into 2+ clusters of alternating chunks (PDVD 39324 evts 339870/339930/
// 339990/340010) -- no recarve seed exists (the pieces are mutually near-
// collinear) so nothing re-assembles them.  Pool the touching band-like
// family members and merge the ones consistent with one band, discriminating
// against the two cases that must stay apart (both measured on 39324):
//  - two distinct parallel bands touching side-by-side (evt 339850: offset
//    means 94.5 cm / 2.27 (sigma_i+sigma_j) apart; hatched single-band pieces
//    read <= 47.4 cm / 1.26),
//  - two genuine bands crossing as an X: the union's residual about a single
//    line GROWS toward the axial ends (outer/inner rms ~2.0 synthetic; the
//    worst single-band case reads 1.48).
static void merge_back_bands(Grouping& live_grouping, const Tree::Scope& scope,
                             std::vector<Cluster*>& family)
{
    const double BAND_PERP_ANG  = 10.;             // member axis within this of perp-to-drift (deg)
    const double BAND_LEN_MIN   = 60 * units::cm;
    const double BAND_WIDTH_MIN = 6 * units::cm;   // rms width: thin track pieces must not pool
    const double TOUCH_MAX      = 2 * units::cm;   // pieces of one hatched band touch
    const double PAR_SEP_ABS    = 85 * units::cm;  // offset-mean gap: above = distinct parallel bands
                                                   // (hatched-piece gaps measure <= 67 cm, genuine
                                                   // parallel-band anchors 130 cm -- 39324)
    const double PAR_SEP_NSIG   = 1.7;             // ... in units of (sigma_i + sigma_j)
    const double CROSS_GROW_MAX = 1.75;            // union outer/inner rms: above = crossing bands
    const double X_SLAB_MAX     = 20 * units::cm;  // pieces of ONE iso band share its x-slab
                                                   // (slab ~12 cm wide; two bands at different
                                                   // drift times sit 50 cm apart -- 39324 339890)
    const double UNION_RMS_MAX  = 50 * units::cm;  // merged group's transverse rms about its own
                                                   // line: genuine re-assembled bands read 17-43,
                                                   // wrongly-fused complexes 57-64 (27409 40900,
                                                   // 39324 339890)

    const geo_point_t drift_dir(1, 0, 0);

    // 1. band-like family members, two tiers: wide members (the recarve band
    // test + its seed width gate) can anchor a band; thin band-aligned DEBRIS
    // (a short narrow chunk shed by the carve) may join a pool and merge into
    // a nearby anchor but can never anchor or merge on its own (the width gate
    // is what keeps thin crossing-track pieces from forming bands).  Debris
    // must also be SHORT: a long thin iso-aligned member is a track, not a
    // chunk (PDHD 27409 evt 40908: 318/503 cm tracks at 9.5 deg, width 2-3.5).
    const double DEBRIS_WIDTH_MIN = 2 * units::cm;
    const double DEBRIS_LEN_MAX   = 100 * units::cm;
    std::vector<size_t> band;       // indices into family
    std::vector<char> wide;         // parallel to `band`
    for (size_t i = 0; i != family.size(); i++) {
        Cluster* c = family.at(i);
        if (!c) continue;
        if (c->get_length() < BAND_LEN_MIN) continue;
        const auto& pca = c->get_pca();
        if (pca.axis.size() < 1 || pca.values.size() < 2 || pca.values.at(0) <= 0) continue;
        geo_point_t axis(pca.axis.at(0).x(), pca.axis.at(0).y(), pca.axis.at(0).z());
        const double perp_ang = fabs(axis.angle(drift_dir) - 3.1415926 / 2.) / 3.1415926 * 180.;
        const int n = c->npoints();
        const double width = (n >= 1) ? std::sqrt(pca.values.at(1) / n) : 0;
        if (sep_debug())
            std::cout << "SEPDBG bandtest m" << i << " len=" << c->get_length() / units::cm
                      << " perp_ang=" << perp_ang << " width=" << width / units::cm
                      << " npts=" << n << std::endl;
        if (perp_ang >= BAND_PERP_ANG) continue;
        if (n < 1 || width < DEBRIS_WIDTH_MIN) continue;
        if (width < BAND_WIDTH_MIN && c->get_length() >= DEBRIS_LEN_MAX) continue;
        band.push_back(i);
        wide.push_back(width >= BAND_WIDTH_MIN);
    }
    if (band.size() < 2) return;

    // 2. pools: connected components of the touch relation (deterministic BFS
    // in family order)
    std::vector<int> pool_id(band.size(), -1);
    int npool = 0;
    for (size_t i = 0; i != band.size(); i++) {
        if (pool_id[i] >= 0) continue;
        pool_id[i] = npool;
        std::vector<size_t> queue = {i};
        while (!queue.empty()) {
            const size_t qi = queue.back();
            queue.pop_back();
            for (size_t j = 0; j != band.size(); j++) {
                if (pool_id[j] >= 0) continue;
                std::tuple<int, int, double> dis =
                    family.at(band[qi])->get_closest_points(*family.at(band[j]));
                if (std::get<2>(dis) >= TOUCH_MAX) continue;
                pool_id[j] = npool;
                queue.push_back(j);
            }
        }
        npool++;
    }

    for (int pid = 0; pid != npool; pid++) {
        std::vector<size_t> pool;  // indices into family, ascending
        std::vector<char> pwide;   // parallel to `pool`
        for (size_t i = 0; i != band.size(); i++)
            if (pool_id[i] == pid) {
                pool.push_back(band[i]);
                pwide.push_back(wide[i]);
            }
        if (pool.size() < 2) continue;

        // 3. single-line fit of the pool union in (y,z), npoints-weighted blob
        // centers (bands are extended transverse to drift)
        double sw = 0, sy = 0, sz = 0, syy = 0, syz = 0, szz = 0;
        for (size_t pi : pool) {
            for (const Blob* blob : family.at(pi)->children()) {
                const geo_point_t bc = blob->center_pos();
                const double w = blob->npoints();
                sw += w;
                sy += w * bc.y();
                sz += w * bc.z();
                syy += w * bc.y() * bc.y();
                syz += w * bc.y() * bc.z();
                szz += w * bc.z() * bc.z();
            }
        }
        if (sw <= 0) continue;
        const double my = sy / sw, mz = sz / sw;
        const double cyy = syy / sw - my * my;
        const double cyz = syz / sw - my * mz;
        const double czz = szz / sw - mz * mz;
        const double theta = 0.5 * std::atan2(2 * cyz, cyy - czz);
        const double uy = std::cos(theta), uz = std::sin(theta);  // along the band
        const double ny = -uz, nz = uy;                           // across the band

        // per-member signed transverse offset mean/sigma + weight + x mean
        // (pieces of one hatched band share its x-slab; the y-z fit is x-blind)
        std::vector<double> mu(pool.size(), 0), sig(pool.size(), 0), w_mem(pool.size(), 0);
        std::vector<double> xmu(pool.size(), 0);
        std::vector<double> ts;  // unweighted blob axial coords, for quartiles
        for (size_t k = 0; k != pool.size(); k++) {
            double w_i = 0, sp = 0, spp = 0, sx = 0;
            for (const Blob* blob : family.at(pool[k])->children()) {
                const geo_point_t bc = blob->center_pos();
                const double w = blob->npoints();
                const double p = (bc.y() - my) * ny + (bc.z() - mz) * nz;
                w_i += w;
                sp += w * p;
                spp += w * p * p;
                sx += w * bc.x();
                ts.push_back((bc.y() - my) * uy + (bc.z() - mz) * uz);
            }
            w_mem[k] = w_i;
            if (w_i <= 0) continue;
            mu[k] = sp / w_i;
            const double var = spp / w_i - mu[k] * mu[k];
            sig[k] = var > 0 ? std::sqrt(var) : 0;
            xmu[k] = sx / w_i;
        }

        // crossing guard: residual about the single line must not grow toward
        // the axial ends (inner = middle two t-quartiles)
        std::sort(ts.begin(), ts.end());
        const double q1 = ts.at(ts.size() / 4);
        const double q3 = ts.at((3 * ts.size()) / 4);
        double w_in = 0, r2_in = 0, w_out = 0, r2_out = 0;
        for (size_t pi : pool) {
            for (const Blob* blob : family.at(pi)->children()) {
                const geo_point_t bc = blob->center_pos();
                const double w = blob->npoints();
                const double p = (bc.y() - my) * ny + (bc.z() - mz) * nz;
                const double t = (bc.y() - my) * uy + (bc.z() - mz) * uz;
                if (t >= q1 && t <= q3) { w_in += w; r2_in += w * p * p; }
                else { w_out += w; r2_out += w * p * p; }
            }
        }
        if (w_in <= 0 || w_out <= 0) continue;
        const double inner_rms = std::sqrt(r2_in / w_in);
        const double outer_rms = std::sqrt(r2_out / w_out);

        if (sep_debug()) {
            std::cout << "SEPDBG mergeback pool";
            for (size_t k = 0; k != pool.size(); k++)
                std::cout << " m" << pool[k] << "(mu=" << mu[k] / units::cm
                          << ",sig=" << sig[k] / units::cm
                          << ",x=" << xmu[k] / units::cm
                          << ",w=" << w_mem[k] << ")";
            std::cout << " inner=" << inner_rms / units::cm
                      << " outer=" << outer_rms / units::cm << std::endl;
        }
        if (inner_rms > 0 && outer_rms / inner_rms > CROSS_GROW_MAX) continue;

        // 4. anchor-based grouping by transverse offset.  Big pieces (the band
        // cores) group via the strict parallel-separation gates; small debris
        // joins the nearest anchor group when plausibly of the same band.
        // (Plain single-linkage over all pieces chains two parallel bands
        // together through their overlap-region debris -- 39324 evt 339850.)
        const double ANCHOR_FRAC = 0.20;  // of pool npoints
        double w_pool = 0;
        int k_big = -1;  // largest WIDE member; thin debris can never anchor
        for (size_t k = 0; k != pool.size(); k++) {
            w_pool += w_mem[k];
            if (pwide[k] && (k_big < 0 || w_mem[k] > w_mem[k_big])) k_big = (int) k;
        }
        if (k_big < 0) continue;  // pool is all debris: nothing to assemble around
        std::vector<bool> is_anchor(pool.size(), false);
        for (size_t k = 0; k != pool.size(); k++)
            is_anchor[k] = pwide[k] && (((int) k == k_big) || (w_mem[k] >= ANCHOR_FRAC * w_pool));
        std::vector<int> grp(pool.size());
        for (size_t k = 0; k != pool.size(); k++) grp[k] = (int) k;
        auto root = [&](int k) {
            while (grp[k] != k) k = grp[k] = grp[grp[k]];
            return k;
        };
        // anchors: single-linkage with both gates (anchors are few, so the
        // debris-chaining failure mode cannot recur).  The x-slab gate keeps
        // bands at different drift times apart even when their y-z projections
        // overlap (39324 evt 339890: two same-slope bands 50 cm apart in x).
        for (size_t a = 0; a + 1 < pool.size(); a++) {
            if (!is_anchor[a]) continue;
            for (size_t b = a + 1; b < pool.size(); b++) {
                if (!is_anchor[b]) continue;
                const double dmu = fabs(mu[a] - mu[b]);
                if (dmu >= PAR_SEP_ABS) continue;
                if (dmu >= PAR_SEP_NSIG * (sig[a] + sig[b])) continue;
                if (fabs(xmu[a] - xmu[b]) >= X_SLAB_MAX) continue;
                const int ra = root(a), rb = root(b);
                if (ra != rb) grp[std::max(ra, rb)] = std::min(ra, rb);
            }
        }
        // debris: nearest same-slab anchor by offset, absolute gate only
        // (debris sigmas are unreliable); too-far debris stays its own cluster
        for (size_t k = 0; k != pool.size(); k++) {
            if (is_anchor[k]) continue;
            double best = PAR_SEP_ABS;
            int ka = -1;
            for (size_t a = 0; a != pool.size(); a++) {
                if (!is_anchor[a]) continue;
                if (fabs(xmu[k] - xmu[a]) >= X_SLAB_MAX) continue;
                const double dmu = fabs(mu[k] - mu[a]);
                if (dmu < best) { best = dmu; ka = (int) a; }
            }
            if (ka < 0) continue;
            const int rk = root((int) k), ra = root(ka);
            if (rk != ra) grp[std::max(rk, ra)] = std::min(rk, ra);
        }

        // 5. merge each multi-member group into its lowest-family-index member
        for (size_t a = 0; a != pool.size(); a++) {
            if (root(a) != (int) a) continue;  // not a group leader
            std::vector<size_t> members;       // pool-local indices, ascending
            for (size_t b = a; b < pool.size(); b++)
                if (root(b) == (int) a) members.push_back(b);
            if (members.size() < 2) continue;

            // union-rms cap: refit the GROUP's own y-z line and require its
            // transverse rms to stay band-sized.  A genuine re-assembled band
            // reads 17-43 cm; two distinct same-slab complexes wrongly grouped
            // read 57-64 cm (27409 evt 40900, 39324 evt 339890).
            {
                double gw = 0, gy = 0, gz = 0, gyy = 0, gyz = 0, gzz = 0;
                for (size_t b : members) {
                    for (const Blob* blob : family.at(pool[b])->children()) {
                        const geo_point_t bc = blob->center_pos();
                        const double w = blob->npoints();
                        gw += w;
                        gy += w * bc.y();
                        gz += w * bc.z();
                        gyy += w * bc.y() * bc.y();
                        gyz += w * bc.y() * bc.z();
                        gzz += w * bc.z() * bc.z();
                    }
                }
                if (gw <= 0) continue;
                const double gmy = gy / gw, gmz = gz / gw;
                const double gcyy = gyy / gw - gmy * gmy;
                const double gcyz = gyz / gw - gmy * gmz;
                const double gczz = gzz / gw - gmz * gmz;
                const double gth = 0.5 * std::atan2(2 * gcyz, gcyy - gczz);
                const double gny = -std::sin(gth), gnz = std::cos(gth);
                double r2 = 0;
                for (size_t b : members) {
                    for (const Blob* blob : family.at(pool[b])->children()) {
                        const geo_point_t bc = blob->center_pos();
                        const double p = (bc.y() - gmy) * gny + (bc.z() - gmz) * gnz;
                        r2 += blob->npoints() * p * p;
                    }
                }
                const double union_rms = std::sqrt(r2 / gw);
                if (sep_debug())
                    std::cout << "SEPDBG mergeback group lead=m" << pool[a]
                              << " n=" << members.size()
                              << " union_rms=" << union_rms / units::cm << std::endl;
                if (union_rms > UNION_RMS_MAX) {
                    std::cout << "Separate band_merge_back: rejected group of "
                              << members.size() << " (union rms "
                              << union_rms / units::cm << " cm)" << std::endl;
                    continue;
                }
            }

            Cluster* target = family.at(pool[members.front()]);
            const auto tscope = target->get_default_scope();
            const auto ttrans = target->get_scope_transform(tscope);
            size_t nmerged = 0;
            for (size_t k = 1; k != members.size(); k++) {
                Cluster* donor = family.at(pool[members[k]]);
                target->take_children(*donor, true);
                family.at(pool[members[k]]) = nullptr;
                live_grouping.destroy_child(donor);
                assert(donor == nullptr);
                nmerged++;
            }
            // full cache clear via the scope round-trip (see collinear_recover)
            target->set_default_scope(tscope);
            target->set_scope_filter(tscope, true);
            target->set_scope_transform(tscope, ttrans);

            std::cout << "Separate band_merge_back: merged " << nmerged + 1
                      << " band pieces, len " << target->get_length() / units::cm
                      << " cm, inner/outer rms " << inner_rms / units::cm << " / "
                      << outer_rms / units::cm << " cm" << std::endl;
        }
    }

    // 6. band interior steal: a chunk of a band can stay FUSED inside a
    // non-band sibling (e.g. a crossing drift-angled track -- Separate_2's
    // relink glues them where they touch; 39324 evt 339930).  The member-level
    // merge above cannot take it without dragging the track into the band.
    // For each band-like member B, steal from each long non-band sibling H the
    // grown run of blobs that sit far off H's own main line yet touch B.
    {
        const double H_LEN_MIN        = 100 * units::cm;
        const double STEAL_OFF_MIN    = 15 * units::cm;  // blob center off H's main line
        const double STEAL_SEED_TOUCH = 5 * units::cm;   // run must touch B
        const double STEAL_GROW       = 15 * units::cm;  // growth radius within the off-line set
        const double STEAL_RUN_MIN    = 15 * units::cm;  // moved run spatial extent
        const double STEAL_X_MAX      = 20 * units::cm;  // run must live in B's x-slab: a fused
                                                         // band chunk is iso (x-narrow), while a
                                                         // drift-track run spans x (39324 339990:
                                                         // a 193 cm run spanning ~37 cm in x was
                                                         // a real track, wrongly stolen)

        // refresh band-like members: the group merges above changed them
        std::vector<size_t> bands2;
        for (size_t i = 0; i != family.size(); i++) {
            Cluster* c = family.at(i);
            if (!c) continue;
            if (c->get_length() < BAND_LEN_MIN) continue;
            const auto& pca = c->get_pca();
            if (pca.axis.size() < 1 || pca.values.size() < 2 || pca.values.at(0) <= 0) continue;
            geo_point_t axis(pca.axis.at(0).x(), pca.axis.at(0).y(), pca.axis.at(0).z());
            if (fabs(axis.angle(drift_dir) - 3.1415926 / 2.) / 3.1415926 * 180. >= BAND_PERP_ANG) continue;
            const int n = c->npoints();
            if (n < 1 || std::sqrt(pca.values.at(1) / n) < BAND_WIDTH_MIN) continue;
            bands2.push_back(i);
        }
        for (size_t bi : bands2) {
            for (size_t hi = 0; hi != family.size(); hi++) {
                if (hi == bi) continue;
                if (std::find(bands2.begin(), bands2.end(), hi) != bands2.end()) continue;
                Cluster* H = family.at(hi);
                if (!H) continue;
                if (H->get_length() < H_LEN_MIN) continue;
                const auto& hpca = H->get_pca();
                if (hpca.axis.empty() || hpca.values.at(0) <= 0) continue;
                geo_point_t haxis(hpca.axis.at(0).x(), hpca.axis.at(0).y(), hpca.axis.at(0).z());
                const geo_point_t hcen = hpca.center;
                Cluster* B = family.at(bi);

                // candidates: H blobs far off H's own main line
                const auto& blobs = H->children();
                std::vector<int> cand;
                std::vector<geo_point_t> cpos;
                for (int k = 0; k != (int) blobs.size(); k++) {
                    const geo_point_t bc = blobs.at(k)->center_pos();
                    const geo_point_t rel = bc - hcen;
                    const geo_point_t perp = rel - haxis * rel.dot(haxis);
                    if (perp.magnitude() < STEAL_OFF_MIN) continue;
                    cand.push_back(k);
                    cpos.push_back(bc);
                }
                if (cand.empty()) {
                    if (sep_debug())
                        std::cout << "SEPDBG bandsteal B=" << bi << " H=" << hi
                                  << " ncand=0" << std::endl;
                    continue;
                }
                // seeds: candidates touching B; grow within the candidate set
                std::vector<char> sel(cand.size(), 0);
                std::vector<size_t> queue;
                for (size_t k = 0; k != cand.size(); k++) {
                    geo_point_t tp(cpos[k].x(), cpos[k].y(), cpos[k].z());
                    if (B->nnearby(tp, STEAL_SEED_TOUCH) > 0) {
                        sel[k] = 1;
                        queue.push_back(k);
                    }
                }
                if (queue.empty()) {
                    if (sep_debug())
                        std::cout << "SEPDBG bandsteal B=" << bi << " H=" << hi
                                  << " ncand=" << cand.size() << " nseed=0" << std::endl;
                    continue;
                }
                while (!queue.empty()) {
                    const size_t qk = queue.back();
                    queue.pop_back();
                    for (size_t j = 0; j != cand.size(); j++) {
                        if (sel[j]) continue;
                        if ((cpos[j] - cpos[qk]).magnitude() < STEAL_GROW) {
                            sel[j] = 1;
                            queue.push_back(j);
                        }
                    }
                }
                // spatial extent of the selected run + its x stats (a genuine
                // fused band chunk lives in B's x-slab; a drift-track run does
                // not and must never be stolen)
                geo_point_t lo(1e18, 1e18, 1e18), hi_p(-1e18, -1e18, -1e18);
                size_t nsel = 0;
                double sel_sx = 0, sel_sw = 0;
                for (size_t k = 0; k != cand.size(); k++) {
                    if (!sel[k]) continue;
                    nsel++;
                    lo = geo_point_t(std::min(lo.x(), cpos[k].x()), std::min(lo.y(), cpos[k].y()),
                                     std::min(lo.z(), cpos[k].z()));
                    hi_p = geo_point_t(std::max(hi_p.x(), cpos[k].x()), std::max(hi_p.y(), cpos[k].y()),
                                       std::max(hi_p.z(), cpos[k].z()));
                    const double w = blobs.at(cand[k])->npoints();
                    sel_sx += w * cpos[k].x();
                    sel_sw += w;
                }
                const double run = (hi_p - lo).magnitude();
                const double run_x_ext = hi_p.x() - lo.x();
                double b_sx = 0, b_sw = 0;
                for (const Blob* bb : B->children()) {
                    const double w = bb->npoints();
                    b_sx += w * bb->center_pos().x();
                    b_sw += w;
                }
                const double dx_band = (sel_sw > 0 && b_sw > 0)
                    ? fabs(sel_sx / sel_sw - b_sx / b_sw) : 1e18;
                if (sep_debug())
                    std::cout << "SEPDBG bandsteal B=" << bi << " H=" << hi
                              << " ncand=" << cand.size() << " nsel=" << nsel
                              << " run=" << run / units::cm
                              << " run_x_ext=" << run_x_ext / units::cm
                              << " dx_band=" << dx_band / units::cm << std::endl;
                if (run < STEAL_RUN_MIN) continue;
                if (run_x_ext > STEAL_X_MAX || dx_band > STEAL_X_MAX) continue;

                const auto bscope = B->get_default_scope();
                const auto btrans = B->get_scope_transform(bscope);
                if (nsel == blobs.size()) {
                    B->take_children(*H, true);
                    family.at(hi) = nullptr;
                    live_grouping.destroy_child(H);
                    assert(H == nullptr);
                }
                else {
                    std::vector<int> b2groupid(H->nchildren(), -1);  // -1 = stay
                    for (size_t k = 0; k != cand.size(); k++)
                        if (sel[k]) b2groupid.at(cand[k]) = 0;
                    auto pieces = live_grouping.separate(H, b2groupid, false);
                    Cluster* piece = pieces.at(0);
                    B->take_children(*piece, true);
                    live_grouping.destroy_child(piece);
                    assert(piece == nullptr);
                }
                B->set_default_scope(bscope);
                B->set_scope_filter(bscope, true);
                B->set_scope_transform(bscope, btrans);

                std::cout << "Separate band_merge_back: band stole " << nsel
                          << " blobs (run " << run / units::cm << " cm) from sibling "
                          << hi << std::endl;
            }
        }
    }
}

// Step B (band_recarve): one cluster per crossing isochronous band.  When two
// wide isochronous bands cross at a shallow angle, the recursion above carves
// them into a path-tube, arbitrary "saved" slices and a leftover bin -- the
// final Separate_2 is pure 5 cm connectivity and the touching bands stay
// arbitrarily mixed.  Pool the touching band-like family members and re-carve
// them with a k=2 line fit in the (y,z) plane (bands are extended transverse
// to drift), assigning each blob to the nearer band axis.  Aborts (keeping the
// existing carve) when no valid seed pair exists or the fit degenerates.
static void recarve_two_bands(Grouping& live_grouping, const Tree::Scope& scope,
                              std::vector<Cluster*>& family)
{
    const double BAND_PERP_ANG = 10.;            // member axis within this of perp-to-drift (deg)
    const double BAND_LEN_MIN  = 60 * units::cm;
    const double TOUCH_MAX     = 2 * units::cm;  // members of one band complex touch
    const double SEED_ANG_MIN  = 15.;            // yz angle between the two seed axes (deg);
                                                 // crossing thin tracks that survive the width
                                                 // gate pair up at ~11 deg (PDHD 27409 evt 40900),
                                                 // genuine two-band seeds at ~25-39 deg (39324)
    const double SEED_ANG_MAX  = 45.;
    const double SEED_WIDTH_MIN = 6 * units::cm; // rms width across the axis: both seeds must be
                                                 // ribbon-like (39324 seed pair: 17.1/6.2 cm), not
                                                 // thin line-track slices (4.6 cm and below)
    const int    NITER         = 10;             // fixed k=2 line-fit iterations
    const double MIN_FRAC      = 0.10;           // degeneracy guard (npoints fraction per side)

    const geo_point_t drift_dir(1, 0, 0);

    // 1. band-like family members
    std::vector<size_t> band;  // indices into family
    for (size_t i = 0; i != family.size(); i++) {
        Cluster* c = family.at(i);
        if (!c) continue;
        if (c->get_length() < BAND_LEN_MIN) continue;
        const auto& pca = c->get_pca();
        if (pca.axis.size() < 1 || pca.values.at(0) <= 0) continue;
        geo_point_t axis(pca.axis.at(0).x(), pca.axis.at(0).y(), pca.axis.at(0).z());
        if (fabs(axis.angle(drift_dir) - 3.1415926 / 2.) / 3.1415926 * 180. >= BAND_PERP_ANG) continue;
        band.push_back(i);
    }
    if (band.size() < 2) return;

    // unit (y,z) projection of a member's main axis
    auto yz_axis = [&](size_t i) -> std::pair<double, double> {
        const auto& pca = family.at(i)->get_pca();
        const double ay = pca.axis.at(0).y(), az = pca.axis.at(0).z();
        const double n = std::sqrt(ay * ay + az * az);
        return {ay / n, az / n};
    };

    // rms width across the main axis (pca values are scatter sums, not
    // per-point variances).  Distinguishes a wide isochronous ribbon from a
    // thin line-track that merely lies near-perpendicular to drift.
    auto rms_width = [&](size_t i) -> double {
        Cluster* c = family.at(i);
        const auto& pca = c->get_pca();
        const int n = c->npoints();
        if (n < 1 || pca.values.size() < 2) return 0;
        return std::sqrt(pca.values.at(1) / n);
    };
    auto yz_angle = [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
        double dot = fabs(a.first * b.first + a.second * b.second);
        if (dot > 1.) dot = 1.;
        return std::acos(dot) / 3.1415926 * 180.;
    };


    // 2. seed pair: touching members with distinct yz axes; widest angle wins,
    // ties resolved by smaller family indices (loop order).
    int s1 = -1, s2 = -1;
    double best_ang = -1;
    for (size_t i = 0; i + 1 < band.size(); i++) {
        if (rms_width(band[i]) < SEED_WIDTH_MIN) continue;
        for (size_t j = i + 1; j < band.size(); j++) {
            if (rms_width(band[j]) < SEED_WIDTH_MIN) continue;
            const double ang = yz_angle(yz_axis(band[i]), yz_axis(band[j]));
            if (ang < SEED_ANG_MIN || ang > SEED_ANG_MAX) continue;
            std::tuple<int, int, double> dis = family.at(band[i])->get_closest_points(*family.at(band[j]));
            if (std::get<2>(dis) >= TOUCH_MAX) continue;
            if (ang > best_ang) { best_ang = ang; s1 = band[i]; s2 = band[j]; }
        }
    }
    if (s1 < 0) return;

    // 3. pool: grow from the seed pair over band-like members by touch
    std::vector<size_t> pool = {(size_t) s1, (size_t) s2};
    bool grown = true;
    while (grown) {
        grown = false;
        for (size_t bi : band) {
            if (std::find(pool.begin(), pool.end(), bi) != pool.end()) continue;
            for (size_t pi : pool) {
                std::tuple<int, int, double> dis = family.at(bi)->get_closest_points(*family.at(pi));
                if (std::get<2>(dis) < TOUCH_MAX) {
                    pool.push_back(bi);
                    grown = true;
                    break;
                }
            }
        }
    }
    std::sort(pool.begin(), pool.end());  // family (creation) order

    // 4. k=2 weighted line fit in (y,z) over pooled blob centers
    struct Line2 { double cy, cz, dy, dz; };
    Line2 line[2];
    for (int k = 0; k != 2; k++) {
        const size_t si = (k == 0) ? s1 : s2;
        const auto& pca = family.at(si)->get_pca();
        const auto a = yz_axis(si);
        line[k] = {pca.center.y(), pca.center.z(), a.first, a.second};
    }
    auto line_dis = [](const Line2& l, double by, double bz) {
        const double vy = by - l.cy, vz = bz - l.cz;
        return fabs(vy * l.dz - vz * l.dy);
    };
    auto nearer = [&](double by, double bz) {
        return line_dis(line[0], by, bz) <= line_dis(line[1], by, bz) ? 0 : 1;
    };

    for (int iter = 0; iter != NITER; iter++) {
        double sw[2] = {0, 0}, sy[2] = {0, 0}, sz[2] = {0, 0};
        double syy[2] = {0, 0}, syz[2] = {0, 0}, szz[2] = {0, 0};
        for (size_t pi : pool) {
            for (const Blob* blob : family.at(pi)->children()) {
                const geo_point_t bc = blob->center_pos();
                const double w = blob->npoints();
                const int k = nearer(bc.y(), bc.z());
                sw[k] += w;
                sy[k] += w * bc.y();
                sz[k] += w * bc.z();
                syy[k] += w * bc.y() * bc.y();
                syz[k] += w * bc.y() * bc.z();
                szz[k] += w * bc.z() * bc.z();
            }
        }
        if (sw[0] <= 0 || sw[1] <= 0) return;  // one side emptied -> no 2-band structure
        for (int k = 0; k != 2; k++) {
            const double my = sy[k] / sw[k], mz = sz[k] / sw[k];
            const double cyy = syy[k] / sw[k] - my * my;
            const double cyz = syz[k] / sw[k] - my * mz;
            const double czz = szz[k] / sw[k] - mz * mz;
            const double theta = 0.5 * std::atan2(2 * cyz, cyy - czz);
            line[k] = {my, mz, std::cos(theta), std::sin(theta)};
        }
    }

    // 5. final assignment + degeneracy guard, BEFORE any mutation
    double w_side[2] = {0, 0};
    for (size_t pi : pool)
        for (const Blob* blob : family.at(pi)->children()) {
            const geo_point_t bc = blob->center_pos();
            w_side[nearer(bc.y(), bc.z())] += blob->npoints();
        }
    const double w_tot = w_side[0] + w_side[1];
    if (w_tot <= 0 || std::min(w_side[0], w_side[1]) / w_tot < MIN_FRAC) return;

    std::cout << "Separate band_recarve: pooled " << pool.size() << " members, seed yz angle "
              << best_ang << " deg, sides " << w_side[0] << " / " << w_side[1] << " npoints" << std::endl;

    // 6. materialize: pool everything into a fresh cluster, then separate by
    // nearer band axis.  Scope triple per the junk-bin pattern below
    // (a fresh make_child cluster is invisible to later passes without it).
    Cluster& pooled = live_grouping.make_child();
    pooled.set_default_scope(scope);
    pooled.set_scope_filter(scope, true);
    pooled.set_scope_transform(scope, family.at(pool.front())->get_scope_transform(scope));
    for (size_t pi : pool) {
        Cluster* member = family.at(pi);
        pooled.take_children(*member, true);
        family.at(pi) = nullptr;
        live_grouping.destroy_child(member);
        assert(member == nullptr);
    }
    std::vector<int> b2groupid;
    b2groupid.reserve(pooled.nchildren());
    for (const Blob* blob : pooled.children()) {
        const geo_point_t bc = blob->center_pos();
        b2groupid.push_back(nearer(bc.y(), bc.z()));
    }
    Cluster* pooled_ptr = &pooled;
    auto out = live_grouping.separate(pooled_ptr, b2groupid, true);
    assert(pooled_ptr == nullptr);
    for (auto& [gid, cl] : out) family.push_back(cl);
}


// Step C (track_recarve): split one member holding two long crossing track
// arms (an "X").  Separate_1's path carve can leave both arms of two crossing
// cosmics in one piece, and the final Separate_2 (pure 5 cm connectivity)
// re-links them at the crossing, so no connectivity step can hold them apart.
// Model the member's blob centers with a k=2 3D line fit and accept the split
// only when:
//   - both fitted arms are long (axial extent) and carry a fair share of points,
//   - both arms are genuinely thin (small rms residual to their line -- a wide
//     isochronous band never fits two thin lines and is left to band_recarve),
//   - the two lines CROSS: closest approach within a few cm, located in the
//     interior of both arms.  A kinked single track (a "V") has its arms meet
//     at an end and is rejected; near-parallel diverging pairs (no crossing)
//     are also rejected and left to the carve.
static void recarve_crossing_tracks(Grouping& live_grouping, const Tree::Scope& scope,
                                    std::vector<Cluster*>& family)
{
    const double LEN_MIN       = 100 * units::cm;  // member length
    const double RATIO_MIN     = 0.01;             // eval1/eval0: an X opens the 2nd axis
    const double RATIO_MAX     = 0.8;              // beyond this it is a blob, not 2 arms
                                                   // (a wide-opening X reaches ~0.56)
    const int    NITER         = 10;               // fixed assign/refit iterations
    const double MIN_FRAC      = 0.20;             // npoints fraction per arm
    const double ARM_LEN_MIN   = 60 * units::cm;   // axial extent of each arm
    const double CROSS_DIS_MAX = 10 * units::cm;   // fitted lines must approach this close
    const double CROSS_FRAC_LO = 0.15;             // crossing interior to both arms
    const double CROSS_FRAC_HI = 0.85;
    const double RESID_MAX     = 10 * units::cm;   // npoints-weighted rms residual per arm
                                                   // (blob centers scatter several cm on real
                                                   // arms; bands are vetoed by the crossing
                                                   // guard, not by this)

    struct Line3 {
        geo_point_t c, d;  // point on line, unit direction
    };
    auto line_dis = [](const Line3& l, const geo_point_t& p) {
        geo_point_t v(p.x() - l.c.x(), p.y() - l.c.y(), p.z() - l.c.z());
        const double t = v.dot(l.d);
        geo_point_t r(v.x() - t * l.d.x(), v.y() - t * l.d.y(), v.z() - t * l.d.z());
        return r.magnitude();
    };
    // principal direction of a weighted 3x3 covariance via fixed-count power
    // iteration (deterministic; `start` seeds the iteration).
    auto principal_dir = [](const double cov[3][3], geo_point_t start) -> geo_point_t {
        double v[3] = {start.x(), start.y(), start.z()};
        double n0 = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        if (n0 <= 0) { v[0] = 1; v[1] = 0; v[2] = 0; n0 = 1; }
        for (int i = 0; i != 3; i++) v[i] /= n0;
        for (int it = 0; it != 30; it++) {
            double w[3] = {0, 0, 0};
            for (int i = 0; i != 3; i++)
                for (int j = 0; j != 3; j++) w[i] += cov[i][j] * v[j];
            const double n = std::sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
            if (n <= 0) break;
            for (int i = 0; i != 3; i++) v[i] = w[i] / n;
        }
        return geo_point_t(v[0], v[1], v[2]);
    };

    const size_t nfam = family.size();  // snapshot: do not revisit our own outputs
    for (size_t fi = 0; fi != nfam; fi++) {
        Cluster* c = family.at(fi);
        if (!c) continue;
        if (c->get_length() < LEN_MIN) continue;
        const auto& pca = c->get_pca();
        if (pca.values.size() < 2 || pca.values.at(0) <= 0) continue;
        const double ratio1 = pca.values.at(1) / pca.values.at(0);
        if (ratio1 < RATIO_MIN || ratio1 > RATIO_MAX) continue;

        geo_point_t axis0(pca.axis.at(0).x(), pca.axis.at(0).y(), pca.axis.at(0).z());
        geo_point_t axis1(pca.axis.at(1).x(), pca.axis.at(1).y(), pca.axis.at(1).z());
        const geo_point_t center = pca.center;

        // seed: split blob centers by the sign of their 2nd-axis component
        // (the X's opening direction), then iterate assign/refit.
        const auto& blobs = c->children();
        std::vector<int> side(blobs.size(), 0);
        for (size_t bi = 0; bi != blobs.size(); bi++) {
            const geo_point_t bc = blobs.at(bi)->center_pos();
            geo_point_t v(bc.x() - center.x(), bc.y() - center.y(), bc.z() - center.z());
            side[bi] = (v.dot(axis1) >= 0) ? 0 : 1;
        }

        Line3 line[2] = {{center, axis0}, {center, axis0}};
        bool degenerate = false;
        for (int iter = 0; iter != NITER + 1; iter++) {
            // refit each side: weighted mean + principal direction
            double sw[2] = {0, 0};
            double mean[2][3] = {{0, 0, 0}, {0, 0, 0}};
            for (size_t bi = 0; bi != blobs.size(); bi++) {
                const geo_point_t bc = blobs.at(bi)->center_pos();
                const double w = blobs.at(bi)->npoints();
                const int k = side[bi];
                sw[k] += w;
                mean[k][0] += w * bc.x();  mean[k][1] += w * bc.y();  mean[k][2] += w * bc.z();
            }
            if (sw[0] <= 0 || sw[1] <= 0) { degenerate = true; break; }
            double cov[2][3][3] = {};
            for (int k = 0; k != 2; k++)
                for (int i = 0; i != 3; i++) mean[k][i] /= sw[k];
            for (size_t bi = 0; bi != blobs.size(); bi++) {
                const geo_point_t bc = blobs.at(bi)->center_pos();
                const double w = blobs.at(bi)->npoints();
                const int k = side[bi];
                const double d[3] = {bc.x() - mean[k][0], bc.y() - mean[k][1], bc.z() - mean[k][2]};
                for (int i = 0; i != 3; i++)
                    for (int j = 0; j != 3; j++) cov[k][i][j] += w * d[i] * d[j];
            }
            for (int k = 0; k != 2; k++) {
                const geo_point_t dir = principal_dir(cov[k], line[k].d);
                line[k] = {geo_point_t(mean[k][0], mean[k][1], mean[k][2]), dir};
            }
            if (iter == NITER) break;
            // reassign
            for (size_t bi = 0; bi != blobs.size(); bi++) {
                const geo_point_t bc = blobs.at(bi)->center_pos();
                side[bi] = (line_dis(line[0], bc) <= line_dis(line[1], bc)) ? 0 : 1;
            }
        }
        if (degenerate) continue;

        // validation, all BEFORE any mutation
        double w_side[2] = {0, 0}, resid2[2] = {0, 0};
        double tmin[2] = {1e18, 1e18}, tmax[2] = {-1e18, -1e18};
        for (size_t bi = 0; bi != blobs.size(); bi++) {
            const geo_point_t bc = blobs.at(bi)->center_pos();
            const double w = blobs.at(bi)->npoints();
            const int k = side[bi];
            w_side[k] += w;
            const double dis = line_dis(line[k], bc);
            resid2[k] += w * dis * dis;
            geo_point_t v(bc.x() - line[k].c.x(), bc.y() - line[k].c.y(), bc.z() - line[k].c.z());
            const double t = v.dot(line[k].d);
            if (t < tmin[k]) tmin[k] = t;
            if (t > tmax[k]) tmax[k] = t;
        }
        const double w_tot = w_side[0] + w_side[1];
        if (w_tot <= 0 || std::min(w_side[0], w_side[1]) / w_tot < MIN_FRAC) continue;
        bool ok = true;
        for (int k = 0; k != 2; k++) {
            if (tmax[k] - tmin[k] < ARM_LEN_MIN) ok = false;
            if (std::sqrt(resid2[k] / w_side[k]) > RESID_MAX) ok = false;
        }
        if (!ok) continue;

        // crossing: closest approach of the two fitted (infinite) lines must be
        // close and interior to both arms.
        const geo_point_t w0(line[0].c.x() - line[1].c.x(), line[0].c.y() - line[1].c.y(),
                             line[0].c.z() - line[1].c.z());
        const double b = line[0].d.dot(line[1].d);
        const double denom = 1. - b * b;
        if (denom < 1e-4) continue;  // near-parallel: no crossing
        const double d0 = line[0].d.dot(w0), d1 = line[1].d.dot(w0);
        const double t0 = (b * d1 - d0) / denom;
        const double t1 = (d1 - b * d0) / denom;
        geo_point_t p0(line[0].c.x() + t0 * line[0].d.x(), line[0].c.y() + t0 * line[0].d.y(),
                       line[0].c.z() + t0 * line[0].d.z());
        geo_point_t p1(line[1].c.x() + t1 * line[1].d.x(), line[1].c.y() + t1 * line[1].d.y(),
                       line[1].c.z() + t1 * line[1].d.z());
        geo_point_t gap(p0.x() - p1.x(), p0.y() - p1.y(), p0.z() - p1.z());
        const double f0 = (t0 - tmin[0]) / (tmax[0] - tmin[0]);
        const double f1 = (t1 - tmin[1]) / (tmax[1] - tmin[1]);
        if (gap.magnitude() > CROSS_DIS_MAX) continue;
        // Kink veto: a kinked single track (a "V") meets at the END of both
        // arms.  Require the crossing to be interior to at least one arm: an
        // "X" (both interior) or a "T" (one track ends on the other, which it
        // crosses mid-length) are two distinct tracks and split; a "V" is not.
        const bool int0 = (f0 >= CROSS_FRAC_LO && f0 <= CROSS_FRAC_HI);
        const bool int1 = (f1 >= CROSS_FRAC_LO && f1 <= CROSS_FRAC_HI);
        if (!int0 && !int1) continue;
        // The non-interior arm's crossing must still lie within the arm (with
        // slack), else the two lines never actually meet inside the member.
        if (f0 < -0.1 || f0 > 1.1 || f1 < -0.1 || f1 > 1.1) continue;

        std::cout << "Separate track_recarve: len " << c->get_length() / units::cm
                  << " cm split into arms " << w_side[0] << " / " << w_side[1]
                  << " npoints, cross frac " << f0 << " / " << f1
                  << ", resid " << std::sqrt(resid2[0] / w_side[0]) / units::cm << " / "
                  << std::sqrt(resid2[1] / w_side[1]) / units::cm << " cm" << std::endl;

        std::vector<int> b2groupid(side.begin(), side.end());
        auto out = live_grouping.separate(c, b2groupid, true);
        assert(c == nullptr);
        family.at(fi) = nullptr;
        for (auto& [gid, cl] : out) family.push_back(cl);
    }
}


// Step A2 (track_repartition): fix an interior segment swap between two
// crossing thin tracks ALREADY split into two family members.  At a crossing
// the carve + Separate_2 relink can fuse a mid-track chunk of track B into
// track A's cluster (PDVD 39324 evt 339990 group0123: a 48 cm chunk of one
// track, 26 cm off the host's axis, embedded in the other's cluster).
// collinear_interior cannot reach it -- it moves whole sibling fragments and
// this chunk is fused into the host.  Pool the pair's blobs, refit a k=2 3D
// line model seeded by the current membership (machinery and validation
// duplicated from recarve_crossing_tracks above) and reassign by nearer line.
// A surgical no-op guard skips the pair when the would-be moves are only
// crossing-region dribble: mutate only if the moved blobs include a run far
// from the host line with substantial axial extent.
static void repartition_crossing_tracks(Grouping& live_grouping, const Tree::Scope& scope,
                                        std::vector<Cluster*>& family)
{
    const double LEN_MIN        = 100 * units::cm;  // each member
    const double THIN_RATIO_MAX = 0.15;             // pca eval1/eval0 per member
    const int    NITER          = 10;
    const double MIN_FRAC       = 0.20;
    const double ARM_LEN_MIN    = 60 * units::cm;
    const double CROSS_DIS_MAX  = 10 * units::cm;
    const double CROSS_FRAC_LO  = 0.15;
    const double CROSS_FRAC_HI  = 0.85;
    const double RESID_MAX      = 10 * units::cm;
    const double MOVE_PERP_MIN  = 15 * units::cm;   // moved run: distance from host line
    const double MOVE_RUN_MIN   = 15 * units::cm;   // ... and axial extent (the 339990
                                                    // chunk reads 48 cm at 20-40 cm)

    struct Line3 {
        geo_point_t c, d;  // point on line, unit direction
    };
    auto line_dis = [](const Line3& l, const geo_point_t& p) {
        geo_point_t v(p.x() - l.c.x(), p.y() - l.c.y(), p.z() - l.c.z());
        const double t = v.dot(l.d);
        geo_point_t r(v.x() - t * l.d.x(), v.y() - t * l.d.y(), v.z() - t * l.d.z());
        return r.magnitude();
    };
    auto principal_dir = [](const double cov[3][3], geo_point_t start) -> geo_point_t {
        double v[3] = {start.x(), start.y(), start.z()};
        double n0 = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        if (n0 <= 0) { v[0] = 1; v[1] = 0; v[2] = 0; n0 = 1; }
        for (int i = 0; i != 3; i++) v[i] /= n0;
        for (int it = 0; it != 30; it++) {
            double w[3] = {0, 0, 0};
            for (int i = 0; i != 3; i++)
                for (int j = 0; j != 3; j++) w[i] += cov[i][j] * v[j];
            const double n = std::sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
            if (n <= 0) break;
            for (int i = 0; i != 3; i++) v[i] = w[i] / n;
        }
        return geo_point_t(v[0], v[1], v[2]);
    };

    const size_t nfam = family.size();  // snapshot: do not revisit our own outputs
    for (size_t fi = 0; fi + 1 < nfam; fi++) {
        Cluster* a = family.at(fi);
        if (!a) continue;
        if (a->get_length() < LEN_MIN) continue;
        const auto& pca_a = a->get_pca();
        if (pca_a.values.size() < 2 || pca_a.values.at(0) <= 0) continue;
        if (pca_a.values.at(1) / pca_a.values.at(0) >= THIN_RATIO_MAX) continue;

        for (size_t fj = fi + 1; fj < nfam; fj++) {
            Cluster* b = family.at(fj);
            if (!b) continue;
            if (b->get_length() < LEN_MIN) continue;
            const auto& pca_b = b->get_pca();
            if (pca_b.values.size() < 2 || pca_b.values.at(0) <= 0) continue;
            if (pca_b.values.at(1) / pca_b.values.at(0) >= THIN_RATIO_MAX) continue;

            // the members' pca main lines must genuinely cross
            Line3 host[2] = {
                {pca_a.center, geo_point_t(pca_a.axis.at(0).x(), pca_a.axis.at(0).y(), pca_a.axis.at(0).z())},
                {pca_b.center, geo_point_t(pca_b.axis.at(0).x(), pca_b.axis.at(0).y(), pca_b.axis.at(0).z())}};
            const geo_point_t w0(host[0].c.x() - host[1].c.x(), host[0].c.y() - host[1].c.y(),
                                 host[0].c.z() - host[1].c.z());
            const double bb = host[0].d.dot(host[1].d);
            const double denom = 1. - bb * bb;
            if (denom < 1e-4) continue;  // near-parallel
            const double d0 = host[0].d.dot(w0), d1 = host[1].d.dot(w0);
            const double tc0 = (bb * d1 - d0) / denom;
            const double tc1 = (d1 - bb * d0) / denom;
            geo_point_t p0(host[0].c.x() + tc0 * host[0].d.x(), host[0].c.y() + tc0 * host[0].d.y(),
                           host[0].c.z() + tc0 * host[0].d.z());
            geo_point_t p1(host[1].c.x() + tc1 * host[1].d.x(), host[1].c.y() + tc1 * host[1].d.y(),
                           host[1].c.z() + tc1 * host[1].d.z());
            geo_point_t gap(p0.x() - p1.x(), p0.y() - p1.y(), p0.z() - p1.z());
            if (gap.magnitude() > CROSS_DIS_MAX) continue;
            // axial extents from blob centers, crossing interior to >=1 member
            double hmin[2] = {1e18, 1e18}, hmax[2] = {-1e18, -1e18};
            Cluster* mem[2] = {a, b};
            for (int k = 0; k != 2; k++)
                for (const Blob* blob : mem[k]->children()) {
                    const geo_point_t bc = blob->center_pos();
                    geo_point_t v(bc.x() - host[k].c.x(), bc.y() - host[k].c.y(), bc.z() - host[k].c.z());
                    const double t = v.dot(host[k].d);
                    if (t < hmin[k]) hmin[k] = t;
                    if (t > hmax[k]) hmax[k] = t;
                }
            const double cf0 = (tc0 - hmin[0]) / (hmax[0] - hmin[0]);
            const double cf1 = (tc1 - hmin[1]) / (hmax[1] - hmin[1]);
            const bool int0 = (cf0 >= CROSS_FRAC_LO && cf0 <= CROSS_FRAC_HI);
            const bool int1 = (cf1 >= CROSS_FRAC_LO && cf1 <= CROSS_FRAC_HI);
            if (!int0 && !int1) continue;
            if (cf0 < -0.1 || cf0 > 1.1 || cf1 < -0.1 || cf1 > 1.1) continue;

            // pool the pair's blobs; side seeded by current membership
            std::vector<const Blob*> blobs;
            std::vector<int> side, member_of;
            for (int k = 0; k != 2; k++)
                for (const Blob* blob : mem[k]->children()) {
                    blobs.push_back(blob);
                    side.push_back(k);
                    member_of.push_back(k);
                }

            Line3 line[2] = {host[0], host[1]};
            bool degenerate = false;
            for (int iter = 0; iter != NITER + 1; iter++) {
                double sw[2] = {0, 0};
                double mean[2][3] = {{0, 0, 0}, {0, 0, 0}};
                for (size_t bi = 0; bi != blobs.size(); bi++) {
                    const geo_point_t bc = blobs.at(bi)->center_pos();
                    const double w = blobs.at(bi)->npoints();
                    const int k = side[bi];
                    sw[k] += w;
                    mean[k][0] += w * bc.x();  mean[k][1] += w * bc.y();  mean[k][2] += w * bc.z();
                }
                if (sw[0] <= 0 || sw[1] <= 0) { degenerate = true; break; }
                double cov[2][3][3] = {};
                for (int k = 0; k != 2; k++)
                    for (int i = 0; i != 3; i++) mean[k][i] /= sw[k];
                for (size_t bi = 0; bi != blobs.size(); bi++) {
                    const geo_point_t bc = blobs.at(bi)->center_pos();
                    const double w = blobs.at(bi)->npoints();
                    const int k = side[bi];
                    const double d[3] = {bc.x() - mean[k][0], bc.y() - mean[k][1], bc.z() - mean[k][2]};
                    for (int i = 0; i != 3; i++)
                        for (int j = 0; j != 3; j++) cov[k][i][j] += w * d[i] * d[j];
                }
                for (int k = 0; k != 2; k++) {
                    const geo_point_t dir = principal_dir(cov[k], line[k].d);
                    line[k] = {geo_point_t(mean[k][0], mean[k][1], mean[k][2]), dir};
                }
                if (iter == NITER) break;
                for (size_t bi = 0; bi != blobs.size(); bi++) {
                    const geo_point_t bc = blobs.at(bi)->center_pos();
                    side[bi] = (line_dis(line[0], bc) <= line_dis(line[1], bc)) ? 0 : 1;
                }
            }
            if (degenerate) continue;

            // final assignment from the final fitted lines -- also recomputed
            // per pooled child at materialization below, so it must be a pure
            // function of the lines (not of blob order or pointer identity)
            for (size_t bi = 0; bi != blobs.size(); bi++) {
                const geo_point_t bc = blobs.at(bi)->center_pos();
                side[bi] = (line_dis(line[0], bc) <= line_dis(line[1], bc)) ? 0 : 1;
            }

            // validation (same gates as recarve_crossing_tracks), before mutation
            double w_side[2] = {0, 0}, resid2[2] = {0, 0};
            double tmin[2] = {1e18, 1e18}, tmax[2] = {-1e18, -1e18};
            for (size_t bi = 0; bi != blobs.size(); bi++) {
                const geo_point_t bc = blobs.at(bi)->center_pos();
                const double w = blobs.at(bi)->npoints();
                const int k = side[bi];
                w_side[k] += w;
                const double dis = line_dis(line[k], bc);
                resid2[k] += w * dis * dis;
                geo_point_t v(bc.x() - line[k].c.x(), bc.y() - line[k].c.y(), bc.z() - line[k].c.z());
                const double t = v.dot(line[k].d);
                if (t < tmin[k]) tmin[k] = t;
                if (t > tmax[k]) tmax[k] = t;
            }
            const double w_tot = w_side[0] + w_side[1];
            if (w_tot <= 0 || std::min(w_side[0], w_side[1]) / w_tot < MIN_FRAC) continue;
            bool ok = true;
            for (int k = 0; k != 2; k++) {
                if (tmax[k] - tmin[k] < ARM_LEN_MIN) ok = false;
                if (std::sqrt(resid2[k] / w_side[k]) > RESID_MAX) ok = false;
            }
            if (!ok) continue;

            // surgical no-op guard: among the blobs that change membership,
            // require a run sitting far from its HOST's line with substantial
            // axial extent along its NEW line.  Crossing-region dribble (close
            // to both lines by construction) never qualifies.
            double run_tmin[2] = {1e18, 1e18}, run_tmax[2] = {-1e18, -1e18};
            size_t nmoved = 0;
            for (size_t bi = 0; bi != blobs.size(); bi++) {
                if (side[bi] == member_of[bi]) continue;
                nmoved++;
                const geo_point_t bc = blobs.at(bi)->center_pos();
                if (line_dis(host[member_of[bi]], bc) < MOVE_PERP_MIN) continue;
                const int k = side[bi];
                geo_point_t v(bc.x() - line[k].c.x(), bc.y() - line[k].c.y(), bc.z() - line[k].c.z());
                const double t = v.dot(line[k].d);
                if (t < run_tmin[k]) run_tmin[k] = t;
                if (t > run_tmax[k]) run_tmax[k] = t;
            }
            const bool substantial = (run_tmax[0] - run_tmin[0] >= MOVE_RUN_MIN) ||
                                     (run_tmax[1] - run_tmin[1] >= MOVE_RUN_MIN);
            if (sep_debug())
                std::cout << "SEPDBG repart pair " << fi << "/" << fj
                          << " cross_dis=" << gap.magnitude() / units::cm
                          << " cf=" << cf0 << "/" << cf1 << " nmoved=" << nmoved
                          << " run=" << (run_tmax[0] - run_tmin[0]) / units::cm << "/"
                          << (run_tmax[1] - run_tmin[1]) / units::cm
                          << " substantial=" << substantial << std::endl;
            if (nmoved == 0 || !substantial) continue;

            std::cout << "Separate track_repartition: pair lens "
                      << a->get_length() / units::cm << " / " << b->get_length() / units::cm
                      << " cm, moved " << nmoved << " blobs, sides "
                      << w_side[0] << " / " << w_side[1] << " npoints" << std::endl;

            // materialize: pool into a fresh cluster, separate by final side
            Cluster& pooled = live_grouping.make_child();
            pooled.set_default_scope(scope);
            pooled.set_scope_filter(scope, true);
            pooled.set_scope_transform(scope, a->get_scope_transform(scope));
            pooled.take_children(*a, true);
            family.at(fi) = nullptr;
            live_grouping.destroy_child(a);
            assert(a == nullptr);
            pooled.take_children(*b, true);
            family.at(fj) = nullptr;
            live_grouping.destroy_child(b);
            assert(b == nullptr);
            std::vector<int> b2groupid;
            b2groupid.reserve(pooled.nchildren());
            for (const Blob* blob : pooled.children()) {
                const geo_point_t bc = blob->center_pos();
                b2groupid.push_back((line_dis(line[0], bc) <= line_dis(line[1], bc)) ? 0 : 1);
            }
            Cluster* pooled_ptr = &pooled;
            auto out = live_grouping.separate(pooled_ptr, b2groupid, true);
            assert(pooled_ptr == nullptr);
            for (auto& [gid, cl] : out) family.push_back(cl);
            break;  // fi consumed; next fi
        }
    }
}

// Step C (iso_slab_split): split a member that mixes isochronous bands with
// drift-direction tracks.  An iso band lives at ONE drift time -- a narrow,
// dense x-slab -- while a drift-direction track spans x and can touch several
// bands; the final Separate_2 (pure 5 cm connectivity) then chains
// band-track-band into one cluster FOREVER: no y-z-projected mechanism can
// cut it (the projection is x-blind, and a drift track nearly vanishes in
// y-z).  PDVD 39324 evt 339890: a 44.7k-pt member holding TWO bands at
// x-slabs 50 cm apart plus three drift tracks; evt 339990: one band plus a
// drift track fused by the carve.  Partition by x-slab occupancy:
//  1. dense narrow x-slabs (the bands) from a weighted blob-center histogram;
//  2. out-of-slab blobs form connected components; components with real x
//     extent (the drift signature -- band shoulders are x-narrow) seed
//     tracks, collinear seeds joined across the slabs they pierce;
//  3. each blob: track component membership, else near a track line (the
//     track's continuation THROUGH a band), else the band of its slab, else
//     the nearest slab in x.
// Fires only when >=1 valid slab AND >=1 track exist: a pure band, a pure
// track, an isochronous track (slab, no crossers) and a band plus dribble
// are all structural no-ops.
static void split_iso_slabs(Grouping& live_grouping, const Tree::Scope& scope,
                            std::vector<Cluster*>& family)
{
    const double BIN_W          = 4 * units::cm;   // x histogram bin
    const double SLAB_DENS_FACT = 5.;              // dense bin: > this x the uniform density
    const double SLAB_W_MAX     = 25 * units::cm;  // a slab is one drift time: narrow
    const double SLAB_FRAC_MIN  = 0.20;            // of member npoints
    const double SLAB_YZ_MIN    = 60 * units::cm;  // band-sized in (y,z)
    const double TRACK_GROW     = 10 * units::cm;  // component growth radius
    const double TRACK_EXT_MIN  = 30 * units::cm;  // component 3D extent
    const double TRACK_XEXT_MIN = 20 * units::cm;  // must span x (band shoulders read <=12)
    const double JOIN_ANG_MAX   = 15.;             // collinear component joining (deg)
    const double JOIN_OFF_MAX   = 10 * units::cm;
    const double TRACK_NEAR     = 6 * units::cm;   // in-slab blob assignment to a track line
    const double MEMBER_LEN_MIN = 100 * units::cm;

    struct Line3 {
        geo_point_t c, d;
    };
    auto line_dis = [](const Line3& l, const geo_point_t& p) {
        geo_point_t v(p.x() - l.c.x(), p.y() - l.c.y(), p.z() - l.c.z());
        const double t = v.dot(l.d);
        geo_point_t r(v.x() - t * l.d.x(), v.y() - t * l.d.y(), v.z() - t * l.d.z());
        return r.magnitude();
    };
    auto power = [](const double m[3][3], const double s[3], double out[3]) {
        double v[3] = {s[0], s[1], s[2]};
        double n0 = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        if (n0 <= 0) { v[0] = 1; v[1] = 0; v[2] = 0; n0 = 1; }
        for (int i = 0; i != 3; i++) v[i] /= n0;
        for (int it = 0; it != 30; it++) {
            double w[3] = {0, 0, 0};
            for (int i = 0; i != 3; i++)
                for (int j = 0; j != 3; j++) w[i] += m[i][j] * v[j];
            const double n = std::sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
            if (n <= 0) break;
            for (int i = 0; i != 3; i++) v[i] = w[i] / n;
        }
        for (int i = 0; i != 3; i++) out[i] = v[i];
    };
    // weighted 3D line fit of a blob subset
    auto fit_line = [&power](const std::vector<geo_point_t>& pos, const std::vector<double>& wt,
                             const std::vector<size_t>& idx) -> Line3 {
        double sw = 0, sm[3] = {0, 0, 0};
        for (size_t i : idx) {
            sw += wt[i];
            sm[0] += wt[i] * pos[i].x();
            sm[1] += wt[i] * pos[i].y();
            sm[2] += wt[i] * pos[i].z();
        }
        for (int i = 0; i != 3; i++) sm[i] /= sw;
        double cov[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        for (size_t i : idx) {
            const double d[3] = {pos[i].x() - sm[0], pos[i].y() - sm[1], pos[i].z() - sm[2]};
            for (int a = 0; a != 3; a++)
                for (int b = 0; b != 3; b++) cov[a][b] += wt[i] * d[a] * d[b];
        }
        const double seed[3] = {1, 0, 0};
        double e0[3];
        power(cov, seed, e0);
        return {geo_point_t(sm[0], sm[1], sm[2]), geo_point_t(e0[0], e0[1], e0[2])};
    };

    const size_t nfam = family.size();  // snapshot: do not revisit our own outputs
    for (size_t fi = 0; fi != nfam; fi++) {
        Cluster* m = family.at(fi);
        if (!m) continue;
        if (m->get_length() < MEMBER_LEN_MIN) continue;
        const auto& blobs = m->children();
        if (blobs.size() < 20) continue;

        std::vector<geo_point_t> pos(blobs.size());
        std::vector<double> wt(blobs.size());
        double xlo = 1e18, xhi = -1e18, w_tot = 0;
        for (size_t i = 0; i != blobs.size(); i++) {
            pos[i] = blobs.at(i)->center_pos();
            wt[i] = blobs.at(i)->npoints();
            w_tot += wt[i];
            xlo = std::min(xlo, pos[i].x());
            xhi = std::max(xhi, pos[i].x());
        }
        const int nbins = std::max(1, (int) std::ceil((xhi - xlo) / BIN_W));
        if (nbins < 8 || w_tot <= 0) continue;  // no x structure to find
        std::vector<double> binw(nbins, 0);
        auto xbin = [&](double x) {
            int b = (int) ((x - xlo) / BIN_W);
            return std::min(std::max(b, 0), nbins - 1);
        };
        for (size_t i = 0; i != blobs.size(); i++) binw[xbin(pos[i].x())] += wt[i];
        const double uniform = w_tot / nbins;

        // 1. slabs: contiguous dense-bin runs, narrow, heavy and band-sized
        struct Slab {
            double xmin, xmax, w;
        };
        std::vector<Slab> slabs;
        for (int b = 0; b != nbins;) {
            if (binw[b] <= SLAB_DENS_FACT * uniform) { b++; continue; }
            int e = b;
            double w = 0;
            while (e != nbins && binw[e] > SLAB_DENS_FACT * uniform) { w += binw[e]; e++; }
            const double sxmin = xlo + b * BIN_W, sxmax = xlo + e * BIN_W;
            bool ok = (sxmax - sxmin <= SLAB_W_MAX) && (w >= SLAB_FRAC_MIN * w_tot);
            if (ok) {  // band-sized in (y,z)?
                double ylo = 1e18, yhi = -1e18, zlo = 1e18, zhi = -1e18;
                for (size_t i = 0; i != blobs.size(); i++) {
                    if (pos[i].x() < sxmin || pos[i].x() >= sxmax) continue;
                    ylo = std::min(ylo, pos[i].y());
                    yhi = std::max(yhi, pos[i].y());
                    zlo = std::min(zlo, pos[i].z());
                    zhi = std::max(zhi, pos[i].z());
                }
                const double yz = std::hypot(yhi - ylo, zhi - zlo);
                ok = yz >= SLAB_YZ_MIN;
            }
            if (ok) slabs.push_back({sxmin, sxmax, w});
            b = e;
        }
        if (sep_debug() && !slabs.empty()) {
            std::cout << "SEPDBG slabsplit m" << fi << " nslab=" << slabs.size();
            for (const auto& s : slabs)
                std::cout << " [" << s.xmin / units::cm << "," << s.xmax / units::cm
                          << ")w=" << s.w;
            std::cout << " w_tot=" << w_tot << std::endl;
        }
        if (slabs.empty()) continue;
        auto in_slab = [&](double x) {
            for (size_t s = 0; s != slabs.size(); s++)
                if (x >= slabs[s].xmin && x < slabs[s].xmax) return (int) s;
            return -1;
        };

        // 2. out-of-slab connected components -> track seeds -> collinear join
        std::vector<size_t> outs;
        for (size_t i = 0; i != blobs.size(); i++)
            if (in_slab(pos[i].x()) < 0) outs.push_back(i);
        std::vector<int> comp(outs.size(), -1);
        int ncomp = 0;
        for (size_t i = 0; i != outs.size(); i++) {
            if (comp[i] >= 0) continue;
            comp[i] = ncomp;
            std::vector<size_t> queue = {i};
            while (!queue.empty()) {
                const size_t q = queue.back();
                queue.pop_back();
                for (size_t j = 0; j != outs.size(); j++) {
                    if (comp[j] >= 0) continue;
                    if ((pos[outs[j]] - pos[outs[q]]).magnitude() < TRACK_GROW) {
                        comp[j] = ncomp;
                        queue.push_back(j);
                    }
                }
            }
            ncomp++;
        }
        // track seeds: components with 3D extent AND x extent (drift signature)
        std::vector<std::vector<size_t>> seed_blobs;  // blob indices per seed
        std::vector<Line3> seed_line;
        for (int c = 0; c != ncomp; c++) {
            std::vector<size_t> idx;
            geo_point_t lo(1e18, 1e18, 1e18), hi(-1e18, -1e18, -1e18);
            for (size_t i = 0; i != outs.size(); i++) {
                if (comp[i] != c) continue;
                idx.push_back(outs[i]);
                const geo_point_t& p = pos[outs[i]];
                lo = geo_point_t(std::min(lo.x(), p.x()), std::min(lo.y(), p.y()),
                                 std::min(lo.z(), p.z()));
                hi = geo_point_t(std::max(hi.x(), p.x()), std::max(hi.y(), p.y()),
                                 std::max(hi.z(), p.z()));
            }
            if (idx.empty()) continue;
            const double ext = (hi - lo).magnitude();
            const double xext = hi.x() - lo.x();
            if (sep_debug())
                std::cout << "SEPDBG slabsplit m" << fi << " comp" << c
                          << " n=" << idx.size() << " ext=" << ext / units::cm
                          << " xext=" << xext / units::cm << std::endl;
            if (ext < TRACK_EXT_MIN || xext < TRACK_XEXT_MIN) continue;
            seed_blobs.push_back(idx);
            seed_line.push_back(fit_line(pos, wt, idx));
        }
        if (seed_blobs.empty()) continue;
        // join collinear seeds (a track pierced by a slab reappears beyond it)
        std::vector<int> sgrp(seed_blobs.size());
        for (size_t k = 0; k != sgrp.size(); k++) sgrp[k] = (int) k;
        auto sroot = [&](int k) {
            while (sgrp[k] != k) k = sgrp[k] = sgrp[sgrp[k]];
            return k;
        };
        for (size_t a = 0; a + 1 < seed_blobs.size(); a++)
            for (size_t b = a + 1; b < seed_blobs.size(); b++) {
                const double dot = fabs(seed_line[a].d.dot(seed_line[b].d));
                const double ang = std::acos(std::min(1., dot)) / 3.1415926 * 180.;
                if (ang > JOIN_ANG_MAX) continue;
                if (line_dis(seed_line[a], seed_line[b].c) > JOIN_OFF_MAX) continue;
                if (line_dis(seed_line[b], seed_line[a].c) > JOIN_OFF_MAX) continue;
                const int ra = sroot(a), rb = sroot(b);
                if (ra != rb) sgrp[std::max(ra, rb)] = std::min(ra, rb);
            }
        // track groups: refit one line per group over its joined seed blobs
        std::map<int, std::vector<size_t>> tgroups;  // root -> blob indices
        for (size_t k = 0; k != seed_blobs.size(); k++)
            for (size_t i : seed_blobs[k]) tgroups[sroot(k)].push_back(i);
        std::vector<std::vector<size_t>> track_blobs;
        std::vector<Line3> track_line;
        for (auto& [r, idx] : tgroups) {
            track_blobs.push_back(idx);
            track_line.push_back(fit_line(pos, wt, idx));
        }

        // 3. final assignment: groups 0..T-1 = tracks, T..T+S-1 = slab bands
        const int T = (int) track_blobs.size();
        std::vector<int> b2g(blobs.size(), -1);
        for (int t = 0; t != T; t++)
            for (size_t i : track_blobs[t]) b2g[i] = t;
        for (size_t i = 0; i != blobs.size(); i++) {
            if (b2g[i] >= 0) continue;
            double dmin = 1e18;
            int tbest = -1;
            for (int t = 0; t != T; t++) {
                const double d = line_dis(track_line[t], pos[i]);
                if (d < dmin) { dmin = d; tbest = t; }
            }
            if (dmin < TRACK_NEAR) { b2g[i] = tbest; continue; }
            int s = in_slab(pos[i].x());
            if (s < 0) {  // leftover debris: nearest slab in x
                double dx = 1e18;
                for (size_t k = 0; k != slabs.size(); k++) {
                    const double c = 0.5 * (slabs[k].xmin + slabs[k].xmax);
                    if (fabs(pos[i].x() - c) < dx) { dx = fabs(pos[i].x() - c); s = (int) k; }
                }
            }
            b2g[i] = T + s;
        }
        // degenerate: everything in one group
        bool multi = false;
        for (size_t i = 1; i != blobs.size(); i++)
            if (b2g[i] != b2g[0]) { multi = true; break; }
        if (!multi) continue;

        std::cout << "Separate iso_slab_split: member len " << m->get_length() / units::cm
                  << " cm -> " << T << " track(s) + " << slabs.size() << " band(s)" << std::endl;

        std::vector<int> b2groupid(b2g.begin(), b2g.end());
        Cluster* m_ptr = m;
        auto out = live_grouping.separate(m_ptr, b2groupid, true);
        assert(m_ptr == nullptr);
        family.at(fi) = nullptr;
        for (auto& [gid, cl] : out) family.push_back(cl);
    }
}


static void clustering_separate(Grouping& live_grouping,
                                IDetectorVolumes::pointer dv,
                                IPCTransformSet::pointer pcts,
                                const Tree::Scope& scope,
                                const bool use_ctpc,
                                const int max_hull_points,
                                const bool sbnd_boundary_tag,
                                const bool collinear_recover,
                                const bool collinear_interior,
                                const bool collinear_member_merge,
                                const bool track_repartition,
                                const bool band_merge_back,
                                const bool band_recarve,
                                const bool drift_side_fv_x,
                                const double far_point_x_cut,
                                const double far_point_mid_dis,
                                const bool track_recarve,
                                const double dec1_guard_main_angle,
                                const bool iso_slab_split,
                                const bool tag_family,
                                const bool collinear_global_merge);

class ClusteringSeparate : public IConfigurable, public Clus::IEnsembleVisitor, private NeedDV, private NeedPCTS, private NeedScope {
public:
    ClusteringSeparate() {}
    virtual ~ClusteringSeparate() {}

    void configure(const WireCell::Configuration& config) {
        NeedDV::configure(config);
        NeedPCTS::configure(config);
        NeedScope::configure(config);
        
        use_ctpc_ = get(config, "use_ctpc", true);
        // Cap on points for the per-cluster convex hull used by the separation
        // decision. <0 (default) falls back to Constants::MaxHullPoints (10000),
        // preserving existing behavior bit-for-bit; raise it (e.g. SBND) to let
        // large full-detector overclusters be considered for separation.
        max_hull_points_ = get(config, "max_hull_points", -1);
        // SBND-only two-track upstream-boundary tag (default OFF => bit-identical
        // to prior behavior).  See JudgeSeparateDec_SBND_boundary.
        sbnd_boundary_tag_ = get(config, "sbnd_boundary_tag", false);
        // Post-separation refinements (default OFF => bit-identical to prior
        // behavior).  See recover_collinear_tips / recarve_two_bands.
        collinear_recover_ = get(config, "collinear_recover", false);
        // Interior-bite reclaim extension of collinear_recover (default OFF =>
        // bit-identical); only effective when collinear_recover is also on.
        // See recover_collinear_tips(interior_reclaim).
        collinear_interior_ = get(config, "collinear_interior", false);
        // Rejoin a single straight track the carve cut into long thin pieces
        // (default OFF => bit-identical).  See merge_collinear_members.
        collinear_member_merge_ = get(config, "collinear_member_merge", false);
        // Pairwise k=2 3D repartition of two crossing thin-track members
        // (default OFF => bit-identical).  See repartition_crossing_tracks.
        track_repartition_ = get(config, "track_repartition", false);
        // Re-assemble a single isochronous band hatched into interleaved
        // pieces by the carve (default OFF => bit-identical).  See
        // merge_back_bands.
        band_merge_back_ = get(config, "band_merge_back", false);
        band_recarve_ = get(config, "band_recarve", false);
        // Drift-side FV x-range for multi-APA common-face scopes (default OFF =>
        // bit-identical to prior behavior).  See select_scope_fv(common_face_x).
        drift_side_fv_x_ = get(config, "drift_side_fv_x", false);
        // Drift-x deviation that promotes a boundary point to a "far" point in
        // JudgeSeparateDec_2's two-endpoint test.  Default 140 cm == the prototype
        // expression (effectively dead); set 14 cm for the evidently intended cut.
        far_point_x_cut_ = get(config, "far_point_x_cut", 140 * units::cm);
        // Midpoint-to-cluster cap in the same far-point test (default 25 cm ==
        // prototype).  Raise it so diverging track pairs keep their far points.
        far_point_mid_dis_ = get(config, "far_point_mid_dis", 25 * units::cm);
        // Post-separation k=2 3D-line self-split of a member holding two long
        // crossing track arms (default OFF => bit-identical).
        track_recarve_ = get(config, "track_recarve", false);
        // Dec_1 drift-aligned guard: require the MAIN axis within this angle
        // (deg) of drift for the guard to apply (default <0 = unconditional,
        // bit-identical).  See JudgeSeparateDec_1.
        dec1_guard_main_angle_ = get(config, "dec1_guard_main_angle", -1.0);
        // x-slab-aware split of a member mixing isochronous bands with
        // drift-direction tracks (default OFF => bit-identical).  See
        // split_iso_slabs.
        iso_slab_split_ = get(config, "iso_slab_split", false);
        // Stamp every final family member with a "sep_family" cluster scalar
        // (the parent cluster's ident) so a later pass can decline to undo the
        // split (see ClusteringConnect1 respect_separate_family).  Default OFF
        // => no scalar is written, bit-identical.
        tag_family_ = get(config, "tag_family", false);
        // Grouping-wide end-to-end stitch of two long thin collinear clusters
        // via merge_collinear_members (default OFF => bit-identical).
        collinear_global_merge_ = get(config, "collinear_global_merge", false);
    }

    void visit(Ensemble& ensemble) const {
        auto& live = *ensemble.with_name("live").at(0);
        clustering_separate(live, m_dv, m_pcts, m_scope, use_ctpc_, max_hull_points_, sbnd_boundary_tag_,
                            collinear_recover_, collinear_interior_, collinear_member_merge_,
                            track_repartition_, band_merge_back_,
                            band_recarve_, drift_side_fv_x_,
                            far_point_x_cut_, far_point_mid_dis_, track_recarve_, dec1_guard_main_angle_,
                            iso_slab_split_, tag_family_, collinear_global_merge_);
    }

private:
    double use_ctpc_{true};
    int max_hull_points_{-1};
    bool sbnd_boundary_tag_{false};
    bool collinear_recover_{false};
    bool collinear_interior_{false};
    bool collinear_member_merge_{false};
    bool track_repartition_{false};
    bool band_merge_back_{false};
    bool band_recarve_{false};
    bool drift_side_fv_x_{false};
    double far_point_x_cut_{140 * units::cm};
    double far_point_mid_dis_{25 * units::cm};
    bool track_recarve_{false};
    double dec1_guard_main_angle_{-1.0};
    bool iso_slab_split_{false};
    bool tag_family_{false};
    bool collinear_global_merge_{false};
};


// The original developers do not care.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"

// this algorithm should be able to handle multiple APA/face now ..
static void clustering_separate(
    Grouping& live_grouping,
    const IDetectorVolumes::pointer dv,                // detector volumes
    const IPCTransformSet::pointer pcts,
    const Tree::Scope& scope,
    const bool use_ctpc,
    const int max_hull_points,
    const bool sbnd_boundary_tag,
    const bool collinear_recover,
    const bool collinear_interior,
    const bool collinear_member_merge,
    const bool track_repartition,
    const bool band_merge_back,
    const bool band_recarve,
    const bool drift_side_fv_x,
    const double far_point_x_cut,
    const double far_point_mid_dis,
    const bool track_recarve,
    const double dec1_guard_main_angle,
    const bool iso_slab_split,
    const bool tag_family,
    const bool collinear_global_merge)
{
    // Check that live_grouping has exactly one wpid
	// if (live_grouping.wpids().size() != 1 ) {
	// 	throw std::runtime_error("Live or Dead grouping must have exactly one wpid");
	// }
    geo_point_t drift_dir_abs(1,0,0);

    std::vector<Cluster *> live_clusters = live_grouping.children();  // copy
    // sort the clusters by length descending; use cluster ident() as tiebreaker
    // to guarantee a deterministic order across runs (pointer address is not stable).
    std::sort(live_clusters.begin(), live_clusters.end(), [](const Cluster *c1, const Cluster *c2) {
        if (c1->get_length() != c2->get_length()) return c1->get_length() > c2->get_length();
        return c1->ident() < c2->ident();
    });

    // Scope-aware fiducial volume: per-APA stages use that drift volume's FV;
    // multi-APA (all-APA) stages use the cryostat envelope, except that with
    // drift_side_fv_x a common-face multi-APA scope (drift group) keeps its
    // drift side's x-range.  See select_scope_fv.
    const ScopeFV fv = select_scope_fv(dv, drift_side_fv_x);
    const double det_FV_ymax = fv.ymax;
    const geo_point_t beam_dir = fv.beam_dir;
    const geo_point_t vertical_dir = fv.vertical_dir;

    for (size_t i = 0; i != live_clusters.size(); i++) {
        Cluster *cluster = live_clusters.at(i);
        if (!cluster->get_scope_filter(scope)) continue;
        if (cluster->get_default_scope().hash() != scope.hash()) {
            cluster->set_default_scope(scope);
            // std::cout << "Test: Set default scope: " << pc_name << " " << coords[0] << " " << coords[1] << " " << coords[2] << " " << cluster->get_default_scope().hash() << " " << scope.hash() << std::endl;
        }

        // Cache the length and PCA once — both are expensive for large clusters.
        const double cluster_length = cluster->get_length();

        if (cluster_length > 100 * units::cm) {
            std::vector<geo_point_t> boundary_points;
            std::vector<geo_point_t> independent_points;

            // JudgeSeparateDec_2 populates boundary_points / independent_points;
            // cache the return value so we don't re-run it below.
            bool flag_dec2 =
                JudgeSeparateDec_2(cluster, drift_dir_abs, boundary_points, independent_points, cluster_length, fv,
                                   max_hull_points, far_point_x_cut, far_point_mid_dis);
            // JudgeSeparateDec_1 is cheap compared to Dec_2 but still calls PCA —
            // cache for the second call inside flag_proceed block.
            bool flag_dec1 = JudgeSeparateDec_1(cluster, drift_dir_abs, cluster_length, dec1_guard_main_angle);

            bool flag_proceed = flag_dec2;

            if (sep_debug()) {
                const auto& pca_dbg = cluster->get_pca();
                geo_point_t md(pca_dbg.axis.at(0).x(), pca_dbg.axis.at(0).y(), pca_dbg.axis.at(0).z());
                std::cout << "SEPDBG gate ident=" << cluster->ident() << " len=" << cluster_length / units::cm
                          << " nblob=" << cluster->nchildren()
                          << " dec1=" << flag_dec1 << " dec2=" << flag_dec2
                          << " nindep=" << independent_points.size()
                          << " r1=" << pca_dbg.values.at(1) / pca_dbg.values.at(0)
                          << " r2=" << pca_dbg.values.at(2) / pca_dbg.values.at(0)
                          << " angbeam=" << fabs(md.angle(beam_dir) - 3.1415926 / 2.) / 3.1415926 * 180.
                          << std::endl;
            }

            if (!flag_proceed && flag_dec1 && independent_points.size() > 0) {
                bool flag_top = false;
                for (size_t j = 0; j != independent_points.size(); j++) {
                    if (independent_points.at(j).y() > det_FV_ymax) {
                        flag_top = true;
                        break;
                    }
                }

                // Cache PCA result to avoid repeated calls in the condition chain.
                const auto& pca = cluster->get_pca();
                geo_point_t main_dir(pca.axis.at(0).x(), pca.axis.at(0).y(), pca.axis.at(0).z());
                const double pca_ratio1 = pca.values.at(1) / pca.values.at(0);

                if (flag_top) {
                    if (fabs(main_dir.angle(beam_dir) - 3.1415926 / 2.) / 3.1415926 * 180. < 16 ||
                        fabs(main_dir.angle(beam_dir) - 3.1415926 / 2.) / 3.1415926 * 180. < 33 &&
                            cluster_length > 160 * units::cm ||
                        fabs(main_dir.angle(beam_dir) - 3.1415926 / 2.) / 3.1415926 * 180. < 40 &&
                            cluster_length > 260 * units::cm ||
                        fabs(main_dir.angle(beam_dir) - 3.1415926 / 2.) / 3.1415926 * 180. < 65 &&
                            cluster_length > 360 * units::cm ||
                        fabs(main_dir.angle(beam_dir) - 3.1415926 / 2.) / 3.1415926 * 180. < 45 &&
                            pca_ratio1 > 0.75 ||
                        fabs(main_dir.angle(beam_dir) - 3.1415926 / 2.) / 3.1415926 * 180. < 40 &&
                            pca_ratio1 > 0.55) {
                        flag_proceed = true;
                    }
                    else {
                        if (fabs(main_dir.angle(beam_dir) - 3.1415926 / 2.) / 3.1415926 * 180. < 40 &&
                            pca_ratio1 > 0.2) {
                            // std::vector<Cluster *> temp_sep_clusters = Separate_2(cluster, 10 * units::cm);
                            const auto b2id = Separate_2(cluster, scope, 10 * units::cm);
                            std::set<int> ids;
                            for (const auto& id : b2id) {
                                ids.insert(id);
                            }
                            int num_clusters = 0;

                            for (const auto id : ids) {
                                double length_1 = get_length(cluster, b2id, id);
                                if (length_1 > 60 * units::cm) num_clusters++;
                            }
                            if (num_clusters > 1) flag_proceed = true;
                        }
                    }
                }
                else {
                    if (fabs(main_dir.angle(beam_dir) - 3.1415926 / 2.) / 3.1415926 * 180. < 4 &&
                            cluster_length > 170 * units::cm ||
                        fabs(main_dir.angle(beam_dir) - 3.1415926 / 2.) / 3.1415926 * 180. < 25 &&
                            cluster_length > 210 * units::cm ||
                        fabs(main_dir.angle(beam_dir) - 3.1415926 / 2.) / 3.1415926 * 180. < 28 &&
                            cluster_length > 270 * units::cm ||
                        fabs(main_dir.angle(beam_dir) - 3.1415926 / 2.) / 3.1415926 * 180. < 35 &&
                            cluster_length > 330 * units::cm ||
                        fabs(main_dir.angle(beam_dir) - 3.1415926 / 2.) / 3.1415926 * 180. < 30 &&
                            pca_ratio1 > 0.55) {
                        flag_proceed = true;
                    }
                }

                // SBND-only: catch two beam-inclined cosmics that merged into one
                // cluster and whose PCA looks single-track (so the angle gates above
                // miss them) but which expose two endpoints on the upstream wall plus
                // a third exit elsewhere.  Default OFF; purely additive (false->true).
                if (!flag_proceed && sbnd_boundary_tag &&
                    JudgeSeparateDec_SBND_boundary(cluster, boundary_points, fv, cluster_length, independent_points)) {
                    flag_proceed = true;
                }
            }

            // Two crossing tracks whose arms end INSIDE the volume expose at
            // most one surface contact, so neither Dec_2 nor the angle ladders
            // can ever fire on them.  Offer such non-proceeding clusters
            // directly to the k=2 3D-line splitter -- its own validation
            // (two long thin arms, fair point shares, interior crossing) does
            // the actual discrimination, so the trigger here is loose.
            if (!flag_proceed && track_recarve) {
                std::vector<Cluster *> self = {cluster};
                recarve_crossing_tracks(live_grouping, scope, self);
                continue;  // on a split `cluster` was consumed (dangling)
            }

            if (flag_proceed) {
                // The parent is consumed by Separate_1 below; capture its
                // ident now for the tag_family stamp on the final members.
                const int parent_ident = cluster->ident();
                // Track every output cluster created while separating this
                // original cluster ("family") for the knob-gated refinement
                // post-pass.  Consumed intermediates are value-erased; vector
                // order is creation order (deterministic).
                std::vector<Cluster *> family;
                auto family_add = [&family](const std::vector<Cluster *>& cs) {
                    for (Cluster* c : cs) family.push_back(c);
                };
                auto family_erase = [&family](Cluster* c) {
                    family.erase(std::remove(family.begin(), family.end(), c), family.end());
                };

                // flag_dec1 was already computed above; reuse it here.
                if (flag_dec1) {
                    //	  std::cerr << em("sep prepare sep") << std::endl;

                    const size_t orig_nchildren = cluster->nchildren();
                    //std::cout << "Separate Cluster with " << orig_nchildren << " blobs (ctpc) length " << cluster_length << std::endl;
                    std::vector<Cluster *> sep_clusters =
                        Separate_1(use_ctpc, cluster, boundary_points, independent_points, cluster_length, vertical_dir, beam_dir, dv, pcts, scope);
                    family_add(sep_clusters);

                    std::cout << "Separate Separate_1 for " << orig_nchildren << " " << " returned " << sep_clusters.size() << " clusters" << std::endl;

                    if (sep_clusters.size() >= 2) {  // 1
                        Cluster *cluster2 = sep_clusters.at(1);
                        double length_1 = cluster2->get_length();

                        Cluster *final_sep_cluster = cluster2;

                        if (length_1 > 100 * units::cm) {
                            boundary_points.clear();
                            independent_points.clear();

                            if (JudgeSeparateDec_1(cluster2, drift_dir_abs, length_1, dec1_guard_main_angle) &&
                                JudgeSeparateDec_2(cluster2, drift_dir_abs, boundary_points, independent_points,
                                                   length_1, fv, max_hull_points, far_point_x_cut, far_point_mid_dis)) {
                                std::vector<Cluster *> sep_clusters =
                                    Separate_1(use_ctpc, cluster2, boundary_points, independent_points, length_1, vertical_dir, beam_dir, dv, pcts, scope);
                                family_erase(cluster2);  // consumed by Separate_1
                                family_add(sep_clusters);

                                if (sep_clusters.size() >= 2) {  // 2
                                    Cluster *cluster4 = sep_clusters.at(1);
                                    final_sep_cluster = cluster4;
                                    length_1 = cluster4->get_length();

                                    if (length_1 > 100 * units::cm) {
                                        boundary_points.clear();
                                        independent_points.clear();
                                        if (JudgeSeparateDec_1(cluster4, drift_dir_abs, length_1, dec1_guard_main_angle) &&
                                            JudgeSeparateDec_2(cluster4, drift_dir_abs, boundary_points, independent_points,
                                                               length_1, fv, max_hull_points, far_point_x_cut, far_point_mid_dis)) {
                                            std::vector<Cluster *> sep_clusters = Separate_1(
                                                use_ctpc, cluster4, boundary_points, independent_points, length_1, vertical_dir, beam_dir, dv, pcts, scope);
                                            family_erase(cluster4);  // consumed by Separate_1
                                            family_add(sep_clusters);

                                            if (sep_clusters.size() >= 2) {  // 3
                                                final_sep_cluster = sep_clusters.at(1);
                                            }
                                            else {
                                                final_sep_cluster = 0;
                                            }
                                        }
                                    }
                                }
                                else {
                                    final_sep_cluster = 0;
                                }
                            }
                        }

                        if (final_sep_cluster != 0) {  // 1
                            length_1 = final_sep_cluster->get_length();

                            if (length_1 > 60 * units::cm) {
                                boundary_points.clear();
                                independent_points.clear();
                                JudgeSeparateDec_1(final_sep_cluster, drift_dir_abs, length_1, dec1_guard_main_angle);
                                JudgeSeparateDec_2(final_sep_cluster, drift_dir_abs, boundary_points, independent_points,
                                                   length_1, fv, max_hull_points, far_point_x_cut, far_point_mid_dis);
                                if (independent_points.size() > 0) {
                                    std::vector<Cluster *> sep_clusters = Separate_1(
                                        use_ctpc, final_sep_cluster, boundary_points, independent_points, length_1, vertical_dir, beam_dir, dv, pcts, scope);
                                    family_erase(final_sep_cluster);  // consumed by Separate_1
                                    family_add(sep_clusters);

                                    if (sep_clusters.size() >= 2) {  // 4
                                        final_sep_cluster = sep_clusters.at(1);
                                    }
                                    else {
                                        final_sep_cluster = 0;
                                    }
                                }
                            }

                            if (final_sep_cluster != 0) {  // 2
                                const auto b2id = Separate_2(final_sep_cluster, scope);
                                family_erase(final_sep_cluster);  // consumed by separate()
                                auto frags = live_grouping.separate(final_sep_cluster, b2id, true);
                                assert(final_sep_cluster == nullptr);
                                for (auto& [gid, cl] : frags) family.push_back(cl);
                            }
                        }

                    }
                }
                else if (cluster_length < 6 * units::m) {
                    std::vector<Cluster *> sep_clusters =
                        Separate_1(use_ctpc, cluster, boundary_points, independent_points, cluster_length, vertical_dir, beam_dir, dv, pcts, scope);
                    family_add(sep_clusters);

                    if (sep_clusters.size() >= 2) {
                        Cluster *final_sep_cluster = sep_clusters.at(1);
                        const auto b2id = Separate_2(final_sep_cluster, scope);
                        family_erase(final_sep_cluster);  // consumed by separate()
                        auto frags = live_grouping.separate(final_sep_cluster, b2id, true);
                        assert(final_sep_cluster == nullptr);
                        for (auto& [gid, cl] : frags) family.push_back(cl);
                    }
                }  // else ...

                // Knob-gated refinement of this cluster's separation outputs
                // (no-ops when the knobs are OFF, the default).  Tip recovery
                // first: it pulls collinear track tips out of the leftover bin
                // before the band re-carve pools band-like members.
                if (collinear_recover && family.size() >= 2)
                    recover_collinear_tips(live_grouping, scope, family, collinear_interior);
                if (collinear_member_merge && family.size() >= 2)
                    merge_collinear_members(live_grouping, scope, family);
                if (track_repartition && family.size() >= 2)
                    repartition_crossing_tracks(live_grouping, scope, family);
                if (band_merge_back && family.size() >= 2)
                    merge_back_bands(live_grouping, scope, family);
                if (band_recarve && family.size() >= 2)
                    recarve_two_bands(live_grouping, scope, family);
                if (track_recarve && !family.empty())
                    recarve_crossing_tracks(live_grouping, scope, family);
                if (iso_slab_split && !family.empty())
                    split_iso_slabs(live_grouping, scope, family);
                // The slab split can leave one drift-spanning track as 2+
                // collinear family pieces (its seed joining is line-offset
                // gated, which over-penalizes long gently-curved tracks);
                // give the member-level rejoin a second look at the final
                // family (PDVD 39252 evt 298637).
                if (collinear_member_merge && iso_slab_split && family.size() >= 2)
                    merge_collinear_members(live_grouping, scope, family);
                // Stamp the surviving family members so a later same-stage
                // pass can decline to re-merge what was just deliberately
                // split (PDHD 27980 evt 24: connect1 re-glued a fat 80 cm
                // branch onto the 491 cm cosmic it was separated from).
                if (tag_family) {
                    int nsurv = 0;
                    for (Cluster* m : family)
                        if (m) nsurv++;
                    if (nsurv >= 2)
                        for (Cluster* m : family)
                            if (m) m->set_scalar<int>("sep_family", parent_ident + 1);
                }
            }
        }
    }

    // Grouping-wide collinear stitch (collinear_global_merge): two long thin
    // clusters that touch end-to-end and form one thin track were never
    // merged when they are NOT siblings of one separation family (the
    // member-level rejoin above only sees a family) and their global axes
    // disagree just enough (~10 deg from curvature) to fail connect1's <=5 deg
    // prolongation tests.  PDVD 39252 evt 298637: 321 cm + 104 cm pieces of
    // one cosmic, touching at 0.32 cm, local end directions 3.7 deg apart.
    // Reuses merge_collinear_members and its measured gates unchanged.
    if (collinear_global_merge) {
        const double GLOBAL_SHORT_MIN = 80 * units::cm;  // anchor >=100 cm enforced inside
        std::vector<Cluster*> candidates;
        for (Cluster* c : live_grouping.children()) {
            if (!c->get_scope_filter(scope)) continue;
            if (c->get_length() < GLOBAL_SHORT_MIN) continue;
            if (c->get_default_scope().hash() != scope.hash()) c->set_default_scope(scope);
            candidates.push_back(c);
        }
        std::sort(candidates.begin(), candidates.end(), [](const Cluster* c1, const Cluster* c2) {
            if (c1->get_length() != c2->get_length()) return c1->get_length() > c2->get_length();
            return c1->ident() < c2->ident();
        });
        if (candidates.size() >= 2)
            merge_collinear_members(live_grouping, scope, candidates, GLOBAL_SHORT_MIN);
    }


//  {
//    auto live_clusters = live_grouping.children(); // copy
//     // Process each cluster
//     for (size_t iclus = 0; iclus < live_clusters.size(); ++iclus) {
//         Cluster* cluster = live_clusters.at(iclus);
//         auto& scope = cluster->get_default_scope();
//         std::cout << "Test: " << iclus << " " << cluster->nchildren() << " " << scope.pcname << " " << scope.coords[0] << " " << scope.coords[1] << " " << scope.coords[2] << " " << cluster->get_scope_filter(scope)<< " " << cluster->get_pca().center << std::endl;
//     }
//   }







}

/// @brief PCA based, drift_dir +x, -x the same ...
bool WireCell::Clus::Facade::JudgeSeparateDec_1(const Cluster* cluster, const geo_point_t& drift_dir_abs, const double length,
                                                double guard_main_angle)
{
    // get the main axis — cache PCA result to avoid repeated calls
    const auto& pca = cluster->get_pca();
    geo_point_t dir1(pca.axis.at(0).x(), pca.axis.at(0).y(), pca.axis.at(0).z());
    geo_point_t dir2(pca.axis.at(1).x(), pca.axis.at(1).y(), pca.axis.at(1).z());
    geo_point_t dir3(pca.axis.at(2).x(), pca.axis.at(2).y(), pca.axis.at(2).z());

    double angle1 = fabs(dir2.angle(drift_dir_abs) - 3.1415926 / 2.) / 3.1415926 * 180.;

    // temp_angle1 uses the x-extent of the earliest/latest 3D points as a proxy
    // for the drift-axis extent.  The prototype used num_time_slices * tick_width
    // instead; the toolkit formulation is equivalent for dense clusters and is
    // naturally multi-APA/face because it operates directly on 3D positions.
    auto points = cluster->get_earliest_latest_points();
    double temp_angle1 = asin(fabs(points.first.x() - points.second.x()) / length) / 3.1415926 * 180.;

    double angle2 = fabs(dir3.angle(drift_dir_abs) - 3.1415926 / 2.) / 3.1415926 * 180.;
    double ratio1 = pca.values.at(1) / pca.values.at(0);
    double ratio2 = pca.values.at(2) / pca.values.at(0);

    // std::cout << ratio1 << " " <<  pow(10, exp(1.38115 - 1.19312 * pow(angle1, 1. / 3.)) - 2.2)  << " "  << ratio1 << " " << pow(10, exp(1.38115 - 1.19312 * pow(temp_angle1, 1. / 3.)) - 2.2) << " " << ratio2 << " " << pow(10, exp(1.38115 - 1.19312 * pow(angle2, 1. / 3.)) - 2.2) << " " << ratio1 << " " << angle1 << " " << angle2 << std::endl;

    // PDHD-specific guard: long (>3 m), nearly drift-aligned (angle1 < 5°), very
    // thin (ratio2 < 5%) clusters are through-going cosmics with excellent 3D
    // reconstruction — separating them would be wrong.  The ratio thresholds
    // below would falsely trigger on such tracks because the PCA is dominated by
    // the long drift direction.
    if (angle1 < 5 && ratio2 < 0.05 && length > 300*units::cm) {
        if (guard_main_angle < 0) return false;  // legacy: unconditional guard
        // Knob-gated: only protect clusters whose MAIN axis is actually
        // drift-aligned (the through-going-cosmic topology the guard is for).
        const double main_ang = acos(std::min(1.0, fabs(dir1.dot(drift_dir_abs)))) / 3.1415926 * 180.;
        if (main_ang < guard_main_angle) return false;
    }

    if (ratio1 > pow(10, exp(1.38115 - 1.19312 * pow(angle1, 1. / 3.)) - 2.2) ||
        ratio1 > pow(10, exp(1.38115 - 1.19312 * pow(temp_angle1, 1. / 3.)) - 2.2) ||
        ratio2 > pow(10, exp(1.38115 - 1.19312 * pow(angle2, 1. / 3.)) - 2.2) || 
        ratio1 > 0.75)
        return true;
    return false;
}

bool WireCell::Clus::Facade::JudgeSeparateDec_2(const Cluster* cluster, const geo_point_t& drift_dir_abs,
                               std::vector<geo_point_t>& boundary_points, std::vector<geo_point_t>& independent_points,
                               const double cluster_length, const ScopeFV& fv, int max_hull_points,
                               double far_point_x_cut, double far_point_mid_dis)
{
    // Scope-aware fiducial volume (see select_scope_fv): in a per-APA pass this is
    // the FV of the drift volume being clustered, so "exiting" means leaving that
    // volume; in an all-APA pass it is the global cryostat envelope (multi-APA scope).
    // Note: FV_xmin/xmax etc. are the *inner fiducial* values; the *_margin terms are
    // added outward to reach the physical boundary.
    const double det_FV_xmin = fv.xmin;
    const double det_FV_xmax = fv.xmax;
    const double det_FV_ymin = fv.ymin;
    const double det_FV_ymax = fv.ymax;
    const double det_FV_zmin = fv.zmin;
    const double det_FV_zmax = fv.zmax;
    const double det_FV_xmin_margin = fv.xmin_margin;
    const double det_FV_xmax_margin = fv.xmax_margin;
    const double det_FV_ymax_margin = fv.ymax_margin;
    const double det_FV_zmin_margin = fv.zmin_margin;
    const double det_FV_zmax_margin = fv.zmax_margin;

    boundary_points = cluster->get_hull(max_hull_points);

    // if get_hull failed, return false
    if (boundary_points.size() == 0) {
        return false;
    }
    std::vector<geo_point_t> hy_points;
    std::vector<geo_point_t> ly_points;
    std::vector<geo_point_t> hz_points;
    std::vector<geo_point_t> lz_points;
    std::vector<geo_point_t> hx_points;
    std::vector<geo_point_t> lx_points;

    std::set<int> independent_surfaces;

    for (size_t j = 0; j != boundary_points.size(); j++) {
        if (j == 0) {
            hy_points.push_back(boundary_points.at(j));
            ly_points.push_back(boundary_points.at(j));
            hz_points.push_back(boundary_points.at(j));
            lz_points.push_back(boundary_points.at(j));
            hx_points.push_back(boundary_points.at(j));
            lx_points.push_back(boundary_points.at(j));
        }
        else {
            geo_point_t test_p(boundary_points.at(j).x(), boundary_points.at(j).y(), boundary_points.at(j).z());
            if (cluster->nnearby(test_p, 15 * units::cm) > 75) {
                if (boundary_points.at(j).y() > hy_points.at(0).y()) hy_points.at(0) = boundary_points.at(j);
                if (boundary_points.at(j).y() < ly_points.at(0).y()) ly_points.at(0) = boundary_points.at(j);
                if (boundary_points.at(j).x() > hx_points.at(0).x()) hx_points.at(0) = boundary_points.at(j);
                if (boundary_points.at(j).x() < lx_points.at(0).x()) lx_points.at(0) = boundary_points.at(j);
                if (boundary_points.at(j).z() > hz_points.at(0).z()) hz_points.at(0) = boundary_points.at(j);
                if (boundary_points.at(j).z() < lz_points.at(0).z()) lz_points.at(0) = boundary_points.at(j);
            }
        }
    }

    bool flag_outx = false;
    if (hx_points.at(0).x() > det_FV_xmax + det_FV_xmax_margin || lx_points.at(0).x() < det_FV_xmin - det_FV_xmin_margin) flag_outx = true;

    if (hy_points.at(0).y() > det_FV_ymax) {
        for (size_t j = 0; j != boundary_points.size(); j++) {
            if (boundary_points.at(j).y() > det_FV_ymax) {
                bool flag_save = true;
                for (size_t k = 0; k != hy_points.size(); k++) {
                    double dis2 = pow(hy_points.at(k).x() - boundary_points.at(j).x(), 2) +
                                  pow(hy_points.at(k).y() - boundary_points.at(j).y(), 2) +
                                  pow(hy_points.at(k).z() - boundary_points.at(j).z(), 2);
                    if (dis2 < 625 * units::cm * units::cm) {
                        if (boundary_points.at(j).y() > hy_points.at(k).y()) hy_points.at(k) = boundary_points.at(j);
                        flag_save = false;
                    }
                }
                if (flag_save) hy_points.push_back(boundary_points.at(j));
            }
        }
    }

    if (ly_points.at(0).y() < det_FV_ymin) {
        for (size_t j = 0; j != boundary_points.size(); j++) {
            if (boundary_points.at(j).y() < det_FV_ymin) {
                bool flag_save = true;
                for (size_t k = 0; k != ly_points.size(); k++) {
                    double dis2 = pow(ly_points.at(k).x() - boundary_points.at(j).x(), 2) +
                                  pow(ly_points.at(k).y() - boundary_points.at(j).y(), 2) +
                                  pow(ly_points.at(k).z() - boundary_points.at(j).z(), 2);
                    if (dis2 < 625 * units::cm * units::cm) {
                        if (boundary_points.at(j).y() < ly_points.at(k).y()) ly_points.at(k) = boundary_points.at(j);
                        flag_save = false;
                    }
                }
                if (flag_save) ly_points.push_back(boundary_points.at(j));
            }
        }
    }
    if (hz_points.at(0).z() > det_FV_zmax) {
        for (size_t j = 0; j != boundary_points.size(); j++) {
            if (boundary_points.at(j).z() > det_FV_zmax) {
                bool flag_save = true;
                for (size_t k = 0; k != hz_points.size(); k++) {
                    double dis2 = pow(hz_points.at(k).x() - boundary_points.at(j).x(), 2) +
                                  pow(hz_points.at(k).y() - boundary_points.at(j).y(), 2) +
                                  pow(hz_points.at(k).z() - boundary_points.at(j).z(), 2);
                    if (dis2 < 625 * units::cm * units::cm) {
                        if (boundary_points.at(j).z() > hz_points.at(k).z()) hz_points.at(k) = boundary_points.at(j);
                        flag_save = false;
                    }
                }
                if (flag_save) hz_points.push_back(boundary_points.at(j));
            }
        }
    }
    if (lz_points.at(0).z() < det_FV_zmin) {
        for (size_t j = 0; j != boundary_points.size(); j++) {
            if (boundary_points.at(j).z() < det_FV_zmin) {
                bool flag_save = true;
                for (size_t k = 0; k != lz_points.size(); k++) {
                    double dis2 = pow(lz_points.at(k).x() - boundary_points.at(j).x(), 2) +
                                  pow(lz_points.at(k).y() - boundary_points.at(j).y(), 2) +
                                  pow(lz_points.at(k).z() - boundary_points.at(j).z(), 2);
                    if (dis2 < 625 * units::cm * units::cm) {
                        if (boundary_points.at(j).z() < lz_points.at(k).z()) lz_points.at(k) = boundary_points.at(j);
                        flag_save = false;
                    }
                }
                if (flag_save) lz_points.push_back(boundary_points.at(j));
            }
        }
    }

    int num_outside_points = 0;
    int num_outx_points = 0;

    for (size_t j = 0; j != hy_points.size(); j++) {
        if (hy_points.at(j).x() >= det_FV_xmin && hy_points.at(j).x() <= det_FV_xmax &&
            hy_points.at(j).y() >= det_FV_ymin && hy_points.at(j).y() <= det_FV_ymax &&
            hy_points.at(j).z() >= det_FV_zmin && hy_points.at(j).z() <= det_FV_zmax && (!flag_outx))
            continue;

        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis2 = pow(hy_points.at(j).x() - independent_points.at(k).x(), 2) +
                          pow(hy_points.at(j).y() - independent_points.at(k).y(), 2) +
                          pow(hy_points.at(j).z() - independent_points.at(k).z(), 2);
            if (dis2 < 225 * units::cm * units::cm) flag_save = false;
        }
        if (flag_save) {
            independent_points.push_back(hy_points.at(j));
            if (hy_points.at(j).y() > det_FV_ymax + det_FV_ymax_margin) {
                independent_surfaces.insert(0);
            }
            else if (hy_points.at(j).y() < det_FV_ymin) {
                independent_surfaces.insert(1);
            }
            else if (hy_points.at(j).z() > det_FV_zmax + det_FV_zmax_margin) {
                independent_surfaces.insert(2);
            }
            else if (hy_points.at(j).z() < det_FV_zmin - det_FV_zmin_margin) {
                independent_surfaces.insert(3);
            }
            else if (hy_points.at(j).x() > det_FV_xmax) {
                independent_surfaces.insert(4);
            }
            else if (hy_points.at(j).x() < det_FV_xmin) {
                independent_surfaces.insert(5);
            }

            if (hy_points.at(j).y() > det_FV_ymax + det_FV_ymax_margin || hy_points.at(j).y() < det_FV_ymin ||
                hy_points.at(j).z() < det_FV_zmin - det_FV_zmin_margin || hy_points.at(j).z() > det_FV_zmax + det_FV_zmax_margin ||
                hy_points.at(j).x() < det_FV_xmin || hy_points.at(j).x() > det_FV_xmax)
                num_outside_points++;
            if (hy_points.at(j).x() < det_FV_xmin - det_FV_xmin_margin || hy_points.at(j).x() > det_FV_xmax + det_FV_xmax_margin) num_outx_points++;
        }
    }
    for (size_t j = 0; j != ly_points.size(); j++) {
        if (ly_points.at(j).x() >= det_FV_xmin && ly_points.at(j).x() <= det_FV_xmax &&
            ly_points.at(j).y() >= det_FV_ymin && ly_points.at(j).y() <= det_FV_ymax &&
            ly_points.at(j).z() >= det_FV_zmin && ly_points.at(j).z() <= det_FV_zmax && (!flag_outx))
            continue;

        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis2 = pow(ly_points.at(j).x() - independent_points.at(k).x(), 2) +
                          pow(ly_points.at(j).y() - independent_points.at(k).y(), 2) +
                          pow(ly_points.at(j).z() - independent_points.at(k).z(), 2);
            if (dis2 < 225 * units::cm * units::cm) flag_save = false;
        }
        if (flag_save) {
            independent_points.push_back(ly_points.at(j));

            if (ly_points.at(j).y() < det_FV_ymin) {
                independent_surfaces.insert(1);
            }
            else if (ly_points.at(j).y() > det_FV_ymax + det_FV_ymax_margin) {
                independent_surfaces.insert(0);
            }
            else if (ly_points.at(j).z() > det_FV_zmax + det_FV_zmax_margin) {
                independent_surfaces.insert(2);
            }
            else if (ly_points.at(j).z() < det_FV_zmin - det_FV_zmin_margin) {
                independent_surfaces.insert(3);
            }
            else if (ly_points.at(j).x() > det_FV_xmax) {
                independent_surfaces.insert(4);
            }
            else if (ly_points.at(j).x() < det_FV_xmin) {
                independent_surfaces.insert(5);
            }

            if (ly_points.at(j).y() > det_FV_ymax + det_FV_ymax_margin || ly_points.at(j).y() < det_FV_ymin ||
                ly_points.at(j).z() < det_FV_zmin - det_FV_zmin_margin || ly_points.at(j).z() > det_FV_zmax + det_FV_zmax_margin ||
                ly_points.at(j).x() < det_FV_xmin || ly_points.at(j).x() > det_FV_xmax)
                num_outside_points++;
            if (ly_points.at(j).x() < det_FV_xmin - det_FV_xmin_margin || ly_points.at(j).x() > det_FV_xmax + det_FV_xmax_margin) num_outx_points++;
        }
    }
    for (size_t j = 0; j != hz_points.size(); j++) {
        if (hz_points.at(j).x() >= det_FV_xmin && hz_points.at(j).x() <= det_FV_xmax &&
            hz_points.at(j).y() >= det_FV_ymin && hz_points.at(j).y() <= det_FV_ymax &&
            hz_points.at(j).z() >= det_FV_zmin && hz_points.at(j).z() <= det_FV_zmax && (!flag_outx))
            continue;

        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis2 = pow(hz_points.at(j).x() - independent_points.at(k).x(), 2) +
                          pow(hz_points.at(j).y() - independent_points.at(k).y(), 2) +
                          pow(hz_points.at(j).z() - independent_points.at(k).z(), 2);
            if (dis2 < 225 * units::cm * units::cm) flag_save = false;
        }
        if (flag_save) {
            independent_points.push_back(hz_points.at(j));

            if (hz_points.at(j).z() > det_FV_zmax + det_FV_zmax_margin) {
                independent_surfaces.insert(2);
            }
            else if (hz_points.at(j).z() < det_FV_zmin - det_FV_zmin_margin) {
                independent_surfaces.insert(3);
            }
            else if (hz_points.at(j).y() > det_FV_ymax + det_FV_ymax_margin) {
                independent_surfaces.insert(0);
            }
            else if (hz_points.at(j).y() < det_FV_ymin) {
                independent_surfaces.insert(1);
            }
            else if (hz_points.at(j).x() > det_FV_xmax) {
                independent_surfaces.insert(4);
            }
            else if (hz_points.at(j).x() < det_FV_xmin) {
                independent_surfaces.insert(5);
            }

            if (hz_points.at(j).y() > det_FV_ymax + det_FV_ymax_margin || hz_points.at(j).y() < det_FV_ymin ||
                hz_points.at(j).z() < det_FV_zmin - det_FV_zmin_margin || hz_points.at(j).z() > det_FV_zmax + det_FV_zmax_margin ||
                hz_points.at(j).x() < det_FV_xmin || hz_points.at(j).x() > det_FV_xmax)
                num_outside_points++;
            if (hz_points.at(j).x() < det_FV_xmin - det_FV_xmin_margin || hz_points.at(j).x() > det_FV_xmax + det_FV_xmax_margin) num_outx_points++;
        }
    }
    for (size_t j = 0; j != lz_points.size(); j++) {
        if (lz_points.at(j).x() >= det_FV_xmin && lz_points.at(j).x() <= det_FV_xmax &&
            lz_points.at(j).y() >= det_FV_ymin && lz_points.at(j).y() <= det_FV_ymax &&
            lz_points.at(j).z() >= det_FV_zmin && lz_points.at(j).z() <= det_FV_zmax && (!flag_outx))
            continue;

        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis2 = pow(lz_points.at(j).x() - independent_points.at(k).x(), 2) +
                          pow(lz_points.at(j).y() - independent_points.at(k).y(), 2) +
                          pow(lz_points.at(j).z() - independent_points.at(k).z(), 2);
            if (dis2 < 225 * units::cm * units::cm) flag_save = false;
        }
        if (flag_save) {
            independent_points.push_back(lz_points.at(j));

            if (lz_points.at(j).z() < det_FV_zmin - det_FV_zmin_margin) {
                independent_surfaces.insert(3);
            }
            else if (lz_points.at(j).z() > det_FV_zmax + det_FV_zmax_margin) {
                independent_surfaces.insert(2);
            }
            else if (lz_points.at(j).y() > det_FV_ymax + det_FV_ymax_margin) {
                independent_surfaces.insert(0);
            }
            else if (lz_points.at(j).y() < det_FV_ymin) {
                independent_surfaces.insert(1);
            }
            else if (lz_points.at(j).x() > det_FV_xmax) {
                independent_surfaces.insert(4);
            }
            else if (lz_points.at(j).x() < det_FV_xmin) {
                independent_surfaces.insert(5);
            }

            if (lz_points.at(j).y() > det_FV_ymax + det_FV_ymax_margin || lz_points.at(j).y() < det_FV_ymin ||
                lz_points.at(j).z() < det_FV_zmin - det_FV_zmin_margin || lz_points.at(j).z() > det_FV_zmax + det_FV_zmax_margin ||
                lz_points.at(j).x() < det_FV_xmin || lz_points.at(j).x() > det_FV_xmax)
                num_outside_points++;
            if (lz_points.at(j).x() < det_FV_xmin - det_FV_xmin_margin || lz_points.at(j).x() > det_FV_xmax + det_FV_xmax_margin) num_outx_points++;
        }
    }
    for (size_t j = 0; j != hx_points.size(); j++) {
        if (hx_points.at(j).x() >= det_FV_xmin && hx_points.at(j).x() <= det_FV_xmax &&
            hx_points.at(j).y() >= det_FV_ymin && hx_points.at(j).y() <= det_FV_ymax &&
            hx_points.at(j).z() >= det_FV_zmin && hx_points.at(j).z() <= det_FV_zmax && (!flag_outx))
            continue;

        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis2 = pow(hx_points.at(j).x() - independent_points.at(k).x(), 2) +
                          pow(hx_points.at(j).y() - independent_points.at(k).y(), 2) +
                          pow(hx_points.at(j).z() - independent_points.at(k).z(), 2);
            if (dis2 < 225 * units::cm * units::cm) flag_save = false;
        }
        if (flag_save) {
            independent_points.push_back(hx_points.at(j));

            if (hx_points.at(j).y() > det_FV_ymax + det_FV_ymax_margin || hx_points.at(j).y() < det_FV_ymin ||
                hx_points.at(j).z() < det_FV_zmin - det_FV_zmin_margin || hx_points.at(j).z() > det_FV_zmax + det_FV_zmax_margin ||
                hx_points.at(j).x() < det_FV_xmin || hx_points.at(j).x() > det_FV_xmax)
                num_outside_points++;
            if (hx_points.at(j).x() < det_FV_xmin - det_FV_xmin_margin || hx_points.at(j).x() > det_FV_xmax + det_FV_xmax_margin) {
                num_outx_points++;
            }

            // Bug fix: prototype mistakenly read lx_points.at(j) inside this hx_points
            // loop, causing an out-of-bounds access when hx_points.size() > lx_points.size().
            // Corrected to hx_points.at(j).
            if (hx_points.at(j).x() > det_FV_xmax) {
                independent_surfaces.insert(4);
            }
            else if (hx_points.at(j).x() < det_FV_xmin) {
                independent_surfaces.insert(5);
            }
            else if (hx_points.at(j).y() > det_FV_ymax + det_FV_ymax_margin) {
                independent_surfaces.insert(0);
            }
            else if (hx_points.at(j).y() < det_FV_ymin) {
                independent_surfaces.insert(1);
            }
            else if (hx_points.at(j).z() > det_FV_zmax + det_FV_zmax_margin) {
                independent_surfaces.insert(2);
            }
            else if (hx_points.at(j).z() < det_FV_zmin - det_FV_zmin_margin) {
                independent_surfaces.insert(3);
            }
        }
    }
    for (size_t j = 0; j != lx_points.size(); j++) {
        if (lx_points.at(j).x() >= det_FV_xmin && lx_points.at(j).x() <= det_FV_xmax &&
            lx_points.at(j).y() >= det_FV_ymin && lx_points.at(j).y() <= det_FV_ymax &&
            lx_points.at(j).z() >= det_FV_zmin && lx_points.at(j).z() <= det_FV_zmax && (!flag_outx))
            continue;

        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis2 = pow(lx_points.at(j).x() - independent_points.at(k).x(), 2) +
                          pow(lx_points.at(j).y() - independent_points.at(k).y(), 2) +
                          pow(lx_points.at(j).z() - independent_points.at(k).z(), 2);
            if (dis2 < 225 * units::cm * units::cm) flag_save = false;
        }
        if (flag_save) {
            independent_points.push_back(lx_points.at(j));

            if (lx_points.at(j).y() > det_FV_ymax + det_FV_ymax_margin || lx_points.at(j).y() < det_FV_ymin ||
                lx_points.at(j).z() < det_FV_zmin - det_FV_zmin_margin || lx_points.at(j).z() > det_FV_zmax + det_FV_zmax_margin ||
                lx_points.at(j).x() < det_FV_xmin || lx_points.at(j).x() > det_FV_xmax)
                num_outside_points++;
            if (lx_points.at(j).x() < det_FV_xmin - det_FV_xmin_margin || lx_points.at(j).x() > det_FV_xmax + det_FV_xmax_margin) {
                num_outx_points++;
            }

            if (lx_points.at(j).x() < det_FV_xmin) {
                independent_surfaces.insert(5);
            }
            else if (lx_points.at(j).x() > det_FV_xmax) {
                independent_surfaces.insert(4);
            }
            else if (lx_points.at(j).y() > det_FV_ymax + det_FV_ymax_margin) {
                independent_surfaces.insert(0);
            }
            else if (lx_points.at(j).y() < det_FV_ymin) {
                independent_surfaces.insert(1);
            }
            else if (lx_points.at(j).z() > det_FV_zmax + det_FV_zmax_margin) {
                independent_surfaces.insert(2);
            }
            else if (lx_points.at(j).z() < det_FV_zmin - det_FV_zmin_margin) {
                independent_surfaces.insert(3);
            }
        }
    }

    int num_far_points = 0;

    if (independent_points.size() == 2 && (independent_surfaces.size() > 1 || flag_outx)) {
        geo_vector_t dir_1(independent_points.at(1).x() - independent_points.at(0).x(),
                           independent_points.at(1).y() - independent_points.at(0).y(),
                           independent_points.at(1).z() - independent_points.at(0).z());
        dir_1 = dir_1.norm();
        for (size_t j = 0; j != boundary_points.size(); j++) {
            geo_vector_t dir_2(boundary_points.at(j).x() - independent_points.at(0).x(),
                               boundary_points.at(j).y() - independent_points.at(0).y(),
                               boundary_points.at(j).z() - independent_points.at(0).z());
            double angle_12 = dir_1.angle(dir_2);
            geo_vector_t dir_3 = dir_2 - dir_1 * dir_2.magnitude() * cos(angle_12);
            double angle_3 = dir_3.angle(drift_dir_abs);
            // std::cout << dir_3.Mag()/units::cm << " " << fabs(angle_3-3.1415926/2.)/3.1415926*180. << " " <<
            // fabs(dir_3.X()/units::cm) << std::endl;
            if (fabs(angle_3 - 3.1415926 / 2.) / 3.1415926 * 180. < 7.5) {
                // Default far_point_x_cut (140 cm) keeps the prototype expression
                // `fabs(dir_3.x()/units::cm) > 14*units::cm` bit-for-bit.
                if (fabs(dir_3.x()) > far_point_x_cut) num_far_points++;
                if (fabs(dir_1.angle(drift_dir_abs) - 3.1415926 / 2.) / 3.1415926 * 180. > 15) {
                    if (dir_3.magnitude() > 20 * units::cm) num_far_points++;
                }
            }
            else {
                if (dir_3.magnitude() > 20 * units::cm) num_far_points++;
            }
        }

        // find the middle points and close distance ...
        geo_point_t middle_point((independent_points.at(1).x() + independent_points.at(0).x()) / 2.,
                                 (independent_points.at(1).y() + independent_points.at(0).y()) / 2.,
                                 (independent_points.at(1).z() + independent_points.at(0).z()) / 2.);
        double middle_dis = cluster->get_closest_dis(middle_point);
        // std::cout << middle_dis/units::cm << " " << num_far_points << std::endl;
        if (middle_dis > far_point_mid_dis) {
            num_far_points = 0;
        }
    }

    double max_x = -1e9, min_x = 1e9;
    double max_y = -1e9, min_y = 1e9;
    double max_z = -1e9, min_z = 1e9;
    for (auto it = independent_points.begin(); it != independent_points.end(); it++) {
        if ((*it).x() > max_x) max_x = (*it).x();
        if ((*it).x() < min_x) min_x = (*it).x();
        if ((*it).y() > max_y) max_y = (*it).y();
        if ((*it).y() < min_y) min_y = (*it).y();
        if ((*it).z() > max_z) max_z = (*it).z();
        if ((*it).z() < min_z) min_z = (*it).z();
        // std::cout << (*it).x()/units::cm << " " << (*it).y()/units::cm << " " << (*it).z()/units::cm << std::endl;
    }
    if (hx_points.size() > 0) {
        if (hx_points.at(0).x() > max_x + 10 * units::cm) max_x = hx_points.at(0).x();
    }
    if (lx_points.size() > 0) {
        if (lx_points.at(0).x() < min_x - 10 * units::cm) min_x = lx_points.at(0).x();
    }

    if (max_x - min_x < 2.5 * units::cm &&
        pow(max_y - min_y, 2) + pow(max_z - min_z, 2) + pow(max_x - min_x, 2) > 22500 * units::cm * units::cm) {
        independent_points.clear();
        return false;
    }
    if (max_x - min_x < 2.5 * units::cm && independent_points.size() == 2 && num_outx_points == 0) {
        independent_points.clear();
        return false;
    }


    if (sep_debug()) {
        std::cout << "SEPDBG dec2 ident=" << cluster->ident() << " len=" << cluster_length / units::cm
                  << " nout=" << num_outside_points << " noutx=" << num_outx_points
                  << " nsurf=" << independent_surfaces.size() << " nindep=" << independent_points.size()
                  << " nfar=" << num_far_points << " outx=" << flag_outx
                  << " hy=" << hy_points.at(0).y() / units::cm << " ly=" << ly_points.at(0).y() / units::cm
                  << " hz=" << hz_points.at(0).z() / units::cm << " lz=" << lz_points.at(0).z() / units::cm
                  << " hx=" << hx_points.at(0).x() / units::cm << " lx=" << lx_points.at(0).x() / units::cm
                  << std::endl;
    }

    if ((num_outside_points > 1 && independent_surfaces.size() > 1 ||
         num_outside_points > 2 && cluster_length > 250 * units::cm || num_outx_points > 0) &&
        (independent_points.size() > 2 || independent_points.size() == 2 && num_far_points > 0))
        return true;

    // about to return false ...
    independent_points.clear();

    for (size_t j = 0; j != hy_points.size(); j++) {
        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis2 = pow(hy_points.at(j).x() - independent_points.at(k).x(), 2) +
                          pow(hy_points.at(j).y() - independent_points.at(k).y(), 2) +
                          pow(hy_points.at(j).z() - independent_points.at(k).z(), 2);
            if (dis2 < 225 * units::cm * units::cm) flag_save = false;
        }
        if (flag_save) independent_points.push_back(hy_points.at(j));
    }

    for (size_t j = 0; j != ly_points.size(); j++) {
        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis2 = pow(ly_points.at(j).x() - independent_points.at(k).x(), 2) +
                          pow(ly_points.at(j).y() - independent_points.at(k).y(), 2) +
                          pow(ly_points.at(j).z() - independent_points.at(k).z(), 2);
            if (dis2 < 225 * units::cm * units::cm) flag_save = false;
        }
        if (flag_save) independent_points.push_back(ly_points.at(j));
    }

    for (size_t j = 0; j != hx_points.size(); j++) {
        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis2 = pow(hx_points.at(j).x() - independent_points.at(k).x(), 2) +
                          pow(hx_points.at(j).y() - independent_points.at(k).y(), 2) +
                          pow(hx_points.at(j).z() - independent_points.at(k).z(), 2);
            if (dis2 < 225 * units::cm * units::cm) flag_save = false;
        }
        if (flag_save) independent_points.push_back(hx_points.at(j));
    }

    for (size_t j = 0; j != lx_points.size(); j++) {
        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis2 = pow(lx_points.at(j).x() - independent_points.at(k).x(), 2) +
                          pow(lx_points.at(j).y() - independent_points.at(k).y(), 2) +
                          pow(lx_points.at(j).z() - independent_points.at(k).z(), 2);
            if (dis2 < 225 * units::cm * units::cm) flag_save = false;
        }
        if (flag_save) independent_points.push_back(lx_points.at(j));
    }

    for (size_t j = 0; j != hz_points.size(); j++) {
        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis2 = pow(hz_points.at(j).x() - independent_points.at(k).x(), 2) +
                          pow(hz_points.at(j).y() - independent_points.at(k).y(), 2) +
                          pow(hz_points.at(j).z() - independent_points.at(k).z(), 2);
            if (dis2 < 225 * units::cm * units::cm) flag_save = false;
        }
        if (flag_save) independent_points.push_back(hz_points.at(j));
    }

    for (size_t j = 0; j != lz_points.size(); j++) {
        bool flag_save = true;
        for (size_t k = 0; k != independent_points.size(); k++) {
            double dis2 = pow(lz_points.at(j).x() - independent_points.at(k).x(), 2) +
                          pow(lz_points.at(j).y() - independent_points.at(k).y(), 2) +
                          pow(lz_points.at(j).z() - independent_points.at(k).z(), 2);
            if (dis2 < 225 * units::cm * units::cm) flag_save = false;
        }
        if (flag_save) independent_points.push_back(lz_points.at(j));
    }

    return false;
}

#define _INDEV_
#ifdef _INDEV_

std::vector<Cluster *> WireCell::Clus::Facade::Separate_1(const bool use_ctpc, Cluster *cluster,
                                                     std::vector<geo_point_t> &boundary_points,
                                                     std::vector<geo_point_t> &independent_points,
                                                     double length, geo_point_t dir_cosmic, geo_point_t dir_beam, const IDetectorVolumes::pointer dv, const IPCTransformSet::pointer pcts, const Tree::Scope& scope)
{
    const std::string graph_flavor = use_ctpc ? "ctpc" : "basic";

    auto* grouping = cluster->grouping();

    auto gwpids = grouping->wpids();

    std::map<int, std::map<int, std::map<int, std::pair<double, double>>>> af_dead_u_index ;
    std::map<int, std::map<int, std::map<int, std::pair<double, double>>>> af_dead_v_index ;
    std::map<int, std::map<int, std::map<int, std::pair<double, double>>>> af_dead_w_index ;
    std::map<int, std::map<int, std::shared_ptr<Multi2DPointCloud>>> af_temp_cloud;
    for (auto wpid : gwpids) {
        int apa = wpid.apa();
        int face = wpid.face();
        af_dead_u_index[apa][face] = grouping->get_dead_winds(apa, face, 0); // raw 
        af_dead_v_index[apa][face] = grouping->get_dead_winds(apa, face, 1); // raw
        af_dead_w_index[apa][face] = grouping->get_dead_winds(apa, face, 2); // raw 

        // Create wpids for all three planes with this APA and face
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

        af_temp_cloud[apa][face] = std::make_shared<Multi2DPointCloud>(angle_u, angle_v, angle_w); // 2D Dynamic Point Cloud
    }

    // std::cout << "Test: " << pc_name << " " << coords[0] << " " << coords[1] << " " << coords[2] << std::endl;

    geo_point_t cluster_center = cluster->get_pca().center;
    geo_point_t main_dir = cluster->get_pca().axis.at(0);
    geo_point_t second_dir = cluster->get_pca().axis.at(1);
    

    // special case, if one of the cosmic is very close to the beam direction
    if (cluster->get_pca().values.at(1) > 0.08 * cluster->get_pca().values.at(0) &&
        fabs(main_dir.angle(dir_beam) - 3.1415926 / 2.) > 75 / 180. * 3.1415926 &&
        fabs(second_dir.angle(dir_cosmic) - 3.1415926 / 2.) > 60 / 180. * 3.1415926) {
        main_dir = second_dir;
    }

    main_dir = main_dir.norm();
    if (main_dir.y() > 0)
        main_dir = main_dir * -1;  // make sure it is pointing down????

    geo_point_t start_wcpoint;
    geo_point_t end_wcpoint;
    geo_point_t drift_dir_abs(1, 0, 0);
    geo_point_t dir;

    double min_dis = 1e9;
    // double max_pca_dis;
    int min_index = 0;
    double max_dis = -1e9;
    // double min_pca_dis;
    int max_index = 0;
    // if (flag_debug_porting) {
    //     std::cout << "independent_points.size(): " << independent_points.size() << std::endl;
    // }
    for (size_t j = 0; j != independent_points.size(); j++) {
        geo_point_t dir(independent_points.at(j).x() - cluster_center.x(),
                        independent_points.at(j).y() - cluster_center.y(),
                        independent_points.at(j).z() - cluster_center.z());
        geo_point_t temp_p(independent_points.at(j).x(), independent_points.at(j).y(), independent_points.at(j).z());
        double dis = dir.dot(main_dir);
        // double dis_to_pca = dir.cross(main_dir).magnitude();
        // std::cout << j << " " << dis << " " << dir.Mag() << " " << sqrt(dir.Mag()*dir.Mag() - dis*dis) << std::endl;
        bool flag_connect = false;
        int num_points = cluster->nnearby(temp_p, 15 * units::cm);
        if (num_points > 100) {
            flag_connect = true;
        }
        else if (num_points > 75) {
            num_points = cluster->nnearby(temp_p, 30 * units::cm);
            if (num_points > 160) flag_connect = true;
        }

        // std::cout << dis / units::cm << " A " << cluster->get_num_points(temp_p,15*units::cm) << " " <<
        // cluster->get_num_points(temp_p,30*units::cm)  << std::endl;
        if (dis < min_dis && flag_connect) {
            min_dis = dis;
            min_index = j;
            // min_pca_dis = dis_to_pca;
        }
        if (dis > max_dis && flag_connect) {
            max_dis = dis;
            max_index = j;
            // max_pca_dis = dis_to_pca;
        }
    }
 


    size_t start_wcpoint_idx = 0;
    size_t end_wcpoint_idx = 0;
    {
        start_wcpoint = independent_points.at(min_index);

        // change direction if certain thing happened ...

        {
            geo_point_t p1(independent_points.at(max_index).x(), independent_points.at(max_index).y(),
                           independent_points.at(max_index).z());
            geo_point_t p2(independent_points.at(min_index).x(), independent_points.at(min_index).y(),
                           independent_points.at(min_index).z());
            geo_point_t temp_dir1 = cluster->vhough_transform(p1, 15 * units::cm);
            geo_point_t temp_dir2 = cluster->vhough_transform(p2, 15 * units::cm);

            bool flag_change = false;

            if (fabs(temp_dir1.angle(main_dir) - 3.1415926 / 2.) > fabs(temp_dir2.angle(main_dir) - 3.1415926 / 2.)) {
                if (fabs(temp_dir2.angle(main_dir) - 3.1415926 / 2.) > 80 / 180. * 3.1415926 &&
                    fabs(temp_dir1.angle(main_dir) - 3.1415926 / 2.) <
                        fabs(temp_dir2.angle(main_dir) - 3.1415926 / 2.) + 2.5 / 180. * 3.1415926) {
                }
                else {
                    flag_change = true;
                    start_wcpoint = independent_points.at(max_index);
                    main_dir = main_dir * -1;
                    max_index = min_index;
                }
            }

            if ((!flag_change) &&
                fabs(temp_dir1.angle(drift_dir_abs) - 3.1415926 / 2.) > fabs(temp_dir2.angle(drift_dir_abs) - 3.1415926 / 2.) &&
                fabs(temp_dir2.angle(drift_dir_abs) - 3.1415926 / 2.) / 3.1415926 * 180. < 10 &&
                fabs(temp_dir2.angle(main_dir) - 3.1415926 / 2.) / 3.1415926 * 180. < 80) {
                start_wcpoint = independent_points.at(max_index);
                main_dir = main_dir * -1;
                max_index = min_index;
            }

            if ((!flag_change) && fabs(temp_dir2.angle(drift_dir_abs) - 3.1415926 / 2.) < 1. / 180. * 3.1415926 &&
                fabs(temp_dir1.angle(drift_dir_abs) - 3.1415926 / 2.) > 3. / 180. * 3.1415926 &&
                fabs(temp_dir1.angle(main_dir) - 3.1415926 / 2.) / 3.1415926 * 180. > 70) {
                start_wcpoint = independent_points.at(max_index);
                main_dir = main_dir * -1;
                max_index = min_index;
            }
        }

        geo_point_t start_point(start_wcpoint.x(), start_wcpoint.y(), start_wcpoint.z());
        {
            // geo_point_t drift_dir_abs(1, 0, 0);
            dir = cluster->vhough_transform(start_point, 100 * units::cm);
            geo_point_t dir1 = cluster->vhough_transform(start_point, 30 * units::cm);
            if (dir.angle(dir1) > 20 * 3.1415926 / 180.) {
                if (fabs(dir.angle(drift_dir_abs) - 3.1415926 / 2.) < 5 * 3.1415926 / 180. ||
                    fabs(dir1.angle(drift_dir_abs) - 3.1415926 / 2.) < 5 * 3.1415926 / 180.) {
                    dir = cluster->vhough_transform(start_point, 200 * units::cm);
                }
                else {
                    dir = dir1;
                }
            }
        }
        dir = dir.norm();

        geo_point_t inv_dir = dir * (-1);
        start_wcpoint = cluster->get_furthest_wcpoint(start_wcpoint, inv_dir, 1 * units::cm, 0);
        end_wcpoint = cluster->get_furthest_wcpoint(start_wcpoint, dir);

   
        geo_point_t test_dir(end_wcpoint.x() - start_wcpoint.x(), end_wcpoint.y() - start_wcpoint.y(),
                             end_wcpoint.z() - start_wcpoint.z());
        start_wcpoint_idx = cluster->get_closest_point_index(start_wcpoint);
        end_wcpoint_idx = cluster->get_closest_point_index(end_wcpoint);
        if (fabs(test_dir.angle(drift_dir_abs) - 3.1415926 / 2.) < 2.5 * 3.1415926 / 180.) {
            cluster->adjust_wcpoints_parallel(start_wcpoint_idx, end_wcpoint_idx);
            start_wcpoint = cluster->point3d(start_wcpoint_idx);
            end_wcpoint = cluster->point3d(end_wcpoint_idx);
  
        }
    }
    if (pow(start_wcpoint.x() - end_wcpoint.x(), 2) + pow(start_wcpoint.y() - end_wcpoint.y(), 2) +
             pow(start_wcpoint.z() - end_wcpoint.z(), 2) < (length / 3.) * (length / 3.)) {
        // reverse the case ...
        start_wcpoint = independent_points.at(max_index);
        geo_point_t start_point(start_wcpoint.x(), start_wcpoint.y(), start_wcpoint.z());
        {
            // geo_point_t drift_dir_abs(1, 0, 0);
            dir = cluster->vhough_transform(start_point, 100 * units::cm);
            geo_point_t dir1 = cluster->vhough_transform(start_point, 30 * units::cm);
            if (dir.angle(dir1) > 20 * 3.1415926 / 180.) {
                if (fabs(dir.angle(drift_dir_abs) - 3.1415926 / 2.) < 5 * 3.1415926 / 180. ||
                    fabs(dir1.angle(drift_dir_abs) - 3.1415926 / 2.) < 5 * 3.1415926 / 180.) {
                    dir = cluster->vhough_transform(start_point, 200 * units::cm);
                }
                else {
                    dir = dir1;
                }
            }
        }
        dir = dir.norm();
        geo_point_t inv_dir = dir * (-1);
        start_wcpoint = cluster->get_furthest_wcpoint(start_wcpoint, inv_dir, 1 * units::cm, 0);
        end_wcpoint = cluster->get_furthest_wcpoint(start_wcpoint, dir);

        

        if (pow(start_wcpoint.x() - end_wcpoint.x(), 2) + pow(start_wcpoint.y() - end_wcpoint.y(), 2) +
                 pow(start_wcpoint.z() - end_wcpoint.z(), 2) < (length / 3.) * (length / 3.)) {  // reverse again ...
            start_wcpoint = end_wcpoint;
            geo_point_t start_point(start_wcpoint.x(), start_wcpoint.y(), start_wcpoint.z());
            {
                dir = cluster->vhough_transform(start_point, 100 * units::cm);
                geo_point_t dir1 = cluster->vhough_transform(start_point, 30 * units::cm);
                if (dir.angle(dir1) > 20 * 3.1415926 / 180.) {
                    if (fabs(dir.angle(drift_dir_abs) - 3.1415926 / 2.) < 5 * 3.1415926 / 180. ||
                        fabs(dir1.angle(drift_dir_abs) - 3.1415926 / 2.) < 5 * 3.1415926 / 180.) {
                        dir = cluster->vhough_transform(start_point, 200 * units::cm);
                    }
                    else {
                        dir = dir1;
                    }
                }
            }
            dir = dir.norm();
            geo_point_t inv_dir = dir * (-1);
            start_wcpoint = cluster->get_furthest_wcpoint(start_wcpoint, inv_dir, 1 * units::cm, 0);
            end_wcpoint = cluster->get_furthest_wcpoint(start_wcpoint, dir);
        }

    
      
        geo_point_t test_dir(end_wcpoint.x() - start_wcpoint.x(), end_wcpoint.y() - start_wcpoint.y(),
                             end_wcpoint.z() - start_wcpoint.z());
        start_wcpoint_idx = cluster->get_closest_point_index(start_wcpoint);
        end_wcpoint_idx = cluster->get_closest_point_index(end_wcpoint);
        if (fabs(test_dir.angle(drift_dir_abs) - 3.1415926 / 2.) < 2.5 * 3.1415926 / 180.) {
            cluster->adjust_wcpoints_parallel(start_wcpoint_idx, end_wcpoint_idx);
            start_wcpoint = cluster->point3d(start_wcpoint_idx);
            end_wcpoint = cluster->point3d(end_wcpoint_idx);
      
        }
    }
  
    const auto& path_wcps = cluster->graph_algorithms(graph_flavor, dv, pcts).shortest_path(start_wcpoint_idx, end_wcpoint_idx);


    std::vector<bool> flag_u_pts, flag_v_pts, flag_w_pts;
    std::vector<bool> flag1_u_pts, flag1_v_pts, flag1_w_pts;
    std::vector<bool> flag2_u_pts, flag2_v_pts, flag2_w_pts;
    flag_u_pts.resize(cluster->npoints(), false);
    flag_v_pts.resize(cluster->npoints(), false);
    flag_w_pts.resize(cluster->npoints(), false);

    flag1_u_pts.resize(cluster->npoints(), false);
    flag1_v_pts.resize(cluster->npoints(), false);
    flag1_w_pts.resize(cluster->npoints(), false);

    flag2_u_pts.resize(cluster->npoints(), false);
    flag2_v_pts.resize(cluster->npoints(), false);
    flag2_w_pts.resize(cluster->npoints(), false);

    std::vector<geo_point_t> pts;

    // double acc_dis = 0;

    size_t prev_wcp_idx = path_wcps.front();
    for (auto it = path_wcps.begin(); it != path_wcps.end(); it++) {
        geo_point_t current_wcp = cluster->point3d(*it);
        geo_point_t prev_wcp = cluster->point3d(prev_wcp_idx);
        double dis = sqrt(pow(current_wcp.x() - prev_wcp.x(), 2) + pow(current_wcp.y() - prev_wcp.y(), 2) +
                          pow(current_wcp.z() - prev_wcp.z(), 2));
        // acc_dis += dis;
        // if (cluster->nchildren()==3449) std::cout << current_wcp << " " << dis/units::cm << " " << acc_dis/units::cm << std::endl;

        if (dis <= 1.0 * units::cm) {
            geo_point_t current_pt(current_wcp.x(), current_wcp.y(), current_wcp.z());
            pts.push_back(current_pt);
        }
        else {
            int num_points = int(dis / (1.0 * units::cm)) + 1;
            // double dis_seg = dis / num_points;
            for (int k = 0; k != num_points; k++) {
                geo_point_t current_pt(prev_wcp.x() + (k + 1.) / num_points * (current_wcp.x() - prev_wcp.x()),
                                       prev_wcp.y() + (k + 1.) / num_points * (current_wcp.y() - prev_wcp.y()),
                                       prev_wcp.z() + (k + 1.) / num_points * (current_wcp.z() - prev_wcp.z()));
                pts.push_back(current_pt);
            }
        }
        prev_wcp_idx = (*it);
    }
    for (const auto &pt : pts) {
        auto test_wpid = cluster->wpid(pt);
        if (test_wpid.apa()!=-1){
            af_temp_cloud.at(test_wpid.apa()).at(test_wpid.face())->add(pt);
        }
    }
    // if (flag_debug_porting) {
    //     std::cout << "temp_cloud->get_num_points() " << temp_cloud->get_num_points() << std::endl;
    // }

    const auto& winds = cluster->wire_indices();

    // For clusters spanning multiple (apa,face) pairs (e.g. PDHD dual-drift), a path
    // point recorded under face A may be the closest representative for a cluster point
    // that lives in face B.  To avoid missing such matches, when is_multi_face is true
    // we take the minimum 2D distance across ALL face projections rather than only the
    // test point's own face.  Dead-channel lookups still use the test point's own wpid.
    const bool is_multi_face = (gwpids.size() > 1);

    // Helper: minimum get_closest_2d_dis across all (apa,face) in af_temp_cloud.
    auto min_2d_dis = [&](const geo_point_t& p, int plane) -> double {
        double best = 1e9;
        for (const auto& [apa, face_map] : af_temp_cloud) {
            for (const auto& [face, cloud] : face_map) {
                double d = cloud->get_closest_2d_dis(p, plane).second;
                if (d < best) best = d;
            }
        }
        return best;
    };

    for (size_t j = 0; j != flag_u_pts.size(); j++) {
        geo_point_t test_p = cluster->point3d(j);
        auto test_wpid = cluster->wire_plane_id(j);
        // Symmetric guard to the one in the path-building loop above:
        // skip points that don't map to any (apa,face) to avoid at(-1) crashes.
        if (test_wpid.apa() == -1) continue;

        // 2D distance from test_p to the path projection.
        // Single-face (common case): use only the test point's own face projection.
        // Multi-face: widen to the minimum across all face projections so that a path
        // crossing an APA boundary is not "invisible" to cluster points near the boundary.
        double dis;
        if (is_multi_face) {
            dis = min_2d_dis(test_p, 0);
        } else {
            dis = af_temp_cloud.at(test_wpid.apa()).at(test_wpid.face())->get_closest_2d_dis(test_p, 0).second;
        }
        if (dis <= 1.5 * units::cm) {
            flag_u_pts.at(j) = true;
        }
        if (dis <= 2.4 * units::cm) {
            flag1_u_pts.at(j) = true;
        }
        else {
            auto& dead_u_index = af_dead_u_index.at(test_wpid.apa()).at(test_wpid.face());
            if (dead_u_index.find(winds[0][j]) != dead_u_index.end()) {
                // dead channels are corresponding to raw points
                if (cluster->point3d_raw(j).x() >= dead_u_index[winds[0][j]].first &&
                    cluster->point3d_raw(j).x() <= dead_u_index[winds[0][j]].second) {
                    if (dis < 10 * units::cm) flag1_u_pts.at(j) = true;
                    flag2_u_pts.at(j) = true;
                }
            }
        }
        if (is_multi_face) {
            dis = min_2d_dis(test_p, 1);
        } else {
            dis = af_temp_cloud.at(test_wpid.apa()).at(test_wpid.face())->get_closest_2d_dis(test_p, 1).second;
        }
        if (dis <= 1.5 * units::cm) {
            flag_v_pts.at(j) = true;
        }
        if (dis <= 2.4 * units::cm) {
            flag1_v_pts.at(j) = true;
        }
        else {
            auto& dead_v_index = af_dead_v_index.at(test_wpid.apa()).at(test_wpid.face());
            // dead channels are corresponding to raw points
            if (dead_v_index.find(winds[1][j]) != dead_v_index.end()) {
                if (cluster->point3d_raw(j).x() >= dead_v_index[winds[1][j]].first &&
                    cluster->point3d_raw(j).x() <= dead_v_index[winds[1][j]].second) {
                    if (dis < 10.0 * units::cm) flag1_v_pts.at(j) = true;
                    flag2_v_pts.at(j) = true;
                }
            }
        }
        if (is_multi_face) {
            dis = min_2d_dis(test_p, 2);
        } else {
            dis = af_temp_cloud.at(test_wpid.apa()).at(test_wpid.face())->get_closest_2d_dis(test_p, 2).second;
        }
        if (dis <= 1.5 * units::cm) {
            flag_w_pts.at(j) = true;
        }
        if (dis <= 2.4 * units::cm) {
            flag1_w_pts.at(j) = true;
        }
        else {
            auto& dead_w_index = af_dead_w_index.at(test_wpid.apa()).at(test_wpid.face());
            // dead channels are corresponding to raw points
            if (dead_w_index.find(winds[2][j]) != dead_w_index.end()) {
                if (cluster->point3d_raw(j).x() >= dead_w_index[winds[2][j]].first &&
                    cluster->point3d_raw(j).x() <= dead_w_index[winds[2][j]].second) {
                    if (dis < 10 * units::cm) flag1_w_pts.at(j) = true;
                    flag2_w_pts.at(j) = true;
                }
            }
        }
    }

    // special treatment of first and last point
    {
        auto wpid_front = cluster->wpid(pts.front());
        auto idx_front =  cluster->get_closest_point_index(pts.front());


        auto wpid_back = cluster->wpid(pts.back());
        auto idx_back =  cluster->get_closest_point_index(pts.back());

        std::vector<size_t> indices = cluster->get_closest_2d_index(cluster->point3d_raw(idx_front), 2.1 * units::cm, wpid_front.apa(), wpid_front.face(), 0);
        for (size_t k = 0; k != indices.size(); k++) {
            flag_u_pts.at(indices.at(k)) = true;
        }
        indices = cluster->get_closest_2d_index(cluster->point3d_raw(idx_front), 2.1 * units::cm, wpid_front.apa(), wpid_front.face(), 1);
        for (size_t k = 0; k != indices.size(); k++) {
            flag_v_pts.at(indices.at(k)) = true;
        }
        indices = cluster->get_closest_2d_index(cluster->point3d_raw(idx_front), 2.1 * units::cm, wpid_front.apa(), wpid_front.face(), 2);
        for (size_t k = 0; k != indices.size(); k++) {
            flag_w_pts.at(indices.at(k)) = true;
        }
        indices = cluster->get_closest_2d_index(cluster->point3d_raw(idx_back), 2.1 * units::cm, wpid_back.apa(), wpid_back.face(), 0);
        for (size_t k = 0; k != indices.size(); k++) {
            flag_u_pts.at(indices.at(k)) = true;
        }
        indices = cluster->get_closest_2d_index(cluster->point3d_raw(idx_back), 2.1 * units::cm, wpid_back.apa(), wpid_back.face(), 1);
        for (size_t k = 0; k != indices.size(); k++) {
            flag_v_pts.at(indices.at(k)) = true;
        }
        indices = cluster->get_closest_2d_index(cluster->point3d_raw(idx_back), 2.1 * units::cm, wpid_back.apa(), wpid_back.face(), 2);
        for (size_t k = 0; k != indices.size(); k++) {
            flag_w_pts.at(indices.at(k)) = true;
        }
    }

    // std::vector<Blob*>
    const auto &mcells = cluster->children();
    // Pointer-keyed maps: iteration order is non-deterministic, but these maps are
    // only ever accessed by direct key lookup (operator[] / find), never iterated.
    // Determinism is preserved because iteration over mcells (a vector) drives all loops.
    std::map<const Blob *, int> mcell_np_map, mcell_np_map1;
    for (auto it = mcells.begin(); it != mcells.end(); it++) {
        mcell_np_map[*it] = 0;
        mcell_np_map1[*it] = 0;
    }
    for (size_t j = 0; j != flag_u_pts.size(); j++) {
        const Blob* mcell = cluster->blob_with_point(j);
        
        if (flag_u_pts.at(j) && flag_v_pts.at(j) && flag1_w_pts.at(j) ||
            flag_u_pts.at(j) && flag_w_pts.at(j) && flag1_v_pts.at(j) ||
            flag_w_pts.at(j) && flag_v_pts.at(j) && flag1_u_pts.at(j)) {
            mcell_np_map[mcell]++;
        }

        if (flag_u_pts.at(j) && flag_v_pts.at(j) && (flag2_w_pts.at(j) || flag1_w_pts.at(j)) ||
            flag_u_pts.at(j) && flag_w_pts.at(j) && (flag2_v_pts.at(j) || flag1_v_pts.at(j)) ||
            flag_w_pts.at(j) && flag_v_pts.at(j) && (flag2_u_pts.at(j) || flag1_u_pts.at(j))) {
            mcell_np_map1[mcell]++;
        }
    }


    std::vector<Cluster *> final_clusters;

    std::vector<int> b2groupid(cluster->nchildren(), 0);
    std::set<int> groupids;

    for (size_t idx=0; idx < mcells.size(); idx++) {  
        Blob *mcell = mcells.at(idx);

        const size_t total_wires = mcell->u_wire_index_max() - mcell->u_wire_index_min() +
                             mcell->v_wire_index_max() - mcell->v_wire_index_min() +
                             mcell->w_wire_index_max() - mcell->w_wire_index_min();

        if (mcell_np_map[mcell] > 0.5 * mcell->nbpoints() ||
            (mcell_np_map[mcell] > 0.25 * mcell->nbpoints() && total_wires < 25)) {
            // cluster1->AddCell(mcell, mcell->GetTimeSlice());
            b2groupid[idx] = 0;
            groupids.insert(0);
        }
        else if (mcell_np_map1[mcell] >= 0.95 * mcell->nbpoints()) {
            // delete mcell;  // ghost cell ...
            b2groupid[idx] = -1; // to be deleted
            groupids.insert(-1);
        }
        else {
            // cluster2->AddCell(mcell, mcell->GetTimeSlice());
            b2groupid[idx] = 1;
            groupids.insert(1);
        }
    }
    auto clusters_step0 = grouping->separate(cluster, b2groupid, true);
    assert(cluster == nullptr);


    std::vector<Cluster*> other_clusters;

    if (clusters_step0.find(1) != clusters_step0.end()) {
        // other_clusters = Separate_2(clusters_step0[1], 5 * units::cm);
        const auto b2id = Separate_2(clusters_step0[1], scope, 5 * units::cm);
        auto other_clusters1 = grouping->separate(clusters_step0[1],b2id, true); // the cluster is now nullptr
        assert(clusters_step0[1] == nullptr);

        for (auto it = other_clusters1.begin(); it != other_clusters1.end(); it++) {
            other_clusters.push_back(it->second);
        }
    }


    if (clusters_step0.find(0) != clusters_step0.end()) {
        // merge some clusters from other_clusters to clusters_step0[0]
        {
            // cluster1->Create_point_cloud();
            // ToyPointCloud *cluster1_cloud = cluster1->get_point_cloud();
            std::vector<Cluster *> temp_merge_clusters;
            // check against other clusters
            for (size_t i = 0; i != other_clusters.size(); i++) {
                // other_clusters.at(i)->Create_point_cloud();
                std::tuple<int, int, double> temp_dis = other_clusters.at(i)->get_closest_points(*clusters_step0[0]);
                if (std::get<2>(temp_dis) < 0.5 * units::cm) {
                    double length_1 = other_clusters.at(i)->get_length();
                    geo_point_t p1(end_wcpoint.x(), end_wcpoint.y(), end_wcpoint.z());
                    double close_dis = other_clusters.at(i)->get_closest_dis(p1);

                    if (close_dis < 10 * units::cm && length_1 < 50 * units::cm) {
                        geo_point_t temp_dir1 = clusters_step0[0]->vhough_transform(p1, 15 * units::cm);
                        geo_point_t temp_dir2 = other_clusters.at(i)->vhough_transform(p1, 15 * units::cm);
                        if (temp_dir1.angle(temp_dir2) / 3.1415926 * 180. > 145 && length_1 < 30 * units::cm &&
                                close_dis < 3 * units::cm ||
                            fabs(temp_dir1.angle(drift_dir_abs) - 3.1415926 / 2.) / 3.1415926 * 180. < 3 &&
                                fabs(temp_dir2.angle(drift_dir_abs) - 3.1415926 / 2.) / 3.1415926 * 180. < 3) {
                            temp_merge_clusters.push_back(other_clusters.at(i));
                        }
                    }
                }
            }

            auto scope = clusters_step0[0]->get_default_scope();
            auto scope_transform = clusters_step0[0]->get_scope_transform(scope);

            // std::cout << "Xin1:  " << clusters_step0[0]->npoints() << " " << clusters_step0[0]->kd3d().npoints()  << " " << clusters_step0[0]->sv().nodes().size() << " " << clusters_step0[0]->sv().npoints() << std::endl;
                        
            for (auto temp_cluster : temp_merge_clusters) {
                other_clusters.erase(find(other_clusters.begin(),other_clusters.end(),temp_cluster));
                clusters_step0[0]->take_children(*temp_cluster, true);
                grouping->destroy_child(temp_cluster);
                assert(temp_cluster == nullptr);
            }
            clusters_step0[0]->set_default_scope(scope);
            clusters_step0[0]->set_scope_filter(scope, true);
            clusters_step0[0]->set_scope_transform(scope, scope_transform);

            // std::cout << "Xin2:  " << clusters_step0[0]->npoints() << " " << clusters_step0[0]->kd3d().npoints()  << " " << clusters_step0[0]->sv().nodes().size() << " " << clusters_step0[0]->sv().npoints() << std::endl;


            final_clusters.push_back(clusters_step0[0]);
        }

        // ToyPointCloud *cluster1_cloud = cluster1->get_point_cloud();
        std::vector<Cluster *> saved_clusters;
        std::vector<Cluster *> to_be_merged_clusters;
        for (size_t i = 0; i != other_clusters.size(); i++) {
            // How to write???
            bool flag_save = false;
            double length_1 = other_clusters.at(i)->get_length();
      
            std::tuple<int, int, double> temp_dis = other_clusters.at(i)->get_closest_points(*clusters_step0[0]);
            if (length_1 < 30 * units::cm && std::get<2>(temp_dis) < 5 * units::cm) {
                int temp_total_points = other_clusters.at(i)->npoints();
                int temp_close_points = 0;
                const int threshold_70 = int(0.7 * temp_total_points) + 1;
                for (int j = 0; j != temp_total_points; j++) {
                    if (clusters_step0[0]->get_closest_dis(other_clusters.at(i)->point3d(j)) < 10 * units::cm)
                        temp_close_points++;
                    if (temp_close_points >= threshold_70) break;               // already qualifies
                    if (temp_close_points + (temp_total_points - j - 1) < threshold_70) break; // can't qualify
                }
                if (temp_close_points >= threshold_70) {
                    saved_clusters.push_back(other_clusters.at(i));
                    flag_save = true;
                }
            }
            else if (std::get<2>(temp_dis) < 2.5 * units::cm && length_1 >= 30 * units::cm) {
                int temp_total_points = other_clusters.at(i)->npoints();
                int temp_close_points = 0;
                const int threshold_85 = int(0.85 * temp_total_points) + 1;
                for (int j = 0; j != temp_total_points; j++) {
                    if (clusters_step0[0]->get_closest_dis(other_clusters.at(i)->point3d(j)) < 10 * units::cm)
                        temp_close_points++;
                    if (temp_close_points >= threshold_85) break;               // already qualifies
                    if (temp_close_points + (temp_total_points - j - 1) < threshold_85) break; // can't qualify
                }
                if (temp_close_points >= threshold_85) {
                    saved_clusters.push_back(other_clusters.at(i));
                    flag_save = true;
                }
            }

            if (!flag_save) to_be_merged_clusters.push_back(other_clusters.at(i));
        }

        // add a protection

        // Pre-compute geometry for qualifying to_be_merged clusters (length >= 10 cm) once,
        // so the inner loop avoids repeated get_pca()/get_length() calls and can apply a
        // cheap center-distance lower bound before calling the expensive get_closest_points().
        struct TBMInfo {
            Cluster*    cluster;
            geo_point_t center;
            double      half_len;
            geo_point_t dir;
        };
        std::vector<TBMInfo> tbm_info;
        tbm_info.reserve(to_be_merged_clusters.size());
        for (Cluster* c : to_be_merged_clusters) {
            if (c->get_length() < 10 * units::cm) continue;
            const auto& pca = c->get_pca();
            tbm_info.push_back({c, pca.center, c->get_length() / 2.0,
                                 geo_point_t(pca.axis.at(0).x(), pca.axis.at(0).y(), pca.axis.at(0).z())});
        }

        std::vector<Cluster *> temp_save_clusters;

        for (size_t i = 0; i != saved_clusters.size(); i++) {
            Cluster *cluster1 = saved_clusters.at(i);
            if (cluster1->get_length() < 5 * units::cm) continue;
            const auto& pca1 = cluster1->get_pca();
            geo_point_t dir1(pca1.axis.at(0).x(), pca1.axis.at(0).y(), pca1.axis.at(0).z());
            geo_point_t center1 = pca1.center;
            double half_len1 = cluster1->get_length() / 2.0;

            for (const auto& ti : tbm_info) {
                // Center-distance lower bound: if center_dis - half_len1 - ti.half_len > 15 cm,
                // the closest points between the two clusters must exceed 15 cm — skip without
                // calling the KDtree.  Pure arithmetic, no change to the decision logic.
                double dx = center1.x() - ti.center.x();
                double dy = center1.y() - ti.center.y();
                double dz = center1.z() - ti.center.z();
                double max_reach = 15 * units::cm + half_len1 + ti.half_len;
                if (dx*dx + dy*dy + dz*dz > max_reach * max_reach) continue;

                std::tuple<int, int, double> temp_dis = cluster1->get_closest_points(*ti.cluster);
                if (std::get<2>(temp_dis) < 15 * units::cm &&
                    fabs(dir1.angle(ti.dir) - 3.1415926 / 2.) / 3.1415926 * 180 > 75) {
                    temp_save_clusters.push_back(cluster1);
                    break;
                }
            }
        }
        for (size_t i = 0; i != temp_save_clusters.size(); i++) {
            Cluster *cluster1 = temp_save_clusters.at(i);
            to_be_merged_clusters.push_back(cluster1);
            saved_clusters.erase(find(saved_clusters.begin(), saved_clusters.end(), cluster1));
        }

        Cluster& cluster2 = grouping->make_child();
        cluster2.set_default_scope(scope);
        cluster2.set_scope_filter(scope, true);
        // Inherit scope_transform from the first to-be-merged cluster when possible;
        // fall back to the path cluster's transform so cluster2 is always valid even
        // when to_be_merged_clusters is empty (nothing left after the saved-cluster split).
        if (to_be_merged_clusters.size() > 0)
            cluster2.set_scope_transform(scope, to_be_merged_clusters[0]->get_scope_transform(scope));
        else
            cluster2.set_scope_transform(scope, clusters_step0[0]->get_scope_transform(scope));
        for (size_t i = 0; i != to_be_merged_clusters.size(); i++) {
            cluster2.take_children(*to_be_merged_clusters[i], true);
            grouping->destroy_child(to_be_merged_clusters[i]);
            assert(to_be_merged_clusters[i] == nullptr);
        }        
        
        to_be_merged_clusters.clear();

        final_clusters.push_back(&cluster2);
        for (size_t i = 0; i != saved_clusters.size(); i++) {
            final_clusters.push_back(saved_clusters.at(i));
        }
        //saved_clusters.clear();
    }
    else {
        for (size_t i = 0; i != other_clusters.size(); i++) {
            final_clusters.push_back(other_clusters.at(i));
        }
    }
  
    return final_clusters;
}

#endif //_INDEV_

/// blob -> cluster_id
std::vector<int> WireCell::Clus::Facade::Separate_2(Cluster *cluster, 
                                                          const Tree::Scope& scope,
                                                          const double dis_cut)
{
    if (cluster->nchildren() == 0) {
        return std::vector<int>();
    }

    // std::cout << "Test: cluster has " << cluster->nchildren() << " blobs" << std::endl;
    auto& time_cells_set_map = cluster->time_blob_map();
    // Safe access to nested maps
    

    // std::cout << "Separate_2 nchildren: " << cluster->nchildren() << std::endl;
    const auto& mcells = cluster->children();

    // create graph for points between connected mcells, need to separate apa, face, and then ...
    std::map<int, std::map<int, std::vector<int> > > af_time_slices; // apa,face --> time slices 
    for (auto it = cluster->time_blob_map().begin(); it != cluster->time_blob_map().end(); it++) {
        int apa = it->first;
        for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++) {
            int face = it1->first;
            std::vector<int> time_slices_vec;
            for (auto it2 = it1->second.begin(); it2 != it1->second.end(); it2++) {
                time_slices_vec.push_back(it2->first);
            }
            af_time_slices[apa][face] = time_slices_vec;
        }
    }

    std::vector<std::pair<const Blob *, const Blob *>> connected_mcells;

    for (auto it = af_time_slices.begin(); it != af_time_slices.end(); it++) {
        int apa = it->first;
        for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++) {
            int face = it1->first;
            std::vector<int>& time_slices = it1->second;
 
            for (size_t i = 0; i != time_slices.size(); i++) {
                const BlobSet &mcells_set = time_cells_set_map.at(apa).at(face).at(time_slices.at(i));
                // std::cout << "time_slices.at(i)" << time_slices.at(i) << " mcells_set.size() " << mcells_set.size() << std::endl;

                // create graph for points in mcell inside the same time slice
                if (mcells_set.size() >= 2) {
                    for (auto it2 = mcells_set.begin(); it2 != mcells_set.end(); it2++) {
                        const Blob *mcell1 = *it2;
                        auto it2p = it2;
                        if (it2p != mcells_set.end()) {
                            it2p++;
                            for (auto it3 = it2p; it3 != mcells_set.end(); it3++) {
                                const Blob *mcell2 = *(it3);
                                if (mcell1->overlap_fast(*mcell2, 5)) connected_mcells.push_back(std::make_pair(mcell1, mcell2));
                            }
                        }
                    }
                }
                // create graph for points between connected mcells in adjacent time slices + 1, if not, + 2
                std::vector<BlobSet> vec_mcells_set;
                if (i + 1 < time_slices.size()) {
                    if (time_slices.at(i + 1) - time_slices.at(i) == (int)(1*cluster->grouping()->get_nticks_per_slice().at(apa).at(face))) {
                        vec_mcells_set.push_back(time_cells_set_map.at(apa).at(face).at(time_slices.at(i + 1)));
                        if (i + 2 < time_slices.size())
                            if (time_slices.at(i + 2) - time_slices.at(i) == (int)(2*cluster->grouping()->get_nticks_per_slice().at(apa).at(face)))
                                vec_mcells_set.push_back(time_cells_set_map.at(apa).at(face).at(time_slices.at(i + 2)));
                    }
                    else if (time_slices.at(i + 1) - time_slices.at(i) == (int)(2*cluster->grouping()->get_nticks_per_slice().at(apa).at(face))) {
                        vec_mcells_set.push_back(time_cells_set_map.at(apa).at(face).at(time_slices.at(i + 1)));
                    }
                }
                // std::cout << "time_slices.at(i)" << time_slices.at(i) << " vec_mcells_set.size() " << vec_mcells_set.size() << std::endl;
                bool flag = false;
                for (size_t j = 0; j != vec_mcells_set.size(); j++) {
                    if (flag) break;
                    BlobSet &next_mcells_set = vec_mcells_set.at(j);
                    for (auto it1 = mcells_set.begin(); it1 != mcells_set.end(); it1++) {
                        const Blob *mcell1 = (*it1);
                        for (auto it2 = next_mcells_set.begin(); it2 != next_mcells_set.end(); it2++) {
                            const Blob *mcell2 = (*it2);
                            if (mcell1->overlap_fast(*mcell2, 2)) {
                                flag = true;
                                connected_mcells.push_back(std::make_pair(mcell1, mcell2));
                            }
                        }
                    }
                }
            }
        }
    }

    // form ...

    const int N = mcells.size();
    Weighted::Graph graph(N);

    // Pointer-keyed map: iteration-order safe — only used for direct key lookups, not iteration.
    std::map<const Blob *, int> mcell_index_map;
    for (size_t i = 0; i != mcells.size(); i++) {
        Blob *curr_mcell = mcells.at(i);
        mcell_index_map[curr_mcell] = i;

        // auto v = vertex(i, graph);  // retrieve vertex descriptor
        // (graph)[v].ident = i;
    }

    for (auto it = connected_mcells.begin(); it != connected_mcells.end(); it++) {
        int index1 = mcell_index_map[it->first];
        int index2 = mcell_index_map[it->second];
        // auto edge = add_edge(index1, index2, graph);
        // if (edge.second) {
        //     (graph)[edge.first].dist = 1;
        // }
        /*auto edge =*/ add_edge(index1, index2, 1.0,graph);
    }

    {
        // std::string hack_pc_name = "3d";
        // std::vector<std::string> hack_coords = {"x", "y", "z"};
        // std::cout << "Separate_2: num_edges: " << num_edges(graph) << std::endl;
        std::vector<int> component(num_vertices(graph));
        const int num = connected_components(graph, &component[0]);

        if (num > 1) {
            std::vector<std::shared_ptr<Simple3DPointCloud>> pt_clouds;
            std::vector<std::vector<int>> vec_vec(num);
            for (int j = 0; j != num; j++) {
                pt_clouds.push_back(std::make_shared<Simple3DPointCloud>());
            }
            std::vector<int>::size_type i;
            for (i = 0; i != component.size(); ++i) {
                vec_vec.at(component[i]).push_back(i);
                Blob *mcell = mcells.at(i);
                for (const auto & pt : mcell->points(scope.pcname, scope.coords)) {
                    const std::vector<double> newpt = {pt.x(), pt.y(), pt.z()};
                    pt_clouds.at(component[i])->add(newpt);
                }
            }

            for (int j = 0; j != num; j++) {
                for (int k = j + 1; k != num; k++) {
                    std::tuple<int, int, double> temp_results = pt_clouds.at(j)->get_closest_points(*(pt_clouds.at(k)));
                    if (std::get<2>(temp_results) < dis_cut) {
                        int index1 = vec_vec[j].front();
                        int index2 = vec_vec[k].front();
                        // auto edge = add_edge(index1, index2, graph);
                        // if (edge.second) {
                        //     (graph)[edge.first].dist = 1;
                        // }
                        /*auto edge =*/ add_edge(index1, index2, 1.0,graph);
                    }
                }
            }
        }

        // std::cout << num << std::endl;
    }

  
    std::vector<int> component(num_vertices(graph));
    /*const int num =*/ connected_components(graph, &component[0]);
    return component;
}


