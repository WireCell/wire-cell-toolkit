// doc pdvd/31 round 3: BlobSampler silently lost the charge of every sampled
// point landing on a wrapped strip's continuation.
//
// Gen::AnodePlane builds a plane's channel vector by walking that plane's wires
// and SKIPPING every wrapped continuation (AnodePlane.cxx:244-247):
//
//     for (auto w : wires) {
//         if (w->segment() > 0) { continue; }
//         ...
//         plane_channels.push_back(ich);
//     }
//
// That is a correct channel LIST -- each channel appears once per anode,
// attached to the plane holding its segment-0 wire.  It is not a wire->channel
// lookup table, but BlobSampler used it as one: `p_chi2i[channel_ident]`, i.e.
// unordered_map::operator[], which INSERTS 0 on a miss.  Points on a wrapped
// continuation therefore read channels[0]'s activity -- normally absent, so
// charge_val AND charge_unc both stayed 0, which Cluster::calc_charge_wcp
// (Facade_Cluster.cxx:1087-1091) reads as "this plane saw nothing" and does not
// hold against the point.  Silent by construction: no warning, no sentinel.
//
// Measured on PDVD 039349/14: ~11% of the event's sampled points lost an
// induction plane, and 98.6% of the points along the charge-starved stretch
// that opened doc 31.
//
// These tests pin the two things the fix rests on:
//
//   1. the GEOMETRIC PREMISE -- that wrapped planes really do carry wires whose
//      channel the plane's own segment-0 list omits, with SBND and uBooNE as
//      negative controls where the count must be exactly zero;
//   2. the DEFAULT-OFF CONTRACT -- `wrapped_channel_charge` defaults false.
//      PDHD is affected worse than PDVD (28.8% of wires vs 11.3%) and its
//      config is silent on this knob, so the C++ default is the only thing
//      holding it byte-identical.  A PDHD flip is a separate owner call.
//
// Repro of the same census offline, all four detectors:
//   python3 wcp-porting-img/pdvd/docs/nf_sp_img_clus/scripts/steiner_orphan_channel_census.py

// NB: reach BlobSampler through the factory, not WireCellClus/BlobSampler.h --
// that header holds a vector<unique_ptr<Sampler>> of an incomplete type and so
// cannot be included outside its own translation unit.
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IBlobSampler.h"

#include "WireCellAux/SimpleBlob.h"
#include "WireCellAux/SimpleSlice.h"

#include "WireCellUtil/RayHelpers.h"
#include "WireCellUtil/RayTiling.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/WireSchema.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/doctest.h"

#include <set>
#include <string>

using namespace WireCell;

namespace {

    struct OrphanCount {
        int wires{0};        // wires in the detector
        int orphans{0};      // segment>0 wires whose channel their plane omits
        int planes{0};       // planes carrying at least one orphan
    };

    // Count wires that Gen::AnodePlane's plane_channels rule cannot resolve:
    // segment > 0 AND no segment-0 wire of the SAME plane carries the channel.
    //
    // A strip that wraps back into its own plane is NOT an orphan -- its
    // channel is listed via the segment-0 half and the lookup succeeds.  That
    // distinction is the whole point: it is what makes the defect a property of
    // the (anode, face, plane), not of the charge.
    OrphanCount count_orphans(const std::string& filename)
    {
        OrphanCount out;
        const auto store = WireSchema::load(filename.c_str());
        for (const auto& anode : store.anodes()) {
            for (const auto& face : store.faces(anode)) {
                for (const auto& plane : store.planes(face)) {
                    const auto wires = store.wires(plane);
                    std::set<int> listed;   // what plane_channels would hold
                    for (const auto& w : wires) {
                        if (w.segment == 0) listed.insert(w.channel);
                    }
                    int here = 0;
                    for (const auto& w : wires) {
                        if (w.segment > 0 && !listed.count(w.channel)) ++here;
                    }
                    out.wires += (int) wires.size();
                    out.orphans += here;
                    if (here) ++out.planes;
                }
            }
        }
        return out;
    }

    // doc pdvd/31 round 4.  The SAME AnodePlane rule, seen from the other side:
    // `plane_channels` is pushed back in wire order while skipping segment>0, so
    //
    //     channels[i]->ident() == wires[i]->channel()
    //
    // holds ONLY while no segment>0 wire has been skipped yet, and
    // channels.size() == wires.size() - n_seg_gt0.
    //
    // Two live sites index that vector by WIRE index and so depend on exactly
    // this invariant: ImproveCluster_1::make_iblobs_improved
    // (improvecluster_1.cxx:833-846), which is the retiler the Steiner stage
    // runs (`cm.improve_cluster_2`), and the identical line in base
    // RetileCluster (retile_cluster.cxx:426-437), which no production pipeline
    // currently reaches.  Where the invariant fails, the retile writes a slice's
    // activity under the WRONG IChannel, and past channels.size() it drops the
    // wire's activity entirely.
    struct PlaneRow {
        int anode{-1}, face{-1}, plane{-1};   // idents, not array positions
        int nwires{0}, seg_gt0{0}, first_bad{-1};
    };
    struct PlaneMap {
        int planes{0};            // planes examined
        int planes_broken{0};     // planes where channels[i] stops naming wire i
        int seg_gt0{0};           // total segment>0 wires
        int first_bad{-1};        // smallest wire index at which it first breaks
        std::vector<PlaneRow> broken;
    };

    // doc pdvd/31 round 5: the PRECONDITION of the retile-mapping fix.
    //
    // The fix resolves a wire's channel through IAnodePlane::channel(ident),
    // which reads Gen::AnodePlane's m_ichannels -- a map filled ONLY from
    // segment-0 wires (AnodePlane.cxx:253, inside the same loop that skips
    // segment>0).  It therefore resolves an orphan of one plane only because
    // that channel's segment-0 wire lives in ANOTHER plane of the SAME anode.
    //
    // Returns the number of channels that break that assumption: seen on a
    // segment>0 wire somewhere in an anode, and on no segment-0 wire anywhere
    // in that same anode.  Must be 0, or IAnodePlane::channel() returns nullptr
    // and the fix silently drops exactly the wires it exists to rescue.
    int unresolvable_channels(const std::string& filename)
    {
        int bad = 0;
        const auto store = WireSchema::load(filename.c_str());
        for (const auto& anode : store.anodes()) {
            std::set<int> seg0, segn;
            for (const auto& face : store.faces(anode)) {
                for (const auto& plane : store.planes(face)) {
                    for (const auto& w : store.wires(plane)) {
                        (w.segment == 0 ? seg0 : segn).insert(w.channel);
                    }
                }
            }
            for (int ch : segn) {
                if (!seg0.count(ch)) ++bad;
            }
        }
        return bad;
    }

    PlaneMap channel_index_invariant(const std::string& filename)
    {
        PlaneMap out;
        const auto store = WireSchema::load(filename.c_str());
        for (const auto& anode : store.anodes()) {
            for (const auto& face : store.faces(anode)) {
                for (const auto& plane : store.planes(face)) {
                    const auto wires = store.wires(plane);
                    ++out.planes;
                    int here = 0, first = -1;
                    for (size_t i = 0; i < wires.size(); ++i) {
                        if (wires[i].segment > 0) {
                            if (first < 0) first = (int) i;
                            ++here;
                        }
                    }
                    out.seg_gt0 += here;
                    if (here) {
                        ++out.planes_broken;
                        if (out.first_bad < 0 || first < out.first_bad) out.first_bad = first;
                        out.broken.push_back(PlaneRow{anode.ident, face.ident, plane.ident,
                                                      (int) wires.size(), here, first});
                    }
                }
            }
        }
        return out;
    }
}

TEST_CASE("pdvd doc31: wrapped planes omit their continuations' channels")
{
    // PDVD production geometry (protodunevd/params.jsonnet:245).  1568 of 12288
    // channels carry two segments, split at the CRU boundary y = +-1685 mm;
    // every anode contributes 98 orphan U and 98 orphan V wires in one of its
    // two faces.  W is never wrapped.
    const auto pdvd = count_orphans("protodunevd-wires-larsoft-v7-uvwfit.json.bz2");
    CHECK(pdvd.wires == 13856);
    CHECK(pdvd.orphans == 1568);
    CHECK(pdvd.planes == 16);

    // PDHD is affected harder (pdhd/params.jsonnet:187).  It is named here so a
    // future reader cannot mistake this for a PDVD-only defect.
    const auto pdhd = count_orphans("protodunehd-wires-larsoft-v1.json.bz2");
    CHECK(pdhd.wires == 22208);
    CHECK(pdhd.orphans == 6400);
    CHECK(pdhd.planes == 16);
}

TEST_CASE("pdvd doc31: unwrapped detectors have no orphans at all")
{
    // The negative control.  These two have no multi-segment channel anywhere,
    // so the lookup never misses, the ident map is never built, and the knob
    // cannot change their output no matter how it is set.  If either of these
    // ever becomes non-zero, the byte-identity argument for it is void.
    for (const std::string fname : {"sbnd-wires-geometry-v0206.json.bz2",
                                    "microboone-celltree-wires-v2.1.json.bz2"}) {
        const auto oc = count_orphans(fname);
        CHECK_MESSAGE(oc.orphans == 0, "unexpected wrapped continuations in ", fname.c_str());
        CHECK_MESSAGE(oc.planes == 0, "unexpected affected planes in ", fname.c_str());
    }
}

TEST_CASE("pdhd doc04: wrapped_channel_charge defaults ON and false is still reachable")
{
    // Flipped 2026-09-06 (doc pdhd/04 sec 9, owner decision: these are bug
    // fixes, so the fixed path is the default).  The old test asserted false
    // and said "this default is the whole of PDHD's protection" -- that reading
    // was exactly the defect: PDHD's clustering config omitted the key, so the
    // default WAS its behaviour, and the behaviour was the bug.
    //
    // SBND and uBooNE are not protected by this default and never were: they
    // have zero segment>0 wires, so the ident-resolved branch is unreachable
    // there.  That is the "unwrapped detectors have no orphans at all" case
    // above, which is now load-bearing for the byte-identity claim.
    PluginManager::instance().add("WireCellClus");
    auto icfg = Factory::lookup<IConfigurable>("BlobSampler", "doc31_wrapped_channel_probe");
    REQUIRE(icfg);

    auto cfg = icfg->default_configuration();
    REQUIRE_MESSAGE(cfg.isMember("wrapped_channel_charge"), "missing knob: wrapped_channel_charge");
    CHECK(cfg["wrapped_channel_charge"].asBool() == true);

    // The escape hatch must survive the flip.  Every config that threads this
    // knob now emits the key unconditionally (the `[if x then 'x']: true`
    // suppression idiom would have made `false` unreachable), so an explicit
    // false has to round-trip.
    cfg["wrapped_channel_charge"] = false;
    icfg->configure(cfg);
    CHECK(icfg->default_configuration()["wrapped_channel_charge"].asBool() == false);

    cfg["wrapped_channel_charge"] = true;
    icfg->configure(cfg);
    CHECK(icfg->default_configuration()["wrapped_channel_charge"].asBool() == true);
}

TEST_CASE("pdvd doc31 round4: channels[wire_index] is not a valid lookup on wrapped planes")
{
    // The invariant two retile sites assume.  This is NOT the same statement as
    // the orphan census above: an orphan is a channel the plane never lists at
    // all, while ANY segment>0 wire -- orphan or a strip that wraps back inside
    // its own plane -- shifts every later wire's position in `plane_channels`.
    // So this is the count that bounds the retile defect, and it is >= orphans.
    const auto pdvd = channel_index_invariant("protodunevd-wires-larsoft-v7-uvwfit.json.bz2");
    const auto pdhd = channel_index_invariant("protodunehd-wires-larsoft-v1.json.bz2");
    const auto sbnd = channel_index_invariant("sbnd-wires-geometry-v0206.json.bz2");
    const auto ubon = channel_index_invariant("microboone-celltree-wires-v2.1.json.bz2");

    MESSAGE("PDVD planes=" << pdvd.planes << " broken=" << pdvd.planes_broken
            << " seg>0=" << pdvd.seg_gt0 << " first_bad_wire_index=" << pdvd.first_bad);
    MESSAGE("PDHD planes=" << pdhd.planes << " broken=" << pdhd.planes_broken
            << " seg>0=" << pdhd.seg_gt0 << " first_bad_wire_index=" << pdhd.first_bad);
    for (const auto& r : pdvd.broken) {
        MESSAGE("PDVD broken plane a" << r.anode << " f" << r.face << " p" << r.plane
                << " nwires=" << r.nwires << " seg>0=" << r.seg_gt0
                << " first_bad=" << r.first_bad);
    }

    // The negative controls, and the reason SBND/uBooNE cannot be affected by
    // the retile mapping any more than they could by the sampler one: with no
    // segment>0 wire anywhere, channels[i] names wire i on every plane.
    CHECK(sbnd.seg_gt0 == 0);
    CHECK(sbnd.planes_broken == 0);
    CHECK(sbnd.planes == 6);      // 2 anodes x 1 face x 3 planes
    CHECK(ubon.seg_gt0 == 0);
    CHECK(ubon.planes_broken == 0);
    CHECK(ubon.planes == 3);      // 1 anode x 1 face x 3 planes

    // PDVD: 16 of 48 planes, 98 wrapped continuations each (= the 1568 above,
    // so on PDVD every segment>0 wire is also an orphan).
    CHECK(pdvd.planes == 48);
    CHECK(pdvd.planes_broken == 16);
    CHECK(pdvd.seg_gt0 == 1568);
    CHECK(pdvd.first_bad >= 0);

    // PDHD: 16 of 24 planes.  Note 11968 segment>0 wires against only 6400
    // orphans -- 5568 PDHD continuations wrap back INSIDE their own plane, so
    // they are invisible to the orphan census yet still shift plane_channels.
    // That is why this count, not the orphan count, is the one that bounds the
    // retile defect.  PDHD does not run the Steiner stage (no CreateSteinerGraph
    // in any pdhd config), so the retile sites are not reached there -- recorded
    // so the geometry fact is not mistaken for an exposure claim.
    CHECK(pdhd.planes == 24);
    CHECK(pdhd.planes_broken == 16);
    CHECK(pdhd.seg_gt0 == 11968);
    CHECK(pdhd.seg_gt0 > 6400);   // strictly more than the orphan count

    // PDVD's structure, pinned because doc 31 section 9 reasons from it.  Every
    // broken plane is 287 wires with 98 continuations, and they sit in ONE
    // contiguous band at one end: first_bad is 0 (band at the BOTTOM, so every
    // wire index is shifted by +98 and the top 98 are dropped) or 189 = 287-98
    // (band at the TOP, so indices 0..188 are correct and only the top 98 are
    // dropped).  Never interleaved -- which is what makes the consequence
    // predictable per plane instead of per wire.
    int n_bottom = 0, n_top = 0;
    for (const auto& r : pdvd.broken) {
        CHECK(r.nwires == 287);
        CHECK(r.seg_gt0 == 98);
        CHECK((r.first_bad == 0 || r.first_bad == 189));
        if (r.first_bad == 0) ++n_bottom; else ++n_top;
    }
    CHECK(n_bottom == 8);
    CHECK(n_top == 8);

    // The face doc 31's flagship track lives on (anode 4, face 0).  U is the
    // fully-shifted one and V is the top-band one, which is exactly the
    // asymmetry section 9.2 explains the terminal starvation with: below the
    // vertex the track's V wire is in the dropped top band while above it is
    // not, and U is wrong in both halves.
    int a4f0_u = -2, a4f0_v = -2;
    for (const auto& r : pdvd.broken) {
        if (r.anode == 4 && r.face == 0 && r.plane == 0) a4f0_u = r.first_bad;
        if (r.anode == 4 && r.face == 0 && r.plane == 1) a4f0_v = r.first_bad;
    }
    CHECK(a4f0_u == 0);
    CHECK(a4f0_v == 189);
}

TEST_CASE("pdhd doc04: wrapped_channel_activity defaults ON")
{
    // ImproveCluster_2 is the retiler the Steiner stage actually runs
    // (cm.improve_cluster_2 on PDVD, SBND and uBooNE alike).  Flipped to a true
    // default 2026-09-06 alongside wrapped_channel_charge (doc pdhd/04 sec 9):
    // same misconception, same call-site family, and a bug fix belongs on by
    // default.  SBND and uBooNE do not need the key omitted to stay unchanged --
    // with no segment>0 wire the resolved-by-ident branch cannot be taken.
    //
    // Only the default is checked: configure() runs NeedDV, which requires a
    // live DetectorVolumes instance this test has no business standing up.  The
    // ON path is exercised end-to-end by the doc-31 round-5 arms instead.
    PluginManager::instance().add("WireCellClus");
    auto icfg = Factory::lookup<IConfigurable>("ImproveCluster_2", "doc31_retile_probe");
    REQUIRE(icfg);

    auto cfg = icfg->default_configuration();
    REQUIRE_MESSAGE(cfg.isMember("wrapped_channel_activity"),
                    "missing knob: wrapped_channel_activity");
    CHECK(cfg["wrapped_channel_activity"].asBool() == true);
}

TEST_CASE("pdvd doc31 round5: every continuation's channel resolves within its own anode")
{
    // The precondition the retile-mapping fix rests on.  If this were non-zero
    // for any detector, IAnodePlane::channel() would hand back nullptr for
    // exactly the wires the fix exists to rescue, and the "fix" would silently
    // reproduce the bug.  Checked on all four production geometries so that a
    // future wires file cannot quietly break it.
    CHECK(unresolvable_channels("protodunevd-wires-larsoft-v7-uvwfit.json.bz2") == 0);
    CHECK(unresolvable_channels("protodunehd-wires-larsoft-v1.json.bz2") == 0);
    // Vacuously true for the unwrapped two, but pin them anyway: it is the
    // statement "no channel of this detector needs the anode-level lookup".
    CHECK(unresolvable_channels("sbnd-wires-geometry-v0206.json.bz2") == 0);
    CHECK(unresolvable_channels("microboone-celltree-wires-v2.1.json.bz2") == 0);
}

namespace {

    // doc pdvd/102 sec 6.  The prototype retile samples with calc_sampling_points
    // (CalcPoints.cxx:75-160), whose toolkit port is the "charge_stepped"
    // strategy.  Unlike "stepped" it DECIDES which wires to sample from their
    // charge (ChargeStepped::get_wire_charge / is_plane_bad), on top of the
    // per-point charge every strategy gets from make_dataset.  Both of those
    // lookups carry their own copy of the wrapped-continuation fix, so the
    // fix has to be pinned for this strategy separately.
    struct SampledCharges {
        size_t npts{0};
        std::set<double> charge[3];   // distinct charge_val per plane over all points
    };

    // PDHD apa1 exactly as the compiled PDHD clustering config builds it: one
    // sensitive face (faces[0] is null), the production wires file.
    IAnodePlane::pointer pdhd_apa1()
    {
        static IAnodePlane::pointer anode;
        if (anode) return anode;
        PluginManager& pm = PluginManager::instance();
        pm.add("WireCellAux");
        pm.add("WireCellGen");
        pm.add("WireCellClus");
        {
            auto icfg = Factory::lookup<IConfigurable>("WireSchemaFile", "doc102_pdhd_wires");
            auto cfg = icfg->default_configuration();
            cfg["filename"] = "protodunehd-wires-larsoft-v1.json.bz2";
            icfg->configure(cfg);
        }
        auto icfg = Factory::lookup<IConfigurable>("AnodePlane", "doc102_pdhd_apa1");
        auto cfg = icfg->default_configuration();
        cfg["ident"] = 1;
        cfg["nimpacts"] = 10;
        cfg["wire_schema"] = "WireSchemaFile:doc102_pdhd_wires";
        cfg["faces"][0] = Json::nullValue;
        cfg["faces"][1]["anode"] = 3520.945;
        cfg["faces"][1]["response"] = 3430.465;
        cfg["faces"][1]["cathode"] = 1.5875;
        icfg->configure(cfg);
        anode = Factory::find<IAnodePlane>("AnodePlane", "doc102_pdhd_apa1");
        return anode;
    }

    SampledCharges sample_charges(const IBlob::pointer& iblob, const std::string& name,
                                  const std::string& strategy, bool wrapped)
    {
        auto icfg = Factory::lookup<IConfigurable>("BlobSampler", name);
        auto cfg = icfg->default_configuration();
        Json::Value one(Json::objectValue);
        one["name"] = strategy;
        if (strategy == "charge_stepped") {
            one["disable_mix_dead_cell"] = false;   // the PR retile's setting (pdhd/clus.jsonnet bs_live_face)
        }
        cfg["strategy"] = Json::Value(Json::arrayValue);
        cfg["strategy"].append(one);
        cfg["extra"] = Json::Value(Json::arrayValue);
        cfg["extra"].append(".*charge_val");
        cfg["wrapped_channel_charge"] = wrapped;
        icfg->configure(cfg);
        auto bs = Factory::find<IBlobSampler>("BlobSampler", name);
        auto [ds, aux] = bs->sample_blob(iblob, 0);
        SampledCharges out;
        out.npts = ds.size_major();
        const std::string letters[3] = {"u", "v", "w"};
        for (int p = 0; p < 3 && out.npts; ++p) {
            auto arr = ds.get(letters[p] + "charge_val");
            REQUIRE(arr);
            for (double q : arr->elements<double>()) out.charge[p].insert(q);
        }
        return out;
    }
}

TEST_CASE("pdvd doc102: charge_stepped resolves wrapped continuations' charge exactly as stepped")
{
    auto anode = pdhd_apa1();
    REQUIRE(anode);
    IAnodeFace::pointer face;
    for (const auto& f : anode->faces()) {
        if (f) face = f;
    }
    REQUIRE(face);
    const auto& coords = face->raygrid();
    const auto planes = face->planes();
    REQUIRE(planes.size() == 3);

    // A 4 cm x 4 cm patch at (y, z) = (286, 115) cm, the middle of apa1's
    // wrapped band: there BOTH induction planes are continuations whose channel
    // their own plane does not list (U wires 400-799, V 348-747), and W is never
    // wrapped.
    std::vector<Point> pts;
    for (double y = 2840; y <= 2880; y += 1.0) {
        for (double z = 1132; z <= 1172; z += 1.0) {
            pts.emplace_back(3520.945 * units::mm, y * units::mm, z * units::mm);
        }
    }
    auto measures = RayGrid::make_measures(coords, pts);
    auto activities = RayGrid::make_activities(coords, measures);
    auto blobs = RayGrid::make_blobs(coords, activities);
    REQUIRE(!blobs.empty());
    auto width = [](const RayGrid::Strip& s) { return s.bounds.second - s.bounds.first; };
    size_t best = 0;
    for (size_t i = 1; i < blobs.size(); ++i) {
        const auto& a = blobs[i].strips();
        const auto& b = blobs[best].strips();
        if (width(a[2]) * width(a[4]) > width(b[2]) * width(b[4])) best = i;
    }
    const auto& shape = blobs[best];
    const auto& strips = shape.strips();
    REQUIRE(strips.size() == 5);

    // The premise, from the anode itself rather than the wires file: the blob's
    // U and V strips hold orphan continuations, its W strip holds none.
    int orphans[3] = {0, 0, 0};
    for (int p = 0; p < 3; ++p) {
        std::set<int> listed;
        for (const auto& ich : planes[p]->channels()) listed.insert(ich->ident());
        const auto& wires = planes[p]->wires();
        for (int wi = strips[2 + p].bounds.first; wi < strips[2 + p].bounds.second; ++wi) {
            if (wires[wi]->segment() > 0 && !listed.count(wires[wi]->channel())) ++orphans[p];
        }
    }
    MESSAGE("strip widths u/v/w " << width(strips[2]) << "/" << width(strips[3]) << "/" << width(strips[4])
            << ", orphans u/v/w " << orphans[0] << "/" << orphans[1] << "/" << orphans[2]);
    REQUIRE(orphans[0] > 0);
    REQUIRE(orphans[1] > 0);
    REQUIRE(orphans[2] == 0);
    // The prototype's all-wires branch (N_max * N_min <= 2500) is the one that
    // reads charge for every non-stepped wire.
    REQUIRE(std::max({width(strips[2]), width(strips[3]), width(strips[4])}) *
            std::min({width(strips[2]), width(strips[3]), width(strips[4])}) <= 2500);

    // Wrapped planes carry charge BELOW the 4000 threshold but non-zero (live,
    // not dead); W carries charge above it.  Keyed by the anode's channel, so a
    // continuation's activity sits under its real channel exactly as imaging
    // leaves it.
    const double qwrap = 2000, qw = 5000;
    ISlice::map_t activity;
    for (int p = 0; p < 3; ++p) {
        const auto& wires = planes[p]->wires();
        const int lo = std::max(0, strips[2 + p].bounds.first - 3);
        const int hi = std::min((int) wires.size(), strips[2 + p].bounds.second + 3);
        for (int wi = lo; wi < hi; ++wi) {
            auto ich = anode->channel(wires[wi]->channel());
            REQUIRE(ich);
            activity[ich] = ISlice::value_t(p == 2 ? qw : qwrap, 1.0f);
        }
    }
    auto slice = std::make_shared<Aux::SimpleSlice>(nullptr, 0, 0.0, 2 * units::ms, activity);
    IBlob::pointer iblob = std::make_shared<Aux::SimpleBlob>(0, 1.0f, 0.0f, shape, slice, face);

    const auto cs_on = sample_charges(iblob, "doc102_cs_on", "charge_stepped", true);
    const auto cs_off = sample_charges(iblob, "doc102_cs_off", "charge_stepped", false);
    const auto st_on = sample_charges(iblob, "doc102_st_on", "stepped", true);
    const auto st_off = sample_charges(iblob, "doc102_st_off", "stepped", false);
    MESSAGE("points: charge_stepped on/off " << cs_on.npts << "/" << cs_off.npts
            << ", stepped on/off " << st_on.npts << "/" << st_off.npts);
    REQUIRE(cs_on.npts > 0);
    REQUIRE(st_on.npts > 0);

    // With the fix, every sampled point of BOTH strategies reads the wrapped
    // planes' real charge: no point sees a wrapped plane as empty.
    const std::set<double> want_wrap{qwrap}, want_w{qw};
    CHECK(cs_on.charge[0] == want_wrap);
    CHECK(cs_on.charge[1] == want_wrap);
    CHECK(cs_on.charge[2] == want_w);
    CHECK(st_on.charge[0] == want_wrap);
    CHECK(st_on.charge[1] == want_wrap);
    CHECK(st_on.charge[2] == want_w);

    // Negative control: the legacy lookup loses it, for both strategies.
    CHECK((cs_off.charge[0].count(0.0) + cs_off.charge[1].count(0.0)) > 0);
    CHECK((st_off.charge[0].count(0.0) + st_off.charge[1].count(0.0)) > 0);

    // And charge_stepped's wire SELECTION depends on it.  A non-stepped wire at
    // 2000 is live-below-threshold and dropped; the 0 the legacy lookup reads is
    // "dead", which disable_mix_dead_cell=false keeps (BlobSampler.cxx
    // ChargeStepped, `charge != 0 || disable_mix_dead_cell`).  So without the fix
    // the retile would sample the wrapped band MORE densely than live charge
    // warrants.  Stepped never reads charge to select, so its count cannot move.
    CHECK(cs_off.npts > cs_on.npts);
    CHECK(st_off.npts == st_on.npts);
}
