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

TEST_CASE("pdvd doc31: wrapped_channel_charge defaults OFF and round-trips")
{
    PluginManager::instance().add("WireCellClus");
    auto icfg = Factory::lookup<IConfigurable>("BlobSampler", "doc31_wrapped_channel_probe");
    REQUIRE(icfg);

    auto cfg = icfg->default_configuration();
    REQUIRE_MESSAGE(cfg.isMember("wrapped_channel_charge"), "missing knob: wrapped_channel_charge");
    // The load-bearing assertion: PDHD's config never mentions this key, so
    // this default is the whole of its protection.
    CHECK(cfg["wrapped_channel_charge"].asBool() == false);

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
