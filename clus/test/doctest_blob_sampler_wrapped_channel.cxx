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
