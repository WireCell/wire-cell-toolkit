// Pins the LEGACY DEFAULT of every configuration knob added to clus/ by the
// SBND pattern-recognition port (sbnd_xin/docs/pr/*).
//
// WHY THIS TEST EXISTS.  Each of those knobs was introduced under one rule:
// with the key absent, the component must reproduce the pre-knob behavior
// byte-for-byte.  That rule is what makes the branch's A/B gates meaningful --
// a silently flipped default would make every "no behavior change" claim false
// while every gate still passed, because the gate compares two runs of the same
// (already-wrong) default.  A failure here is therefore a FEATURE: it means a
// default moved, which is a stop-and-ask, not something to update the test for.
//
// WHAT THIS TEST DOES NOT SAY.  It pins the C++ default only.  It says nothing
// about what any detector actually runs: SBND deliberately flips many of these
// ON in cfg/pgrapher/experiment/sbnd/wct-pr-perevt.jsonnet (doc 68, the single
// source of the SBND operating point).  That operating point is gated by the
// compiled-config proof (wcsonnet + abtest/cmp_cfg.sh), not by this file.  A
// green run here does NOT mean production is on the legacy path.
//
// Components are reached through the factory, which also pins that each one is
// still registered under the name the jsonnet uses.

#include "WireCellUtil/doctest.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/Units.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/PRSegmentFunctions.h"

using namespace WireCell;

// Fetch a component's default_configuration() by factory type name.
//
// The PluginManager::add is required even though wcdoctest-clus already links
// libWireCellClus: the WIRECELL_FACTORY registrars live in translation units
// that nothing in the test references, so the linker never pulls them in and
// Factory::lookup throws "No factory for class".  (This is why
// doctest_pattern_recognition.cxx's configure_components() wraps every lookup
// in a try/catch that logs "Skipping".)  add() resolves to the same
// build/clus/libWireCellClus.so that is already mapped, so this cannot pick up
// a stale installed copy.
//
// The instance NAME is load-bearing, not cosmetic.  Factory::lookup caches one
// instance per (type, name) for the whole process, and the DEFAULT-named
// instances of two of the components below are shared: uboone-mabc_config.json
// instantiates CreateSteinerGraph and TaggerCheckNeutrino with no "name", and
// doctest_pattern_recognition.cxx's configure_components() calls configure() on
// every entry it can resolve.  Asking for the unnamed instance here would read
// back whatever THAT test configured, in whatever order doctest happened to
// register the cases -- i.e. this file would silently stop testing defaults.
// A private name is virgin by construction.
static const char* kProbe = "knobdefaults_probe";

static Configuration defaults_of(const std::string& type)
{
    PluginManager::instance().add("WireCellClus");
    auto icfg = Factory::lookup<IConfigurable>(type, kProbe);
    REQUIRE(icfg);
    return icfg->default_configuration();
}

// Assert the key exists AND holds the expected value.  Checking the value
// without the isMember() guard would silently pass on a typo'd key, since an
// absent Json::Value converts to 0/false.
#define CHECK_KNOB_BOOL(cfg, key, want)                     \
    do {                                                    \
        REQUIRE_MESSAGE((cfg).isMember(key), "missing knob: " key); \
        CHECK((cfg)[key].asBool() == (want));               \
    } while (0)

#define CHECK_KNOB_NUM(cfg, key, want)                      \
    do {                                                    \
        REQUIRE_MESSAGE((cfg).isMember(key), "missing knob: " key); \
        CHECK((cfg)[key].asDouble() == doctest::Approx(want)); \
    } while (0)

// ---------------------------------------------------------------------------
// TaggerCheckNeutrino -- the knob hub.  Almost every doc pr/* fix threads its
// switch through this component.
// ---------------------------------------------------------------------------

TEST_CASE("clus knob defaults: TaggerCheckNeutrino switches are all OFF")
{
    auto cfg = defaults_of("TaggerCheckNeutrino");

    // Every boolean added by the port defaults to the legacy path.
    CHECK_KNOB_BOOL(cfg, "dir_weak_use_score", false);       // legacy raw dir_weak reads
    CHECK_KNOB_BOOL(cfg, "proton_dir_vote", false);          // doc pr/8
    CHECK_KNOB_BOOL(cfg, "endpoint_trim_retry", false);      // doc pr/9 F1
    CHECK_KNOB_BOOL(cfg, "iso_endpoint", false);             // doc pr/25 isochronous trunk
    CHECK_KNOB_BOOL(cfg, "sp_dedx_use_recomb_model", false); // doc pr/2 2e(i)
    CHECK_KNOB_BOOL(cfg, "nu_skip_cosmic", false);           // doc pr/3 sec 8
    CHECK_KNOB_BOOL(cfg, "nu_skip_cosmic_bundle", false);
    CHECK_KNOB_BOOL(cfg, "skip_cosmic_companions", false);   // doc pr/20 I P4
    // doc pr/26 sec 9: store singlephoton_tagger()'s verdict in
    // TaggerInfo::photon_flag (prototype NeutrinoID.cxx:271).  The C++ default
    // stays the legacy gap even though SBND's cfg now flips it on, so every
    // other detector's tagger ntuple keeps a constant-0 photon_flag branch.
    CHECK_KNOB_BOOL(cfg, "sp_photon_flag", false);

    // Numeric knobs whose legacy value is the INERT one: 0 disables the guard,
    // so an absent key leaves the code path untouched.
    CHECK_KNOB_NUM(cfg, "fit_vertex_min_seg_length", 0.0);   // 0 = all segments enter the fit
    CHECK_KNOB_NUM(cfg, "cathode_x", 0.0);
    CHECK_KNOB_NUM(cfg, "cathode_kink_xcut", 0.0);           // 0 = kink search sees every fit point
    CHECK_KNOB_NUM(cfg, "shower_topo_demote_len", 0.0);      // 0 = long segments stay shower-eligible
    CHECK_KNOB_NUM(cfg, "nu_skip_cosmic_bundle_min_length", 0.0);
    CHECK_KNOB_NUM(cfg, "cosmic_companion_min_length", 0.0);
    // low >= high disables the beam gate; both 0 is the disabled state.
    CHECK_KNOB_NUM(cfg, "beam_window_low", 0.0);
    CHECK_KNOB_NUM(cfg, "beam_window_high", 0.0);
}

TEST_CASE("clus knob defaults: TaggerCheckNeutrino literals are the uBooNE values")
{
    // These knobs did not switch behavior on or off -- they lifted a hard-coded
    // uBooNE constant into config (doc pr/2 sec 2e).  The default must remain
    // the literal that was in the code, or every detector that does not
    // override it silently changes.
    auto cfg = defaults_of("TaggerCheckNeutrino");

    // The two distinct MIP dQ/dx scales that commit 1628328e untangled.  They
    // are NOT interchangeable: 50e3 is the flat-template amplitude handed to
    // do_track_comp, 43e3 the scale for median-dQ/dx thresholds.
    CHECK_KNOB_NUM(cfg, "mip_dqdx", 50000.0);
    CHECK_KNOB_NUM(cfg, "mip_dqdx_median", 43000.0);

    CHECK_KNOB_NUM(cfg, "proton_dir_score_max", 0.25);
    CHECK_KNOB_NUM(cfg, "proton_dir_asym_min", 1.3);

    // Isochronous-endpoint geometry (doc pr/25 round 3), cm where dimensioned.
    CHECK_KNOB_NUM(cfg, "iso_endpoint_min_length", 40.0);
    CHECK_KNOB_NUM(cfg, "iso_endpoint_max_xext", 25.0);
    CHECK_KNOB_NUM(cfg, "iso_endpoint_xext_frac", 0.35);
    CHECK_KNOB_NUM(cfg, "iso_endpoint_xext_quantile", 0.02);
    CHECK_KNOB_NUM(cfg, "iso_endpoint_tube_radius", 4.0);
    CHECK_KNOB_NUM(cfg, "iso_endpoint_min_aspect", 0.12);

    // Detector-extent literals (doc pr/2 sec 2e(iv)); uBooNE y=+117 top.
    CHECK_KNOB_NUM(cfg, "cosmic_y_top_main", 100.0);
    CHECK_KNOB_NUM(cfg, "cosmic_y_top_strict", 102.0);
    CHECK_KNOB_NUM(cfg, "cosmic_y_top_loose", 80.0);
    CHECK_KNOB_NUM(cfg, "cosmic_y_small_piece", 50.0);
    CHECK_KNOB_NUM(cfg, "vertex_z_prior_scale", 200.0);  // uBooNE 1037 cm detector

    // Charge -> energy calibration (doc pr/2 sec 2e(iii)); uBooNE.
    CHECK_KNOB_NUM(cfg, "kine_fudge_factor", 0.95);
    CHECK_KNOB_NUM(cfg, "kine_recom_factor", 0.7);
    CHECK_KNOB_NUM(cfg, "kine_shower_fudge_factor", 0.8);
    CHECK_KNOB_NUM(cfg, "kine_shower_recom_factor", 0.5);
    CHECK_KNOB_NUM(cfg, "kine_proton_recom_factor", 0.35);
    CHECK_KNOB_NUM(cfg, "kine_plane_asym_switch", 0.04);
    CHECK_KNOB_NUM(cfg, "kine_w_value", 23.6);           // eV per electron-ion pair
    CHECK_KNOB_NUM(cfg, "sp_mean_dedx_cut", 2.3);        // MeV/cm
}

TEST_CASE("clus knob defaults: TaggerCheckNeutrino vector knobs round-trip")
{
    // Vector-valued knobs are the easy ones to break: an append() onto a
    // pre-populated array yields a key that looks set but is rejected
    // downstream (get_dir3 wants exactly 3 elements).  Pin length AND content.
    auto cfg = defaults_of("TaggerCheckNeutrino");

    REQUIRE(cfg["kine_plane_weights"].isArray());
    REQUIRE(cfg["kine_plane_weights"].size() == 3);
    CHECK(cfg["kine_plane_weights"][0].asDouble() == doctest::Approx(0.25));  // U
    CHECK(cfg["kine_plane_weights"][1].asDouble() == doctest::Approx(0.25));  // V
    CHECK(cfg["kine_plane_weights"][2].asDouble() == doctest::Approx(1.0));   // W

    // SSM beam-line references (doc pr/2 sec 2e(i-a)): uBooNE BNB target and
    // NuMI absorber unit vectors.
    REQUIRE(cfg["ssm_target_dir"].isArray());
    REQUIRE(cfg["ssm_target_dir"].size() == 3);
    CHECK(cfg["ssm_target_dir"][0].asDouble() == doctest::Approx(0.46));
    CHECK(cfg["ssm_target_dir"][1].asDouble() == doctest::Approx(0.05));
    CHECK(cfg["ssm_target_dir"][2].asDouble() == doctest::Approx(0.885));

    REQUIRE(cfg["ssm_absorber_dir"].isArray());
    REQUIRE(cfg["ssm_absorber_dir"].size() == 3);
    CHECK(cfg["ssm_absorber_dir"][0].asDouble() == doctest::Approx(0.33));
    CHECK(cfg["ssm_absorber_dir"][1].asDouble() == doctest::Approx(0.75));
    CHECK(cfg["ssm_absorber_dir"][2].asDouble() == doctest::Approx(-0.59));

    // Muon median-dQ/dx-vs-length envelope {c0, c1, pivot_cm, power}.  pivot is
    // in CM here (18.0), which is what makes PatternAlgorithms::muon_dqdx_cut()
    // and muon_dqdx_cut_cm() agree; see the check further below.
    REQUIRE(cfg["muon_dqdx_curve"].isArray());
    REQUIRE(cfg["muon_dqdx_curve"].size() == 4);
    CHECK(cfg["muon_dqdx_curve"][0].asDouble() == doctest::Approx(0.8866));
    CHECK(cfg["muon_dqdx_curve"][1].asDouble() == doctest::Approx(0.9533));
    CHECK(cfg["muon_dqdx_curve"][2].asDouble() == doctest::Approx(18.0));
    CHECK(cfg["muon_dqdx_curve"][3].asDouble() == doctest::Approx(0.4234));
}

// ---------------------------------------------------------------------------
// CreateSteinerGraph -- the one knob on this branch whose LEGACY default is
// true.  Worth its own case precisely because it inverts the usual rule: a
// well-meaning "make all new knobs default false" sweep would break it.
// ---------------------------------------------------------------------------

TEST_CASE("clus knob defaults: CreateSteinerGraph replace defaults TRUE")
{
    auto cfg = defaults_of("CreateSteinerGraph");
    // true = the historical always-rebuild path, so configs without the key are
    // byte-identical.  replace=false is the opt-in that keeps an existing
    // steiner_graph so cached GraphAlgorithms refs cannot dangle
    // (SBND evts 54095/163543).
    CHECK_KNOB_BOOL(cfg, "replace", true);
}

// ---------------------------------------------------------------------------
// The clustering visitors added or extended by the port.
// ---------------------------------------------------------------------------

TEST_CASE("clus knob defaults: ClusteringUnmergeBundle")
{
    auto cfg = defaults_of("ClusteringUnmergeBundle");
    CHECK_KNOB_BOOL(cfg, "restore_demoted_mains", false);  // doc pr/20 Part I
    // require_provenance=true RAISES on a legacy pctree rather than silently
    // reproducing pre-pr/20 behavior, so it must stay opt-in.
    CHECK_KNOB_BOOL(cfg, "require_provenance", false);
    REQUIRE(cfg.isMember("wasmain_aname"));
    CHECK(cfg["wasmain_aname"].asString() == "real_cluster_was_main");
}

TEST_CASE("clus knob defaults: ClusteringProtectBundle")
{
    auto cfg = defaults_of("ClusteringProtectBundle");
    // The cathode-rejoin pass is disabled at its default: xcut <= 0 turns the
    // whole pass off, and perp <= 0 turns off the direction fallback within it.
    CHECK_KNOB_NUM(cfg, "cathode_rejoin_xcut", 0.0);
    CHECK_KNOB_NUM(cfg, "cathode_rejoin_perp", 0.0);
    CHECK_KNOB_NUM(cfg, "cathode_x", 0.0);
    // Geometry that only matters once the pass is enabled; pinned so enabling
    // it in one detector cannot silently retune another.
    CHECK_KNOB_NUM(cfg, "cathode_rejoin_dyz", 4 * units::cm);
    CHECK_KNOB_NUM(cfg, "cathode_rejoin_dis", 8 * units::cm);
    CHECK_KNOB_NUM(cfg, "cathode_rejoin_angle", 20.0);  // degrees
    CHECK_KNOB_NUM(cfg, "cathode_rejoin_dir_radius", 15 * units::cm);
    CHECK_KNOB_NUM(cfg, "cathode_rejoin_dir_npts", 20);
    CHECK_KNOB_BOOL(cfg, "skip_convicted", true);
}

TEST_CASE("clus knob defaults: ClusteringExamineBundles")
{
    auto cfg = defaults_of("ClusteringExamineBundles");

    // Off => the provenance array is never created => the "perblob" key set is
    // unchanged and every detector's output is byte-identical
    // (doc pr/20 Part I P1).
    CHECK_KNOB_BOOL(cfg, "save_bundle_main_provenance", false);

    // The four pre-existing knobs, pinned in the same pass.  This component had
    // no default_configuration() at all until doc 70 added one, so NONE of them
    // round-tripped and all five were invisible to a config dump.
    CHECK_KNOB_BOOL(cfg, "use_ctpc", true);
    CHECK_KNOB_BOOL(cfg, "use_flash_t0", false);
    CHECK_KNOB_BOOL(cfg, "flags_from_longest", false);  // historical arbitrary-donor merge
    CHECK_KNOB_NUM(cfg, "flash_t0_window", 80 * units::ns);
    REQUIRE(cfg.isMember("graph_name"));
    CHECK(cfg["graph_name"].asString() == "relaxed");
}

// ---------------------------------------------------------------------------
// The three tagger-check components: demoted-main participation (doc pr/20).
// ---------------------------------------------------------------------------

TEST_CASE("clus knob defaults: tagger checks do not evaluate demoted mains")
{
    for (const std::string tn : {"TaggerCheckFC", "TaggerCheckSTM", "TaggerCheckTGM"}) {
        CAPTURE(tn);
        auto cfg = defaults_of(tn);
        CHECK_KNOB_BOOL(cfg, "evaluate_demoted_mains", false);
    }
    // TGM alone can also exempt demoted mains from its pair veto.
    auto tgm = defaults_of("TaggerCheckTGM");
    CHECK_KNOB_BOOL(tgm, "exempt_demoted_main_pairs", false);
}

// ---------------------------------------------------------------------------
// The POD option structs.  These carry defaults that never pass through
// default_configuration(), so a factory round-trip cannot see them: the
// in-class initializer IS the legacy contract at every call site that takes
// the struct by default argument.
// ---------------------------------------------------------------------------

TEST_CASE("clus knob defaults: TrackPidOptions in-class initializers")
{
    Clus::PR::TrackPidOptions o;  // == the default argument at every call site
    CHECK(o.mip_dqdx == doctest::Approx(50000 / units::cm));
    CHECK(o.proton_dir_vote == false);
    CHECK(o.proton_dir_score_max == doctest::Approx(0.25));
    CHECK(o.proton_dir_asym_min == doctest::Approx(1.3));
    CHECK(o.endpoint_trim_retry == false);
    CHECK(o.start_n == 1);
    CHECK(o.end_n == 1);
}

TEST_CASE("clus knob defaults: KineChargeOptions in-class initializers")
{
    Clus::PR::KineChargeOptions o;
    CHECK(o.fudge_factor == doctest::Approx(0.95));
    CHECK(o.recom_factor == doctest::Approx(0.7));
    CHECK(o.shower_fudge_factor == doctest::Approx(0.8));
    CHECK(o.shower_recom_factor == doctest::Approx(0.5));
    CHECK(o.proton_recom_factor == doctest::Approx(0.35));
    CHECK(o.plane_weights[0] == doctest::Approx(0.25));
    CHECK(o.plane_weights[1] == doctest::Approx(0.25));
    CHECK(o.plane_weights[2] == doctest::Approx(1.0));
    CHECK(o.plane_asym_switch == doctest::Approx(0.04));
    CHECK(o.w_value == doctest::Approx(23.6));
}

// ---------------------------------------------------------------------------
// muon_dqdx_cut vs muon_dqdx_cut_cm.
//
// PatternAlgorithms exposes the same envelope in two unit conventions because
// the ssm_tagger divides its sg_length by units::cm at the source.  The header
// claims the two are bit-identical, resting on pivot/units::cm == 18.0 being
// exact.  Six-plus call sites depend on that; a retune of the curve to a
// non-decade pivot would silently break it, so pin it here rather than trust
// the comment.
// ---------------------------------------------------------------------------

TEST_CASE("clus knob defaults: muon dQ/dx envelope agrees in both unit conventions")
{
    Clus::PR::PatternAlgorithms algo;  // service-free; the default curve is in force

    for (double length_cm : {1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 18.0, 25.0, 37.0,
                             50.0, 100.0, 123.4, 250.0, 500.0}) {
        CAPTURE(length_cm);
        // Exact equality, not Approx: the claim in the header is bit-identity.
        CHECK(algo.muon_dqdx_cut(length_cm * units::cm) == algo.muon_dqdx_cut_cm(length_cm));
    }

    // The envelope is a falling function of length (a long muon's median dQ/dx
    // sits closer to MIP), which is the property every caller relies on.
    CHECK(algo.muon_dqdx_cut(10 * units::cm) > algo.muon_dqdx_cut(100 * units::cm));
}

// ---------------------------------------------------------------------------
// PrDisplayDump -- the PR event-display calib dump (sbnd_xin/docs/pr/26).
//
// This component is READ-ONLY: it walks the PR graph, the steiner point clouds
// and the fitted 2-D charge and writes one JSON.  Its safety claim is therefore
// not about a behavior default but about REACHABILITY: it exists under the name
// the SBND cm_by_name entry uses, and it is inert unless that name appears in
// pipeline_names.  If the factory registration were lost, the compiled config
// would still be valid and the job would still run -- it would just die at
// component lookup, at the end of a 12-second-per-event PR chain.  So pin the
// registration.
//
// The dQdx_scale/offset defaults are pinned too.  They are not used to TRANSFORM
// anything in the dump (charges go out raw, unlike the Bee layers), only to be
// recorded in "meta" so a viewer can reproduce the Bee colouring.  A silent
// change would therefore not fail any gate -- it would just make the display
// disagree with Bee for no visible reason.
// ---------------------------------------------------------------------------

TEST_CASE("clus knob defaults: PrDisplayDump is registered and inert by default")
{
    auto cfg = defaults_of("PrDisplayDump");

    // Same convention as the Bee track_fit layer and SbndPrMagnifyTrackingVisitor.
    CHECK_KNOB_NUM(cfg, "dQdx_scale", 0.1);
    CHECK_KNOB_NUM(cfg, "dQdx_offset", -1000);
    CHECK_KNOB_NUM(cfg, "nticks", 3427);

    // 0 = keep every fitted 2-D cell.  A non-zero default would silently thin
    // the six projection panels.
    CHECK_KNOB_NUM(cfg, "proj_charge_min", 0);

    // Compact JSON: these dumps are multi-MB per event and are read by a
    // machine, not a human.
    CHECK_KNOB_BOOL(cfg, "pretty", false);

    // No anodes and no detector volumes by default: the component must be
    // constructible without geometry, which is what makes it safe to leave
    // defined in cm_by_name for every job whether or not it is named.
    REQUIRE(cfg.isMember("anodes"));
    CHECK(cfg["anodes"].size() == 0);
    REQUIRE(cfg.isMember("detector_volumes"));
    CHECK(cfg["detector_volumes"].asString() == "");

    // The grouping it reads.  "live" is the only grouping that carries a
    // TrackFitting slot; a typo here would produce a silently empty dump.
    REQUIRE(cfg.isMember("grouping"));
    CHECK(cfg["grouping"].asString() == "live");
}
