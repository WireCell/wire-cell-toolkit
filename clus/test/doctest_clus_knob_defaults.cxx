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
#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/TrackFittingPresets.h"

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
    CHECK_KNOB_BOOL(cfg, "v3_extension_guard", false);       // doc pr/24 sec 18 round 5
    CHECK_KNOB_BOOL(cfg, "sp_dedx_use_recomb_model", false); // doc pr/2 2e(i)
    CHECK_KNOB_BOOL(cfg, "nu_skip_cosmic", false);           // doc pr/3 sec 8
    CHECK_KNOB_BOOL(cfg, "nu_skip_cosmic_bundle", false);
    CHECK_KNOB_BOOL(cfg, "skip_cosmic_companions", false);   // doc pr/20 I P4
    CHECK_KNOB_BOOL(cfg, "nu_fallback_demoted_mains", false); // docs/73 sec 12 round 3
    // doc pr/26 sec 9: store singlephoton_tagger()'s verdict in
    // TaggerInfo::photon_flag (prototype NeutrinoID.cxx:271).  The C++ default
    // stays the legacy gap even though SBND's cfg now flips it on, so every
    // other detector's tagger ntuple keeps a constant-0 photon_flag branch.
    CHECK_KNOB_BOOL(cfg, "sp_photon_flag", false);
    // sbnd_xin/docs/75 -- nue/single-photon tagger FV consistency + the
    // DL-swap main-flag leak fix.  Both false = legacy path, byte-identical.
    CHECK_KNOB_BOOL(cfg, "nue_sp_consistent_fv", false);
    CHECK_KNOB_BOOL(cfg, "nu_selected_as_main_snapshot_all", false);
    // doc pr/40: track (proton/pion/muon) mis-identified as electron.
    CHECK_KNOB_BOOL(cfg, "track_pid_persist_dqdx", false);      // F1
    CHECK_KNOB_BOOL(cfg, "shower_reclass_dqdx_guard", false);   // F2
    CHECK_KNOB_BOOL(cfg, "shower_topo_dqdx_guard", false);      // F3
    // doc pr/40 round 2: two follow-on defects from the pr/40 fix round.
    CHECK_KNOB_BOOL(cfg, "track_pid_persist_4mom", false);            // F4
    CHECK_KNOB_BOOL(cfg, "shower_proton_daughter_pion", false);       // F5
    CHECK_KNOB_BOOL(cfg, "reclass_never_computed_ke_floor", false);   // F6
    // doc pr/40 round 4: two follow-on defects from round 2/3's F5 fix.
    CHECK_KNOB_BOOL(cfg, "shower_proton_daughter_pion_dissolve", false); // F7
    CHECK_KNOB_BOOL(cfg, "muon_multi_proton_pion", false);               // F8
    // doc pr/40 round 5: muon mis-identified as electron, three mechanisms.
    CHECK_KNOB_BOOL(cfg, "track_pid_persist_dqdx_electron_guard", false);     // F9
    CHECK_KNOB_BOOL(cfg, "shower_connect_main_vertex_straight_guard", false); // F10
    CHECK_KNOB_BOOL(cfg, "shower_traj_straight_guard", false);                // F11
    // doc pr/40 round 6: boundary-level fixes the round-5 measurement demanded.
    CHECK_KNOB_BOOL(cfg, "shower_absorb_track_guard", false);                 // F12
    CHECK_KNOB_BOOL(cfg, "michel_stem_muon_rescue", false);                   // F14
    // doc pr/74 round 2 -- P1 cascade guard + P2 Michel-terminal check.
    CHECK_KNOB_BOOL(cfg, "shower_in_cascade_guard", false);                   // P1
    CHECK_KNOB_NUM(cfg, "shower_in_max_len", 40.0);                           // cm; inert while P1 off
    CHECK_KNOB_NUM(cfg, "shower_in_mip_hi", 1.3);                             // ratio; inert while P1 off
    CHECK_KNOB_BOOL(cfg, "michel_stem_michel_check", false);                  // P2
    CHECK_KNOB_NUM(cfg, "michel_stem_max_far_len", 40.0);                     // cm; inert while P2 off
    CHECK_KNOB_BOOL(cfg, "shower_stem_backfill", false);                      // K4
    CHECK_KNOB_NUM(cfg, "stem_backfill_max_len", 30.0);                       // cm; inert while K4 off
    CHECK_KNOB_NUM(cfg, "stem_backfill_mip_lo", 0.75);                        // ratio; inert while K4 off
    CHECK_KNOB_NUM(cfg, "stem_backfill_mip_hi", 3.5);                         // ratio; inert while K4 off
    CHECK_KNOB_NUM(cfg, "stem_backfill_min_shower_len", 40.0);                // cm; inert while K4 off
    CHECK_KNOB_BOOL(cfg, "shower_conn3_unreachable", false);                  // K5 (pr/65 rung 2)
    CHECK_KNOB_NUM(cfg, "conn3_unreachable_min_len", 10.0);                   // cm; inert while K5 off
    CHECK_KNOB_NUM(cfg, "conn3_stitch_max", 0.0);                             // cm; doc pr/84 r2 F3; 0 = off
    CHECK_KNOB_BOOL(cfg, "shower_dedup_start_seg", false);                    // doc pr/84 r3 S1; off = legacy
    CHECK_KNOB_BOOL(cfg, "shower_endpoint_skip_orphan_vtx", false);           // doc pr/91 r1 F1; off = legacy end_point search
    CHECK_KNOB_BOOL(cfg, "shower_walk_visited_parity", false);               // doc pr/91 round 3; off = legacy has_node()-gated frontier
    // doc pr/40 round 9 -- straight-track PID guard family + B2 bridge.
    CHECK_KNOB_BOOL(cfg, "shower_connect_from_vertices_straight_guard", false);  // r8 Part A
    CHECK_KNOB_BOOL(cfg, "shower_connect_start_seg_straight_guard", false);      // r7 c2c (D1 re-target)
    CHECK_KNOB_BOOL(cfg, "examine_direction_dirsign_shower_in_guard", false);    // r7 c2a (D2 re-scope)
    CHECK_KNOB_BOOL(cfg, "daughter_shower_angle_reclass_straight_guard", false); // r7 c2b
    CHECK_KNOB_BOOL(cfg, "shower_topo_reexam_straight_guard", false);            // r7 c1 safety net
    CHECK_KNOB_NUM(cfg, "sfv_kink_max", 25.0);                                   // deg; inert while guards off
    CHECK_KNOB_BOOL(cfg, "shower_nv_bridge_track", false);                       // B2
    CHECK_KNOB_NUM(cfg, "shower_nv_bridge_max_gap", 1.8);                        // cm; inert while B2 off
    // doc pr/92 -- stray-satellite drop from kine/PF.
    CHECK_KNOB_BOOL(cfg, "kine_drop_stray_satellites", false);                   // off = legacy Enu sum
    CHECK_KNOB_NUM(cfg, "kine_sat_min_energy", 20.0);                            // MeV; inert while pr/92 off
    CHECK_KNOB_NUM(cfg, "kine_sat_prox_max", 8.0);                               // cm; inert while pr/92 off
    CHECK_KNOB_NUM(cfg, "kine_sat_angle_bad", 60.0);                             // deg; inert while pr/92 off
    CHECK_KNOB_NUM(cfg, "kine_sat_angle_main", 45.0);                            // deg; inert while pr/92 off
    CHECK_KNOB_NUM(cfg, "kine_sat_far_dis", 90.0);                               // cm; inert while pr/92 off
    CHECK_KNOB_NUM(cfg, "kine_sat_axis_dis_cut", 30.0);                          // cm; inert while pr/92 off
    CHECK_KNOB_NUM(cfg, "kine_sat_cont_kink", 25.0);                             // deg; inert while pr/92 off
    CHECK_KNOB_NUM(cfg, "kine_sat_track_max_nseg", 3.0);                         // count; pr/92 r2 topology split
    CHECK_KNOB_NUM(cfg, "kine_sat_em_far_dis", 150.0);                           // cm; pr/92 r2 EM far-drop
    CHECK_KNOB_BOOL(cfg, "shower_traj_michel_stem", false);                   // doc pr/74 round 4 K6
    CHECK_KNOB_NUM(cfg, "michel_stem_traj_min_len", 15.0);                    // cm; inert while K6 off
    CHECK_KNOB_NUM(cfg, "michel_stem_traj_max_len", 45.0);                    // cm; inert while K6 off
    CHECK_KNOB_NUM(cfg, "michel_stem_traj_mip_lo", 1.3);                      // x MIP median; inert while K6 off
    CHECK_KNOB_NUM(cfg, "michel_stem_traj_max_far_len", 40.0);                // cm; inert while K6 off
    CHECK_KNOB_NUM(cfg, "michel_stem_traj_min_kink_deg", 40.0);               // deg; inert while K6 off
    // doc pr/44: long-muon pseudo-shower keeps its muon start segment.
    CHECK_KNOB_BOOL(cfg, "shower_long_muon_keep_type", false);
    // doc pr/40 round 10: Bragg-PID-confident muon/proton start segment kept.
    CHECK_KNOB_BOOL(cfg, "shower_bragg_protect_start_segment", false);
    // doc pr/93 round 3 -- "electron" is really tracks / hadronic-pi0 shower.
    CHECK_KNOB_BOOL(cfg, "shower_reclass_case_b_dqdx_guard", false);          // Cause A (55595)
    CHECK_KNOB_BOOL(cfg, "shower_accept_pid_guard", false);                   // Cause B (348471, 69314)
    CHECK_KNOB_NUM(cfg, "shower_pid_guard_min_len", 50.0);                    // cm; inert while Cause A/B off
    CHECK_KNOB_BOOL(cfg, "shower_vote_track_pid_counts", false);              // Cause C (292643)
    CHECK_KNOB_BOOL(cfg, "shower_cone_absorb_guard", false);               // Cause D (315167)
    // doc pr/93 round 4 -- PF-hierarchy fine-tunes + 137238 cross-cluster muon.
    CHECK_KNOB_BOOL(cfg, "shower_detach_track_stem", false);                  // r4 (348471, 292643)
    // doc pr/99 round 2 -- shower_ghost_member_drop family.
    CHECK_KNOB_BOOL(cfg, "shower_ghost_member_drop", false);                  // pr/99 r2 (395148); false = byte-identical
    CHECK_KNOB_NUM(cfg, "shower_ghost_overlap_frac", 0.7);                    // inert while drop off
    CHECK_KNOB_NUM(cfg, "shower_ghost_dqdx_ratio", 0.25);                     // inert while drop off
    CHECK_KNOB_NUM(cfg, "shower_ghost_min_len", 10.0);                        // cm; inert while drop off

    // doc pr/99 round 3 -- kine-charge ownership + A5 hadronic tag.
    CHECK_KNOB_BOOL(cfg, "kine_charge_dedup", false);                         // pr/99 r3 C1 (168596); false = byte-identical
    CHECK_KNOB_BOOL(cfg, "kine_charge_rebuild", false);                       // pr/99 r3 C1b; false = byte-identical
    // doc pr/101 -- Enu accounting round (K1-K5); all defaults = legacy.
    CHECK_KNOB_BOOL(cfg, "kine_charge_track_ctx", false);                     // pr/101 K1 (37112); false = byte-identical
    CHECK_KNOB_BOOL(cfg, "kine_mass_rules", false);                           // pr/101 K2; false = byte-identical
    CHECK_KNOB_BOOL(cfg, "kine_hadronic_dqdx", false);                        // pr/101 K3; false = byte-identical
    CHECK_KNOB_NUM(cfg, "kine_long_muon_mode", 0);                            // pr/101 K4; 0 = legacy dQdx
    CHECK_KNOB_NUM(cfg, "kine_long_muon_ratio_lo", 0.3);                      // inert unless mode 2
    CHECK_KNOB_NUM(cfg, "kine_long_muon_ratio_hi", 0.5);                      // inert unless mode 2
    CHECK_KNOB_BOOL(cfg, "long_muon_range_empty_chain_fallback", false);      // doc 84 r1 P1; false = byte-identical
    CHECK_KNOB_BOOL(cfg, "long_muon_members_geometry", false);                // doc 84 r2; false = byte-identical
    CHECK_KNOB_BOOL(cfg, "long_muon_cathode_bridge", false);                  // doc 84 r2; false = byte-identical
    CHECK_KNOB_NUM(cfg, "long_muon_cathode_bridge_x", 0.0);                   // cm; inert unless bridge on
    CHECK_KNOB_NUM(cfg, "long_muon_cathode_bridge_xcut", 6.0);                // cm; inert unless bridge on
    CHECK_KNOB_NUM(cfg, "long_muon_cathode_bridge_gap", 20.0);                // cm; inert unless bridge on
    CHECK_KNOB_NUM(cfg, "long_muon_cathode_bridge_angle", 25.0);              // deg; inert unless bridge on
    CHECK_KNOB_NUM(cfg, "long_muon_cathode_bridge_lever", 5.0);               // doc 84 r4 G1; cm; 5.0 == the round-2 hardcode
    CHECK_KNOB_BOOL(cfg, "long_muon_cathode_bridge_track_partner", false);    // doc 84 r4 G2; false = legacy |13|-only partner showers
    CHECK_KNOB_NUM(cfg, "long_muon_cathode_bridge_short_gap", 0.0);           // doc 84 r4 G3; cm; 0 = waiver off
    CHECK_KNOB_NUM(cfg, "long_muon_cathode_bridge_short_gap_angle", 10.0);    // deg; inert while short_gap 0
    CHECK_KNOB_NUM(cfg, "long_muon_cathode_bridge_short_gap_len", 50.0);      // cm; inert while short_gap 0
    CHECK_KNOB_BOOL(cfg, "kine_mainvtx_used_guard", false);                   // pr/101 K5; false = byte-identical
    CHECK_KNOB_BOOL(cfg, "shower_hadronic_tag", false);                       // pr/99 r3 A5 (315167/395148); false = byte-identical
    CHECK_KNOB_NUM(cfg, "shower_hadronic_min_len", 10.0);                     // cm; inert while tag off
    CHECK_KNOB_NUM(cfg, "shower_hadronic_scan_len", 30.0);                    // cm; inert while tag off
    CHECK_KNOB_NUM(cfg, "shower_hadronic_bin", 3.0);                          // cm; inert while tag off
    CHECK_KNOB_NUM(cfg, "shower_hadronic_r_cyl", 8.0);                        // cm; inert while tag off
    CHECK_KNOB_NUM(cfg, "shower_hadronic_r_core", 1.2);                       // cm; inert while tag off
    CHECK_KNOB_NUM(cfg, "shower_hadronic_growth_max", 0.8);                   // ratio; inert while tag off
    CHECK_KNOB_NUM(cfg, "shower_hadronic_growth_bragg", 1.2);                 // ratio; inert while tag off
    CHECK_KNOB_NUM(cfg, "shower_hadronic_bragg_ratio", 3.0);                  // ratio; inert while tag off
    CHECK_KNOB_NUM(cfg, "shower_hadronic_stem_ratio", 0.0);                  // MIP units; 0 = branch off
    CHECK_KNOB_BOOL(cfg, "kine_count_orphan_tracks", false);                  // r4 (315167)
    CHECK_KNOB_NUM(cfg, "kine_orphan_track_min", 50.0);                       // cm; read only when on
    // doc pr/117 round 1 -- EM clustering knob families.
    CHECK_KNOB_BOOL(cfg, "shower_pass4_best_owner", false);                   // pr/117 r1; false = legacy greedy owner
    CHECK_KNOB_BOOL(cfg, "shower_merge_relax", false);                        // pr/117 r1; false = no pass
    CHECK_KNOB_NUM(cfg, "shower_merge_relax_dis", 6.0);                       // cm; inert while off
    CHECK_KNOB_NUM(cfg, "shower_merge_relax_angle", 15.0);                    // deg; inert while off
    CHECK_KNOB_NUM(cfg, "shower_merge_relax_min_len", 5.0);                   // cm; fragment length floor; inert while off
    CHECK_KNOB_BOOL(cfg, "shower_flank_absorb", false);                       // pr/117 r1; false = no pass
    CHECK_KNOB_NUM(cfg, "shower_flank_absorb_max_dis", 6.0);                  // cm; inert while off
    CHECK_KNOB_NUM(cfg, "shower_flank_absorb_max_len", 25.0);                 // cm; inert while off
    // doc pr/118 round 1 -- P2 body-distance + charge-continuity merge.
    CHECK_KNOB_BOOL(cfg, "shower_ex1_conn3_body_dis", false);                 // pr/118 r1; false = start-segment gate
    CHECK_KNOB_BOOL(cfg, "shower_merge_relax_continuity", false);             // pr/118 r1; false = legacy merge_relax only
    CHECK_KNOB_NUM(cfg, "shower_merge_relax_cont_frac", 1.0);                 // T2 charge-presence fraction; inert while off
    CHECK_KNOB_NUM(cfg, "shower_merge_relax_cont_gap", 8.0);                  // cm, T2 stub gap; inert while off
    CHECK_KNOB_NUM(cfg, "shower_merge_relax_cont_qmed", 5000.0);              // T2 line-charge floor; inert while off
    CHECK_KNOB_NUM(cfg, "shower_merge_relax_cont_axis", 7.5);                 // deg, T1+T2 axis cone; inert while off
    CHECK_KNOB_NUM(cfg, "shower_merge_relax_cont_dmax", 120.0);               // cm, T2 junction reach; inert while off
    CHECK_KNOB_NUM(cfg, "shower_merge_relax_cont_t1_gap", 1.0);               // cm, T1 touching gap; inert while off
    CHECK_KNOB_NUM(cfg, "shower_merge_relax_cont_t1_fold", 30.0);             // deg, T1 fold; inert while off
    // doc pr/120 round 1 -- backward-stem + em-straight-track admission guards.
    CHECK_KNOB_BOOL(cfg, "stem_backfill_back_guard", false);                  // pr/120 r1; false = legacy stem_backfill
    CHECK_KNOB_NUM(cfg, "stem_backfill_back_ang", 110.0);                     // deg, backward ceiling; inert while off
    CHECK_KNOB_BOOL(cfg, "shower_ex1_walk_em_track_guard", false);            // pr/120 r1; false = legacy ex1 walk
    CHECK_KNOB_NUM(cfg, "shower_ex1_walk_em_track_len", 20.0);                // cm, e- straight-long floor; inert while off
    // doc pr/121 round 1 -- examine_shower_1 dedup re-homes multi-seg victims.
    CHECK_KNOB_BOOL(cfg, "shower_ex1_dedup_rehome", false);                   // pr/121 r1; false = legacy dedup drop
    // doc pr/123 round 1 -- pass4_angle over-reach (owner line 2026-08-28).
    CHECK_KNOB_BOOL(cfg, "shower_pass4_prune_detached", false);               // pr/123 r1; false = no prune pass
    CHECK_KNOB_NUM(cfg, "shower_pass4_prune_gap", 40.0);                      // cm; inert while prune off
    CHECK_KNOB_NUM(cfg, "shower_pass4_track_guard_len", 0.0);                 // cm; 0 = no length guard
    CHECK_KNOB_NUM(cfg, "shower_pass4_prox_guard_len", 0.0);                  // cm; 0 = pass4_proximity unguarded (doc pr/130 item 1b)
    CHECK_KNOB_NUM(cfg, "shower_pass3_backfill_guard_len", 0.0);              // cm; 0 = pass3 sibling backfill unguarded (doc pr/130 item 1b)
    CHECK_KNOB_NUM(cfg, "stem_backfill_back_dvtx", 0.0);                      // cm; 0 = back guard ignores vertex distance (doc pr/130 item B)
    CHECK_KNOB_BOOL(cfg, "shower_pass4_prefilter_v1_escape", false);          // false = pass-4 angle_v2>30 pre-filter ignores angle_v1 (doc pr/136 r2)
    CHECK_KNOB_NUM(cfg, "shower_pass4_prefilter_v1_max_v2", 0.0);             // deg; 0 = the angle_v1 escape has no angle_v2 ceiling (doc pr/136 r2)
    CHECK_KNOB_NUM(cfg, "shower_pass4_prefilter_v1_max_dis", 0.0);            // cm; 0 = the angle_v1 escape has no proximity bound (doc pr/136 r3)
    // doc pr/132 -- pi0-finder knobs; defaults = the legacy hard-coded constants.
    CHECK_KNOB_NUM(cfg, "pi0_mass_offset", 10.0);                             // MeV; the finders' "+10" offset (doc pr/132 K1)
    CHECK_KNOB_NUM(cfg, "pi0_assoc_angle_deg", 30.0);                         // deg; association cut (doc pr/132 K2)
    CHECK_KNOB_NUM(cfg, "pi0_attached_partner_min_mev", 0.0);                 // MeV; 0 = no nueCC-fake guard (doc pr/132 K3)
    CHECK_KNOB_BOOL(cfg, "pi0_nv_allow_type2", false);                        // doc pr/132 K4; false = conn_type-3-only pool
    CHECK_KNOB_NUM(cfg, "pi0_nv_max_prongs", 2.0);                            // doc pr/132 K5; without-vertex prong cap
    CHECK_KNOB_BOOL(cfg, "pi0_readmit_retyped", false);                       // doc pr/132 K7; retyped showers stay out of pi0 pairing
    CHECK_KNOB_BOOL(cfg, "pi0_admit_type3", false);                           // doc pr/132 K8; with-vertex pool conn_type 2 only
    CHECK_KNOB_NUM(cfg, "pi0_crumb_assoc_mev", 0.0);                          // doc pr/132 K9; 0 = crumbs keep the angle test
    CHECK_KNOB_NUM(cfg, "pi0_collinear_merge_deg", 0.0);                      // doc pr/132 K12; 0 = no virtual collinear merge
    CHECK_KNOB_NUM(cfg, "pi0_nv_partner_min_mev", 0.0);                       // doc pr/132 K13; 0 = no path-2 partner floor
    CHECK_KNOB_BOOL(cfg, "pi0_nv_retry_paired", false);                       // doc pr/132 K14; legacy path-2 early return
    CHECK_KNOB_BOOL(cfg, "pi0_reseat_start_assoc", false);                    // doc pr/132 K15; fit-cloud conn-2 starts
    CHECK_KNOB_NUM(cfg, "shower_em_collinear_deg", 0.0);                      // doc pr/132 K16; 0 = no build-time EM collinear merge
    CHECK_KNOB_NUM(cfg, "shower_em_collinear_dis_cm", 60.0);                  // doc pr/132 K16; inert while deg = 0
    CHECK_KNOB_NUM(cfg, "shower_em_collinear_host_mev", 20.0);                // doc pr/132 K16; inert while deg = 0
    CHECK_KNOB_NUM(cfg, "shower_em_backext_perp_cm", 0.0);                    // doc pr/132 K17; 0 = no start back-extension
    CHECK_KNOB_NUM(cfg, "shower_em_backext_len_cm", 40.0);                    // doc pr/132 K17; inert while perp = 0
    CHECK_KNOB_NUM(cfg, "pi0_accept_merge_dis_cm", 0.0);                      // doc pr/132 K18; 0 = no acceptance-aware merge
    CHECK_KNOB_NUM(cfg, "pi0_bp_vertex_miss_cm", 0.0);                        // doc pr/132 K19; 0 = no NC vertex proposer
    CHECK_KNOB_BOOL(cfg, "pi0_admit_muon_showers", false);                    // doc pr/133 K20; mu-typed showers stay out of the pi0 pools
    CHECK_KNOB_NUM(cfg, "pi0_nc_sig_angle_deg", 0.0);                         // doc pr/133 K21; 0 = v3 bp gate
    CHECK_KNOB_NUM(cfg, "pi0_nc_floor_mev", 0.0);                             // doc pr/133 K21 v2; 0 = legacy 20 MeV partner floor
    CHECK_KNOB_NUM(cfg, "pi0_nc_pf_assoc_deg", 0.0);                          // doc pr/133 K21 v2.2; 0 = no post-fire PF update
    CHECK_KNOB_BOOL(cfg, "pi0_nc_frag_merge", false);                         // doc pr/134 K22; false = fragment-level bp pairing
    CHECK_KNOB_NUM(cfg, "pi0_pf_assoc_deg", 0.0);                             // doc pr/134 K23; 0 = no P1 post-accept PF absorb
    CHECK_KNOB_BOOL(cfg, "pi0_prefer_main_vertex", false);                    // doc pr/134 K24; false = legacy association + ranking
    CHECK_KNOB_NUM(cfg, "pi0_nv_max_vtx_shift_cm", 0.0);                      // doc pr/132 K10; 0 = no decay-point shift cap
    CHECK_KNOB_NUM(cfg, "pi0_nv_mass_window_mev", 60.0);                      // doc pr/132 K11; legacy without-vertex half-window
    // doc pr/124 front A -- 25-40 cm gap-band tier-2 prune.
    CHECK_KNOB_NUM(cfg, "shower_pass4_prune_gap2", 0.0);                      // cm; 0 = no tier-2 band prune
    CHECK_KNOB_NUM(cfg, "shower_pass4_prune2_ang", 40.0);                     // deg; inert while gap2 = 0
    CHECK_KNOB_NUM(cfg, "shower_pass4_prune2_mdqdx", 2.5);                    // x MIP; inert while gap2 = 0
    CHECK_KNOB_NUM(cfg, "shower_pass3_cone_guard_len", 0.0);                  // cm; 0 = no pass3 track-pdg decline
    CHECK_KNOB_BOOL(cfg, "shower_samevtx_track_absorb", false);               // doc pr/125; false = no pass
    CHECK_KNOB_NUM(cfg, "shower_samevtx_absorb_gap", 6.0);                    // cm; inert while pass off
    CHECK_KNOB_NUM(cfg, "shower_samevtx_absorb_max_len", 50.0);               // cm; inert while pass off
    CHECK_KNOB_NUM(cfg, "shower_samevtx_absorb_min_len", 5.0);                // cm; inert while pass off
    CHECK_KNOB_BOOL(cfg, "shower_satellite_absorb", false);                   // doc pr/125; false = no pass
    CHECK_KNOB_NUM(cfg, "shower_satellite_absorb_max_mev", 10.0);             // MeV; inert while pass off
    CHECK_KNOB_NUM(cfg, "shower_satellite_absorb_host_mev", 20.0);            // MeV; inert while pass off
    // doc pr/138 Phase B -- the EM shower splitter, the one pass that CUTS.
    // shower_split false is what keeps the whole chain byte-identical; the rest
    // are Phase A measurements against 172 owner hand labels and are inert
    // while the pass is off.  max_parts 2 is load-bearing: it is the ONLY
    // setting whose boundary was measured exact (median charge agreement with
    // the owner 1.000 over 27 fired two-way splits), and 3+ is an experiment
    // that MISSED its pre-registered target (doc sec B3).
    CHECK_KNOB_BOOL(cfg, "shower_split", false);                              // doc pr/138 B2; false = no pass
    CHECK_KNOB_NUM(cfg, "shower_split_max_valley", 0.95);                     // sec A5.4 knee; inert while pass off
    CHECK_KNOB_NUM(cfg, "shower_split_min_frac", 0.03);                       // per-seed charge share floor; inert while pass off
    CHECK_KNOB_NUM(cfg, "shower_split_max_parts", 2.0);                       // 2 = the measured-exact kernel; inert while pass off
    CHECK_KNOB_NUM(cfg, "shower_split_min_charge", 1e6);                      // raw Fit::dQ; inert while pass off
    CHECK_KNOB_NUM(cfg, "shower_split_min_nseg", 3.0);                        // inert while pass off
    CHECK_KNOB_NUM(cfg, "shower_split_bundle_gap", 4.0);                      // cm; inert while pass off
    CHECK_KNOB_NUM(cfg, "shower_split_snap", 0.80);                           // k>=3 only; inert while pass off
    // doc pr/139 phase 1 -- four follow-ups to the shipped splitter.  Every
    // one defaults to the shipped behaviour, so production stays byte-identical.
    CHECK_KNOB_BOOL(cfg, "shower_split_skip_shared", false);                  // doc pr/139 P1.1; false = shipped (shared peels allowed)
    CHECK_KNOB_BOOL(cfg, "shower_split_shed_shared", false);                  // doc pr/139 sec 15; false = shipped (an all-shared component is refused, not shed)
    CHECK_KNOB_NUM(cfg, "shower_split_max_seeds", 4);                         // doc pr/139 sec 17; 4 = shipped angular-maxima cap
    CHECK_KNOB_NUM(cfg, "shower_split_em_type_max_len", 0);                   // doc pr/139 sec 25; cm; 0 = off (a short daughter keeps its default track typing)
    CHECK_KNOB_NUM(cfg, "shower_split_max_impact", 0.0);                      // doc pr/139 P1.2; cm; 0 = no bound = shipped
    CHECK_KNOB_BOOL(cfg, "shower_split_em_start", false);                     // doc pr/139 P1.3; false = nearest-overall seed = shipped
    CHECK_KNOB_BOOL(cfg, "shower_split_rehome", false);                       // doc pr/139 P1.4; false = no re-home pass
    CHECK_KNOB_NUM(cfg, "shower_split_rehome_gap", 4.0);                      // doc pr/139 P1.4; cm; inert while the pass is off
    CHECK_KNOB_BOOL(cfg, "kine_count_guard_freed", false);                    // pr/123 r2; kine twin of pf_orphan_guard_freed
    CHECK_KNOB_NUM(cfg, "kine_guard_freed_impact", 0.0);                     // pr/129; 0 = pointing test off
    CHECK_KNOB_NUM(cfg, "kine_guard_freed_miss_deg", 90.0);                  // pr/129
    CHECK_KNOB_BOOL(cfg, "kine_count_near_cross_cluster", false);             // pr/128; kine twin of pf_orphan_near_cross_cluster
    CHECK_KNOB_BOOL(cfg, "kine_count_conn4_near", false);                     // pr/128; kine twin of pf_conn4_near_candidate
    CHECK_KNOB_BOOL(cfg, "straight_cont_cross_cluster", false);               // r4 (137238)
    CHECK_KNOB_BOOL(cfg, "sccc_bridge_body", false);                          // r4 second rung
    CHECK_KNOB_NUM(cfg, "sccc_max_gap", 5.0);                                 // cm; base tier
    CHECK_KNOB_NUM(cfg, "sccc_kink_max", 15.0);                               // deg; base tier
    CHECK_KNOB_NUM(cfg, "sccc_gap_aligned", 12.0);                            // cm; aligned tier
    CHECK_KNOB_NUM(cfg, "sccc_kink_tight", 7.5);                              // deg; aligned tier
    // doc pr/43 round 2 -- three PID-consistency knobs (K1/K2/K3).
    CHECK_KNOB_BOOL(cfg, "single_muon_proton_chain_veto", false);
    CHECK_KNOB_BOOL(cfg, "single_muon_long_muon_claim", false);
    CHECK_KNOB_BOOL(cfg, "pid_flag_reconcile", false);
    // doc pr/45: find_other_segments empty-2D-tree (-1) sentinel guard.
    CHECK_KNOB_BOOL(cfg, "other_seg_empty_2d_guard", false);
    // doc pr/46: long-muon stub bridge in find_cont_muon_segment.
    CHECK_KNOB_BOOL(cfg, "long_muon_stub_bridge", false);
    // doc 84 round 1: P3 stub length (6.0 = the legacy literal) + P2 angle relax.
    CHECK_KNOB_NUM(cfg, "long_muon_stub_bridge_len", 6.0);                    // cm; 6.0 = byte-identical
    CHECK_KNOB_BOOL(cfg, "long_muon_angle_relax_long", false);                // false = legacy 10 deg walk
    CHECK_KNOB_NUM(cfg, "long_muon_angle_relax_deg", 16.0);                   // inert unless relax on
    // doc 84 round 1: MCS defaults stay legacy (enable false, single-segment
    // pf_muon source, no chain-comparator log line).
    CHECK_KNOB_BOOL(cfg, "mcs_enable", false);
    REQUIRE_MESSAGE(cfg.isMember("mcs_muon_source"), "missing knob: mcs_muon_source");
    CHECK(cfg["mcs_muon_source"].asString() == "pf_muon");
    CHECK_KNOB_BOOL(cfg, "mcs_range_comparator_chain", false);                // P5, log-only
    CHECK_KNOB_BOOL(cfg, "mcs_bridged_members", false);                       // doc 84 round 3
    // doc pr/48: back-to-back track fixes.  All three ship OFF; the teb_*
    // operating point is inert until two_end_break opens the pass.
    CHECK_KNOB_BOOL(cfg, "two_end_break", false);
    CHECK_KNOB_BOOL(cfg, "kink_walk_dqdx_stop", false);
    CHECK_KNOB_BOOL(cfg, "kink_break_protect", false);
    // doc pr/50: 172230-class near-vertex robustness -- both passes OFF.
    CHECK_KNOB_BOOL(cfg, "fit_blob_coverage_defer", false);
    CHECK_KNOB_BOOL(cfg, "vertex_kink_snap", false);
    // doc pr/51: main-vertex graph audit + DL rerank cross-cluster swap
    // guard (506746) -- both OFF.
    CHECK_KNOB_BOOL(cfg, "main_vertex_graph_audit", false);
    // doc pr/106 sec 10: exclusion-free charge cloud for the DL vertex net -- OFF.
    CHECK_KNOB_BOOL(cfg, "dl_vtx_cloud_no_exclusion", false);
    // doc pr/107: dQ/dx fit keeps every trajectory point (prototype parity) -- OFF.
    CHECK_KNOB_BOOL(cfg, "dqdx_fit_keep_all_points", false);
    // docs/73 sec 12 round 3: empty-2D-index sentinel guard in
    // eliminate_short_vertex_activities case 5 -- OFF.
    CHECK_KNOB_BOOL(cfg, "esva_ignore_empty_2d", false);
    // doc pr/51 round 3: op3 satellite-anchor extension + traditional-path
    // swap-apply -- both OFF.
    CHECK_KNOB_BOOL(cfg, "main_vertex_swap_apply", false);
    // doc pr/51 round 4: diagnostic-only rough-path probe -- OFF.
    CHECK_KNOB_BOOL(cfg, "rough_path_probe", false);
    CHECK_KNOB_BOOL(cfg, "sgp_edge_probe", false);   // doc pr/73: per-edge sentinel, log-only
    CHECK_KNOB_BOOL(cfg, "vertex_scoreboard", false); // doc pr/75: vertex scoreboard, recording-only
    CHECK_KNOB_BOOL(cfg, "dl_vtx_harvest", false);    // doc pr/79 sec 10: live-feature harvest, recording-only
    // doc pr/51 round 5: steiner gap penalty -- scale 0 = the flavor is
    // never built, do_rough_path stays on the unpenalized "steiner_graph";
    // the four sub-knobs are inert while the scale is 0.
    CHECK_KNOB_NUM(cfg, "steiner_gap_penalty", 0.0);
    CHECK_KNOB_NUM(cfg, "sgp_dead_alpha", 0.25);
    CHECK_KNOB_NUM(cfg, "sgp_min_edge", 0.5);
    CHECK_KNOB_NUM(cfg, "sgp_sample_step", 0.3);
    CHECK_KNOB_NUM(cfg, "sgp_point_radius", 0.2);
    // doc pr/51 round 6: weak-charge deficit term, OFF by default.
    CHECK_KNOB_NUM(cfg, "sgp_weak_scale", 0.0);
    CHECK_KNOB_NUM(cfg, "sgp_weak_qref", 2000.0);
    // doc pr/73 round 2 F3a: do_rough_path route excursion cap, cm.  NEGATIVE
    // is the off value, deliberately -- 0 is a meaningful cap (reject any
    // excursion), so unlike the sgp scale knobs above this one cannot use 0 to
    // mean "off" and its off-test is `< 0`, not `<= 0`.
    CHECK_KNOB_NUM(cfg, "sgp_max_sep", -1.0);
    // doc pr/83: oriented break_segment splits -- OFF (legacy slice by boost
    // source/target; on a reversed edge each child gets the wrong half).
    CHECK_KNOB_BOOL(cfg, "break_seg_orient", false);
    // doc pr/51 round 7: robust vertex fit -- master OFF (AddSegment
    // epilogue never runs); the ten satellites are inert while it is off.
    CHECK_KNOB_BOOL(cfg, "mvfit_robust", false);
    CHECK_KNOB_BOOL(cfg, "mvfit_main_only", true);
    CHECK_KNOB_NUM(cfg, "mvfit_min_len", 10.0);
    CHECK_KNOB_NUM(cfg, "mvfit_rin_margin", 2.0);
    CHECK_KNOB_NUM(cfg, "mvfit_rout_frac", 0.5);
    CHECK_KNOB_NUM(cfg, "mvfit_rout_min", 9.0);
    CHECK_KNOB_NUM(cfg, "mvfit_rout_max", 18.0);
    CHECK_KNOB_NUM(cfg, "mvfit_angle", 20.0);
    CHECK_KNOB_NUM(cfg, "mvfit_min_pts", 5);
    CHECK_KNOB_NUM(cfg, "mvfit_min_aniso", 3.0);
    CHECK_KNOB_NUM(cfg, "mvfit_prior_range", 1.0);
    // doc pr/64 round 7: reassign, instead of drop, an association point
    // that loses (or never enters) clustering_points_segments' Stage-C ghost
    // removal to a same-cluster segment that actually wins the global 2D
    // projection contest -- OFF.
    CHECK_KNOB_BOOL(cfg, "assoc_reassign_orphans", false);
    // doc pr/64 round 8: clear a merge survivor's associate_points when
    // examine_structure_final_1/_1p/_3 deletes a segment that had non-empty
    // associate_points -- OFF.
    CHECK_KNOB_BOOL(cfg, "assoc_clear_on_merge", false);
    // doc pr/72 round 2: guard examine_structure_3 against merging a
    // genuine near-vertex track stub into an unrelated shower/track trunk
    // -- OFF; es3sg_require_terminal's own default is true (it is only
    // read/inert while es3_stub_guard is false).
    CHECK_KNOB_BOOL(cfg, "es3_stub_guard", false);
    CHECK_KNOB_BOOL(cfg, "es3sg_require_terminal", true);

    // Numeric knobs whose legacy value is the INERT one: 0 disables the guard,
    // so an absent key leaves the code path untouched.
    CHECK_KNOB_NUM(cfg, "fit_vertex_min_seg_length", 0.0);   // 0 = all segments enter the fit
    CHECK_KNOB_NUM(cfg, "cathode_x", 0.0);
    CHECK_KNOB_NUM(cfg, "cathode_kink_xcut", 0.0);           // 0 = kink search sees every fit point
    // doc pr/47: wide-baseline cathode kink accept.  angle 0 = path never
    // evaluated; skirt/baseline are inert until the angle knob opens it.
    CHECK_KNOB_NUM(cfg, "cathode_wide_kink_angle", 0.0);
    CHECK_KNOB_NUM(cfg, "cathode_wide_kink_skirt", 3.0);
    CHECK_KNOB_NUM(cfg, "cathode_wide_kink_baseline", 15.0);
    // doc pr/48: the teb_* operating point (cm/deg/dimensionless), all inert
    // while two_end_break is false.
    CHECK_KNOB_NUM(cfg, "teb_min_len", 10.0);
    CHECK_KNOB_NUM(cfg, "teb_min_arm", 1.8);
    CHECK_KNOB_NUM(cfg, "teb_min_arm_pts", 4.0);
    CHECK_KNOB_NUM(cfg, "teb_stub_max", 4.0);
    CHECK_KNOB_NUM(cfg, "teb_accept_range", 15.0);
    CHECK_KNOB_NUM(cfg, "teb_rise_r1", 1.3);
    CHECK_KNOB_NUM(cfg, "teb_rise_r2", 1.15);
    CHECK_KNOB_NUM(cfg, "teb_abs_end_min", 1.7);
    CHECK_KNOB_NUM(cfg, "teb_dip_floor", 0.6);
    CHECK_KNOB_NUM(cfg, "kink_dqdx_hot_ratio", 1.7);
    CHECK_KNOB_NUM(cfg, "teb_score_cap_r1", 0.6);
    CHECK_KNOB_NUM(cfg, "teb_score_cap_r2", 0.9);
    CHECK_KNOB_NUM(cfg, "teb_turn_angle", 25.0);
    CHECK_KNOB_NUM(cfg, "teb_turn_baseline", 35.0);
    CHECK_KNOB_NUM(cfg, "teb_turn_skirt", 3.0);
    // doc pr/90 round 2: R2 argmax arm-fill guard + second-prong gate cap;
    // both 0 = legacy even with two_end_break on.
    CHECK_KNOB_NUM(cfg, "teb_turn_min_arm_frac", 0.0);
    // doc pr/90 round 4: chain-topology gate admission (D1) + route R3
    // turn/activity thresholds (D3) + the R2 bragg veto (D4); all OFF =>
    // legacy even with two_end_break on.
    CHECK_KNOB_BOOL(cfg, "teb_chain_topology", false);
    CHECK_KNOB_NUM(cfg, "teb_r3_turn", 0.0);
    CHECK_KNOB_NUM(cfg, "teb_r3_hot", 0.0);
    CHECK_KNOB_NUM(cfg, "teb_bragg_veto_turn", 0.0);
    // doc pr/50: the vks_* operating point (cm/deg/dimensionless), all inert
    // while vertex_kink_snap is false.  vks_hot_ratio 0 = the optional
    // Bragg-hot veto is OFF (it misfires on the failure class itself:
    // misassigned charge near a misplaced vertex reads hot).
    CHECK_KNOB_NUM(cfg, "vks_radius", 5.0);
    CHECK_KNOB_NUM(cfg, "vks_min_dis", 0.5);
    CHECK_KNOB_NUM(cfg, "vks_angle", 25.0);
    CHECK_KNOB_NUM(cfg, "vks_margin", 10.0);
    CHECK_KNOB_NUM(cfg, "vks_collinear", 30.0);
    CHECK_KNOB_NUM(cfg, "vks_skirt", 0.3);
    CHECK_KNOB_NUM(cfg, "vks_baseline", 2.0);
    CHECK_KNOB_NUM(cfg, "vks_min_arm", 1.5);
    CHECK_KNOB_NUM(cfg, "vks_fit_miss", 0.35);
    CHECK_KNOB_NUM(cfg, "vks_hot_ratio", 0.0);
    CHECK_KNOB_NUM(cfg, "vks_carry_prong", 0.0);  // doc pr/85: 0 = no carry, byte-identical
    // doc pr/104: main-vertex junction snap, inert while vertex_junction_snap is false.
    CHECK_KNOB_BOOL(cfg, "vertex_junction_snap", false);
    CHECK_KNOB_NUM(cfg, "vjs_radius", 5.0);
    CHECK_KNOB_NUM(cfg, "vjs_min_arm", 3.0);
    CHECK_KNOB_NUM(cfg, "vjs_min_prongs", 2);
    CHECK_KNOB_NUM(cfg, "vjs_collinear", 150.0);
    CHECK_KNOB_NUM(cfg, "vjs_fit_margin", 0.5);
    CHECK_KNOB_NUM(cfg, "vjs_fit_rms", 1.0);
    CHECK_KNOB_BOOL(cfg, "vjs_override_kink_snap", false);
    CHECK_KNOB_NUM(cfg, "vjs_min_move", 1.0);
    // doc pr/51: the mvga_* operating point (cm/deg/dimensionless), all
    // inert while main_vertex_graph_audit is false.
    CHECK_KNOB_NUM(cfg, "mvga_radius", 15.0);
    CHECK_KNOB_NUM(cfg, "mvga_dup_tol", 1.4);
    CHECK_KNOB_NUM(cfg, "mvga_dup_frac", 0.7);
    CHECK_KNOB_NUM(cfg, "mvga_dup_angle", 20.0);
    CHECK_KNOB_NUM(cfg, "mvga_bridge_mip", 0.5);
    CHECK_KNOB_NUM(cfg, "mvga_reconnect", 5.0);
    CHECK_KNOB_NUM(cfg, "mvga_stub", 2.0);
    CHECK_KNOB_NUM(cfg, "mvga_stub_pts", 4.0);
    CHECK_KNOB_NUM(cfg, "mvga_reseat_angle", 150.0);
    CHECK_KNOB_NUM(cfg, "mvga_satellite", 0.0);  // round 3: 0 = main-vertex-only op3 scope
    CHECK_KNOB_BOOL(cfg, "mvga_interposed", false);  // doc pr/85: terminal-only op3, byte-identical
    CHECK_KNOB_NUM(cfg, "mvga_interposed_angle", 150.0);
    CHECK_KNOB_NUM(cfg, "mvga_interposed_len", 0.0);  // doc pr/86: 0 = use mvga_stub, byte-identical
    CHECK_KNOB_NUM(cfg, "mvga_sat_dup_frac", 0.0);  // doc pr/86 P4: 0 = use mvga_dup_frac, byte-identical
    CHECK_KNOB_BOOL(cfg, "mvga_interposed_deg1", false);  // doc pr/86 P1b: deg-1 anchors out of reach, byte-identical
    CHECK_KNOB_NUM(cfg, "mvga_splice_straighten", 0.0);  // doc pr/86 round 2 R1: 0 = concatenation verbatim, byte-identical
    CHECK_KNOB_NUM(cfg, "mvga_approach_collapse", 0.0);  // doc pr/86 round 2 R2: 0 = op3.5 skipped, byte-identical
    CHECK_KNOB_NUM(cfg, "mvga_straighten_radius", 0.0);  // doc pr/86 round 2: 0 = prototype 0.2 cm veto radius
    CHECK_KNOB_NUM(cfg, "mvga_op1_radius", 0.0);   // doc pr/83 r3: 0 = use mvga_radius (-1 = unscoped), byte-identical
    CHECK_KNOB_NUM(cfg, "mvga_op1_dup_frac", 0.0); // doc pr/83 r3: 0 = use mvga_dup_frac, byte-identical
    CHECK_KNOB_BOOL(cfg, "mvga_op1_post", false);  // doc pr/83 r3 class A: post-op3 dup pass skipped, byte-identical
    CHECK_KNOB_BOOL(cfg, "swap_orphan_dup_audit", false);  // doc pr/83 r3 Mechanism C: abandoned cluster unaudited, byte-identical
    CHECK_KNOB_NUM(cfg, "mvga_proj_dup_frac", 0.0);  // doc pr/83 r4: projective dup collapse disabled, byte-identical
    CHECK_KNOB_NUM(cfg, "mvga_proj_dqdx_ratio", 0.4);  // doc pr/83 r4: stem dQ/dx gate default; inert while frac == 0
    CHECK_KNOB_NUM(cfg, "mvga_proj_angle", 0.0);  // doc pr/83 r4b: 0 = use mvga_dup_angle, byte-identical
    CHECK_KNOB_NUM(cfg, "mvga_ac_veto_radius", 0.0);  // doc pr/99 r2: 0 = legacy straighten_radius rule, byte-identical
    CHECK_KNOB_NUM(cfg, "mvga_ac_chord_max", 0.0);    // doc pr/99 r2: 0 = no chord cap, byte-identical
    CHECK_KNOB_BOOL(cfg, "mvga_ac_no_cascade", false); // doc pr/99 r2: created products stay collapsible, byte-identical
    CHECK_KNOB_NUM(cfg, "mvga_passthru", 0.0);        // doc pr/103: 0 = op0 pass-through split off, byte-identical
    CHECK_KNOB_NUM(cfg, "mvga_passthru_tol", 1.0);    // doc pr/103: miss tolerance (cm); inert while mvga_passthru == 0
    CHECK_KNOB_BOOL(cfg, "mvga_interposed_fallback", false); // doc pr/103: far-angle decline stands, byte-identical
    CHECK_KNOB_NUM(cfg, "mvga_interposed_fallback_min_angle", 0.0); // doc pr/103: no floor; inert while the fallback is off
    CHECK_KNOB_NUM(cfg, "mvga_dup_starved_asym", 0.0); // doc pr/99 r2: 0 = op1-post angle decline stands, byte-identical
    CHECK_KNOB_NUM(cfg, "mvga_dup_starved_mip", 0.0); // doc pr/99 r2: 0 = op1-post angle decline stands, byte-identical
    CHECK_KNOB_NUM(cfg, "mvga_dup_starved_span", 0.0); // doc pr/99 r2: 0 = no span-comparability test, byte-identical
    CHECK_KNOB_NUM(cfg, "shower_topo_demote_len", 0.0);      // 0 = long segments stay shower-eligible
    CHECK_KNOB_NUM(cfg, "nu_skip_cosmic_bundle_min_length", 0.0);
    CHECK_KNOB_NUM(cfg, "cosmic_companion_min_length", 0.0);
    // low >= high disables the beam gate; both 0 is the disabled state.
    CHECK_KNOB_NUM(cfg, "beam_window_low", 0.0);
    CHECK_KNOB_NUM(cfg, "beam_window_high", 0.0);
    // doc pr/72 round 2: the es3sg_* operating point (cm/deg/dimensionless),
    // fitted from a 117-event merge census, all inert while es3_stub_guard
    // is false.
    CHECK_KNOB_NUM(cfg, "es3sg_stub_max", 7.0);
    CHECK_KNOB_NUM(cfg, "es3sg_len_ratio", 2.0);
    CHECK_KNOB_NUM(cfg, "es3sg_ang3_min", 15.0);
    CHECK_KNOB_NUM(cfg, "es3sg_ang_ratio", 1.0);
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
    CHECK_KNOB_NUM(cfg, "v3_extension_min_gain", -1.0);      // doc pr/24 sec 18 round 5, cm

    // doc pr/67 round 3 (S2), cm.  10.0 is the legacy literal that used to be
    // hard-coded in the first clause of the isochronous-snap gate
    // (NeutrinoOtherSegments.cxx find_other_segments).  If this default ever
    // drifts, the knob-off arm stops being byte-identical -- and because the
    // component member is in cm and is scaled by units::cm at the copy into
    // PatternAlgorithms, a value of 1000 here would mean someone declared the
    // component member in internal units and double-scaled it.
    CHECK_KNOB_NUM(cfg, "iso_snap_min_dir_mag", 10.0);

    // doc pr/102 P1+P2 -- all 0 = off; any drift here silently changes the
    // find_other_segments admission/tagging behaviour on every detector.
    // len_admit and uncover_3d are in cm at the component and scaled by
    // units::cm at the copy into PatternAlgorithms.
    CHECK_KNOB_NUM(cfg, "other_seg_keep_isolated_min_nnf", 0.0);
    CHECK_KNOB_NUM(cfg, "other_seg_keep_isolated_len_admit", 0.0);

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
    // doc pr/53 round 6: the C++ default stays "relaxed"; the SBND operating
    // point selects "relaxed_strict" in cfg only (doc 68 single-source rule).
    REQUIRE(cfg.isMember("graph_name"));
    CHECK(cfg["graph_name"].asString() == "relaxed");
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

TEST_CASE("clus knob defaults: ClusteringDeghost")
{
    auto cfg = defaults_of("ClusteringDeghost");

    // doc 79 round 2: the skeleton shortest-path graph flavor.  Default
    // "ctpc" == byte-identical legacy behavior on every detector;
    // "ctpc_fast" (SBND dg_fast cfg knob) arms the busy-cluster lazy walk.
    REQUIRE(cfg.isMember("graph_name"));
    CHECK(cfg["graph_name"].asString() == "ctpc");

    // The four pre-existing keys, pinned in the same pass (this component
    // had no default_configuration() until doc 79 added one, so none of
    // them round-tripped before).
    CHECK_KNOB_BOOL(cfg, "use_ctpc", true);
    CHECK_KNOB_NUM(cfg, "length_cut", 0);
    CHECK_KNOB_BOOL(cfg, "allow_mixed_faces", false);
    CHECK_KNOB_BOOL(cfg, "empty_view_unique", false);
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
    CHECK(o.track_pid_persist_dqdx == false);  // doc pr/40 F1
    CHECK(o.track_pid_persist_4mom == false);  // doc pr/40 round 2 F4
}

// ---------------------------------------------------------------------------
// doc sbnd_xin/docs/pr/40 -- track (proton/pion/muon) mis-identified as
// electron.  segment_dqdx_spares_electron_reclass is the shared evidence
// check F2/F3 both use; pin its behavior directly since it has no factory
// config surface of its own.
// ---------------------------------------------------------------------------

TEST_CASE("clus knob defaults: segment_dqdx_spares_electron_reclass")
{
    using namespace WireCell::Clus::PR;
    const double MIP = 43000 / units::cm;

    auto make_seg_with_dqdx = [&](double dqdx_per_cm) {
        auto seg = std::make_shared<Segment>();
        std::vector<Fit> fits;
        for (int i = 0; i < 5; ++i) {
            Fit f;
            f.point = WireCell::Point(i * units::cm, 0, 0);
            f.dx = 1 * units::cm;
            f.dQ = dqdx_per_cm * units::cm;
            f.index = i;
            fits.push_back(f);
        }
        seg->fits(fits);
        return seg;
    };

    // No valid fits at all: "no evidence", never spares -- NOT the same as
    // MIP-like.  A segment with zero charge information must not be treated
    // as if it looked like a clean muon.
    {
        auto seg = std::make_shared<Segment>();
        CHECK(segment_dqdx_spares_electron_reclass(seg, MIP) == false);
    }

    // Proton-like: comfortably above the 1.75x MIP threshold.
    CHECK(segment_dqdx_spares_electron_reclass(make_seg_with_dqdx(3.0 * MIP), MIP) == true);
    // Muon/MIP-like: comfortably below the 1.2x MIP threshold.
    CHECK(segment_dqdx_spares_electron_reclass(make_seg_with_dqdx(0.9 * MIP), MIP) == true);
    // Ambiguous middle: between 1.2x and 1.75x, spares nothing -- an electron
    // read is still plausible here.
    CHECK(segment_dqdx_spares_electron_reclass(make_seg_with_dqdx(1.4 * MIP), MIP) == false);
    // Degenerate scale: never spares (guards the division).
    CHECK(segment_dqdx_spares_electron_reclass(make_seg_with_dqdx(3.0 * MIP), 0.0) == false);
}

// ---------------------------------------------------------------------------
// doc sbnd_xin/docs/pr/93 Cause B -- segment_confident_nonelectron_pid is the
// spare-test for the two forced set_pdg(11) shower-acceptance sites.  Pin the
// predicate: real particle_info, pdg not in {0, +-11}, AND a real (<1.0)
// template-PID score.  The 100 "unscored" sentinel (median-fallback pdg-13
// stamps) must NOT spare.
// ---------------------------------------------------------------------------

TEST_CASE("clus knob defaults: segment_confident_nonelectron_pid")
{
    using namespace WireCell::Clus::PR;

    auto make_seg_with_pid = [&](int pdg, double score) {
        auto seg = std::make_shared<Segment>();
        auto pinfo = std::make_shared<Aux::ParticleInfo>(
            pdg, 100.0 * units::MeV, "test",
            WireCell::D4Vector<double>(100.0 * units::MeV, 0, 0, 0));
        seg->particle_info(pinfo);
        seg->particle_score(score);
        return seg;
    };

    // Null / no-particle-info: never spares.
    CHECK(segment_confident_nonelectron_pid(nullptr) == false);
    CHECK(segment_confident_nonelectron_pid(std::make_shared<Segment>()) == false);

    // The two motivating tapes: confident proton (348471, score 0.23) and
    // confident pion (69314) both spare.
    CHECK(segment_confident_nonelectron_pid(make_seg_with_pid(2212, 0.23)) == true);
    CHECK(segment_confident_nonelectron_pid(make_seg_with_pid(211, 0.5)) == true);
    CHECK(segment_confident_nonelectron_pid(make_seg_with_pid(-13, 0.9)) == true);

    // Median-fallback muon: pdg 13 but score 100 sentinel -- NOT confident,
    // does not spare (this population is exactly what the vote/acceptance
    // sites legitimately relabel).
    CHECK(segment_confident_nonelectron_pid(make_seg_with_pid(13, 100.0)) == false);

    // Already-electron or unset pdg: never spares regardless of score.
    CHECK(segment_confident_nonelectron_pid(make_seg_with_pid(11, 0.1)) == false);
    CHECK(segment_confident_nonelectron_pid(make_seg_with_pid(-11, 0.1)) == false);
    CHECK(segment_confident_nonelectron_pid(make_seg_with_pid(0, 0.1)) == false);
}

// ---------------------------------------------------------------------------
// doc sbnd_xin/docs/pr/40 round 2 F5 -- an electron cannot father a proton.
// segment_has_proton_daughter is the shared topology check
// shower_proton_daughter_pion uses; pin its behavior directly on a small
// hand-built graph, mirroring doctest_pr_graph_order.cxx's PR::add_segment
// idiom (vertices need no fit point -- find_vertices only reads graph
// structure).
// ---------------------------------------------------------------------------

TEST_CASE("clus knob defaults: segment_has_proton_daughter")
{
    using namespace WireCell::Clus::PR;
    const double MIP = 43000 / units::cm;

    auto make_seg_with_dqdx = [&](double dqdx_per_cm) {
        auto seg = std::make_shared<Segment>();
        std::vector<Fit> fits;
        for (int i = 0; i < 5; ++i) {
            Fit f;
            f.point = WireCell::Point(i * units::cm, 0, 0);
            f.dx = 1 * units::cm;
            f.dQ = dqdx_per_cm * units::cm;
            f.index = i;
            fits.push_back(f);
        }
        seg->fits(fits);
        return seg;
    };
    auto make_proton = [&](double dqdx_per_cm) {
        auto seg = make_seg_with_dqdx(dqdx_per_cm);
        auto pinfo = std::make_shared<Aux::ParticleInfo>(
            2212, 938.272 * units::MeV, "proton",
            WireCell::D4Vector<double>(938.272 * units::MeV, 0, 0, 0));
        seg->particle_info(pinfo);
        return seg;
    };

    // Bare electron-candidate segment, no neighbours at all: never fires.
    {
        Graph g;
        auto main_vertex = std::make_shared<Vertex>();
        auto far_vertex = std::make_shared<Vertex>();
        auto seg = std::make_shared<Segment>();
        add_segment(g, seg, main_vertex, far_vertex);
        CHECK(segment_has_proton_daughter(g, seg, main_vertex, MIP) == false);
    }

    // Emanates from the neutrino vertex, far end abuts a PID'd AND
    // charge-confirmed proton (>1.75x MIP): fires.
    {
        Graph g;
        auto main_vertex = std::make_shared<Vertex>();
        auto far_vertex = std::make_shared<Vertex>();
        auto far2 = std::make_shared<Vertex>();
        auto seg = std::make_shared<Segment>();
        add_segment(g, seg, main_vertex, far_vertex);
        auto proton = make_proton(3.0 * MIP);
        add_segment(g, proton, far_vertex, far2);
        CHECK(segment_has_proton_daughter(g, seg, main_vertex, MIP) == true);
    }

    // Far-end neighbour IS pdg 2212 but its own charge is only ambiguous
    // (between the 1.2x/1.75x MIP band): not independently confirmed, does
    // not fire -- a raw PID label alone is not enough.
    {
        Graph g;
        auto main_vertex = std::make_shared<Vertex>();
        auto far_vertex = std::make_shared<Vertex>();
        auto far2 = std::make_shared<Vertex>();
        auto seg = std::make_shared<Segment>();
        add_segment(g, seg, main_vertex, far_vertex);
        auto proton = make_proton(1.4 * MIP);
        add_segment(g, proton, far_vertex, far2);
        CHECK(segment_has_proton_daughter(g, seg, main_vertex, MIP) == false);
    }

    // Proton hangs off the main_vertex end, not the far end -- a SIBLING,
    // not a daughter (the ordinary, correct nueCC topology).  Never fires.
    {
        Graph g;
        auto main_vertex = std::make_shared<Vertex>();
        auto far_vertex = std::make_shared<Vertex>();
        auto far3 = std::make_shared<Vertex>();
        auto seg = std::make_shared<Segment>();
        add_segment(g, seg, main_vertex, far_vertex);
        auto proton = make_proton(3.0 * MIP);
        add_segment(g, proton, main_vertex, far3);
        CHECK(segment_has_proton_daughter(g, seg, main_vertex, MIP) == false);
    }

    // No main_vertex yet (e.g. stage 3, before determine_main_vertex):
    // never fires, regardless of what the far end looks like.
    {
        Graph g;
        auto main_vertex = std::make_shared<Vertex>();
        auto far_vertex = std::make_shared<Vertex>();
        auto far2 = std::make_shared<Vertex>();
        auto seg = std::make_shared<Segment>();
        add_segment(g, seg, main_vertex, far_vertex);
        auto proton = make_proton(3.0 * MIP);
        add_segment(g, proton, far_vertex, far2);
        CHECK(segment_has_proton_daughter(g, seg, nullptr, MIP) == false);
    }

    // Degenerate scale: never fires (guards the division, same contract as
    // segment_dqdx_spares_electron_reclass).
    {
        Graph g;
        auto main_vertex = std::make_shared<Vertex>();
        auto far_vertex = std::make_shared<Vertex>();
        auto far2 = std::make_shared<Vertex>();
        auto seg = std::make_shared<Segment>();
        add_segment(g, seg, main_vertex, far_vertex);
        auto proton = make_proton(3.0 * MIP);
        add_segment(g, proton, far_vertex, far2);
        CHECK(segment_has_proton_daughter(g, seg, main_vertex, 0.0) == false);
    }
}

// ---------------------------------------------------------------------------
// doc sbnd_xin/docs/pr/40 round 4 F8 -- a muon cannot terminate in a
// multi-proton hadronic vertex.  segment_at_multi_proton_vertex is F8's
// shared topology check, generalizing segment_has_proton_daughter above.
// ---------------------------------------------------------------------------

TEST_CASE("clus knob defaults: segment_at_multi_proton_vertex")
{
    using namespace WireCell::Clus::PR;
    const double MIP = 43000 / units::cm;

    auto make_seg_with_dqdx = [&](double dqdx_per_cm) {
        auto seg = std::make_shared<Segment>();
        std::vector<Fit> fits;
        for (int i = 0; i < 5; ++i) {
            Fit f;
            f.point = WireCell::Point(i * units::cm, 0, 0);
            f.dx = 1 * units::cm;
            f.dQ = dqdx_per_cm * units::cm;
            f.index = i;
            fits.push_back(f);
        }
        seg->fits(fits);
        return seg;
    };
    auto make_proton = [&](double dqdx_per_cm) {
        auto seg = make_seg_with_dqdx(dqdx_per_cm);
        auto pinfo = std::make_shared<Aux::ParticleInfo>(
            2212, 938.272 * units::MeV, "proton",
            WireCell::D4Vector<double>(938.272 * units::MeV, 0, 0, 0));
        seg->particle_info(pinfo);
        return seg;
    };

    // Muon segment's far end (not main_vertex) has TWO charge-confirmed
    // protons attached: fires.
    {
        Graph g;
        auto main_vertex = std::make_shared<Vertex>();
        auto far_vertex  = std::make_shared<Vertex>();
        auto p1_end = std::make_shared<Vertex>();
        auto p2_end = std::make_shared<Vertex>();
        auto muon = std::make_shared<Segment>();
        add_segment(g, muon, main_vertex, far_vertex);
        add_segment(g, make_proton(3.0 * MIP), far_vertex, p1_end);
        add_segment(g, make_proton(3.5 * MIP), far_vertex, p2_end);
        CHECK(segment_at_multi_proton_vertex(g, muon, main_vertex, MIP) == true);
    }

    // Only ONE charge-confirmed proton at the far end: does not fire --
    // the owner's "two protons" wording is read literally (min_protons=2).
    {
        Graph g;
        auto main_vertex = std::make_shared<Vertex>();
        auto far_vertex  = std::make_shared<Vertex>();
        auto p1_end = std::make_shared<Vertex>();
        auto muon = std::make_shared<Segment>();
        add_segment(g, muon, main_vertex, far_vertex);
        add_segment(g, make_proton(3.0 * MIP), far_vertex, p1_end);
        CHECK(segment_at_multi_proton_vertex(g, muon, main_vertex, MIP) == false);
    }

    // The two-proton vertex IS the neutrino vertex: does not fire -- an
    // ordinary numuCC muon-plus-two-protons topology at the neutrino vertex
    // is correct and must be left alone.
    {
        Graph g;
        auto main_vertex = std::make_shared<Vertex>();
        auto far_vertex  = std::make_shared<Vertex>();
        auto p1_end = std::make_shared<Vertex>();
        auto p2_end = std::make_shared<Vertex>();
        auto muon = std::make_shared<Segment>();
        add_segment(g, muon, main_vertex, far_vertex);
        add_segment(g, make_proton(3.0 * MIP), main_vertex, p1_end);
        add_segment(g, make_proton(3.5 * MIP), main_vertex, p2_end);
        CHECK(segment_at_multi_proton_vertex(g, muon, main_vertex, MIP) == false);
    }

    // Two pdg-2212 neighbours but both only ambiguous (below 1.75x MIP):
    // not independently confirmed, does not fire.
    {
        Graph g;
        auto main_vertex = std::make_shared<Vertex>();
        auto far_vertex  = std::make_shared<Vertex>();
        auto p1_end = std::make_shared<Vertex>();
        auto p2_end = std::make_shared<Vertex>();
        auto muon = std::make_shared<Segment>();
        add_segment(g, muon, main_vertex, far_vertex);
        add_segment(g, make_proton(1.4 * MIP), far_vertex, p1_end);
        add_segment(g, make_proton(1.3 * MIP), far_vertex, p2_end);
        CHECK(segment_at_multi_proton_vertex(g, muon, main_vertex, MIP) == false);
    }

    // No main_vertex yet: never fires, regardless of the far-end topology.
    {
        Graph g;
        auto main_vertex = std::make_shared<Vertex>();
        auto far_vertex  = std::make_shared<Vertex>();
        auto p1_end = std::make_shared<Vertex>();
        auto p2_end = std::make_shared<Vertex>();
        auto muon = std::make_shared<Segment>();
        add_segment(g, muon, main_vertex, far_vertex);
        add_segment(g, make_proton(3.0 * MIP), far_vertex, p1_end);
        add_segment(g, make_proton(3.5 * MIP), far_vertex, p2_end);
        CHECK(segment_at_multi_proton_vertex(g, muon, nullptr, MIP) == false);
    }
}

// ---------------------------------------------------------------------------
// doc sbnd_xin/docs/pr/40 round 5 F10/F11 -- shared straightness veto used by
// shower_clustering_connecting_to_main_vertex and segment_is_shower_
// trajectory.  Same shape as the toolkit's own straightness demotion
// (NeutrinoVertexFinder.cxx:1432-1447): length > 10 cm, and either
// direct_length >= 34 cm or direct/arc ratio > 0.93.
// ---------------------------------------------------------------------------

TEST_CASE("clus knob defaults: segment_is_straight_long_track")
{
    using namespace WireCell::Clus::PR;

    auto make_straight_seg = [&](double length_cm) {
        auto seg = std::make_shared<Segment>();
        std::vector<Fit> fits;
        const int n = 20;
        for (int i = 0; i < n; ++i) {
            Fit f;
            f.point = WireCell::Point(i * (length_cm / (n - 1)) * units::cm, 0, 0);
            f.dx = (length_cm / (n - 1)) * units::cm;
            f.dQ = 1.0;
            f.index = i;
            fits.push_back(f);
        }
        seg->fits(fits);
        return seg;
    };
    // A zigzag whose arc length is ~1.43x its straight-line span (each leg is
    // a 45-degree jog forward+sideways of equal size), well under the 0.93
    // ratio cut regardless of overall length.
    auto make_wiggly_seg = [&](double length_cm) {
        auto seg = std::make_shared<Segment>();
        std::vector<Fit> fits;
        const int n = 20;
        const double step = (length_cm / (n - 1)) * units::cm;
        for (int i = 0; i < n; ++i) {
            Fit f;
            f.point = WireCell::Point(i * step * 0.7071, (i % 2 == 0 ? 0.0 : step * 0.7071), 0);
            f.dx = step;
            f.dQ = 1.0;
            f.index = i;
            fits.push_back(f);
        }
        seg->fits(fits);
        return seg;
    };

    // Long (>10cm) and straight: fires.
    CHECK(segment_is_straight_long_track(make_straight_seg(20.0)) == true);

    // Long and wiggly (direct/arc well under 0.93, direct well under 34cm):
    // does not fire.
    CHECK(segment_is_straight_long_track(make_wiggly_seg(20.0)) == false);

    // Short (<=10cm) even though perfectly straight: does not fire -- this is
    // the boundary that leaves SBND evt 55715 seg 15005 (6.1 cm) untouched
    // per the owner's decision not to extend the relabel across it.
    CHECK(segment_is_straight_long_track(make_straight_seg(6.0)) == false);

    // direct_length >= 34cm alone is sufficient regardless of the ratio path
    // (a long straight track well past the ratio-only branch).
    CHECK(segment_is_straight_long_track(make_straight_seg(40.0)) == true);

    // Null segment: never fires.
    CHECK(segment_is_straight_long_track(nullptr) == false);
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
// TrackFitting -- NOT an IConfigurable.  Its Parameters struct is reached by
// TaggerCheckNeutrino/TaggerCheckSTM's own trackfitting_config_file loader
// (a plain JSON -> set_parameter(name, double) walk, doc pr/28 S17), not
// through default_configuration(), so it is pinned by direct construction
// instead of defaults_of().
// ---------------------------------------------------------------------------

TEST_CASE("clus knob defaults: TrackFitting skip_revert_iso_xext_cut is off")
{
    // doc pr/28 S17: abstain guard on the multi-track charge veto's revert
    // (skip_trajectory_point's `p = ps_point`, revived unconditionally by
    // 23bd6783 T1/T2).  C++ default -1 = off, i.e. the revert stays
    // unconditional -- byte-identical to pre-S17 behaviour (gate
    // work-pr28r10-off48b vs work-pr28r10-cleanref48, 48/48 events, 96/96
    // archives).  Round-tripped through set_parameter/get_parameter, the same
    // path the trackfitting_config JSON loader uses.
    Clus::TrackFitting tf;
    CHECK(tf.get_parameter("skip_revert_iso_xext_cut") == doctest::Approx(-1.0));

    tf.set_parameter("skip_revert_iso_xext_cut", 20 * units::cm);
    CHECK(tf.get_parameter("skip_revert_iso_xext_cut") == doctest::Approx(20 * units::cm));

    // TrackFittingPresets::create_with_current_values() (used by the SBND/uBooNE
    // component wiring) sets every Parameters field explicitly, including this
    // one -- pin it separately so a preset edit that drops the line, or sets a
    // non-off value, is caught even if the in-class initializer above stays -1.
    auto preset = Clus::TrackFittingPresets::create_with_current_values();
    CHECK(preset.get_parameters().skip_revert_iso_xext_cut == doctest::Approx(-1.0));
}

TEST_CASE("clus knob defaults: TrackFitting dqdx_fit_keep_all_points is off")
{
    // doc pr/107: the pre-dQ/dx form_map_graph pass drops zero-quantity
    // trajectory points (legacy, toolkit-only); > 0 keeps every point as the
    // prototype does.  Default 0 = legacy = byte-identical; same
    // set_parameter round-trip contract as the other double-sentinel knobs.
    Clus::TrackFitting tf;
    CHECK(tf.get_parameter("dqdx_fit_keep_all_points") == doctest::Approx(0.0));
    tf.set_parameter("dqdx_fit_keep_all_points", 1.0);
    CHECK(tf.get_parameter("dqdx_fit_keep_all_points") == doctest::Approx(1.0));
    auto preset = Clus::TrackFittingPresets::create_with_current_values();
    CHECK(preset.get_parameters().dqdx_fit_keep_all_points == doctest::Approx(0.0));
}

TEST_CASE("clus knob defaults: TrackFitting fit_blob_coverage is off")
{
    // doc pr/49: cross-cluster projection-ghost deweighting of the 2D charge
    // association (18255-57441 V-plane ghost).  C++ default -1 = off --
    // examine_point_association's charge_cut loops take the legacy path and
    // the compiled tree is byte-identical.  >= 0 = on, value = wire/slice
    // tolerance in cells.  Same double-sentinel/set_parameter round-trip
    // contract as skip_revert_iso_xext_cut above.  The two companion
    // numerics ride the C++ defaults and are inert while the main knob is
    // off.  Round 3 (owner decision 2026-08-08): "foreign" is scope-aware
    // (a cluster with a segment in the current fit never counts), and the
    // 3D far-gate default is 0 = disabled -- scope membership replaced the
    // 15 cm distance criterion.  The weight stays 0.1 (deweight, not drop:
    // dead-channel single-view charge must stay usable).
    Clus::TrackFitting tf;
    CHECK(tf.get_parameter("fit_blob_coverage") == doctest::Approx(-1.0));
    CHECK(tf.get_parameter("fit_blob_coverage_ghost_dis") == doctest::Approx(0.0));
    CHECK(tf.get_parameter("fit_blob_coverage_weight") == doctest::Approx(0.1));

    tf.set_parameter("fit_blob_coverage", 0.0);
    CHECK(tf.get_parameter("fit_blob_coverage") == doctest::Approx(0.0));
    tf.set_parameter("fit_blob_coverage_weight", 0.2);
    CHECK(tf.get_parameter("fit_blob_coverage_weight") == doctest::Approx(0.2));

    auto preset = Clus::TrackFittingPresets::create_with_current_values();
    CHECK(preset.get_parameters().fit_blob_coverage == doctest::Approx(-1.0));
    CHECK(preset.get_parameters().fit_blob_coverage_ghost_dis == doctest::Approx(0.0));
    CHECK(preset.get_parameters().fit_blob_coverage_weight == doctest::Approx(0.1));
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

// ---------------------------------------------------------------------------
// MultiAlgBlobClustering BeePFConfig -- the particle-flow (Bee "mc" tree)
// display knobs.  These are in-class initializers on a nested struct rather
// than default_configuration() keys, so pin the struct directly.  Closes the
// gap where none of the pf_* switches were pinned (doc pr/38 Round 4).
// ---------------------------------------------------------------------------

#include "WireCellClus/MultiAlgBlobClustering.h"

TEST_CASE("clus knob defaults: MultiAlgBlobClustering BeePFConfig pf switches are all OFF")
{
    WireCell::Clus::MultiAlgBlobClustering::BeePFConfig pfc;
    CHECK(pfc.prototype_names == false);
    CHECK(pfc.em_ke_min == 0.0);
    CHECK(pfc.np_ke_min == 0.0);
    CHECK(pfc.pf_track_main_cluster_only == false);      // doc pr/34 F1
    CHECK(pfc.pf_shower_vertex_barrier == false);        // doc pr/34 F2 (+ pr/38)
    CHECK(pfc.pf_shower_parent_precedence == false);     // doc pr/34 F3
    CHECK(pfc.pf_pi0_node_per_id == false);              // doc pr/34 F4
    CHECK(pfc.pf_pdg_name_prototype_fallback == false);  // doc pr/34 F5
    CHECK(pfc.pf_orphan_track_parentage == false);       // doc pr/38 Round 4
    CHECK(pfc.pf_orphan_audit_only == false);            // doc pr/65 round 3 (pin was missing; added pr/84 r2)
    CHECK(pfc.pf_direct_when_touching == false);         // doc pr/84 r2 F1
    CHECK(pfc.pf_touch_max == 3.0 * WireCell::units::cm);        // read only when F1 on
    CHECK(pfc.pf_pseudo_gap_from_main == false);         // doc pr/84 r2 F2
    CHECK(pfc.pf_unique_node_ids == false);              // doc pr/84 r3 G1
    CHECK(pfc.pf_drop_stray_satellites == false);        // doc pr/92
    CHECK(pfc.pf_orphan_confident_track == false);       // doc pr/93 r4 (315167)
    CHECK(pfc.pf_orphan_track_min == 50.0 * WireCell::units::cm);  // read only when on
    CHECK(pfc.pf_track_owns_loose_vertex == false);      // doc pr/93 r4 (69314)
}
