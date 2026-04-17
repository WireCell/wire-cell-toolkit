// UbooneNumuBDTScorer.cxx
//
// IEnsembleVisitor that runs TMVA BDT scoring for the numu CC tagger.
// Ports the following functions from prototype_pid/src/NeutrinoID_numu_bdts.h:
//   cal_numu_bdts_xgboost()  — top-level scorer, writes ti.numu_score
//   cal_numu_1_bdt()         — sub-BDT on flag-1 (direct-muon) features
//   cal_numu_2_bdt()         — sub-BDT on flag-2 (long-muon shower) features
//   cal_numu_3_bdt()         — sub-BDT on flag-3 (indirect-muon) features
//   cal_cosmict_10_bdt()     — sub-BDT on upstream-dirt features (used by xgboost)
//
// NOT ported: cal_numu_bdts() (old TMVA variant), cal_cosmict_{2_4,3_5,6,7,8}_bdt()
//   (only needed by the old variant).
//
// Translation conventions vs. prototype:
//   tagger_info.xxx  →  ti.xxx   (Clus::PR::TaggerInfo& passed by reference)
//   kine_info.xxx    →  ki.xxx   (const Clus::PR::KineInfo& passed by reference)
//   match_isFC       →  ti.match_isFC  (placeholder; filled by TaggerCheckNeutrino)
//   "input_data_files/weights/foo.xml"  →  m_*_xml  (configured via wc.resolve)

#include "WireCellRoot/UbooneNumuBDTScorer.h"
#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Persist.h"

#include <cmath>

WIRECELL_FACTORY(UbooneNumuBDTScorer, WireCell::Root::UbooneNumuBDTScorer,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Root;

UbooneNumuBDTScorer::UbooneNumuBDTScorer()
    : log(Log::logger("root.UbooneNumuBDTScorer"))
{
}

void UbooneNumuBDTScorer::configure(const WireCell::Configuration& cfg)
{
    auto resolve = [&](const std::string& p) {
        if (p.empty()) return p;
        std::string r = Persist::resolve(p);
        if (r.empty()) log->error("UbooneNumuBDTScorer: weight file not found: {}", p);
        return r;
    };
    m_grouping_name    = get<std::string>(cfg, "grouping",             "live");
    m_numu1_xml        = resolve(get<std::string>(cfg, "numu1_weights_xml",     ""));
    m_numu2_xml        = resolve(get<std::string>(cfg, "numu2_weights_xml",     ""));
    m_numu3_xml        = resolve(get<std::string>(cfg, "numu3_weights_xml",     ""));
    m_cosmict10_xml    = resolve(get<std::string>(cfg, "cosmict10_weights_xml", ""));
    m_numu_xgboost_xml = resolve(get<std::string>(cfg, "numu_xgboost_xml",      ""));

    init_readers();
}

Configuration UbooneNumuBDTScorer::default_configuration() const
{
    Configuration cfg;
    cfg["grouping"]             = "live";
    cfg["numu1_weights_xml"]    = "";   // e.g. wc.resolve("uboone/weights/numu_tagger1.weights.xml")
    cfg["numu2_weights_xml"]    = "";
    cfg["numu3_weights_xml"]    = "";
    cfg["cosmict10_weights_xml"] = "";  // e.g. wc.resolve("uboone/weights/cos_tagger_10.weights.xml")
    cfg["numu_xgboost_xml"]     = "";   // e.g. wc.resolve("uboone/weights/numu_scalars_scores_0923.xml")
    return cfg;
}

// ===========================================================================
// init_readers — create all TMVA readers and BookMVA once.
// Called from configure() after XML paths are resolved.
// ===========================================================================
void UbooneNumuBDTScorer::init_readers()
{
    // --- cosmict_10 reader (5 variables) ---
    if (!m_cosmict10_xml.empty()) {
        m_reader_cosmict10 = std::make_unique<TMVA::Reader>("!V:Silent");
        m_reader_cosmict10->AddVariable("cosmict_10_vtx_z",        &m_cosmict10_vtx_z);
        m_reader_cosmict10->AddVariable("cosmict_10_flag_shower",  &m_cosmict10_flag_shower);
        m_reader_cosmict10->AddVariable("cosmict_10_flag_dir_weak",&m_cosmict10_flag_dir_weak);
        m_reader_cosmict10->AddVariable("cosmict_10_angle_beam",   &m_cosmict10_angle_beam);
        m_reader_cosmict10->AddVariable("cosmict_10_length",       &m_cosmict10_length);
        m_reader_cosmict10->BookMVA("MyBDT", m_cosmict10_xml);
    }

    // --- numu_1 reader (7 variables) ---
    if (!m_numu1_xml.empty()) {
        m_reader_numu1 = std::make_unique<TMVA::Reader>("!V:Silent");
        m_reader_numu1->AddVariable("numu_cc_1_particle_type",    &m_numu1_particle_type);
        m_reader_numu1->AddVariable("numu_cc_1_length",           &m_numu1_length);
        m_reader_numu1->AddVariable("numu_cc_1_medium_dQ_dx",     &m_numu1_medium_dQ_dx);
        m_reader_numu1->AddVariable("numu_cc_1_dQ_dx_cut",        &m_numu1_dQ_dx_cut);
        m_reader_numu1->AddVariable("numu_cc_1_direct_length",    &m_numu1_direct_length);
        m_reader_numu1->AddVariable("numu_cc_1_n_daughter_tracks",&m_numu1_n_daughter_tracks);
        m_reader_numu1->AddVariable("numu_cc_1_n_daughter_all",   &m_numu1_n_daughter_all);
        m_reader_numu1->BookMVA("MyBDT", m_numu1_xml);
    }

    // --- numu_2 reader (4 variables) ---
    if (!m_numu2_xml.empty()) {
        m_reader_numu2 = std::make_unique<TMVA::Reader>("!V:Silent");
        m_reader_numu2->AddVariable("numu_cc_2_length",           &m_numu2_length);
        m_reader_numu2->AddVariable("numu_cc_2_total_length",     &m_numu2_total_length);
        m_reader_numu2->AddVariable("numu_cc_2_n_daughter_tracks",&m_numu2_n_daughter_tracks);
        m_reader_numu2->AddVariable("numu_cc_2_n_daughter_all",   &m_numu2_n_daughter_all);
        m_reader_numu2->BookMVA("MyBDT", m_numu2_xml);
    }

    // --- numu_3 reader (7 variables) ---
    if (!m_numu3_xml.empty()) {
        m_reader_numu3 = std::make_unique<TMVA::Reader>("!V:Silent");
        m_reader_numu3->AddVariable("numu_cc_3_particle_type",   &m_numu3_particle_type);
        m_reader_numu3->AddVariable("numu_cc_3_max_length",      &m_numu3_max_length);
        m_reader_numu3->AddVariable("numu_cc_3_track_length",   &m_numu3_track_length);
        m_reader_numu3->AddVariable("numu_cc_3_max_length_all",  &m_numu3_max_length_all);
        m_reader_numu3->AddVariable("numu_cc_3_max_muon_length", &m_numu3_max_muon_length);
        m_reader_numu3->AddVariable("numu_cc_3_n_daughter_tracks",&m_numu3_n_daughter_tracks);
        m_reader_numu3->AddVariable("numu_cc_3_n_daughter_all",  &m_numu3_n_daughter_all);
        m_reader_numu3->BookMVA("MyBDT", m_numu3_xml);
    }

    // --- xgboost final reader (~72 variables) ---
    if (!m_numu_xgboost_xml.empty()) {
        // Indices into m_xgb_vars — must match the order below exactly.
        // 0..7:   numu flag-3 features (8)
        // 8..20:  cosmict flag-2 features (13)
        // 21..25: cosmict flag-4 features (5)
        // 26..35: cosmict flag-3 features (10)
        // 36..40: cosmict flag-5 features (5)
        // 41..45: cosmict flag-6 features (5)
        // 46..57: cosmict flag-7 features (12)
        // 58..62: cosmict flag-8 features (5)
        // 63:     cosmict flag-9 (1)
        // 64..68: top-level flags (5)
        // 69:     kine_reco_Enu (1)
        // 70:     match_isFC (1)
        // 71..73: sub-BDT scores (3)
        // Total = 74
        m_xgb_vars.resize(74, 0.0f);

        m_reader_xgboost = std::make_unique<TMVA::Reader>("!V:Silent");
        int idx = 0;

        // --- numu flag-3 features ---
        m_reader_xgboost->AddVariable("numu_cc_flag_3",         &m_xgb_vars[idx++]); // 0
        m_reader_xgboost->AddVariable("numu_cc_3_particle_type",&m_xgb_vars[idx++]); // 1
        m_reader_xgboost->AddVariable("numu_cc_3_max_length",   &m_xgb_vars[idx++]); // 2
        m_reader_xgboost->AddVariable("numu_cc_3_track_length", &m_xgb_vars[idx++]); // 3
        m_reader_xgboost->AddVariable("numu_cc_3_max_length_all",   &m_xgb_vars[idx++]); // 4
        m_reader_xgboost->AddVariable("numu_cc_3_max_muon_length",  &m_xgb_vars[idx++]); // 5
        m_reader_xgboost->AddVariable("numu_cc_3_n_daughter_tracks",&m_xgb_vars[idx++]); // 6
        m_reader_xgboost->AddVariable("numu_cc_3_n_daughter_all",   &m_xgb_vars[idx++]); // 7

        // --- cosmic tagger flag-2 features ---
        m_reader_xgboost->AddVariable("cosmict_flag_2",               &m_xgb_vars[idx++]); // 8
        m_reader_xgboost->AddVariable("cosmict_2_filled",             &m_xgb_vars[idx++]); // 9
        m_reader_xgboost->AddVariable("cosmict_2_particle_type",      &m_xgb_vars[idx++]); // 10
        m_reader_xgboost->AddVariable("cosmict_2_n_muon_tracks",      &m_xgb_vars[idx++]); // 11
        m_reader_xgboost->AddVariable("cosmict_2_total_shower_length",&m_xgb_vars[idx++]); // 12
        m_reader_xgboost->AddVariable("cosmict_2_flag_inside",        &m_xgb_vars[idx++]); // 13
        m_reader_xgboost->AddVariable("cosmict_2_angle_beam",         &m_xgb_vars[idx++]); // 14
        m_reader_xgboost->AddVariable("cosmict_2_flag_dir_weak",      &m_xgb_vars[idx++]); // 15
        m_reader_xgboost->AddVariable("cosmict_2_dQ_dx_end",          &m_xgb_vars[idx++]); // 16
        m_reader_xgboost->AddVariable("cosmict_2_dQ_dx_front",        &m_xgb_vars[idx++]); // 17
        m_reader_xgboost->AddVariable("cosmict_2_theta",              &m_xgb_vars[idx++]); // 18
        m_reader_xgboost->AddVariable("cosmict_2_phi",                &m_xgb_vars[idx++]); // 19
        m_reader_xgboost->AddVariable("cosmict_2_valid_tracks",       &m_xgb_vars[idx++]); // 20

        // --- cosmic tagger flag-4 features ---
        m_reader_xgboost->AddVariable("cosmict_flag_4",            &m_xgb_vars[idx++]); // 21
        m_reader_xgboost->AddVariable("cosmict_4_filled",          &m_xgb_vars[idx++]); // 22
        m_reader_xgboost->AddVariable("cosmict_4_flag_inside",     &m_xgb_vars[idx++]); // 23
        m_reader_xgboost->AddVariable("cosmict_4_angle_beam",      &m_xgb_vars[idx++]); // 24
        m_reader_xgboost->AddVariable("cosmict_4_connected_showers",&m_xgb_vars[idx++]); // 25

        // --- cosmic tagger flag-3 features ---
        m_reader_xgboost->AddVariable("cosmict_flag_3",        &m_xgb_vars[idx++]); // 26
        m_reader_xgboost->AddVariable("cosmict_3_filled",      &m_xgb_vars[idx++]); // 27
        m_reader_xgboost->AddVariable("cosmict_3_flag_inside", &m_xgb_vars[idx++]); // 28
        m_reader_xgboost->AddVariable("cosmict_3_angle_beam",  &m_xgb_vars[idx++]); // 29
        m_reader_xgboost->AddVariable("cosmict_3_flag_dir_weak",&m_xgb_vars[idx++]); // 30
        m_reader_xgboost->AddVariable("cosmict_3_dQ_dx_end",   &m_xgb_vars[idx++]); // 31
        m_reader_xgboost->AddVariable("cosmict_3_dQ_dx_front", &m_xgb_vars[idx++]); // 32
        m_reader_xgboost->AddVariable("cosmict_3_theta",       &m_xgb_vars[idx++]); // 33
        m_reader_xgboost->AddVariable("cosmict_3_phi",         &m_xgb_vars[idx++]); // 34
        m_reader_xgboost->AddVariable("cosmict_3_valid_tracks",&m_xgb_vars[idx++]); // 35

        // --- cosmic tagger flag-5 features ---
        m_reader_xgboost->AddVariable("cosmict_flag_5",            &m_xgb_vars[idx++]); // 36
        m_reader_xgboost->AddVariable("cosmict_5_filled",          &m_xgb_vars[idx++]); // 37
        m_reader_xgboost->AddVariable("cosmict_5_flag_inside",     &m_xgb_vars[idx++]); // 38
        m_reader_xgboost->AddVariable("cosmict_5_angle_beam",      &m_xgb_vars[idx++]); // 39
        m_reader_xgboost->AddVariable("cosmict_5_connected_showers",&m_xgb_vars[idx++]); // 40

        // --- cosmic tagger flag-6 features ---
        m_reader_xgboost->AddVariable("cosmict_flag_6",       &m_xgb_vars[idx++]); // 41
        m_reader_xgboost->AddVariable("cosmict_6_filled",     &m_xgb_vars[idx++]); // 42
        m_reader_xgboost->AddVariable("cosmict_6_flag_dir_weak",&m_xgb_vars[idx++]); // 43
        m_reader_xgboost->AddVariable("cosmict_6_flag_inside",&m_xgb_vars[idx++]); // 44
        m_reader_xgboost->AddVariable("cosmict_6_angle",      &m_xgb_vars[idx++]); // 45

        // --- cosmic tagger flag-7 features ---
        m_reader_xgboost->AddVariable("cosmict_flag_7",              &m_xgb_vars[idx++]); // 46
        m_reader_xgboost->AddVariable("cosmict_7_filled",            &m_xgb_vars[idx++]); // 47
        m_reader_xgboost->AddVariable("cosmict_7_flag_sec",          &m_xgb_vars[idx++]); // 48
        m_reader_xgboost->AddVariable("cosmict_7_n_muon_tracks",     &m_xgb_vars[idx++]); // 49
        m_reader_xgboost->AddVariable("cosmict_7_total_shower_length",&m_xgb_vars[idx++]); // 50
        m_reader_xgboost->AddVariable("cosmict_7_flag_inside",       &m_xgb_vars[idx++]); // 51
        m_reader_xgboost->AddVariable("cosmict_7_angle_beam",        &m_xgb_vars[idx++]); // 52
        m_reader_xgboost->AddVariable("cosmict_7_flag_dir_weak",     &m_xgb_vars[idx++]); // 53
        m_reader_xgboost->AddVariable("cosmict_7_dQ_dx_end",         &m_xgb_vars[idx++]); // 54
        m_reader_xgboost->AddVariable("cosmict_7_dQ_dx_front",       &m_xgb_vars[idx++]); // 55
        m_reader_xgboost->AddVariable("cosmict_7_theta",             &m_xgb_vars[idx++]); // 56
        m_reader_xgboost->AddVariable("cosmict_7_phi",               &m_xgb_vars[idx++]); // 57

        // --- cosmic tagger flag-8 features ---
        m_reader_xgboost->AddVariable("cosmict_flag_8",      &m_xgb_vars[idx++]); // 58
        m_reader_xgboost->AddVariable("cosmict_8_filled",    &m_xgb_vars[idx++]); // 59
        m_reader_xgboost->AddVariable("cosmict_8_flag_out",  &m_xgb_vars[idx++]); // 60
        m_reader_xgboost->AddVariable("cosmict_8_muon_length",&m_xgb_vars[idx++]); // 61
        m_reader_xgboost->AddVariable("cosmict_8_acc_length",&m_xgb_vars[idx++]); // 62

        // --- cosmic tagger flag-9 features ---
        m_reader_xgboost->AddVariable("cosmict_flag_9",  &m_xgb_vars[idx++]); // 63

        // --- top-level cosmic / numu flags ---
        m_reader_xgboost->AddVariable("cosmic_flag",    &m_xgb_vars[idx++]); // 64
        m_reader_xgboost->AddVariable("cosmic_filled",  &m_xgb_vars[idx++]); // 65
        m_reader_xgboost->AddVariable("cosmict_flag",   &m_xgb_vars[idx++]); // 66
        m_reader_xgboost->AddVariable("numu_cc_flag",   &m_xgb_vars[idx++]); // 67
        m_reader_xgboost->AddVariable("cosmict_flag_1", &m_xgb_vars[idx++]); // 68

        // --- kinematics + fiducial flag ---
        m_reader_xgboost->AddVariable("kine_reco_Enu", &m_xgb_vars[idx++]); // 69
        m_reader_xgboost->AddVariable("match_isFC",    &m_xgb_vars[idx++]); // 70

        // --- sub-BDT scores (filled above) ---
        m_reader_xgboost->AddVariable("cosmict_10_score", &m_xgb_vars[idx++]); // 71
        m_reader_xgboost->AddVariable("numu_1_score",     &m_xgb_vars[idx++]); // 72
        m_reader_xgboost->AddVariable("numu_2_score",     &m_xgb_vars[idx++]); // 73

        m_reader_xgboost->BookMVA("MyBDT", m_numu_xgboost_xml);
    }
}

void UbooneNumuBDTScorer::visit(Clus::Facade::Ensemble& ensemble) const
{
    if (!m_reader_xgboost) {
        log->warn("UbooneNumuBDTScorer: numu_xgboost_xml not set or file not found — skipping numu BDT scoring");
        return;
    }

    auto groupings = ensemble.with_name(m_grouping_name);
    if (groupings.empty()) {
        log->debug("UbooneNumuBDTScorer: no grouping '{}'", m_grouping_name);
        return;
    }

    auto& grouping = *groupings.at(0);
    auto tf = grouping.get_track_fitting();
    if (!tf) {
        log->warn("UbooneNumuBDTScorer: no TrackFitting in grouping '{}'", m_grouping_name);
        return;
    }

    Clus::PR::TaggerInfo& ti = tf->get_tagger_info_mutable();
    const Clus::PR::KineInfo& ki = tf->get_kine_info();

    cal_numu_bdts_xgboost(ti, ki);
}

// ===========================================================================
// cal_cosmict_10_bdt
//
// Scores upstream-dirt clusters (per-cluster, vector features).
// Iterates over cosmict_10_* vectors; returns the minimum BDT score
// (most cosmic-like element wins).
//
// Prototype: NeutrinoID_numu_bdts.h::cal_cosmict_10_bdt()
// ===========================================================================
float UbooneNumuBDTScorer::cal_cosmict_10_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    if (!m_reader_cosmict10) return default_val;

    float val = default_val;

    if (!ti.cosmict_10_length.empty()) {
        val = 1e9f;
        for (size_t i = 0; i < ti.cosmict_10_length.size(); ++i) {
            m_cosmict10_vtx_z        = ti.cosmict_10_vtx_z.at(i);
            m_cosmict10_flag_shower  = ti.cosmict_10_flag_shower.at(i);
            m_cosmict10_flag_dir_weak= ti.cosmict_10_flag_dir_weak.at(i);
            m_cosmict10_angle_beam   = ti.cosmict_10_angle_beam.at(i);
            m_cosmict10_length       = ti.cosmict_10_length.at(i);

            if (std::isnan(m_cosmict10_angle_beam)) m_cosmict10_angle_beam = 0;

            float tmp_bdt = m_reader_cosmict10->EvaluateMVA("MyBDT");
            if (tmp_bdt < val) val = tmp_bdt;
        }
    }

    return val;
}

// ===========================================================================
// cal_numu_1_bdt
//
// Scores segments directly at main_vertex (flag-1 features, per-segment vector).
// Returns the maximum BDT score over all flag-1 candidates.
//
// Prototype: NeutrinoID_numu_bdts.h::cal_numu_1_bdt()
// ===========================================================================
float UbooneNumuBDTScorer::cal_numu_1_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    if (!m_reader_numu1) return default_val;

    float val = default_val;

    if (!ti.numu_cc_1_particle_type.empty()) {
        val = -1e9f;
        for (size_t i = 0; i < ti.numu_cc_1_particle_type.size(); ++i) {
            (void)ti.numu_cc_flag_1.at(i); // loaded but not used as BDT input variable
            m_numu1_particle_type    = ti.numu_cc_1_particle_type.at(i);
            m_numu1_length           = ti.numu_cc_1_length.at(i);
            m_numu1_medium_dQ_dx     = ti.numu_cc_1_medium_dQ_dx.at(i);
            m_numu1_dQ_dx_cut        = ti.numu_cc_1_dQ_dx_cut.at(i);
            m_numu1_direct_length    = ti.numu_cc_1_direct_length.at(i);
            m_numu1_n_daughter_tracks = ti.numu_cc_1_n_daughter_tracks.at(i);
            m_numu1_n_daughter_all   = ti.numu_cc_1_n_daughter_all.at(i);

            if (std::isinf(m_numu1_dQ_dx_cut)) m_numu1_dQ_dx_cut = 10;

            float tmp_bdt = m_reader_numu1->EvaluateMVA("MyBDT");
            if (tmp_bdt > val) val = tmp_bdt;
        }
    }

    return val;
}

// ===========================================================================
// cal_numu_2_bdt
//
// Scores long-muon showers (flag-2 features, per-shower vector).
// Returns the maximum BDT score.
//
// Prototype: NeutrinoID_numu_bdts.h::cal_numu_2_bdt()
// ===========================================================================
float UbooneNumuBDTScorer::cal_numu_2_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    if (!m_reader_numu2) return default_val;

    float val = default_val;

    if (!ti.numu_cc_2_length.empty()) {
        val = -1e9f;
        for (size_t i = 0; i < ti.numu_cc_2_length.size(); ++i) {
            m_numu2_length            = ti.numu_cc_2_length.at(i);
            m_numu2_total_length      = ti.numu_cc_2_total_length.at(i);
            m_numu2_n_daughter_tracks = ti.numu_cc_2_n_daughter_tracks.at(i);
            m_numu2_n_daughter_all    = ti.numu_cc_2_n_daughter_all.at(i);

            float tmp_bdt = m_reader_numu2->EvaluateMVA("MyBDT");
            if (tmp_bdt > val) val = tmp_bdt;
        }
    }

    return val;
}

// ===========================================================================
// cal_numu_3_bdt
//
// Scores the flag-3 (indirect muon) check using scalar features.
// Single EvaluateMVA call — no iteration.
//
// Prototype: NeutrinoID_numu_bdts.h::cal_numu_3_bdt()
// ===========================================================================
float UbooneNumuBDTScorer::cal_numu_3_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    if (!m_reader_numu3) return default_val;
    if (ti.numu_cc_flag_3 == 0) return default_val;

    m_numu3_particle_type    = ti.numu_cc_3_particle_type;
    m_numu3_max_length       = ti.numu_cc_3_max_length;
    m_numu3_track_length     = ti.numu_cc_3_track_length;
    m_numu3_max_length_all   = ti.numu_cc_3_max_length_all;
    m_numu3_max_muon_length  = ti.numu_cc_3_max_muon_length;
    m_numu3_n_daughter_tracks = ti.numu_cc_3_n_daughter_tracks;
    m_numu3_n_daughter_all   = ti.numu_cc_3_n_daughter_all;

    return m_reader_numu3->EvaluateMVA("MyBDT");
}

// ===========================================================================
// cal_numu_bdts_xgboost
//
// Top-level scorer for numu CC identification.
// 1. Evaluates four sub-BDTs, storing intermediate scores in ti.
// 2. Runs a final xgboost TMVA model (~70 input features) → ti.numu_score.
//
// Prototype: NeutrinoID_numu_bdts.h::cal_numu_bdts_xgboost()
// ===========================================================================
void UbooneNumuBDTScorer::cal_numu_bdts_xgboost(Clus::PR::TaggerInfo& ti,
                                                  const Clus::PR::KineInfo& ki) const
{
    // Fill sub-BDT scores first; they feed into the main xgboost model.
    ti.numu_1_score     = cal_numu_1_bdt   (ti, -0.4f);
    ti.numu_2_score     = cal_numu_2_bdt   (ti, -0.1f);
    ti.numu_3_score     = cal_numu_3_bdt   (ti, -0.2f);
    ti.cosmict_10_score = cal_cosmict_10_bdt(ti,  0.7f);

    // Copy TaggerInfo/KineInfo fields into the persistent float buffer.
    // Index order must match AddVariable order in init_readers().

    // --- numu flag-3 features ---
    m_xgb_vars[ 0] = ti.numu_cc_flag_3;
    m_xgb_vars[ 1] = ti.numu_cc_3_particle_type;
    m_xgb_vars[ 2] = ti.numu_cc_3_max_length;
    m_xgb_vars[ 3] = ti.numu_cc_3_track_length;
    m_xgb_vars[ 4] = ti.numu_cc_3_max_length_all;
    m_xgb_vars[ 5] = ti.numu_cc_3_max_muon_length;
    m_xgb_vars[ 6] = ti.numu_cc_3_n_daughter_tracks;
    m_xgb_vars[ 7] = ti.numu_cc_3_n_daughter_all;

    // --- cosmic tagger flag-2 features ---
    m_xgb_vars[ 8] = ti.cosmict_flag_2;
    m_xgb_vars[ 9] = ti.cosmict_2_filled;
    m_xgb_vars[10] = ti.cosmict_2_particle_type;
    m_xgb_vars[11] = ti.cosmict_2_n_muon_tracks;
    m_xgb_vars[12] = ti.cosmict_2_total_shower_length;
    m_xgb_vars[13] = ti.cosmict_2_flag_inside;
    m_xgb_vars[14] = ti.cosmict_2_angle_beam;
    m_xgb_vars[15] = ti.cosmict_2_flag_dir_weak;
    m_xgb_vars[16] = ti.cosmict_2_dQ_dx_end;
    m_xgb_vars[17] = ti.cosmict_2_dQ_dx_front;
    m_xgb_vars[18] = ti.cosmict_2_theta;
    m_xgb_vars[19] = ti.cosmict_2_phi;
    m_xgb_vars[20] = ti.cosmict_2_valid_tracks;

    // --- cosmic tagger flag-4 features ---
    m_xgb_vars[21] = ti.cosmict_flag_4;
    m_xgb_vars[22] = ti.cosmict_4_filled;
    m_xgb_vars[23] = ti.cosmict_4_flag_inside;
    m_xgb_vars[24] = ti.cosmict_4_angle_beam;
    m_xgb_vars[25] = ti.cosmict_4_connected_showers;

    // --- cosmic tagger flag-3 features ---
    m_xgb_vars[26] = ti.cosmict_flag_3;
    m_xgb_vars[27] = ti.cosmict_3_filled;
    m_xgb_vars[28] = ti.cosmict_3_flag_inside;
    m_xgb_vars[29] = ti.cosmict_3_angle_beam;
    m_xgb_vars[30] = ti.cosmict_3_flag_dir_weak;
    m_xgb_vars[31] = ti.cosmict_3_dQ_dx_end;
    m_xgb_vars[32] = ti.cosmict_3_dQ_dx_front;
    m_xgb_vars[33] = ti.cosmict_3_theta;
    m_xgb_vars[34] = ti.cosmict_3_phi;
    m_xgb_vars[35] = ti.cosmict_3_valid_tracks;

    // --- cosmic tagger flag-5 features ---
    m_xgb_vars[36] = ti.cosmict_flag_5;
    m_xgb_vars[37] = ti.cosmict_5_filled;
    m_xgb_vars[38] = ti.cosmict_5_flag_inside;
    m_xgb_vars[39] = ti.cosmict_5_angle_beam;
    m_xgb_vars[40] = ti.cosmict_5_connected_showers;

    // --- cosmic tagger flag-6 features ---
    m_xgb_vars[41] = ti.cosmict_flag_6;
    m_xgb_vars[42] = ti.cosmict_6_filled;
    m_xgb_vars[43] = ti.cosmict_6_flag_dir_weak;
    m_xgb_vars[44] = ti.cosmict_6_flag_inside;
    m_xgb_vars[45] = ti.cosmict_6_angle;

    // --- cosmic tagger flag-7 features ---
    m_xgb_vars[46] = ti.cosmict_flag_7;
    m_xgb_vars[47] = ti.cosmict_7_filled;
    m_xgb_vars[48] = ti.cosmict_7_flag_sec;
    m_xgb_vars[49] = ti.cosmict_7_n_muon_tracks;
    m_xgb_vars[50] = ti.cosmict_7_total_shower_length;
    m_xgb_vars[51] = ti.cosmict_7_flag_inside;
    m_xgb_vars[52] = ti.cosmict_7_angle_beam;
    m_xgb_vars[53] = ti.cosmict_7_flag_dir_weak;
    m_xgb_vars[54] = ti.cosmict_7_dQ_dx_end;
    m_xgb_vars[55] = ti.cosmict_7_dQ_dx_front;
    m_xgb_vars[56] = ti.cosmict_7_theta;
    m_xgb_vars[57] = ti.cosmict_7_phi;

    // --- cosmic tagger flag-8 features ---
    m_xgb_vars[58] = ti.cosmict_flag_8;
    m_xgb_vars[59] = ti.cosmict_8_filled;
    m_xgb_vars[60] = ti.cosmict_8_flag_out;
    m_xgb_vars[61] = ti.cosmict_8_muon_length;
    m_xgb_vars[62] = ti.cosmict_8_acc_length;

    // --- cosmic tagger flag-9 features ---
    m_xgb_vars[63] = ti.cosmict_flag_9;

    // --- top-level cosmic / numu flags ---
    m_xgb_vars[64] = ti.cosmic_flag;
    m_xgb_vars[65] = ti.cosmic_filled;
    m_xgb_vars[66] = ti.cosmict_flag;
    m_xgb_vars[67] = ti.numu_cc_flag;
    m_xgb_vars[68] = ti.cosmict_flag_1;

    // --- kinematics + fiducial flag ---
    m_xgb_vars[69] = static_cast<float>(ki.kine_reco_Enu);
    m_xgb_vars[70] = ti.match_isFC;

    // --- sub-BDT scores (filled above) ---
    m_xgb_vars[71] = ti.cosmict_10_score;
    m_xgb_vars[72] = ti.numu_1_score;
    m_xgb_vars[73] = ti.numu_2_score;

    // Guard against NaN inputs that can cause TMVA to crash.
    if (std::isnan(m_xgb_vars[24])) m_xgb_vars[24] = 0; // cosmict_4_angle_beam
    if (std::isnan(m_xgb_vars[52])) m_xgb_vars[52] = 0; // cosmict_7_angle_beam
    if (std::isnan(m_xgb_vars[56])) m_xgb_vars[56] = 0; // cosmict_7_theta
    if (std::isnan(m_xgb_vars[57])) m_xgb_vars[57] = 0; // cosmict_7_phi

    double val1 = m_reader_xgboost->EvaluateMVA("MyBDT");

    // Clamp to avoid division by zero in log-odds transformation
    val1 = std::max(-0.9999, std::min(0.9999, val1));
    // Convert raw TMVA output to log-likelihood ratio (matches prototype).
    ti.numu_score = static_cast<float>(std::log10((1.0 + val1) / (1.0 - val1)));
}
