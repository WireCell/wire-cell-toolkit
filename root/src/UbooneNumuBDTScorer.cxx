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

#include "TMVA/Reader.h"

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
    m_grouping_name  = get<std::string>(cfg, "grouping",        "live");
    m_numu1_xml      = get<std::string>(cfg, "numu1_weights_xml",    "");
    m_numu2_xml      = get<std::string>(cfg, "numu2_weights_xml",    "");
    m_numu3_xml      = get<std::string>(cfg, "numu3_weights_xml",    "");
    m_cosmict10_xml  = get<std::string>(cfg, "cosmict10_weights_xml","");
    m_numu_xgboost_xml = get<std::string>(cfg, "numu_xgboost_xml",  "");
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

void UbooneNumuBDTScorer::visit(Clus::Facade::Ensemble& ensemble) const
{
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
    float val = default_val;

    // Local scalar vars are required: TMVA::Reader needs float* for per-entry variables.
    float cosmict_10_vtx_z;
    float cosmict_10_flag_shower;
    float cosmict_10_flag_dir_weak;
    float cosmict_10_angle_beam;
    float cosmict_10_length;

    TMVA::Reader reader_cosmict_10;
    reader_cosmict_10.AddVariable("cosmict_10_vtx_z",        &cosmict_10_vtx_z);
    reader_cosmict_10.AddVariable("cosmict_10_flag_shower",  &cosmict_10_flag_shower);
    reader_cosmict_10.AddVariable("cosmict_10_flag_dir_weak",&cosmict_10_flag_dir_weak);
    reader_cosmict_10.AddVariable("cosmict_10_angle_beam",   &cosmict_10_angle_beam);
    reader_cosmict_10.AddVariable("cosmict_10_length",       &cosmict_10_length);

    reader_cosmict_10.BookMVA("MyBDT", m_cosmict10_xml);

    if (!ti.cosmict_10_length.empty()) {
        val = 1e9f;
        for (size_t i = 0; i < ti.cosmict_10_length.size(); ++i) {
            cosmict_10_vtx_z        = ti.cosmict_10_vtx_z.at(i);
            cosmict_10_flag_shower  = ti.cosmict_10_flag_shower.at(i);
            cosmict_10_flag_dir_weak= ti.cosmict_10_flag_dir_weak.at(i);
            cosmict_10_angle_beam   = ti.cosmict_10_angle_beam.at(i);
            cosmict_10_length       = ti.cosmict_10_length.at(i);

            if (std::isnan(cosmict_10_angle_beam)) cosmict_10_angle_beam = 0;

            float tmp_bdt = reader_cosmict_10.EvaluateMVA("MyBDT");
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
    float val = default_val;

    float numu_cc_1_particle_type;
    float numu_cc_1_length;
    float numu_cc_1_medium_dQ_dx;
    float numu_cc_1_dQ_dx_cut;
    float numu_cc_1_direct_length;
    float numu_cc_1_n_daughter_tracks;
    float numu_cc_1_n_daughter_all;

    TMVA::Reader reader_numu_1;
    reader_numu_1.AddVariable("numu_cc_1_particle_type",    &numu_cc_1_particle_type);
    reader_numu_1.AddVariable("numu_cc_1_length",           &numu_cc_1_length);
    reader_numu_1.AddVariable("numu_cc_1_medium_dQ_dx",     &numu_cc_1_medium_dQ_dx);
    reader_numu_1.AddVariable("numu_cc_1_dQ_dx_cut",        &numu_cc_1_dQ_dx_cut);
    reader_numu_1.AddVariable("numu_cc_1_direct_length",    &numu_cc_1_direct_length);
    reader_numu_1.AddVariable("numu_cc_1_n_daughter_tracks",&numu_cc_1_n_daughter_tracks);
    reader_numu_1.AddVariable("numu_cc_1_n_daughter_all",   &numu_cc_1_n_daughter_all);

    reader_numu_1.BookMVA("MyBDT", m_numu1_xml);

    if (!ti.numu_cc_1_particle_type.empty()) {
        val = -1e9f;
        for (size_t i = 0; i < ti.numu_cc_1_particle_type.size(); ++i) {
            (void)ti.numu_cc_flag_1.at(i); // loaded but not used as BDT input variable
            numu_cc_1_particle_type   = ti.numu_cc_1_particle_type.at(i);
            numu_cc_1_length          = ti.numu_cc_1_length.at(i);
            numu_cc_1_medium_dQ_dx    = ti.numu_cc_1_medium_dQ_dx.at(i);
            numu_cc_1_dQ_dx_cut       = ti.numu_cc_1_dQ_dx_cut.at(i);
            numu_cc_1_direct_length   = ti.numu_cc_1_direct_length.at(i);
            numu_cc_1_n_daughter_tracks = ti.numu_cc_1_n_daughter_tracks.at(i);
            numu_cc_1_n_daughter_all  = ti.numu_cc_1_n_daughter_all.at(i);

            if (std::isinf(numu_cc_1_dQ_dx_cut)) numu_cc_1_dQ_dx_cut = 10;

            float tmp_bdt = reader_numu_1.EvaluateMVA("MyBDT");
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
    float val = default_val;

    float numu_cc_2_length;
    float numu_cc_2_total_length;
    float numu_cc_2_n_daughter_tracks;
    float numu_cc_2_n_daughter_all;

    TMVA::Reader reader_numu_2;
    reader_numu_2.AddVariable("numu_cc_2_length",           &numu_cc_2_length);
    reader_numu_2.AddVariable("numu_cc_2_total_length",     &numu_cc_2_total_length);
    reader_numu_2.AddVariable("numu_cc_2_n_daughter_tracks",&numu_cc_2_n_daughter_tracks);
    reader_numu_2.AddVariable("numu_cc_2_n_daughter_all",   &numu_cc_2_n_daughter_all);

    reader_numu_2.BookMVA("MyBDT", m_numu2_xml);

    if (!ti.numu_cc_2_length.empty()) {
        val = -1e9f;
        for (size_t i = 0; i < ti.numu_cc_2_length.size(); ++i) {
            numu_cc_2_length           = ti.numu_cc_2_length.at(i);
            numu_cc_2_total_length     = ti.numu_cc_2_total_length.at(i);
            numu_cc_2_n_daughter_tracks= ti.numu_cc_2_n_daughter_tracks.at(i);
            numu_cc_2_n_daughter_all   = ti.numu_cc_2_n_daughter_all.at(i);

            float tmp_bdt = reader_numu_2.EvaluateMVA("MyBDT");
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
    float val = default_val;

    TMVA::Reader reader_numu_3;
    // Variable name in XML is "numu_cc_3_acc_track_length" — must match training.
    reader_numu_3.AddVariable("numu_cc_3_particle_type",   &ti.numu_cc_3_particle_type);
    reader_numu_3.AddVariable("numu_cc_3_max_length",      &ti.numu_cc_3_max_length);
    reader_numu_3.AddVariable("numu_cc_3_acc_track_length",&ti.numu_cc_3_acc_track_length);
    reader_numu_3.AddVariable("numu_cc_3_max_length_all",  &ti.numu_cc_3_max_length_all);
    reader_numu_3.AddVariable("numu_cc_3_max_muon_length", &ti.numu_cc_3_max_muon_length);
    reader_numu_3.AddVariable("numu_cc_3_n_daughter_tracks",&ti.numu_cc_3_n_daughter_tracks);
    reader_numu_3.AddVariable("numu_cc_3_n_daughter_all",  &ti.numu_cc_3_n_daughter_all);

    reader_numu_3.BookMVA("MyBDT", m_numu3_xml);

    val = reader_numu_3.EvaluateMVA("MyBDT");

    return val;
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

    // The main xgboost model uses scalar TaggerInfo fields plus kine and match_isFC.
    // Need a local mutable copy of kine_reco_Enu for TMVA (needs float*).
    float kine_reco_Enu = static_cast<float>(ki.kine_reco_Enu);

    TMVA::Reader reader;

    // --- numu flag-3 features ---
    reader.AddVariable("numu_cc_flag_3",         &ti.numu_cc_flag_3);
    reader.AddVariable("numu_cc_3_particle_type",&ti.numu_cc_3_particle_type);
    reader.AddVariable("numu_cc_3_max_length",   &ti.numu_cc_3_max_length);
    // NOTE: training variable name "numu_cc_3_track_length" maps to acc_track_length member
    reader.AddVariable("numu_cc_3_track_length", &ti.numu_cc_3_acc_track_length);
    reader.AddVariable("numu_cc_3_max_length_all",   &ti.numu_cc_3_max_length_all);
    reader.AddVariable("numu_cc_3_max_muon_length",  &ti.numu_cc_3_max_muon_length);
    reader.AddVariable("numu_cc_3_n_daughter_tracks",&ti.numu_cc_3_n_daughter_tracks);
    reader.AddVariable("numu_cc_3_n_daughter_all",   &ti.numu_cc_3_n_daughter_all);

    // --- cosmic tagger flag-2 features ---
    reader.AddVariable("cosmict_flag_2",               &ti.cosmict_flag_2);
    reader.AddVariable("cosmict_2_filled",             &ti.cosmict_2_filled);
    reader.AddVariable("cosmict_2_particle_type",      &ti.cosmict_2_particle_type);
    reader.AddVariable("cosmict_2_n_muon_tracks",      &ti.cosmict_2_n_muon_tracks);
    reader.AddVariable("cosmict_2_total_shower_length",&ti.cosmict_2_total_shower_length);
    reader.AddVariable("cosmict_2_flag_inside",        &ti.cosmict_2_flag_inside);
    reader.AddVariable("cosmict_2_angle_beam",         &ti.cosmict_2_angle_beam);
    reader.AddVariable("cosmict_2_flag_dir_weak",      &ti.cosmict_2_flag_dir_weak);
    reader.AddVariable("cosmict_2_dQ_dx_end",          &ti.cosmict_2_dQ_dx_end);
    reader.AddVariable("cosmict_2_dQ_dx_front",        &ti.cosmict_2_dQ_dx_front);
    reader.AddVariable("cosmict_2_theta",              &ti.cosmict_2_theta);
    reader.AddVariable("cosmict_2_phi",                &ti.cosmict_2_phi);
    reader.AddVariable("cosmict_2_valid_tracks",       &ti.cosmict_2_valid_tracks);

    // --- cosmic tagger flag-4 features ---
    reader.AddVariable("cosmict_flag_4",            &ti.cosmict_flag_4);
    reader.AddVariable("cosmict_4_filled",          &ti.cosmict_4_filled);
    reader.AddVariable("cosmict_4_flag_inside",     &ti.cosmict_4_flag_inside);
    reader.AddVariable("cosmict_4_angle_beam",      &ti.cosmict_4_angle_beam);
    reader.AddVariable("cosmict_4_connected_showers",&ti.cosmict_4_connected_showers);

    // --- cosmic tagger flag-3 features ---
    reader.AddVariable("cosmict_flag_3",        &ti.cosmict_flag_3);
    reader.AddVariable("cosmict_3_filled",      &ti.cosmict_3_filled);
    reader.AddVariable("cosmict_3_flag_inside", &ti.cosmict_3_flag_inside);
    reader.AddVariable("cosmict_3_angle_beam",  &ti.cosmict_3_angle_beam);
    reader.AddVariable("cosmict_3_flag_dir_weak",&ti.cosmict_3_flag_dir_weak);
    reader.AddVariable("cosmict_3_dQ_dx_end",   &ti.cosmict_3_dQ_dx_end);
    reader.AddVariable("cosmict_3_dQ_dx_front", &ti.cosmict_3_dQ_dx_front);
    reader.AddVariable("cosmict_3_theta",       &ti.cosmict_3_theta);
    reader.AddVariable("cosmict_3_phi",         &ti.cosmict_3_phi);
    reader.AddVariable("cosmict_3_valid_tracks",&ti.cosmict_3_valid_tracks);

    // --- cosmic tagger flag-5 features ---
    reader.AddVariable("cosmict_flag_5",            &ti.cosmict_flag_5);
    reader.AddVariable("cosmict_5_filled",          &ti.cosmict_5_filled);
    reader.AddVariable("cosmict_5_flag_inside",     &ti.cosmict_5_flag_inside);
    reader.AddVariable("cosmict_5_angle_beam",      &ti.cosmict_5_angle_beam);
    reader.AddVariable("cosmict_5_connected_showers",&ti.cosmict_5_connected_showers);

    // --- cosmic tagger flag-6 features ---
    reader.AddVariable("cosmict_flag_6",       &ti.cosmict_flag_6);
    reader.AddVariable("cosmict_6_filled",     &ti.cosmict_6_filled);
    reader.AddVariable("cosmict_6_flag_dir_weak",&ti.cosmict_6_flag_dir_weak);
    reader.AddVariable("cosmict_6_flag_inside",&ti.cosmict_6_flag_inside);
    reader.AddVariable("cosmict_6_angle",      &ti.cosmict_6_angle);

    // --- cosmic tagger flag-7 features ---
    reader.AddVariable("cosmict_flag_7",              &ti.cosmict_flag_7);
    reader.AddVariable("cosmict_7_filled",            &ti.cosmict_7_filled);
    reader.AddVariable("cosmict_7_flag_sec",          &ti.cosmict_7_flag_sec);
    reader.AddVariable("cosmict_7_n_muon_tracks",     &ti.cosmict_7_n_muon_tracks);
    reader.AddVariable("cosmict_7_total_shower_length",&ti.cosmict_7_total_shower_length);
    reader.AddVariable("cosmict_7_flag_inside",       &ti.cosmict_7_flag_inside);
    reader.AddVariable("cosmict_7_angle_beam",        &ti.cosmict_7_angle_beam);
    reader.AddVariable("cosmict_7_flag_dir_weak",     &ti.cosmict_7_flag_dir_weak);
    reader.AddVariable("cosmict_7_dQ_dx_end",         &ti.cosmict_7_dQ_dx_end);
    reader.AddVariable("cosmict_7_dQ_dx_front",       &ti.cosmict_7_dQ_dx_front);
    reader.AddVariable("cosmict_7_theta",             &ti.cosmict_7_theta);
    reader.AddVariable("cosmict_7_phi",               &ti.cosmict_7_phi);

    // --- cosmic tagger flag-8 features ---
    reader.AddVariable("cosmict_flag_8",      &ti.cosmict_flag_8);
    reader.AddVariable("cosmict_8_filled",    &ti.cosmict_8_filled);
    reader.AddVariable("cosmict_8_flag_out",  &ti.cosmict_8_flag_out);
    reader.AddVariable("cosmict_8_muon_length",&ti.cosmict_8_muon_length);
    reader.AddVariable("cosmict_8_acc_length",&ti.cosmict_8_acc_length);

    // --- cosmic tagger flag-9 features ---
    reader.AddVariable("cosmict_flag_9",  &ti.cosmict_flag_9);

    // --- top-level cosmic / numu flags ---
    reader.AddVariable("cosmic_flag",    &ti.cosmic_flag);
    reader.AddVariable("cosmic_filled",  &ti.cosmic_filled);
    reader.AddVariable("cosmict_flag",   &ti.cosmict_flag);
    reader.AddVariable("numu_cc_flag",   &ti.numu_cc_flag);
    reader.AddVariable("cosmict_flag_1", &ti.cosmict_flag_1);

    // --- kinematics + fiducial flag ---
    reader.AddVariable("kine_reco_Enu", &kine_reco_Enu);
    reader.AddVariable("match_isFC",    &ti.match_isFC);

    // --- sub-BDT scores (filled above) ---
    reader.AddVariable("cosmict_10_score", &ti.cosmict_10_score);
    reader.AddVariable("numu_1_score",     &ti.numu_1_score);
    reader.AddVariable("numu_2_score",     &ti.numu_2_score);

    reader.BookMVA("MyBDT", m_numu_xgboost_xml);

    // Guard against NaN inputs that can cause TMVA to crash.
    if (std::isnan(ti.cosmict_4_angle_beam)) ti.cosmict_4_angle_beam = 0;
    if (std::isnan(ti.cosmict_7_angle_beam)) ti.cosmict_7_angle_beam = 0;
    if (std::isnan(ti.cosmict_7_theta))      ti.cosmict_7_theta = 0;
    if (std::isnan(ti.cosmict_7_phi))        ti.cosmict_7_phi = 0;

    double val1 = reader.EvaluateMVA("MyBDT");

    // Convert raw TMVA output to log-likelihood ratio (matches prototype).
    ti.numu_score = static_cast<float>(std::log10((1.0 + val1) / (1.0 - val1)));
}
