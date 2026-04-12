// UbooneNueBDTScorer.cxx
//
// IEnsembleVisitor that runs TMVA BDT scoring for the nueCC tagger.
// Ports the following functions from prototype_pid/src/NeutrinoID_nue_bdts.h:
//   cal_bdts_xgboost()     — top-level scorer, writes ti.nue_score
//   cal_mipid_bdt()        — sub-BDT on mip_id features
//   cal_gap_bdt()          — sub-BDT on gap features
//   cal_hol_lol_bdt()      — sub-BDT on hol+lol_3 features
//   cal_cme_anc_bdt()      — sub-BDT on cme+anc features
//   cal_mgo_mgt_bdt()      — sub-BDT on mgo+mgt features
//   cal_br1_bdt()          — sub-BDT on br1 features
//   cal_br3_bdt()          — sub-BDT on br3 (non-vector) features
//   cal_br3_3_bdt()        — sub-BDT on br3_3 per-segment vector
//   cal_br3_5_bdt()        — sub-BDT on br3_5 per-segment vector
//   cal_br3_6_bdt()        — sub-BDT on br3_6 per-segment vector
//   cal_stemdir_br2_bdt()  — sub-BDT on stem_dir+br2 features
//   cal_trimuon_bdt()      — sub-BDT on stem_len+brm+lem features
//   cal_br4_tro_bdt()      — sub-BDT on br4+tro_3 features
//   cal_mipquality_bdt()   — sub-BDT on mip_quality features
//   cal_pio_1_bdt()        — sub-BDT on pio_1 features
//   cal_pio_2_bdt()        — sub-BDT on pio_2 per-pi0 vector
//   cal_stw_spt_bdt()      — sub-BDT on stw_1+spt features
//   cal_vis_1_bdt()        — sub-BDT on vis_1 features
//   cal_vis_2_bdt()        — sub-BDT on vis_2 features
//   cal_stw_2_bdt()        — sub-BDT on stw_2 per-segment vector
//   cal_stw_3_bdt()        — sub-BDT on stw_3 per-segment vector
//   cal_stw_4_bdt()        — sub-BDT on stw_4 per-segment vector
//   cal_sig_1_bdt()        — sub-BDT on sig_1 per-segment vector
//   cal_sig_2_bdt()        — sub-BDT on sig_2 per-segment vector
//   cal_lol_1_bdt()        — sub-BDT on lol_1 per-segment vector
//   cal_lol_2_bdt()        — sub-BDT on lol_2 per-segment vector
//   cal_tro_1_bdt()        — sub-BDT on tro_1 per-segment vector
//   cal_tro_2_bdt()        — sub-BDT on tro_2 per-segment vector
//   cal_tro_4_bdt()        — sub-BDT on tro_4 per-segment vector
//   cal_tro_5_bdt()        — sub-BDT on tro_5 per-segment vector
//
// NOT ported: cal_bdts() (old TMVA combination variant).
//
// Translation conventions vs. prototype:
//   tagger_info.xxx                    →  ti.xxx
//   kine_info.kine_reco_Enu            →  ki.kine_reco_Enu
//   match_isFC (NeutrinoID member var) →  ti.match_isFC
//   "input_data_files/weights/foo.xml" →  m_*_xml  (configured via wc.resolve)
//   TMath::Log10((1+v)/(1-v))          →  std::log10((1.0+v)/(1.0-v))

#include "WireCellRoot/UbooneNueBDTScorer.h"
#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Persist.h"

#include "TMVA/Reader.h"

#include <cmath>

WIRECELL_FACTORY(UbooneNueBDTScorer, WireCell::Root::UbooneNueBDTScorer,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)

using namespace WireCell;
using namespace WireCell::Root;

UbooneNueBDTScorer::UbooneNueBDTScorer()
    : log(Log::logger("root.UbooneNueBDTScorer"))
{
}

void UbooneNueBDTScorer::configure(const WireCell::Configuration& cfg)
{
    auto resolve = [&](const std::string& p) {
        if (p.empty()) return p;
        std::string r = Persist::resolve(p);
        if (r.empty()) log->error("UbooneNueBDTScorer: weight file not found: {}", p);
        return r;
    };
    m_grouping_name   = get<std::string>(cfg, "grouping", "live");

    m_mipid_xml       = resolve(get<std::string>(cfg, "mipid_weights_xml",       ""));
    m_gap_xml         = resolve(get<std::string>(cfg, "gap_weights_xml",         ""));
    m_hol_lol_xml     = resolve(get<std::string>(cfg, "hol_lol_weights_xml",     ""));
    m_cme_anc_xml     = resolve(get<std::string>(cfg, "cme_anc_weights_xml",     ""));
    m_mgo_mgt_xml     = resolve(get<std::string>(cfg, "mgo_mgt_weights_xml",     ""));
    m_br1_xml         = resolve(get<std::string>(cfg, "br1_weights_xml",         ""));
    m_br3_xml         = resolve(get<std::string>(cfg, "br3_weights_xml",         ""));
    m_br3_3_xml       = resolve(get<std::string>(cfg, "br3_3_weights_xml",       ""));
    m_br3_5_xml       = resolve(get<std::string>(cfg, "br3_5_weights_xml",       ""));
    m_br3_6_xml       = resolve(get<std::string>(cfg, "br3_6_weights_xml",       ""));
    m_stemdir_br2_xml = resolve(get<std::string>(cfg, "stemdir_br2_weights_xml", ""));
    m_trimuon_xml     = resolve(get<std::string>(cfg, "trimuon_weights_xml",     ""));
    m_br4_tro_xml     = resolve(get<std::string>(cfg, "br4_tro_weights_xml",     ""));
    m_mipquality_xml  = resolve(get<std::string>(cfg, "mipquality_weights_xml",  ""));
    m_pio_1_xml       = resolve(get<std::string>(cfg, "pio_1_weights_xml",       ""));
    m_pio_2_xml       = resolve(get<std::string>(cfg, "pio_2_weights_xml",       ""));
    m_stw_spt_xml     = resolve(get<std::string>(cfg, "stw_spt_weights_xml",     ""));
    m_vis_1_xml       = resolve(get<std::string>(cfg, "vis_1_weights_xml",       ""));
    m_vis_2_xml       = resolve(get<std::string>(cfg, "vis_2_weights_xml",       ""));
    m_stw_2_xml       = resolve(get<std::string>(cfg, "stw_2_weights_xml",       ""));
    m_stw_3_xml       = resolve(get<std::string>(cfg, "stw_3_weights_xml",       ""));
    m_stw_4_xml       = resolve(get<std::string>(cfg, "stw_4_weights_xml",       ""));
    m_sig_1_xml       = resolve(get<std::string>(cfg, "sig_1_weights_xml",       ""));
    m_sig_2_xml       = resolve(get<std::string>(cfg, "sig_2_weights_xml",       ""));
    m_lol_1_xml       = resolve(get<std::string>(cfg, "lol_1_weights_xml",       ""));
    m_lol_2_xml       = resolve(get<std::string>(cfg, "lol_2_weights_xml",       ""));
    m_tro_1_xml       = resolve(get<std::string>(cfg, "tro_1_weights_xml",       ""));
    m_tro_2_xml       = resolve(get<std::string>(cfg, "tro_2_weights_xml",       ""));
    m_tro_4_xml       = resolve(get<std::string>(cfg, "tro_4_weights_xml",       ""));
    m_tro_5_xml       = resolve(get<std::string>(cfg, "tro_5_weights_xml",       ""));
    m_nue_xgboost_xml = resolve(get<std::string>(cfg, "nue_xgboost_xml",         ""));

    init_readers();
}

Configuration UbooneNueBDTScorer::default_configuration() const
{
    Configuration cfg;
    cfg["grouping"]              = "live";
    cfg["mipid_weights_xml"]       = "";  // e.g. wc.resolve("uboone/weights/mipid_BDT.weights.xml")
    cfg["gap_weights_xml"]         = "";
    cfg["hol_lol_weights_xml"]     = "";
    cfg["cme_anc_weights_xml"]     = "";
    cfg["mgo_mgt_weights_xml"]     = "";
    cfg["br1_weights_xml"]         = "";
    cfg["br3_weights_xml"]         = "";
    cfg["br3_3_weights_xml"]       = "";
    cfg["br3_5_weights_xml"]       = "";
    cfg["br3_6_weights_xml"]       = "";
    cfg["stemdir_br2_weights_xml"] = "";
    cfg["trimuon_weights_xml"]     = "";
    cfg["br4_tro_weights_xml"]     = "";
    cfg["mipquality_weights_xml"]  = "";
    cfg["pio_1_weights_xml"]       = "";
    cfg["pio_2_weights_xml"]       = "";
    cfg["stw_spt_weights_xml"]     = "";
    cfg["vis_1_weights_xml"]       = "";
    cfg["vis_2_weights_xml"]       = "";
    cfg["stw_2_weights_xml"]       = "";
    cfg["stw_3_weights_xml"]       = "";
    cfg["stw_4_weights_xml"]       = "";
    cfg["sig_1_weights_xml"]       = "";
    cfg["sig_2_weights_xml"]       = "";
    cfg["lol_1_weights_xml"]       = "";
    cfg["lol_2_weights_xml"]       = "";
    cfg["tro_1_weights_xml"]       = "";
    cfg["tro_2_weights_xml"]       = "";
    cfg["tro_4_weights_xml"]       = "";
    cfg["tro_5_weights_xml"]       = "";
    cfg["nue_xgboost_xml"]         = "";  // e.g. wc.resolve("uboone/weights/XGB_nue_seed2_0923.xml")
    return cfg;
}

// ===========================================================================
// init_readers
//
// Create all TMVA readers once at configure time. Each reader's variables
// are bound to mutable member floats so that evaluate-time code only needs
// to copy values into those members and call EvaluateMVA.
// ===========================================================================
void UbooneNueBDTScorer::init_readers()
{
    // ---- Vector sub-BDT readers ----

    if (!m_br3_3_xml.empty()) {
        m_rdr_br3_3 = std::make_unique<TMVA::Reader>();
        m_rdr_br3_3->AddVariable("br3_3_v_energy",     &m_br3_3_v_energy);
        m_rdr_br3_3->AddVariable("br3_3_v_angle",      &m_br3_3_v_angle);
        m_rdr_br3_3->AddVariable("br3_3_v_dir_length", &m_br3_3_v_dir_length);
        m_rdr_br3_3->AddVariable("br3_3_v_length",     &m_br3_3_v_length);
        m_rdr_br3_3->BookMVA("MyBDT", m_br3_3_xml);
    }

    if (!m_br3_5_xml.empty()) {
        m_rdr_br3_5 = std::make_unique<TMVA::Reader>();
        m_rdr_br3_5->AddVariable("br3_5_v_dir_length",            &m_br3_5_v_dir_length);
        m_rdr_br3_5->AddVariable("br3_5_v_total_length",          &m_br3_5_v_total_length);
        m_rdr_br3_5->AddVariable("br3_5_v_flag_avoid_muon_check", &m_br3_5_v_flag_avoid_muon_check);
        m_rdr_br3_5->AddVariable("br3_5_v_n_seg",                 &m_br3_5_v_n_seg);
        m_rdr_br3_5->AddVariable("br3_5_v_angle",                 &m_br3_5_v_angle);
        m_rdr_br3_5->AddVariable("br3_5_v_sg_length",             &m_br3_5_v_sg_length);
        m_rdr_br3_5->AddVariable("br3_5_v_energy",                &m_br3_5_v_energy);
        m_rdr_br3_5->AddVariable("br3_5_v_n_segs",                &m_br3_5_v_n_segs);
        m_rdr_br3_5->AddVariable("br3_5_v_shower_main_length",    &m_br3_5_v_shower_main_length);
        m_rdr_br3_5->AddVariable("br3_5_v_shower_total_length",   &m_br3_5_v_shower_total_length);
        m_rdr_br3_5->BookMVA("MyBDT", m_br3_5_xml);
    }

    if (!m_br3_6_xml.empty()) {
        m_rdr_br3_6 = std::make_unique<TMVA::Reader>();
        m_rdr_br3_6->AddVariable("br3_6_v_angle",                  &m_br3_6_v_angle);
        m_rdr_br3_6->AddVariable("br3_6_v_angle1",                 &m_br3_6_v_angle1);
        m_rdr_br3_6->AddVariable("br3_6_v_flag_shower_trajectory", &m_br3_6_v_flag_shower_trajectory);
        m_rdr_br3_6->AddVariable("br3_6_v_direct_length",          &m_br3_6_v_direct_length);
        m_rdr_br3_6->AddVariable("br3_6_v_length",                 &m_br3_6_v_length);
        m_rdr_br3_6->AddVariable("br3_6_v_n_other_vtx_segs",       &m_br3_6_v_n_other_vtx_segs);
        m_rdr_br3_6->AddVariable("br3_6_v_energy",                 &m_br3_6_v_energy);
        m_rdr_br3_6->BookMVA("MyBDT", m_br3_6_xml);
    }

    if (!m_pio_2_xml.empty()) {
        m_rdr_pio_2 = std::make_unique<TMVA::Reader>();
        m_rdr_pio_2->AddVariable("pio_2_v_dis2",       &m_pio_2_v_dis2);
        m_rdr_pio_2->AddVariable("pio_2_v_angle2",     &m_pio_2_v_angle2);
        m_rdr_pio_2->AddVariable("pio_2_v_acc_length", &m_pio_2_v_acc_length);
        m_rdr_pio_2->AddVariable("pio_mip_id",         &m_pio_2_v_pio_mip_id);
        m_rdr_pio_2->BookMVA("MyBDT", m_pio_2_xml);
    }

    if (!m_stw_2_xml.empty()) {
        m_rdr_stw_2 = std::make_unique<TMVA::Reader>();
        m_rdr_stw_2->AddVariable("stw_2_v_medium_dQ_dx", &m_stw_2_v_medium_dQ_dx);
        m_rdr_stw_2->AddVariable("stw_2_v_energy",        &m_stw_2_v_energy);
        m_rdr_stw_2->AddVariable("stw_2_v_angle",         &m_stw_2_v_angle);
        m_rdr_stw_2->AddVariable("stw_2_v_dir_length",    &m_stw_2_v_dir_length);
        m_rdr_stw_2->AddVariable("stw_2_v_max_dQ_dx",     &m_stw_2_v_max_dQ_dx);
        m_rdr_stw_2->BookMVA("MyBDT", m_stw_2_xml);
    }

    if (!m_stw_3_xml.empty()) {
        m_rdr_stw_3 = std::make_unique<TMVA::Reader>();
        m_rdr_stw_3->AddVariable("stw_3_v_angle",          &m_stw_3_v_angle);
        m_rdr_stw_3->AddVariable("stw_3_v_dir_length",     &m_stw_3_v_dir_length);
        m_rdr_stw_3->AddVariable("stw_3_v_energy",         &m_stw_3_v_energy);
        m_rdr_stw_3->AddVariable("stw_3_v_medium_dQ_dx",   &m_stw_3_v_medium_dQ_dx);
        m_rdr_stw_3->BookMVA("MyBDT", m_stw_3_xml);
    }

    if (!m_stw_4_xml.empty()) {
        m_rdr_stw_4 = std::make_unique<TMVA::Reader>();
        m_rdr_stw_4->AddVariable("stw_4_v_angle",  &m_stw_4_v_angle);
        m_rdr_stw_4->AddVariable("stw_4_v_dis",    &m_stw_4_v_dis);
        m_rdr_stw_4->AddVariable("stw_4_v_energy", &m_stw_4_v_energy);
        m_rdr_stw_4->BookMVA("MyBDT", m_stw_4_xml);
    }

    if (!m_sig_1_xml.empty()) {
        m_rdr_sig_1 = std::make_unique<TMVA::Reader>();
        m_rdr_sig_1->AddVariable("sig_1_v_angle",               &m_sig_1_v_angle);
        m_rdr_sig_1->AddVariable("sig_1_v_flag_single_shower",  &m_sig_1_v_flag_single_shower);
        m_rdr_sig_1->AddVariable("sig_1_v_energy",              &m_sig_1_v_energy);
        m_rdr_sig_1->AddVariable("sig_1_v_energy_1",            &m_sig_1_v_energy_1);
        m_rdr_sig_1->BookMVA("MyBDT", m_sig_1_xml);
    }

    if (!m_sig_2_xml.empty()) {
        m_rdr_sig_2 = std::make_unique<TMVA::Reader>();
        m_rdr_sig_2->AddVariable("sig_2_v_energy",              &m_sig_2_v_energy);
        m_rdr_sig_2->AddVariable("sig_2_v_shower_angle",        &m_sig_2_v_shower_angle);
        m_rdr_sig_2->AddVariable("sig_2_v_flag_single_shower",  &m_sig_2_v_flag_single_shower);
        m_rdr_sig_2->AddVariable("sig_2_v_medium_dQ_dx",        &m_sig_2_v_medium_dQ_dx);
        m_rdr_sig_2->AddVariable("sig_2_v_start_dQ_dx",         &m_sig_2_v_start_dQ_dx);
        m_rdr_sig_2->BookMVA("MyBDT", m_sig_2_xml);
    }

    if (!m_lol_1_xml.empty()) {
        m_rdr_lol_1 = std::make_unique<TMVA::Reader>();
        m_rdr_lol_1->AddVariable("lol_1_v_energy",      &m_lol_1_v_energy);
        m_rdr_lol_1->AddVariable("lol_1_v_vtx_n_segs",  &m_lol_1_v_vtx_n_segs);
        m_rdr_lol_1->AddVariable("lol_1_v_nseg",        &m_lol_1_v_nseg);
        m_rdr_lol_1->AddVariable("lol_1_v_angle",       &m_lol_1_v_angle);
        m_rdr_lol_1->BookMVA("MyBDT", m_lol_1_xml);
    }

    if (!m_lol_2_xml.empty()) {
        m_rdr_lol_2 = std::make_unique<TMVA::Reader>();
        m_rdr_lol_2->AddVariable("lol_2_v_length",             &m_lol_2_v_length);
        m_rdr_lol_2->AddVariable("lol_2_v_angle",              &m_lol_2_v_angle);
        m_rdr_lol_2->AddVariable("lol_2_v_type",               &m_lol_2_v_type);
        m_rdr_lol_2->AddVariable("lol_2_v_vtx_n_segs",         &m_lol_2_v_vtx_n_segs);
        m_rdr_lol_2->AddVariable("lol_2_v_energy",             &m_lol_2_v_energy);
        m_rdr_lol_2->AddVariable("lol_2_v_shower_main_length", &m_lol_2_v_shower_main_length);
        m_rdr_lol_2->AddVariable("lol_2_v_flag_dir_weak",      &m_lol_2_v_flag_dir_weak);
        m_rdr_lol_2->BookMVA("MyBDT", m_lol_2_xml);
    }

    if (!m_tro_1_xml.empty()) {
        m_rdr_tro_1 = std::make_unique<TMVA::Reader>();
        m_rdr_tro_1->AddVariable("tro_1_v_particle_type",       &m_tro_1_v_particle_type);
        m_rdr_tro_1->AddVariable("tro_1_v_flag_dir_weak",       &m_tro_1_v_flag_dir_weak);
        m_rdr_tro_1->AddVariable("tro_1_v_min_dis",             &m_tro_1_v_min_dis);
        m_rdr_tro_1->AddVariable("tro_1_v_sg1_length",          &m_tro_1_v_sg1_length);
        m_rdr_tro_1->AddVariable("tro_1_v_shower_main_length",  &m_tro_1_v_shower_main_length);
        m_rdr_tro_1->AddVariable("tro_1_v_max_n_vtx_segs",      &m_tro_1_v_max_n_vtx_segs);
        m_rdr_tro_1->AddVariable("tro_1_v_tmp_length",          &m_tro_1_v_tmp_length);
        m_rdr_tro_1->AddVariable("tro_1_v_medium_dQ_dx",        &m_tro_1_v_medium_dQ_dx);
        m_rdr_tro_1->AddVariable("tro_1_v_dQ_dx_cut",           &m_tro_1_v_dQ_dx_cut);
        m_rdr_tro_1->AddVariable("tro_1_v_flag_shower_topology",&m_tro_1_v_flag_shower_topology);
        m_rdr_tro_1->BookMVA("MyBDT", m_tro_1_xml);
    }

    if (!m_tro_2_xml.empty()) {
        m_rdr_tro_2 = std::make_unique<TMVA::Reader>();
        m_rdr_tro_2->AddVariable("tro_2_v_energy",      &m_tro_2_v_energy);
        m_rdr_tro_2->AddVariable("tro_2_v_stem_length", &m_tro_2_v_stem_length);
        m_rdr_tro_2->AddVariable("tro_2_v_iso_angle",   &m_tro_2_v_iso_angle);
        m_rdr_tro_2->AddVariable("tro_2_v_max_length",  &m_tro_2_v_max_length);
        m_rdr_tro_2->AddVariable("tro_2_v_angle",       &m_tro_2_v_angle);
        m_rdr_tro_2->BookMVA("MyBDT", m_tro_2_xml);
    }

    if (!m_tro_4_xml.empty()) {
        m_rdr_tro_4 = std::make_unique<TMVA::Reader>();
        m_rdr_tro_4->AddVariable("tro_4_v_dir2_mag",              &m_tro_4_v_dir2_mag);
        m_rdr_tro_4->AddVariable("tro_4_v_angle",                 &m_tro_4_v_angle);
        m_rdr_tro_4->AddVariable("tro_4_v_angle1",                &m_tro_4_v_angle1);
        m_rdr_tro_4->AddVariable("tro_4_v_angle2",                &m_tro_4_v_angle2);
        m_rdr_tro_4->AddVariable("tro_4_v_length",                &m_tro_4_v_length);
        m_rdr_tro_4->AddVariable("tro_4_v_length1",               &m_tro_4_v_length1);
        m_rdr_tro_4->AddVariable("tro_4_v_medium_dQ_dx",          &m_tro_4_v_medium_dQ_dx);
        m_rdr_tro_4->AddVariable("tro_4_v_end_dQ_dx",             &m_tro_4_v_end_dQ_dx);
        m_rdr_tro_4->AddVariable("tro_4_v_energy",                &m_tro_4_v_energy);
        m_rdr_tro_4->AddVariable("tro_4_v_shower_main_length",    &m_tro_4_v_shower_main_length);
        m_rdr_tro_4->AddVariable("tro_4_v_flag_shower_trajectory",&m_tro_4_v_flag_shower_trajectory);
        m_rdr_tro_4->BookMVA("MyBDT", m_tro_4_xml);
    }

    if (!m_tro_5_xml.empty()) {
        m_rdr_tro_5 = std::make_unique<TMVA::Reader>();
        m_rdr_tro_5->AddVariable("tro_5_v_max_angle",   &m_tro_5_v_max_angle);
        m_rdr_tro_5->AddVariable("tro_5_v_min_angle",   &m_tro_5_v_min_angle);
        m_rdr_tro_5->AddVariable("tro_5_v_max_length",  &m_tro_5_v_max_length);
        m_rdr_tro_5->AddVariable("tro_5_v_iso_angle",   &m_tro_5_v_iso_angle);
        m_rdr_tro_5->AddVariable("tro_5_v_n_vtx_segs", &m_tro_5_v_n_vtx_segs);
        m_rdr_tro_5->AddVariable("tro_5_v_min_count",   &m_tro_5_v_min_count);
        m_rdr_tro_5->AddVariable("tro_5_v_max_count",   &m_tro_5_v_max_count);
        m_rdr_tro_5->AddVariable("tro_5_v_energy",      &m_tro_5_v_energy);
        m_rdr_tro_5->BookMVA("MyBDT", m_tro_5_xml);
    }

    // ---- Final XGBoost reader ----

    if (!m_nue_xgboost_xml.empty()) {
        m_rdr_xgboost = std::make_unique<TMVA::Reader>();
        m_rdr_xgboost->AddVariable("match_isFC",                   &m_xgb_match_isFC);
        m_rdr_xgboost->AddVariable("kine_reco_Enu",                &m_xgb_kine_reco_Enu);
        m_rdr_xgboost->AddVariable("cme_mu_energy",                &m_xgb_cme_mu_energy);
        m_rdr_xgboost->AddVariable("cme_energy",                   &m_xgb_cme_energy);
        m_rdr_xgboost->AddVariable("cme_mu_length",                &m_xgb_cme_mu_length);
        m_rdr_xgboost->AddVariable("cme_length",                   &m_xgb_cme_length);
        m_rdr_xgboost->AddVariable("cme_angle_beam",               &m_xgb_cme_angle_beam);
        m_rdr_xgboost->AddVariable("anc_angle",                    &m_xgb_anc_angle);
        m_rdr_xgboost->AddVariable("anc_max_angle",                &m_xgb_anc_max_angle);
        m_rdr_xgboost->AddVariable("anc_max_length",               &m_xgb_anc_max_length);
        m_rdr_xgboost->AddVariable("anc_acc_forward_length",       &m_xgb_anc_acc_forward_length);
        m_rdr_xgboost->AddVariable("anc_acc_backward_length",      &m_xgb_anc_acc_backward_length);
        m_rdr_xgboost->AddVariable("anc_acc_forward_length1",      &m_xgb_anc_acc_forward_length1);
        m_rdr_xgboost->AddVariable("anc_shower_main_length",       &m_xgb_anc_shower_main_length);
        m_rdr_xgboost->AddVariable("anc_shower_total_length",      &m_xgb_anc_shower_total_length);
        m_rdr_xgboost->AddVariable("anc_flag_main_outside",        &m_xgb_anc_flag_main_outside);
        m_rdr_xgboost->AddVariable("gap_flag_prolong_u",           &m_xgb_gap_flag_prolong_u);
        m_rdr_xgboost->AddVariable("gap_flag_prolong_v",           &m_xgb_gap_flag_prolong_v);
        m_rdr_xgboost->AddVariable("gap_flag_prolong_w",           &m_xgb_gap_flag_prolong_w);
        m_rdr_xgboost->AddVariable("gap_flag_parallel",            &m_xgb_gap_flag_parallel);
        m_rdr_xgboost->AddVariable("gap_n_points",                 &m_xgb_gap_n_points);
        m_rdr_xgboost->AddVariable("gap_n_bad",                    &m_xgb_gap_n_bad);
        m_rdr_xgboost->AddVariable("gap_energy",                   &m_xgb_gap_energy);
        m_rdr_xgboost->AddVariable("gap_num_valid_tracks",         &m_xgb_gap_num_valid_tracks);
        m_rdr_xgboost->AddVariable("gap_flag_single_shower",       &m_xgb_gap_flag_single_shower);
        m_rdr_xgboost->AddVariable("hol_1_n_valid_tracks",         &m_xgb_hol_1_n_valid_tracks);
        m_rdr_xgboost->AddVariable("hol_1_min_angle",              &m_xgb_hol_1_min_angle);
        m_rdr_xgboost->AddVariable("hol_1_energy",                 &m_xgb_hol_1_energy);
        m_rdr_xgboost->AddVariable("hol_1_min_length",             &m_xgb_hol_1_min_length);
        m_rdr_xgboost->AddVariable("hol_2_min_angle",              &m_xgb_hol_2_min_angle);
        m_rdr_xgboost->AddVariable("hol_2_medium_dQ_dx",           &m_xgb_hol_2_medium_dQ_dx);
        m_rdr_xgboost->AddVariable("hol_2_ncount",                 &m_xgb_hol_2_ncount);
        m_rdr_xgboost->AddVariable("lol_3_angle_beam",             &m_xgb_lol_3_angle_beam);
        m_rdr_xgboost->AddVariable("lol_3_n_valid_tracks",         &m_xgb_lol_3_n_valid_tracks);
        m_rdr_xgboost->AddVariable("lol_3_min_angle",              &m_xgb_lol_3_min_angle);
        m_rdr_xgboost->AddVariable("lol_3_vtx_n_segs",             &m_xgb_lol_3_vtx_n_segs);
        m_rdr_xgboost->AddVariable("lol_3_shower_main_length",     &m_xgb_lol_3_shower_main_length);
        m_rdr_xgboost->AddVariable("lol_3_n_out",                  &m_xgb_lol_3_n_out);
        m_rdr_xgboost->AddVariable("lol_3_n_sum",                  &m_xgb_lol_3_n_sum);
        m_rdr_xgboost->AddVariable("hol_1_flag_all_shower",        &m_xgb_hol_1_flag_all_shower);
        m_rdr_xgboost->AddVariable("mgo_energy",                   &m_xgb_mgo_energy);
        m_rdr_xgboost->AddVariable("mgo_max_energy",               &m_xgb_mgo_max_energy);
        m_rdr_xgboost->AddVariable("mgo_total_energy",             &m_xgb_mgo_total_energy);
        m_rdr_xgboost->AddVariable("mgo_n_showers",                &m_xgb_mgo_n_showers);
        m_rdr_xgboost->AddVariable("mgo_max_energy_1",             &m_xgb_mgo_max_energy_1);
        m_rdr_xgboost->AddVariable("mgo_max_energy_2",             &m_xgb_mgo_max_energy_2);
        m_rdr_xgboost->AddVariable("mgo_total_other_energy",       &m_xgb_mgo_total_other_energy);
        m_rdr_xgboost->AddVariable("mgo_n_total_showers",          &m_xgb_mgo_n_total_showers);
        m_rdr_xgboost->AddVariable("mgo_total_other_energy_1",     &m_xgb_mgo_total_other_energy_1);
        m_rdr_xgboost->AddVariable("mgt_flag_single_shower",       &m_xgb_mgt_flag_single_shower);
        m_rdr_xgboost->AddVariable("mgt_max_energy",               &m_xgb_mgt_max_energy);
        m_rdr_xgboost->AddVariable("mgt_total_other_energy",       &m_xgb_mgt_total_other_energy);
        m_rdr_xgboost->AddVariable("mgt_max_energy_1",             &m_xgb_mgt_max_energy_1);
        m_rdr_xgboost->AddVariable("mgt_e_indirect_max_energy",    &m_xgb_mgt_e_indirect_max_energy);
        m_rdr_xgboost->AddVariable("mgt_e_direct_max_energy",      &m_xgb_mgt_e_direct_max_energy);
        m_rdr_xgboost->AddVariable("mgt_n_direct_showers",         &m_xgb_mgt_n_direct_showers);
        m_rdr_xgboost->AddVariable("mgt_e_direct_total_energy",    &m_xgb_mgt_e_direct_total_energy);
        m_rdr_xgboost->AddVariable("mgt_flag_indirect_max_pio",    &m_xgb_mgt_flag_indirect_max_pio);
        m_rdr_xgboost->AddVariable("mgt_e_indirect_total_energy",  &m_xgb_mgt_e_indirect_total_energy);
        m_rdr_xgboost->AddVariable("mip_quality_energy",           &m_xgb_mip_quality_energy);
        m_rdr_xgboost->AddVariable("mip_quality_overlap",          &m_xgb_mip_quality_overlap);
        m_rdr_xgboost->AddVariable("mip_quality_n_showers",        &m_xgb_mip_quality_n_showers);
        m_rdr_xgboost->AddVariable("mip_quality_n_tracks",         &m_xgb_mip_quality_n_tracks);
        m_rdr_xgboost->AddVariable("mip_quality_flag_inside_pi0",  &m_xgb_mip_quality_flag_inside_pi0);
        m_rdr_xgboost->AddVariable("mip_quality_n_pi0_showers",    &m_xgb_mip_quality_n_pi0_showers);
        m_rdr_xgboost->AddVariable("mip_quality_shortest_length",  &m_xgb_mip_quality_shortest_length);
        m_rdr_xgboost->AddVariable("mip_quality_acc_length",       &m_xgb_mip_quality_acc_length);
        m_rdr_xgboost->AddVariable("mip_quality_shortest_angle",   &m_xgb_mip_quality_shortest_angle);
        m_rdr_xgboost->AddVariable("mip_quality_flag_proton",      &m_xgb_mip_quality_flag_proton);
        m_rdr_xgboost->AddVariable("br1_1_shower_type",            &m_xgb_br1_1_shower_type);
        m_rdr_xgboost->AddVariable("br1_1_vtx_n_segs",             &m_xgb_br1_1_vtx_n_segs);
        m_rdr_xgboost->AddVariable("br1_1_energy",                 &m_xgb_br1_1_energy);
        m_rdr_xgboost->AddVariable("br1_1_n_segs",                 &m_xgb_br1_1_n_segs);
        m_rdr_xgboost->AddVariable("br1_1_flag_sg_topology",       &m_xgb_br1_1_flag_sg_topology);
        m_rdr_xgboost->AddVariable("br1_1_flag_sg_trajectory",     &m_xgb_br1_1_flag_sg_trajectory);
        m_rdr_xgboost->AddVariable("br1_1_sg_length",              &m_xgb_br1_1_sg_length);
        m_rdr_xgboost->AddVariable("br1_2_n_connected",            &m_xgb_br1_2_n_connected);
        m_rdr_xgboost->AddVariable("br1_2_max_length",             &m_xgb_br1_2_max_length);
        m_rdr_xgboost->AddVariable("br1_2_n_connected_1",          &m_xgb_br1_2_n_connected_1);
        m_rdr_xgboost->AddVariable("br1_2_n_shower_segs",          &m_xgb_br1_2_n_shower_segs);
        m_rdr_xgboost->AddVariable("br1_2_max_length_ratio",       &m_xgb_br1_2_max_length_ratio);
        m_rdr_xgboost->AddVariable("br1_2_shower_length",          &m_xgb_br1_2_shower_length);
        m_rdr_xgboost->AddVariable("br1_3_n_connected_p",          &m_xgb_br1_3_n_connected_p);
        m_rdr_xgboost->AddVariable("br1_3_max_length_p",           &m_xgb_br1_3_max_length_p);
        m_rdr_xgboost->AddVariable("br1_3_n_shower_main_segs",     &m_xgb_br1_3_n_shower_main_segs);
        m_rdr_xgboost->AddVariable("br3_1_energy",                 &m_xgb_br3_1_energy);
        m_rdr_xgboost->AddVariable("br3_1_n_shower_segments",      &m_xgb_br3_1_n_shower_segments);
        m_rdr_xgboost->AddVariable("br3_1_sg_flag_trajectory",     &m_xgb_br3_1_sg_flag_trajectory);
        m_rdr_xgboost->AddVariable("br3_1_sg_direct_length",       &m_xgb_br3_1_sg_direct_length);
        m_rdr_xgboost->AddVariable("br3_1_sg_length",              &m_xgb_br3_1_sg_length);
        m_rdr_xgboost->AddVariable("br3_1_total_main_length",      &m_xgb_br3_1_total_main_length);
        m_rdr_xgboost->AddVariable("br3_1_total_length",           &m_xgb_br3_1_total_length);
        m_rdr_xgboost->AddVariable("br3_1_iso_angle",              &m_xgb_br3_1_iso_angle);
        m_rdr_xgboost->AddVariable("br3_1_sg_flag_topology",       &m_xgb_br3_1_sg_flag_topology);
        m_rdr_xgboost->AddVariable("br3_2_n_ele",                  &m_xgb_br3_2_n_ele);
        m_rdr_xgboost->AddVariable("br3_2_n_other",                &m_xgb_br3_2_n_other);
        m_rdr_xgboost->AddVariable("br3_2_other_fid",              &m_xgb_br3_2_other_fid);
        m_rdr_xgboost->AddVariable("br3_4_acc_length",             &m_xgb_br3_4_acc_length);
        m_rdr_xgboost->AddVariable("br3_4_total_length",           &m_xgb_br3_4_total_length);
        m_rdr_xgboost->AddVariable("br3_7_min_angle",              &m_xgb_br3_7_min_angle);
        m_rdr_xgboost->AddVariable("br3_8_max_dQ_dx",              &m_xgb_br3_8_max_dQ_dx);
        m_rdr_xgboost->AddVariable("br3_8_n_main_segs",            &m_xgb_br3_8_n_main_segs);
        m_rdr_xgboost->AddVariable("vis_1_n_vtx_segs",             &m_xgb_vis_1_n_vtx_segs);
        m_rdr_xgboost->AddVariable("vis_1_energy",                 &m_xgb_vis_1_energy);
        m_rdr_xgboost->AddVariable("vis_1_num_good_tracks",        &m_xgb_vis_1_num_good_tracks);
        m_rdr_xgboost->AddVariable("vis_1_max_angle",              &m_xgb_vis_1_max_angle);
        m_rdr_xgboost->AddVariable("vis_1_max_shower_angle",       &m_xgb_vis_1_max_shower_angle);
        m_rdr_xgboost->AddVariable("vis_1_tmp_length1",            &m_xgb_vis_1_tmp_length1);
        m_rdr_xgboost->AddVariable("vis_1_tmp_length2",            &m_xgb_vis_1_tmp_length2);
        m_rdr_xgboost->AddVariable("vis_2_n_vtx_segs",             &m_xgb_vis_2_n_vtx_segs);
        m_rdr_xgboost->AddVariable("vis_2_min_angle",              &m_xgb_vis_2_min_angle);
        m_rdr_xgboost->AddVariable("vis_2_min_weak_track",         &m_xgb_vis_2_min_weak_track);
        m_rdr_xgboost->AddVariable("vis_2_angle_beam",             &m_xgb_vis_2_angle_beam);
        m_rdr_xgboost->AddVariable("vis_2_min_angle1",             &m_xgb_vis_2_min_angle1);
        m_rdr_xgboost->AddVariable("vis_2_iso_angle1",             &m_xgb_vis_2_iso_angle1);
        m_rdr_xgboost->AddVariable("vis_2_min_medium_dQ_dx",       &m_xgb_vis_2_min_medium_dQ_dx);
        m_rdr_xgboost->AddVariable("vis_2_min_length",             &m_xgb_vis_2_min_length);
        m_rdr_xgboost->AddVariable("vis_2_sg_length",              &m_xgb_vis_2_sg_length);
        m_rdr_xgboost->AddVariable("vis_2_max_angle",              &m_xgb_vis_2_max_angle);
        m_rdr_xgboost->AddVariable("vis_2_max_weak_track",         &m_xgb_vis_2_max_weak_track);
        m_rdr_xgboost->AddVariable("pio_1_mass",                   &m_xgb_pio_1_mass);
        m_rdr_xgboost->AddVariable("pio_1_pio_type",               &m_xgb_pio_1_pio_type);
        m_rdr_xgboost->AddVariable("pio_1_energy_1",               &m_xgb_pio_1_energy_1);
        m_rdr_xgboost->AddVariable("pio_1_energy_2",               &m_xgb_pio_1_energy_2);
        m_rdr_xgboost->AddVariable("pio_1_dis_1",                  &m_xgb_pio_1_dis_1);
        m_rdr_xgboost->AddVariable("pio_1_dis_2",                  &m_xgb_pio_1_dis_2);
        m_rdr_xgboost->AddVariable("pio_mip_id",                   &m_xgb_pio_mip_id);
        m_rdr_xgboost->AddVariable("stem_dir_flag_single_shower",  &m_xgb_stem_dir_flag_single_shower);
        m_rdr_xgboost->AddVariable("stem_dir_angle",               &m_xgb_stem_dir_angle);
        m_rdr_xgboost->AddVariable("stem_dir_energy",              &m_xgb_stem_dir_energy);
        m_rdr_xgboost->AddVariable("stem_dir_angle1",              &m_xgb_stem_dir_angle1);
        m_rdr_xgboost->AddVariable("stem_dir_angle2",              &m_xgb_stem_dir_angle2);
        m_rdr_xgboost->AddVariable("stem_dir_angle3",              &m_xgb_stem_dir_angle3);
        m_rdr_xgboost->AddVariable("stem_dir_ratio",               &m_xgb_stem_dir_ratio);
        m_rdr_xgboost->AddVariable("br2_num_valid_tracks",         &m_xgb_br2_num_valid_tracks);
        m_rdr_xgboost->AddVariable("br2_n_shower_main_segs",       &m_xgb_br2_n_shower_main_segs);
        m_rdr_xgboost->AddVariable("br2_max_angle",                &m_xgb_br2_max_angle);
        m_rdr_xgboost->AddVariable("br2_sg_length",                &m_xgb_br2_sg_length);
        m_rdr_xgboost->AddVariable("br2_flag_sg_trajectory",       &m_xgb_br2_flag_sg_trajectory);
        m_rdr_xgboost->AddVariable("stem_len_energy",              &m_xgb_stem_len_energy);
        m_rdr_xgboost->AddVariable("stem_len_length",              &m_xgb_stem_len_length);
        m_rdr_xgboost->AddVariable("stem_len_flag_avoid_muon_check",&m_xgb_stem_len_flag_avoid_muon_check);
        m_rdr_xgboost->AddVariable("stem_len_num_daughters",       &m_xgb_stem_len_num_daughters);
        m_rdr_xgboost->AddVariable("stem_len_daughter_length",     &m_xgb_stem_len_daughter_length);
        m_rdr_xgboost->AddVariable("brm_n_mu_segs",                &m_xgb_brm_n_mu_segs);
        m_rdr_xgboost->AddVariable("brm_Ep",                       &m_xgb_brm_Ep);
        m_rdr_xgboost->AddVariable("brm_acc_length",               &m_xgb_brm_acc_length);
        m_rdr_xgboost->AddVariable("brm_shower_total_length",      &m_xgb_brm_shower_total_length);
        m_rdr_xgboost->AddVariable("brm_connected_length",         &m_xgb_brm_connected_length);
        m_rdr_xgboost->AddVariable("brm_n_size",                   &m_xgb_brm_n_size);
        m_rdr_xgboost->AddVariable("brm_n_shower_main_segs",       &m_xgb_brm_n_shower_main_segs);
        m_rdr_xgboost->AddVariable("brm_n_mu_main",                &m_xgb_brm_n_mu_main);
        m_rdr_xgboost->AddVariable("lem_shower_main_length",       &m_xgb_lem_shower_main_length);
        m_rdr_xgboost->AddVariable("lem_n_3seg",                   &m_xgb_lem_n_3seg);
        m_rdr_xgboost->AddVariable("lem_e_charge",                 &m_xgb_lem_e_charge);
        m_rdr_xgboost->AddVariable("lem_e_dQdx",                   &m_xgb_lem_e_dQdx);
        m_rdr_xgboost->AddVariable("lem_shower_num_main_segs",     &m_xgb_lem_shower_num_main_segs);
        m_rdr_xgboost->AddVariable("brm_acc_direct_length",        &m_xgb_brm_acc_direct_length);
        m_rdr_xgboost->AddVariable("stw_1_energy",                 &m_xgb_stw_1_energy);
        m_rdr_xgboost->AddVariable("stw_1_dis",                    &m_xgb_stw_1_dis);
        m_rdr_xgboost->AddVariable("stw_1_dQ_dx",                  &m_xgb_stw_1_dQ_dx);
        m_rdr_xgboost->AddVariable("stw_1_flag_single_shower",     &m_xgb_stw_1_flag_single_shower);
        m_rdr_xgboost->AddVariable("stw_1_n_pi0",                  &m_xgb_stw_1_n_pi0);
        m_rdr_xgboost->AddVariable("stw_1_num_valid_tracks",       &m_xgb_stw_1_num_valid_tracks);
        m_rdr_xgboost->AddVariable("spt_shower_main_length",       &m_xgb_spt_shower_main_length);
        m_rdr_xgboost->AddVariable("spt_shower_total_length",      &m_xgb_spt_shower_total_length);
        m_rdr_xgboost->AddVariable("spt_angle_beam",               &m_xgb_spt_angle_beam);
        m_rdr_xgboost->AddVariable("spt_angle_vertical",           &m_xgb_spt_angle_vertical);
        m_rdr_xgboost->AddVariable("spt_max_dQ_dx",                &m_xgb_spt_max_dQ_dx);
        m_rdr_xgboost->AddVariable("spt_angle_beam_1",             &m_xgb_spt_angle_beam_1);
        m_rdr_xgboost->AddVariable("spt_angle_drift",              &m_xgb_spt_angle_drift);
        m_rdr_xgboost->AddVariable("spt_angle_drift_1",            &m_xgb_spt_angle_drift_1);
        m_rdr_xgboost->AddVariable("spt_num_valid_tracks",         &m_xgb_spt_num_valid_tracks);
        m_rdr_xgboost->AddVariable("spt_n_vtx_segs",               &m_xgb_spt_n_vtx_segs);
        m_rdr_xgboost->AddVariable("spt_max_length",               &m_xgb_spt_max_length);
        m_rdr_xgboost->AddVariable("mip_energy",                   &m_xgb_mip_energy);
        m_rdr_xgboost->AddVariable("mip_n_end_reduction",          &m_xgb_mip_n_end_reduction);
        m_rdr_xgboost->AddVariable("mip_n_first_mip",              &m_xgb_mip_n_first_mip);
        m_rdr_xgboost->AddVariable("mip_n_first_non_mip",          &m_xgb_mip_n_first_non_mip);
        m_rdr_xgboost->AddVariable("mip_n_first_non_mip_1",        &m_xgb_mip_n_first_non_mip_1);
        m_rdr_xgboost->AddVariable("mip_n_first_non_mip_2",        &m_xgb_mip_n_first_non_mip_2);
        m_rdr_xgboost->AddVariable("mip_vec_dQ_dx_0",              &m_xgb_mip_vec_dQ_dx_0);
        m_rdr_xgboost->AddVariable("mip_vec_dQ_dx_1",              &m_xgb_mip_vec_dQ_dx_1);
        m_rdr_xgboost->AddVariable("mip_max_dQ_dx_sample",         &m_xgb_mip_max_dQ_dx_sample);
        m_rdr_xgboost->AddVariable("mip_n_below_threshold",        &m_xgb_mip_n_below_threshold);
        m_rdr_xgboost->AddVariable("mip_n_below_zero",             &m_xgb_mip_n_below_zero);
        m_rdr_xgboost->AddVariable("mip_n_lowest",                 &m_xgb_mip_n_lowest);
        m_rdr_xgboost->AddVariable("mip_n_highest",                &m_xgb_mip_n_highest);
        m_rdr_xgboost->AddVariable("mip_lowest_dQ_dx",             &m_xgb_mip_lowest_dQ_dx);
        m_rdr_xgboost->AddVariable("mip_highest_dQ_dx",            &m_xgb_mip_highest_dQ_dx);
        m_rdr_xgboost->AddVariable("mip_medium_dQ_dx",             &m_xgb_mip_medium_dQ_dx);
        m_rdr_xgboost->AddVariable("mip_stem_length",              &m_xgb_mip_stem_length);
        m_rdr_xgboost->AddVariable("mip_length_main",              &m_xgb_mip_length_main);
        m_rdr_xgboost->AddVariable("mip_length_total",             &m_xgb_mip_length_total);
        m_rdr_xgboost->AddVariable("mip_angle_beam",               &m_xgb_mip_angle_beam);
        m_rdr_xgboost->AddVariable("mip_iso_angle",                &m_xgb_mip_iso_angle);
        m_rdr_xgboost->AddVariable("mip_n_vertex",                 &m_xgb_mip_n_vertex);
        m_rdr_xgboost->AddVariable("mip_n_good_tracks",            &m_xgb_mip_n_good_tracks);
        m_rdr_xgboost->AddVariable("mip_E_indirect_max_energy",    &m_xgb_mip_E_indirect_max_energy);
        m_rdr_xgboost->AddVariable("mip_flag_all_above",           &m_xgb_mip_flag_all_above);
        m_rdr_xgboost->AddVariable("mip_min_dQ_dx_5",              &m_xgb_mip_min_dQ_dx_5);
        m_rdr_xgboost->AddVariable("mip_n_other_vertex",           &m_xgb_mip_n_other_vertex);
        m_rdr_xgboost->AddVariable("mip_n_stem_size",              &m_xgb_mip_n_stem_size);
        m_rdr_xgboost->AddVariable("mip_flag_stem_trajectory",     &m_xgb_mip_flag_stem_trajectory);
        m_rdr_xgboost->AddVariable("mip_min_dis",                  &m_xgb_mip_min_dis);
        m_rdr_xgboost->AddVariable("mip_vec_dQ_dx_2",              &m_xgb_mip_vec_dQ_dx_2);
        m_rdr_xgboost->AddVariable("mip_vec_dQ_dx_3",              &m_xgb_mip_vec_dQ_dx_3);
        m_rdr_xgboost->AddVariable("mip_vec_dQ_dx_4",              &m_xgb_mip_vec_dQ_dx_4);
        m_rdr_xgboost->AddVariable("mip_vec_dQ_dx_5",              &m_xgb_mip_vec_dQ_dx_5);
        m_rdr_xgboost->AddVariable("mip_vec_dQ_dx_6",              &m_xgb_mip_vec_dQ_dx_6);
        m_rdr_xgboost->AddVariable("mip_vec_dQ_dx_7",              &m_xgb_mip_vec_dQ_dx_7);
        m_rdr_xgboost->AddVariable("mip_vec_dQ_dx_8",              &m_xgb_mip_vec_dQ_dx_8);
        m_rdr_xgboost->AddVariable("mip_vec_dQ_dx_9",              &m_xgb_mip_vec_dQ_dx_9);
        m_rdr_xgboost->AddVariable("mip_vec_dQ_dx_10",             &m_xgb_mip_vec_dQ_dx_10);
        m_rdr_xgboost->AddVariable("mip_vec_dQ_dx_11",             &m_xgb_mip_vec_dQ_dx_11);
        m_rdr_xgboost->AddVariable("mip_vec_dQ_dx_12",             &m_xgb_mip_vec_dQ_dx_12);
        m_rdr_xgboost->AddVariable("mip_vec_dQ_dx_13",             &m_xgb_mip_vec_dQ_dx_13);
        m_rdr_xgboost->AddVariable("mip_vec_dQ_dx_14",             &m_xgb_mip_vec_dQ_dx_14);
        m_rdr_xgboost->AddVariable("mip_vec_dQ_dx_15",             &m_xgb_mip_vec_dQ_dx_15);
        m_rdr_xgboost->AddVariable("mip_vec_dQ_dx_16",             &m_xgb_mip_vec_dQ_dx_16);
        m_rdr_xgboost->AddVariable("mip_vec_dQ_dx_17",             &m_xgb_mip_vec_dQ_dx_17);
        m_rdr_xgboost->AddVariable("mip_vec_dQ_dx_18",             &m_xgb_mip_vec_dQ_dx_18);
        m_rdr_xgboost->AddVariable("mip_vec_dQ_dx_19",             &m_xgb_mip_vec_dQ_dx_19);
        m_rdr_xgboost->AddVariable("br3_3_score",                  &m_xgb_br3_3_score);
        m_rdr_xgboost->AddVariable("br3_5_score",                  &m_xgb_br3_5_score);
        m_rdr_xgboost->AddVariable("br3_6_score",                  &m_xgb_br3_6_score);
        m_rdr_xgboost->AddVariable("pio_2_score",                  &m_xgb_pio_2_score);
        m_rdr_xgboost->AddVariable("stw_2_score",                  &m_xgb_stw_2_score);
        m_rdr_xgboost->AddVariable("stw_3_score",                  &m_xgb_stw_3_score);
        m_rdr_xgboost->AddVariable("stw_4_score",                  &m_xgb_stw_4_score);
        m_rdr_xgboost->AddVariable("sig_1_score",                  &m_xgb_sig_1_score);
        m_rdr_xgboost->AddVariable("sig_2_score",                  &m_xgb_sig_2_score);
        m_rdr_xgboost->AddVariable("lol_1_score",                  &m_xgb_lol_1_score);
        m_rdr_xgboost->AddVariable("lol_2_score",                  &m_xgb_lol_2_score);
        m_rdr_xgboost->AddVariable("tro_1_score",                  &m_xgb_tro_1_score);
        m_rdr_xgboost->AddVariable("tro_2_score",                  &m_xgb_tro_2_score);
        m_rdr_xgboost->AddVariable("tro_4_score",                  &m_xgb_tro_4_score);
        m_rdr_xgboost->AddVariable("tro_5_score",                  &m_xgb_tro_5_score);
        m_rdr_xgboost->AddVariable("br4_1_shower_main_length",     &m_xgb_br4_1_shower_main_length);
        m_rdr_xgboost->AddVariable("br4_1_shower_total_length",    &m_xgb_br4_1_shower_total_length);
        m_rdr_xgboost->AddVariable("br4_1_min_dis",                &m_xgb_br4_1_min_dis);
        m_rdr_xgboost->AddVariable("br4_1_energy",                 &m_xgb_br4_1_energy);
        m_rdr_xgboost->AddVariable("br4_1_flag_avoid_muon_check",  &m_xgb_br4_1_flag_avoid_muon_check);
        m_rdr_xgboost->AddVariable("br4_1_n_vtx_segs",             &m_xgb_br4_1_n_vtx_segs);
        m_rdr_xgboost->AddVariable("br4_2_ratio_45",               &m_xgb_br4_2_ratio_45);
        m_rdr_xgboost->AddVariable("br4_2_ratio_35",               &m_xgb_br4_2_ratio_35);
        m_rdr_xgboost->AddVariable("br4_2_ratio_25",               &m_xgb_br4_2_ratio_25);
        m_rdr_xgboost->AddVariable("br4_2_ratio_15",               &m_xgb_br4_2_ratio_15);
        m_rdr_xgboost->AddVariable("br4_2_ratio1_45",              &m_xgb_br4_2_ratio1_45);
        m_rdr_xgboost->AddVariable("br4_2_ratio1_35",              &m_xgb_br4_2_ratio1_35);
        m_rdr_xgboost->AddVariable("br4_2_ratio1_25",              &m_xgb_br4_2_ratio1_25);
        m_rdr_xgboost->AddVariable("br4_2_ratio1_15",              &m_xgb_br4_2_ratio1_15);
        m_rdr_xgboost->AddVariable("br4_2_iso_angle",              &m_xgb_br4_2_iso_angle);
        m_rdr_xgboost->AddVariable("br4_2_iso_angle1",             &m_xgb_br4_2_iso_angle1);
        m_rdr_xgboost->AddVariable("br4_2_angle",                  &m_xgb_br4_2_angle);
        m_rdr_xgboost->AddVariable("tro_3_stem_length",            &m_xgb_tro_3_stem_length);
        m_rdr_xgboost->AddVariable("tro_3_n_muon_segs",            &m_xgb_tro_3_n_muon_segs);
        m_rdr_xgboost->AddVariable("br4_1_n_main_segs",            &m_xgb_br4_1_n_main_segs);
        m_rdr_xgboost->BookMVA("MyBDT", m_nue_xgboost_xml);
    }
}

void UbooneNueBDTScorer::visit(Clus::Facade::Ensemble& ensemble) const
{
    if (!m_rdr_xgboost) {
        log->warn("UbooneNueBDTScorer: nue_xgboost_xml not set or file not found — skipping nue BDT scoring");
        return;
    }

    auto groupings = ensemble.with_name(m_grouping_name);
    if (groupings.empty()) {
        log->debug("UbooneNueBDTScorer: no grouping '{}'", m_grouping_name);
        return;
    }

    auto& grouping = *groupings.at(0);
    auto tf = grouping.get_track_fitting();
    if (!tf) {
        log->warn("UbooneNueBDTScorer: no TrackFitting in grouping '{}'", m_grouping_name);
        return;
    }

    Clus::PR::TaggerInfo& ti  = tf->get_tagger_info_mutable();
    const Clus::PR::KineInfo& ki = tf->get_kine_info();

    cal_bdts_xgboost(ti, ki);
}

// ===========================================================================
// cal_mipid_bdt
//
// Scores MIP identification features (scalar, per-event).
// Gate: ti.mip_filled == 1.
//
// Prototype: NeutrinoID_nue_bdts.h::cal_mipid_bdt()
// ===========================================================================
float UbooneNueBDTScorer::cal_mipid_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    float val = default_val;

    TMVA::Reader reader;
    reader.AddVariable("mip_energy",             &ti.mip_energy);
    reader.AddVariable("mip_n_end_reduction",    &ti.mip_n_end_reduction);
    reader.AddVariable("mip_n_first_mip",        &ti.mip_n_first_mip);
    reader.AddVariable("mip_n_first_non_mip",    &ti.mip_n_first_non_mip);
    reader.AddVariable("mip_n_first_non_mip_1",  &ti.mip_n_first_non_mip_1);
    reader.AddVariable("mip_n_first_non_mip_2",  &ti.mip_n_first_non_mip_2);
    reader.AddVariable("mip_vec_dQ_dx_0",        &ti.mip_vec_dQ_dx_0);
    reader.AddVariable("mip_vec_dQ_dx_1",        &ti.mip_vec_dQ_dx_1);
    reader.AddVariable("mip_max_dQ_dx_sample",   &ti.mip_max_dQ_dx_sample);
    reader.AddVariable("mip_n_below_threshold",  &ti.mip_n_below_threshold);
    reader.AddVariable("mip_n_below_zero",       &ti.mip_n_below_zero);
    reader.AddVariable("mip_n_lowest",           &ti.mip_n_lowest);
    reader.AddVariable("mip_n_highest",          &ti.mip_n_highest);
    reader.AddVariable("mip_lowest_dQ_dx",       &ti.mip_lowest_dQ_dx);
    reader.AddVariable("mip_highest_dQ_dx",      &ti.mip_highest_dQ_dx);
    reader.AddVariable("mip_medium_dQ_dx",       &ti.mip_medium_dQ_dx);
    reader.AddVariable("mip_stem_length",        &ti.mip_stem_length);
    reader.AddVariable("mip_length_main",        &ti.mip_length_main);
    reader.AddVariable("mip_length_total",       &ti.mip_length_total);
    reader.AddVariable("mip_angle_beam",         &ti.mip_angle_beam);
    reader.AddVariable("mip_iso_angle",          &ti.mip_iso_angle);
    reader.AddVariable("mip_n_vertex",           &ti.mip_n_vertex);
    reader.AddVariable("mip_n_good_tracks",      &ti.mip_n_good_tracks);
    reader.AddVariable("mip_E_indirect_max_energy", &ti.mip_E_indirect_max_energy);
    reader.AddVariable("mip_flag_all_above",     &ti.mip_flag_all_above);
    reader.AddVariable("mip_min_dQ_dx_5",        &ti.mip_min_dQ_dx_5);
    reader.AddVariable("mip_n_other_vertex",     &ti.mip_n_other_vertex);
    reader.AddVariable("mip_n_stem_size",        &ti.mip_n_stem_size);
    reader.AddVariable("mip_flag_stem_trajectory",&ti.mip_flag_stem_trajectory);
    reader.AddVariable("mip_min_dis",            &ti.mip_min_dis);

    reader.BookMVA("MyBDT", m_mipid_xml);

    if (ti.mip_filled == 1)
        val = reader.EvaluateMVA("MyBDT");

    return val;
}

// ===========================================================================
// cal_gap_bdt
//
// Scores gap-identification features (scalar, per-event).
// Gate: ti.gap_filled == 1.
//
// Prototype: NeutrinoID_nue_bdts.h::cal_gap_bdt()
// ===========================================================================
float UbooneNueBDTScorer::cal_gap_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    float val = default_val;

    TMVA::Reader reader;
    reader.AddVariable("gap_flag_prolong_u",     &ti.gap_flag_prolong_u);
    reader.AddVariable("gap_flag_prolong_v",     &ti.gap_flag_prolong_v);
    reader.AddVariable("gap_flag_prolong_w",     &ti.gap_flag_prolong_w);
    reader.AddVariable("gap_flag_parallel",      &ti.gap_flag_parallel);
    reader.AddVariable("gap_n_points",           &ti.gap_n_points);
    reader.AddVariable("gap_n_bad",              &ti.gap_n_bad);
    reader.AddVariable("gap_energy",             &ti.gap_energy);
    reader.AddVariable("gap_num_valid_tracks",   &ti.gap_num_valid_tracks);
    reader.AddVariable("gap_flag_single_shower", &ti.gap_flag_single_shower);

    reader.BookMVA("MyBDT", m_gap_xml);

    if (ti.gap_filled == 1)
        val = reader.EvaluateMVA("MyBDT");

    return val;
}

// ===========================================================================
// cal_hol_lol_bdt
//
// Scores high-overlap-lol_3 combined features (scalar, per-event).
// Gate: ti.br_filled == 1.
//
// Prototype: NeutrinoID_nue_bdts.h::cal_hol_lol_bdt()
// ===========================================================================
float UbooneNueBDTScorer::cal_hol_lol_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    float val = default_val;

    TMVA::Reader reader;
    reader.AddVariable("hol_1_n_valid_tracks",     &ti.hol_1_n_valid_tracks);
    reader.AddVariable("hol_1_min_angle",          &ti.hol_1_min_angle);
    reader.AddVariable("hol_1_energy",             &ti.hol_1_energy);
    reader.AddVariable("hol_1_flag_all_shower",    &ti.hol_1_flag_all_shower);
    reader.AddVariable("hol_1_min_length",         &ti.hol_1_min_length);
    reader.AddVariable("hol_2_min_angle",          &ti.hol_2_min_angle);
    reader.AddVariable("hol_2_medium_dQ_dx",       &ti.hol_2_medium_dQ_dx);
    reader.AddVariable("hol_2_ncount",             &ti.hol_2_ncount);
    reader.AddVariable("lol_3_angle_beam",         &ti.lol_3_angle_beam);
    reader.AddVariable("lol_3_n_valid_tracks",     &ti.lol_3_n_valid_tracks);
    reader.AddVariable("lol_3_min_angle",          &ti.lol_3_min_angle);
    reader.AddVariable("lol_3_vtx_n_segs",         &ti.lol_3_vtx_n_segs);
    reader.AddVariable("lol_3_shower_main_length", &ti.lol_3_shower_main_length);
    reader.AddVariable("lol_3_n_out",              &ti.lol_3_n_out);
    reader.AddVariable("lol_3_n_sum",              &ti.lol_3_n_sum);

    reader.BookMVA("MyBDT", m_hol_lol_xml);

    if (ti.br_filled == 1)
        val = reader.EvaluateMVA("MyBDT");

    return val;
}

// ===========================================================================
// cal_cme_anc_bdt
//
// Scores compare-muon-energy + angular-cut combined features (scalar).
// Gate: ti.br_filled == 1.
//
// Prototype: NeutrinoID_nue_bdts.h::cal_cme_anc_bdt()
// ===========================================================================
float UbooneNueBDTScorer::cal_cme_anc_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    float val = default_val;

    TMVA::Reader reader;
    reader.AddVariable("cme_mu_energy",              &ti.cme_mu_energy);
    reader.AddVariable("cme_energy",                 &ti.cme_energy);
    reader.AddVariable("cme_mu_length",              &ti.cme_mu_length);
    reader.AddVariable("cme_length",                 &ti.cme_length);
    reader.AddVariable("cme_angle_beam",             &ti.cme_angle_beam);
    reader.AddVariable("anc_angle",                  &ti.anc_angle);
    reader.AddVariable("anc_max_angle",              &ti.anc_max_angle);
    reader.AddVariable("anc_max_length",             &ti.anc_max_length);
    reader.AddVariable("anc_acc_forward_length",     &ti.anc_acc_forward_length);
    reader.AddVariable("anc_acc_backward_length",    &ti.anc_acc_backward_length);
    reader.AddVariable("anc_acc_forward_length1",    &ti.anc_acc_forward_length1);
    reader.AddVariable("anc_shower_main_length",     &ti.anc_shower_main_length);
    reader.AddVariable("anc_shower_total_length",    &ti.anc_shower_total_length);
    reader.AddVariable("anc_flag_main_outside",      &ti.anc_flag_main_outside);

    reader.BookMVA("MyBDT", m_cme_anc_xml);

    if (ti.br_filled == 1)
        val = reader.EvaluateMVA("MyBDT");

    return val;
}

// ===========================================================================
// cal_mgo_mgt_bdt
//
// Scores multiple-gamma-other + multiple-gamma-track combined features.
// Gate: ti.br_filled == 1.
//
// Prototype: NeutrinoID_nue_bdts.h::cal_mgo_mgt_bdt()
// ===========================================================================
float UbooneNueBDTScorer::cal_mgo_mgt_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    float val = default_val;

    TMVA::Reader reader;
    reader.AddVariable("mgo_energy",                  &ti.mgo_energy);
    reader.AddVariable("mgo_max_energy",              &ti.mgo_max_energy);
    reader.AddVariable("mgo_total_energy",            &ti.mgo_total_energy);
    reader.AddVariable("mgo_n_showers",               &ti.mgo_n_showers);
    reader.AddVariable("mgo_max_energy_1",            &ti.mgo_max_energy_1);
    reader.AddVariable("mgo_max_energy_2",            &ti.mgo_max_energy_2);
    reader.AddVariable("mgo_total_other_energy",      &ti.mgo_total_other_energy);
    reader.AddVariable("mgo_n_total_showers",         &ti.mgo_n_total_showers);
    reader.AddVariable("mgo_total_other_energy_1",    &ti.mgo_total_other_energy_1);
    reader.AddVariable("mgt_flag_single_shower",      &ti.mgt_flag_single_shower);
    reader.AddVariable("mgt_max_energy",              &ti.mgt_max_energy);
    reader.AddVariable("mgt_total_other_energy",      &ti.mgt_total_other_energy);
    reader.AddVariable("mgt_max_energy_1",            &ti.mgt_max_energy_1);
    reader.AddVariable("mgt_e_indirect_max_energy",   &ti.mgt_e_indirect_max_energy);
    reader.AddVariable("mgt_e_direct_max_energy",     &ti.mgt_e_direct_max_energy);
    reader.AddVariable("mgt_n_direct_showers",        &ti.mgt_n_direct_showers);
    reader.AddVariable("mgt_e_direct_total_energy",   &ti.mgt_e_direct_total_energy);
    reader.AddVariable("mgt_flag_indirect_max_pio",   &ti.mgt_flag_indirect_max_pio);
    reader.AddVariable("mgt_e_indirect_total_energy", &ti.mgt_e_indirect_total_energy);

    reader.BookMVA("MyBDT", m_mgo_mgt_xml);

    if (ti.br_filled == 1)
        val = reader.EvaluateMVA("MyBDT");

    return val;
}

// ===========================================================================
// cal_br1_bdt
//
// Scores bad-reconstruction-1 features (3 sub-checks: br1_1, br1_2, br1_3).
// Gate: ti.br_filled == 1.
//
// Prototype: NeutrinoID_nue_bdts.h::cal_br1_bdt()
// ===========================================================================
float UbooneNueBDTScorer::cal_br1_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    float val = default_val;

    TMVA::Reader reader;
    reader.AddVariable("br1_1_shower_type",        &ti.br1_1_shower_type);
    reader.AddVariable("br1_1_vtx_n_segs",         &ti.br1_1_vtx_n_segs);
    reader.AddVariable("br1_1_energy",             &ti.br1_1_energy);
    reader.AddVariable("br1_1_n_segs",             &ti.br1_1_n_segs);
    reader.AddVariable("br1_1_flag_sg_topology",   &ti.br1_1_flag_sg_topology);
    reader.AddVariable("br1_1_flag_sg_trajectory", &ti.br1_1_flag_sg_trajectory);
    reader.AddVariable("br1_1_sg_length",          &ti.br1_1_sg_length);
    reader.AddVariable("br1_2_n_connected",        &ti.br1_2_n_connected);
    reader.AddVariable("br1_2_max_length",         &ti.br1_2_max_length);
    reader.AddVariable("br1_2_n_connected_1",      &ti.br1_2_n_connected_1);
    reader.AddVariable("br1_2_n_shower_segs",      &ti.br1_2_n_shower_segs);
    reader.AddVariable("br1_2_max_length_ratio",   &ti.br1_2_max_length_ratio);
    reader.AddVariable("br1_2_shower_length",      &ti.br1_2_shower_length);
    reader.AddVariable("br1_3_n_connected_p",      &ti.br1_3_n_connected_p);
    reader.AddVariable("br1_3_max_length_p",       &ti.br1_3_max_length_p);
    reader.AddVariable("br1_3_n_shower_main_segs", &ti.br1_3_n_shower_main_segs);

    reader.BookMVA("MyBDT", m_br1_xml);

    if (ti.br_filled == 1)
        val = reader.EvaluateMVA("MyBDT");

    return val;
}

// ===========================================================================
// cal_br3_bdt
//
// Scores bad-reconstruction-2 scalar features (br3_1,2,4,7,8 sub-checks).
// Gate: ti.br_filled == 1.
//
// Prototype: NeutrinoID_nue_bdts.h::cal_br3_bdt()
// ===========================================================================
float UbooneNueBDTScorer::cal_br3_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    float val = default_val;

    TMVA::Reader reader;
    reader.AddVariable("br3_1_energy",             &ti.br3_1_energy);
    reader.AddVariable("br3_1_n_shower_segments",  &ti.br3_1_n_shower_segments);
    reader.AddVariable("br3_1_sg_flag_trajectory", &ti.br3_1_sg_flag_trajectory);
    reader.AddVariable("br3_1_sg_direct_length",   &ti.br3_1_sg_direct_length);
    reader.AddVariable("br3_1_sg_length",          &ti.br3_1_sg_length);
    reader.AddVariable("br3_1_total_main_length",  &ti.br3_1_total_main_length);
    reader.AddVariable("br3_1_total_length",       &ti.br3_1_total_length);
    reader.AddVariable("br3_1_iso_angle",          &ti.br3_1_iso_angle);
    reader.AddVariable("br3_1_sg_flag_topology",   &ti.br3_1_sg_flag_topology);
    reader.AddVariable("br3_2_n_ele",              &ti.br3_2_n_ele);
    reader.AddVariable("br3_2_n_other",            &ti.br3_2_n_other);
    reader.AddVariable("br3_2_other_fid",          &ti.br3_2_other_fid);
    reader.AddVariable("br3_4_acc_length",         &ti.br3_4_acc_length);
    reader.AddVariable("br3_4_total_length",       &ti.br3_4_total_length);
    reader.AddVariable("br3_7_min_angle",          &ti.br3_7_min_angle);
    reader.AddVariable("br3_8_max_dQ_dx",          &ti.br3_8_max_dQ_dx);
    reader.AddVariable("br3_8_n_main_segs",        &ti.br3_8_n_main_segs);

    reader.BookMVA("MyBDT", m_br3_xml);

    if (ti.br_filled == 1)
        val = reader.EvaluateMVA("MyBDT");

    return val;
}

// ===========================================================================
// cal_br3_3_bdt
//
// Per-segment vector BDT (br3_3 sub-check). Returns minimum score.
// No fill gate — returns default_val if vector is empty.
//
// Prototype: NeutrinoID_nue_bdts.h::cal_br3_3_bdt()
// ===========================================================================
float UbooneNueBDTScorer::cal_br3_3_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    if (!m_rdr_br3_3) return default_val;

    float val = 1e9f;

    for (size_t i = 0; i != ti.br3_3_v_energy.size(); ++i) {
        m_br3_3_v_energy     = ti.br3_3_v_energy.at(i);
        m_br3_3_v_angle      = ti.br3_3_v_angle.at(i);
        m_br3_3_v_dir_length = ti.br3_3_v_dir_length.at(i);
        m_br3_3_v_length     = ti.br3_3_v_length.at(i);

        float tmp_val = m_rdr_br3_3->EvaluateMVA("MyBDT");
        if (tmp_val < val) val = tmp_val;
    }
    if (val > 1e8f) val = default_val;

    return val;
}

// ===========================================================================
// cal_br3_5_bdt
//
// Per-segment vector BDT (br3_5 sub-check). Returns minimum score.
// Note: br3_5_v_n_main_segs is NOT added (commented out in prototype).
//
// Prototype: NeutrinoID_nue_bdts.h::cal_br3_5_bdt()
// ===========================================================================
float UbooneNueBDTScorer::cal_br3_5_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    if (!m_rdr_br3_5) return default_val;

    float val = 1e9f;

    for (size_t i = 0; i != ti.br3_5_v_dir_length.size(); ++i) {
        m_br3_5_v_dir_length            = ti.br3_5_v_dir_length.at(i);
        m_br3_5_v_total_length          = ti.br3_5_v_total_length.at(i);
        m_br3_5_v_flag_avoid_muon_check = ti.br3_5_v_flag_avoid_muon_check.at(i);
        m_br3_5_v_n_seg                 = ti.br3_5_v_n_seg.at(i);
        m_br3_5_v_angle                 = ti.br3_5_v_angle.at(i);
        m_br3_5_v_sg_length             = ti.br3_5_v_sg_length.at(i);
        m_br3_5_v_energy                = ti.br3_5_v_energy.at(i);
        m_br3_5_v_n_segs                = ti.br3_5_v_n_segs.at(i);
        m_br3_5_v_shower_main_length    = ti.br3_5_v_shower_main_length.at(i);
        m_br3_5_v_shower_total_length   = ti.br3_5_v_shower_total_length.at(i);

        float tmp_val = m_rdr_br3_5->EvaluateMVA("MyBDT");
        if (tmp_val < val) val = tmp_val;
    }
    if (val > 1e8f) val = default_val;

    return val;
}

// ===========================================================================
// cal_br3_6_bdt
//
// Per-segment vector BDT (br3_6 sub-check). Returns minimum score.
//
// Prototype: NeutrinoID_nue_bdts.h::cal_br3_6_bdt()
// ===========================================================================
float UbooneNueBDTScorer::cal_br3_6_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    if (!m_rdr_br3_6) return default_val;

    float val = 1e9f;

    for (size_t i = 0; i != ti.br3_6_v_angle.size(); ++i) {
        m_br3_6_v_angle                  = ti.br3_6_v_angle.at(i);
        m_br3_6_v_angle1                 = ti.br3_6_v_angle1.at(i);
        m_br3_6_v_flag_shower_trajectory = ti.br3_6_v_flag_shower_trajectory.at(i);
        m_br3_6_v_direct_length          = ti.br3_6_v_direct_length.at(i);
        m_br3_6_v_length                 = ti.br3_6_v_length.at(i);
        m_br3_6_v_n_other_vtx_segs       = ti.br3_6_v_n_other_vtx_segs.at(i);
        m_br3_6_v_energy                 = ti.br3_6_v_energy.at(i);

        float tmp_val = m_rdr_br3_6->EvaluateMVA("MyBDT");
        if (tmp_val < val) val = tmp_val;
    }
    if (val > 1e8f) val = default_val;

    return val;
}

// ===========================================================================
// cal_stemdir_br2_bdt
//
// Scores stem-direction + bad-reconstruction-2 features (scalar).
// Gate: ti.br_filled == 1.
//
// Prototype: NeutrinoID_nue_bdts.h::cal_stemdir_br2_bdt()
// ===========================================================================
float UbooneNueBDTScorer::cal_stemdir_br2_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    float val = default_val;

    TMVA::Reader reader;
    reader.AddVariable("stem_dir_flag_single_shower", &ti.stem_dir_flag_single_shower);
    reader.AddVariable("stem_dir_angle",              &ti.stem_dir_angle);
    reader.AddVariable("stem_dir_energy",             &ti.stem_dir_energy);
    reader.AddVariable("stem_dir_angle1",             &ti.stem_dir_angle1);
    reader.AddVariable("stem_dir_angle2",             &ti.stem_dir_angle2);
    reader.AddVariable("stem_dir_angle3",             &ti.stem_dir_angle3);
    reader.AddVariable("stem_dir_ratio",              &ti.stem_dir_ratio);
    reader.AddVariable("br2_num_valid_tracks",        &ti.br2_num_valid_tracks);
    reader.AddVariable("br2_n_shower_main_segs",      &ti.br2_n_shower_main_segs);
    reader.AddVariable("br2_max_angle",               &ti.br2_max_angle);
    reader.AddVariable("br2_sg_length",               &ti.br2_sg_length);
    reader.AddVariable("br2_flag_sg_trajectory",      &ti.br2_flag_sg_trajectory);

    reader.BookMVA("MyBDT", m_stemdir_br2_xml);

    if (ti.br_filled == 1)
        val = reader.EvaluateMVA("MyBDT");

    return val;
}

// ===========================================================================
// cal_trimuon_bdt
//
// Scores stem-length + broken-muon + low-energy-michel combined features.
// Gate: ti.br_filled == 1.
//
// Prototype: NeutrinoID_nue_bdts.h::cal_trimuon_bdt()
// ===========================================================================
float UbooneNueBDTScorer::cal_trimuon_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    float val = default_val;

    TMVA::Reader reader;
    reader.AddVariable("stem_len_energy",                 &ti.stem_len_energy);
    reader.AddVariable("stem_len_length",                 &ti.stem_len_length);
    reader.AddVariable("stem_len_flag_avoid_muon_check",  &ti.stem_len_flag_avoid_muon_check);
    reader.AddVariable("stem_len_num_daughters",          &ti.stem_len_num_daughters);
    reader.AddVariable("stem_len_daughter_length",        &ti.stem_len_daughter_length);
    reader.AddVariable("brm_n_mu_segs",                   &ti.brm_n_mu_segs);
    reader.AddVariable("brm_Ep",                          &ti.brm_Ep);
    reader.AddVariable("brm_acc_length",                  &ti.brm_acc_length);
    reader.AddVariable("brm_shower_total_length",         &ti.brm_shower_total_length);
    reader.AddVariable("brm_connected_length",            &ti.brm_connected_length);
    reader.AddVariable("brm_n_size",                      &ti.brm_n_size);
    reader.AddVariable("brm_acc_direct_length",           &ti.brm_acc_direct_length);
    reader.AddVariable("brm_n_shower_main_segs",          &ti.brm_n_shower_main_segs);
    reader.AddVariable("brm_n_mu_main",                   &ti.brm_n_mu_main);
    reader.AddVariable("lem_shower_main_length",          &ti.lem_shower_main_length);
    reader.AddVariable("lem_n_3seg",                      &ti.lem_n_3seg);
    reader.AddVariable("lem_e_charge",                    &ti.lem_e_charge);
    reader.AddVariable("lem_e_dQdx",                      &ti.lem_e_dQdx);
    reader.AddVariable("lem_shower_num_main_segs",        &ti.lem_shower_num_main_segs);

    reader.BookMVA("MyBDT", m_trimuon_xml);

    if (ti.br_filled == 1)
        val = reader.EvaluateMVA("MyBDT");

    return val;
}

// ===========================================================================
// cal_br4_tro_bdt
//
// Scores bad-reconstruction-3 + track-overclustering-3 features (scalar).
// Gate: ti.br_filled == 1.
//
// Prototype: NeutrinoID_nue_bdts.h::cal_br4_tro_bdt()
// ===========================================================================
float UbooneNueBDTScorer::cal_br4_tro_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    float val = default_val;

    TMVA::Reader reader;
    reader.AddVariable("br4_1_shower_main_length",    &ti.br4_1_shower_main_length);
    reader.AddVariable("br4_1_shower_total_length",   &ti.br4_1_shower_total_length);
    reader.AddVariable("br4_1_min_dis",               &ti.br4_1_min_dis);
    reader.AddVariable("br4_1_energy",                &ti.br4_1_energy);
    reader.AddVariable("br4_1_flag_avoid_muon_check", &ti.br4_1_flag_avoid_muon_check);
    reader.AddVariable("br4_1_n_vtx_segs",            &ti.br4_1_n_vtx_segs);
    reader.AddVariable("br4_1_n_main_segs",           &ti.br4_1_n_main_segs);
    reader.AddVariable("br4_2_ratio_45",              &ti.br4_2_ratio_45);
    reader.AddVariable("br4_2_ratio_35",              &ti.br4_2_ratio_35);
    reader.AddVariable("br4_2_ratio_25",              &ti.br4_2_ratio_25);
    reader.AddVariable("br4_2_ratio_15",              &ti.br4_2_ratio_15);
    reader.AddVariable("br4_2_ratio1_45",             &ti.br4_2_ratio1_45);
    reader.AddVariable("br4_2_ratio1_35",             &ti.br4_2_ratio1_35);
    reader.AddVariable("br4_2_ratio1_25",             &ti.br4_2_ratio1_25);
    reader.AddVariable("br4_2_ratio1_15",             &ti.br4_2_ratio1_15);
    reader.AddVariable("br4_2_iso_angle",             &ti.br4_2_iso_angle);
    reader.AddVariable("br4_2_iso_angle1",            &ti.br4_2_iso_angle1);
    reader.AddVariable("br4_2_angle",                 &ti.br4_2_angle);
    reader.AddVariable("tro_3_stem_length",           &ti.tro_3_stem_length);
    reader.AddVariable("tro_3_n_muon_segs",           &ti.tro_3_n_muon_segs);

    reader.BookMVA("MyBDT", m_br4_tro_xml);

    if (ti.br_filled == 1)
        val = reader.EvaluateMVA("MyBDT");

    return val;
}

// ===========================================================================
// cal_mipquality_bdt
//
// Scores MIP quality features (scalar, per-event).
// Gate: ti.mip_quality_filled == 1.
//
// Prototype: NeutrinoID_nue_bdts.h::cal_mipquality_bdt()
// ===========================================================================
float UbooneNueBDTScorer::cal_mipquality_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    float val = default_val;

    TMVA::Reader reader;
    reader.AddVariable("mip_quality_energy",            &ti.mip_quality_energy);
    reader.AddVariable("mip_quality_overlap",           &ti.mip_quality_overlap);
    reader.AddVariable("mip_quality_n_showers",         &ti.mip_quality_n_showers);
    reader.AddVariable("mip_quality_n_tracks",          &ti.mip_quality_n_tracks);
    reader.AddVariable("mip_quality_flag_inside_pi0",   &ti.mip_quality_flag_inside_pi0);
    reader.AddVariable("mip_quality_n_pi0_showers",     &ti.mip_quality_n_pi0_showers);
    reader.AddVariable("mip_quality_shortest_length",   &ti.mip_quality_shortest_length);
    reader.AddVariable("mip_quality_acc_length",        &ti.mip_quality_acc_length);
    reader.AddVariable("mip_quality_shortest_angle",    &ti.mip_quality_shortest_angle);
    reader.AddVariable("mip_quality_flag_proton",       &ti.mip_quality_flag_proton);

    reader.BookMVA("MyBDT", m_mipquality_xml);

    if (ti.mip_quality_filled == 1)
        val = reader.EvaluateMVA("MyBDT");

    return val;
}

// ===========================================================================
// cal_pio_1_bdt
//
// Scores pi0 type-1 (vertex-attached pi0) features (scalar).
// Gate: ti.pio_filled == 1 && ti.pio_flag_pio == 1.
//
// Prototype: NeutrinoID_nue_bdts.h::cal_pio_1_bdt()
// ===========================================================================
float UbooneNueBDTScorer::cal_pio_1_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    float val = default_val;

    TMVA::Reader reader;
    reader.AddVariable("pio_1_mass",     &ti.pio_1_mass);
    reader.AddVariable("pio_1_pio_type", &ti.pio_1_pio_type);
    reader.AddVariable("pio_1_energy_1", &ti.pio_1_energy_1);
    reader.AddVariable("pio_1_energy_2", &ti.pio_1_energy_2);
    reader.AddVariable("pio_1_dis_1",    &ti.pio_1_dis_1);
    reader.AddVariable("pio_1_dis_2",    &ti.pio_1_dis_2);
    reader.AddVariable("pio_mip_id",     &ti.pio_mip_id);

    reader.BookMVA("MyBDT", m_pio_1_xml);

    if (ti.pio_filled == 1 && ti.pio_flag_pio == 1)
        val = reader.EvaluateMVA("MyBDT");

    return val;
}

// ===========================================================================
// cal_pio_2_bdt
//
// Per-pi0 vector BDT for type-2 (non-vertex-attached) pi0s.
// Gate: ti.pio_filled == 1 && ti.pio_flag_pio == 0.
// Returns minimum score; default_val if condition not met or vector empty.
//
// Prototype: NeutrinoID_nue_bdts.h::cal_pio_2_bdt()
// ===========================================================================
float UbooneNueBDTScorer::cal_pio_2_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    if (!m_rdr_pio_2) return default_val;

    float val = 1e9f;

    if (ti.pio_filled == 1 && ti.pio_flag_pio == 0) {
        m_pio_2_v_pio_mip_id = ti.pio_mip_id;
        for (size_t i = 0; i != ti.pio_2_v_dis2.size(); ++i) {
            m_pio_2_v_dis2       = ti.pio_2_v_dis2.at(i);
            m_pio_2_v_angle2     = ti.pio_2_v_angle2.at(i);
            m_pio_2_v_acc_length = ti.pio_2_v_acc_length.at(i);

            float tmp_val = m_rdr_pio_2->EvaluateMVA("MyBDT");
            if (tmp_val < val) val = tmp_val;
        }
    }
    if (val > 1e8f) val = default_val;

    return val;
}

// ===========================================================================
// cal_stw_spt_bdt
//
// Scores shower-to-wall-1 + single-point features (scalar).
// Gate: ti.br_filled == 1.
//
// Prototype: NeutrinoID_nue_bdts.h::cal_stw_spt_bdt()
// ===========================================================================
float UbooneNueBDTScorer::cal_stw_spt_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    float val = default_val;

    TMVA::Reader reader;
    reader.AddVariable("stw_1_energy",            &ti.stw_1_energy);
    reader.AddVariable("stw_1_dis",               &ti.stw_1_dis);
    reader.AddVariable("stw_1_dQ_dx",             &ti.stw_1_dQ_dx);
    reader.AddVariable("stw_1_flag_single_shower", &ti.stw_1_flag_single_shower);
    reader.AddVariable("stw_1_n_pi0",             &ti.stw_1_n_pi0);
    reader.AddVariable("stw_1_num_valid_tracks",  &ti.stw_1_num_valid_tracks);
    reader.AddVariable("spt_shower_main_length",  &ti.spt_shower_main_length);
    reader.AddVariable("spt_shower_total_length", &ti.spt_shower_total_length);
    reader.AddVariable("spt_angle_beam",          &ti.spt_angle_beam);
    reader.AddVariable("spt_angle_vertical",      &ti.spt_angle_vertical);
    reader.AddVariable("spt_max_dQ_dx",           &ti.spt_max_dQ_dx);
    reader.AddVariable("spt_angle_beam_1",        &ti.spt_angle_beam_1);
    reader.AddVariable("spt_angle_drift",         &ti.spt_angle_drift);
    reader.AddVariable("spt_angle_drift_1",       &ti.spt_angle_drift_1);
    reader.AddVariable("spt_num_valid_tracks",    &ti.spt_num_valid_tracks);
    reader.AddVariable("spt_n_vtx_segs",          &ti.spt_n_vtx_segs);
    reader.AddVariable("spt_max_length",          &ti.spt_max_length);

    reader.BookMVA("MyBDT", m_stw_spt_xml);

    if (ti.br_filled == 1)
        val = reader.EvaluateMVA("MyBDT");

    return val;
}

// ===========================================================================
// cal_vis_1_bdt
//
// Scores vertex-inside-shower type-1 features (scalar).
// Gate: ti.vis_1_filled == 1.
//
// Prototype: NeutrinoID_nue_bdts.h::cal_vis_1_bdt()
// ===========================================================================
float UbooneNueBDTScorer::cal_vis_1_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    float val = default_val;

    TMVA::Reader reader;
    reader.AddVariable("vis_1_n_vtx_segs",        &ti.vis_1_n_vtx_segs);
    reader.AddVariable("vis_1_energy",             &ti.vis_1_energy);
    reader.AddVariable("vis_1_num_good_tracks",    &ti.vis_1_num_good_tracks);
    reader.AddVariable("vis_1_max_angle",          &ti.vis_1_max_angle);
    reader.AddVariable("vis_1_max_shower_angle",   &ti.vis_1_max_shower_angle);
    reader.AddVariable("vis_1_tmp_length1",        &ti.vis_1_tmp_length1);
    reader.AddVariable("vis_1_tmp_length2",        &ti.vis_1_tmp_length2);

    reader.BookMVA("MyBDT", m_vis_1_xml);

    if (ti.vis_1_filled == 1)
        val = reader.EvaluateMVA("MyBDT");

    return val;
}

// ===========================================================================
// cal_vis_2_bdt
//
// Scores vertex-inside-shower type-2 features (scalar).
// Gate: ti.vis_2_filled == 1.
//
// Prototype: NeutrinoID_nue_bdts.h::cal_vis_2_bdt()
// ===========================================================================
float UbooneNueBDTScorer::cal_vis_2_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    float val = default_val;

    TMVA::Reader reader;
    reader.AddVariable("vis_2_n_vtx_segs",        &ti.vis_2_n_vtx_segs);
    reader.AddVariable("vis_2_min_angle",          &ti.vis_2_min_angle);
    reader.AddVariable("vis_2_min_weak_track",     &ti.vis_2_min_weak_track);
    reader.AddVariable("vis_2_angle_beam",         &ti.vis_2_angle_beam);
    reader.AddVariable("vis_2_min_angle1",         &ti.vis_2_min_angle1);
    reader.AddVariable("vis_2_iso_angle1",         &ti.vis_2_iso_angle1);
    reader.AddVariable("vis_2_min_medium_dQ_dx",   &ti.vis_2_min_medium_dQ_dx);
    reader.AddVariable("vis_2_min_length",         &ti.vis_2_min_length);
    reader.AddVariable("vis_2_sg_length",          &ti.vis_2_sg_length);
    reader.AddVariable("vis_2_max_angle",          &ti.vis_2_max_angle);
    reader.AddVariable("vis_2_max_weak_track",     &ti.vis_2_max_weak_track);

    reader.BookMVA("MyBDT", m_vis_2_xml);

    if (ti.vis_2_filled == 1)
        val = reader.EvaluateMVA("MyBDT");

    return val;
}

// ===========================================================================
// cal_stw_2_bdt
//
// Per-segment vector BDT (shower-to-wall type 2). Returns minimum score.
//
// Prototype: NeutrinoID_nue_bdts.h::cal_stw_2_bdt()
// ===========================================================================
float UbooneNueBDTScorer::cal_stw_2_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    if (!m_rdr_stw_2) return default_val;

    float val = 1e9f;

    for (size_t i = 0; i != ti.stw_2_v_medium_dQ_dx.size(); ++i) {
        m_stw_2_v_medium_dQ_dx = ti.stw_2_v_medium_dQ_dx.at(i);
        m_stw_2_v_energy        = ti.stw_2_v_energy.at(i);
        m_stw_2_v_angle         = ti.stw_2_v_angle.at(i);
        m_stw_2_v_dir_length    = ti.stw_2_v_dir_length.at(i);
        m_stw_2_v_max_dQ_dx     = ti.stw_2_v_max_dQ_dx.at(i);

        float tmp_val = m_rdr_stw_2->EvaluateMVA("MyBDT");
        if (tmp_val < val) val = tmp_val;
    }
    if (val > 1e8f) val = default_val;

    return val;
}

// ===========================================================================
// cal_stw_3_bdt
//
// Per-segment vector BDT (shower-to-wall type 3). Returns minimum score.
//
// Prototype: NeutrinoID_nue_bdts.h::cal_stw_3_bdt()
// ===========================================================================
float UbooneNueBDTScorer::cal_stw_3_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    if (!m_rdr_stw_3) return default_val;

    float val = 1e9f;

    for (size_t i = 0; i != ti.stw_3_v_angle.size(); ++i) {
        m_stw_3_v_angle        = ti.stw_3_v_angle.at(i);
        m_stw_3_v_dir_length   = ti.stw_3_v_dir_length.at(i);
        m_stw_3_v_energy       = ti.stw_3_v_energy.at(i);
        m_stw_3_v_medium_dQ_dx = ti.stw_3_v_medium_dQ_dx.at(i);

        float tmp_val = m_rdr_stw_3->EvaluateMVA("MyBDT");
        if (tmp_val < val) val = tmp_val;
    }
    if (val > 1e8f) val = default_val;

    return val;
}

// ===========================================================================
// cal_stw_4_bdt
//
// Per-segment vector BDT (shower-to-wall type 4). Returns minimum score.
//
// Prototype: NeutrinoID_nue_bdts.h::cal_stw_4_bdt()
// ===========================================================================
float UbooneNueBDTScorer::cal_stw_4_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    if (!m_rdr_stw_4) return default_val;

    float val = 1e9f;

    for (size_t i = 0; i != ti.stw_4_v_angle.size(); ++i) {
        m_stw_4_v_angle  = ti.stw_4_v_angle.at(i);
        m_stw_4_v_dis    = ti.stw_4_v_dis.at(i);
        m_stw_4_v_energy = ti.stw_4_v_energy.at(i);

        float tmp_val = m_rdr_stw_4->EvaluateMVA("MyBDT");
        if (tmp_val < val) val = tmp_val;
    }
    if (val > 1e8f) val = default_val;

    return val;
}

// ===========================================================================
// cal_sig_1_bdt
//
// Per-segment vector BDT (single-shower pi0 tagger type 1). Returns minimum score.
//
// Prototype: NeutrinoID_nue_bdts.h::cal_sig_1_bdt()
// ===========================================================================
float UbooneNueBDTScorer::cal_sig_1_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    if (!m_rdr_sig_1) return default_val;

    float val = 1e9f;

    for (size_t i = 0; i != ti.sig_1_v_angle.size(); ++i) {
        m_sig_1_v_angle              = ti.sig_1_v_angle.at(i);
        m_sig_1_v_flag_single_shower = ti.sig_1_v_flag_single_shower.at(i);
        m_sig_1_v_energy             = ti.sig_1_v_energy.at(i);
        m_sig_1_v_energy_1           = ti.sig_1_v_energy_1.at(i);

        float tmp_val = m_rdr_sig_1->EvaluateMVA("MyBDT");
        if (tmp_val < val) val = tmp_val;
    }
    if (val > 1e8f) val = default_val;

    return val;
}

// ===========================================================================
// cal_sig_2_bdt
//
// Per-segment vector BDT (single-shower pi0 tagger type 2). Returns minimum score.
//
// Prototype: NeutrinoID_nue_bdts.h::cal_sig_2_bdt()
// ===========================================================================
float UbooneNueBDTScorer::cal_sig_2_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    if (!m_rdr_sig_2) return default_val;

    float val = 1e9f;

    for (size_t i = 0; i != ti.sig_2_v_energy.size(); ++i) {
        m_sig_2_v_energy             = ti.sig_2_v_energy.at(i);
        m_sig_2_v_shower_angle       = ti.sig_2_v_shower_angle.at(i);
        m_sig_2_v_flag_single_shower = ti.sig_2_v_flag_single_shower.at(i);
        m_sig_2_v_medium_dQ_dx       = ti.sig_2_v_medium_dQ_dx.at(i);
        m_sig_2_v_start_dQ_dx        = ti.sig_2_v_start_dQ_dx.at(i);

        float tmp_val = m_rdr_sig_2->EvaluateMVA("MyBDT");
        if (tmp_val < val) val = tmp_val;
    }
    if (val > 1e8f) val = default_val;

    return val;
}

// ===========================================================================
// cal_lol_1_bdt
//
// Per-segment vector BDT (low-energy-overlapping type 1). Returns minimum score.
//
// Prototype: NeutrinoID_nue_bdts.h::cal_lol_1_bdt()
// ===========================================================================
float UbooneNueBDTScorer::cal_lol_1_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    if (!m_rdr_lol_1) return default_val;

    float val = 1e9f;

    for (size_t i = 0; i != ti.lol_1_v_energy.size(); ++i) {
        m_lol_1_v_energy     = ti.lol_1_v_energy.at(i);
        m_lol_1_v_vtx_n_segs = ti.lol_1_v_vtx_n_segs.at(i);
        m_lol_1_v_nseg       = ti.lol_1_v_nseg.at(i);
        m_lol_1_v_angle      = ti.lol_1_v_angle.at(i);

        float tmp_val = m_rdr_lol_1->EvaluateMVA("MyBDT");
        if (tmp_val < val) val = tmp_val;
    }
    if (val > 1e8f) val = default_val;

    return val;
}

// ===========================================================================
// cal_lol_2_bdt
//
// Per-segment vector BDT (low-energy-overlapping type 2). Returns minimum score.
//
// Prototype: NeutrinoID_nue_bdts.h::cal_lol_2_bdt()
// ===========================================================================
float UbooneNueBDTScorer::cal_lol_2_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    if (!m_rdr_lol_2) return default_val;

    float val = 1e9f;

    for (size_t i = 0; i != ti.lol_2_v_length.size(); ++i) {
        m_lol_2_v_length             = ti.lol_2_v_length.at(i);
        m_lol_2_v_angle              = ti.lol_2_v_angle.at(i);
        m_lol_2_v_type               = ti.lol_2_v_type.at(i);
        m_lol_2_v_vtx_n_segs         = ti.lol_2_v_vtx_n_segs.at(i);
        m_lol_2_v_energy             = ti.lol_2_v_energy.at(i);
        m_lol_2_v_shower_main_length = ti.lol_2_v_shower_main_length.at(i);
        m_lol_2_v_flag_dir_weak      = ti.lol_2_v_flag_dir_weak.at(i);

        float tmp_val = m_rdr_lol_2->EvaluateMVA("MyBDT");
        if (tmp_val < val) val = tmp_val;
    }
    if (val > 1e8f) val = default_val;

    return val;
}

// ===========================================================================
// cal_tro_1_bdt
//
// Per-segment vector BDT (track-overclustering type 1). Returns minimum score.
//
// Prototype: NeutrinoID_nue_bdts.h::cal_tro_1_bdt()
// ===========================================================================
float UbooneNueBDTScorer::cal_tro_1_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    if (!m_rdr_tro_1) return default_val;

    float val = 1e9f;

    for (size_t i = 0; i != ti.tro_1_v_particle_type.size(); ++i) {
        m_tro_1_v_particle_type      = ti.tro_1_v_particle_type.at(i);
        m_tro_1_v_flag_dir_weak      = ti.tro_1_v_flag_dir_weak.at(i);
        m_tro_1_v_min_dis            = ti.tro_1_v_min_dis.at(i);
        m_tro_1_v_sg1_length         = ti.tro_1_v_sg1_length.at(i);
        m_tro_1_v_shower_main_length = ti.tro_1_v_shower_main_length.at(i);
        m_tro_1_v_max_n_vtx_segs     = ti.tro_1_v_max_n_vtx_segs.at(i);
        m_tro_1_v_tmp_length         = ti.tro_1_v_tmp_length.at(i);
        m_tro_1_v_medium_dQ_dx       = ti.tro_1_v_medium_dQ_dx.at(i);
        m_tro_1_v_dQ_dx_cut          = ti.tro_1_v_dQ_dx_cut.at(i);
        m_tro_1_v_flag_shower_topology = ti.tro_1_v_flag_shower_topology.at(i);

        float tmp_val = m_rdr_tro_1->EvaluateMVA("MyBDT");
        if (tmp_val < val) val = tmp_val;
    }
    if (val > 1e8f) val = default_val;

    return val;
}

// ===========================================================================
// cal_tro_2_bdt
//
// Per-segment vector BDT (track-overclustering type 2). Returns minimum score.
//
// Prototype: NeutrinoID_nue_bdts.h::cal_tro_2_bdt()
// ===========================================================================
float UbooneNueBDTScorer::cal_tro_2_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    if (!m_rdr_tro_2) return default_val;

    float val = 1e9f;

    for (size_t i = 0; i != ti.tro_2_v_energy.size(); ++i) {
        m_tro_2_v_energy      = ti.tro_2_v_energy.at(i);
        m_tro_2_v_stem_length = ti.tro_2_v_stem_length.at(i);
        m_tro_2_v_iso_angle   = ti.tro_2_v_iso_angle.at(i);
        m_tro_2_v_max_length  = ti.tro_2_v_max_length.at(i);
        m_tro_2_v_angle       = ti.tro_2_v_angle.at(i);

        float tmp_val = m_rdr_tro_2->EvaluateMVA("MyBDT");
        if (tmp_val < val) val = tmp_val;
    }
    if (val > 1e8f) val = default_val;

    return val;
}

// ===========================================================================
// cal_tro_4_bdt
//
// Per-segment vector BDT (track-overclustering type 4). Returns minimum score.
//
// Prototype: NeutrinoID_nue_bdts.h::cal_tro_4_bdt()
// ===========================================================================
float UbooneNueBDTScorer::cal_tro_4_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    if (!m_rdr_tro_4) return default_val;

    float val = 1e9f;

    for (size_t i = 0; i != ti.tro_4_v_dir2_mag.size(); ++i) {
        m_tro_4_v_dir2_mag              = ti.tro_4_v_dir2_mag.at(i);
        m_tro_4_v_angle                 = ti.tro_4_v_angle.at(i);
        m_tro_4_v_angle1                = ti.tro_4_v_angle1.at(i);
        m_tro_4_v_angle2                = ti.tro_4_v_angle2.at(i);
        m_tro_4_v_length                = ti.tro_4_v_length.at(i);
        m_tro_4_v_length1               = ti.tro_4_v_length1.at(i);
        m_tro_4_v_medium_dQ_dx          = ti.tro_4_v_medium_dQ_dx.at(i);
        m_tro_4_v_end_dQ_dx             = ti.tro_4_v_end_dQ_dx.at(i);
        m_tro_4_v_energy                = ti.tro_4_v_energy.at(i);
        m_tro_4_v_shower_main_length    = ti.tro_4_v_shower_main_length.at(i);
        m_tro_4_v_flag_shower_trajectory= ti.tro_4_v_flag_shower_trajectory.at(i);

        float tmp_val = m_rdr_tro_4->EvaluateMVA("MyBDT");
        if (tmp_val < val) val = tmp_val;
    }
    if (val > 1e8f) val = default_val;

    return val;
}

// ===========================================================================
// cal_tro_5_bdt
//
// Per-segment vector BDT (track-overclustering type 5). Returns minimum score.
//
// Prototype: NeutrinoID_nue_bdts.h::cal_tro_5_bdt()
// ===========================================================================
float UbooneNueBDTScorer::cal_tro_5_bdt(Clus::PR::TaggerInfo& ti, float default_val) const
{
    if (!m_rdr_tro_5) return default_val;

    float val = 1e9f;

    for (size_t i = 0; i != ti.tro_5_v_max_angle.size(); ++i) {
        m_tro_5_v_max_angle   = ti.tro_5_v_max_angle.at(i);
        m_tro_5_v_min_angle   = ti.tro_5_v_min_angle.at(i);
        m_tro_5_v_max_length  = ti.tro_5_v_max_length.at(i);
        m_tro_5_v_iso_angle   = ti.tro_5_v_iso_angle.at(i);
        m_tro_5_v_n_vtx_segs  = ti.tro_5_v_n_vtx_segs.at(i);
        m_tro_5_v_min_count   = ti.tro_5_v_min_count.at(i);
        m_tro_5_v_max_count   = ti.tro_5_v_max_count.at(i);
        m_tro_5_v_energy      = ti.tro_5_v_energy.at(i);

        float tmp_val = m_rdr_tro_5->EvaluateMVA("MyBDT");
        if (tmp_val < val) val = tmp_val;
    }
    if (val > 1e8f) val = default_val;

    return val;
}

// ===========================================================================
// cal_bdts_xgboost
//
// Top-level nueCC XGBoost scorer.
//
// Step 1: Compute the 15 vector sub-BDT scores (these appear as features in the
//         final model in addition to the raw scalar features).
// Step 2: Apply variable protection (clamping and NaN guards).
// Step 3: Copy ti fields to xgboost member float buffer, then evaluate.
// Step 4: Transform via log-odds: nue_score = log10((1+val)/(1-val)).
//
// Prototype: NeutrinoID_nue_bdts.h::cal_bdts_xgboost()
// ===========================================================================
void UbooneNueBDTScorer::cal_bdts_xgboost(Clus::PR::TaggerInfo& ti,
                                           const Clus::PR::KineInfo& ki) const
{
    // ---- Step 1: evaluate 15 vector sub-BDTs and store their scores ----------
    ti.br3_3_score = cal_br3_3_bdt(ti, 0.3f);
    ti.br3_5_score = cal_br3_5_bdt(ti, 0.42f);
    ti.br3_6_score = cal_br3_6_bdt(ti, 0.75f);
    ti.pio_2_score = cal_pio_2_bdt(ti, 0.2f);
    ti.stw_2_score = cal_stw_2_bdt(ti, 0.7f);
    ti.stw_3_score = cal_stw_3_bdt(ti, 0.5f);
    ti.stw_4_score = cal_stw_4_bdt(ti, 0.7f);
    ti.sig_1_score = cal_sig_1_bdt(ti, 0.59f);
    ti.sig_2_score = cal_sig_2_bdt(ti, 0.55f);
    ti.lol_1_score = cal_lol_1_bdt(ti, 0.85f);
    ti.lol_2_score = cal_lol_2_bdt(ti, 0.7f);
    ti.tro_1_score = cal_tro_1_bdt(ti, 0.28f);
    ti.tro_2_score = cal_tro_2_bdt(ti, 0.35f);
    ti.tro_4_score = cal_tro_4_bdt(ti, 0.33f);
    ti.tro_5_score = cal_tro_5_bdt(ti, 0.5f);

    // ---- Step 2: variable protection ----------------------------------------
    if (ti.mip_min_dis > 1000.f) ti.mip_min_dis = 1000.f;
    if (ti.mip_quality_shortest_length > 1000.f) ti.mip_quality_shortest_length = 1000.f;
    if (std::isnan(ti.mip_quality_shortest_angle)) ti.mip_quality_shortest_angle = 0.f;
    if (std::isnan(ti.stem_dir_ratio)) ti.stem_dir_ratio = 1.0f;

    // ---- Step 3: copy ti fields to xgboost member float buffer, evaluate ----
    if (!m_rdr_xgboost) {
        ti.nue_score = -15.f;
        return;
    }

    m_xgb_match_isFC                   = ti.match_isFC;
    m_xgb_kine_reco_Enu                = static_cast<float>(ki.kine_reco_Enu);
    m_xgb_cme_mu_energy                = ti.cme_mu_energy;
    m_xgb_cme_energy                   = ti.cme_energy;
    m_xgb_cme_mu_length                = ti.cme_mu_length;
    m_xgb_cme_length                   = ti.cme_length;
    m_xgb_cme_angle_beam               = ti.cme_angle_beam;
    m_xgb_anc_angle                    = ti.anc_angle;
    m_xgb_anc_max_angle                = ti.anc_max_angle;
    m_xgb_anc_max_length               = ti.anc_max_length;
    m_xgb_anc_acc_forward_length       = ti.anc_acc_forward_length;
    m_xgb_anc_acc_backward_length      = ti.anc_acc_backward_length;
    m_xgb_anc_acc_forward_length1      = ti.anc_acc_forward_length1;
    m_xgb_anc_shower_main_length       = ti.anc_shower_main_length;
    m_xgb_anc_shower_total_length      = ti.anc_shower_total_length;
    m_xgb_anc_flag_main_outside        = ti.anc_flag_main_outside;
    m_xgb_gap_flag_prolong_u           = ti.gap_flag_prolong_u;
    m_xgb_gap_flag_prolong_v           = ti.gap_flag_prolong_v;
    m_xgb_gap_flag_prolong_w           = ti.gap_flag_prolong_w;
    m_xgb_gap_flag_parallel            = ti.gap_flag_parallel;
    m_xgb_gap_n_points                 = ti.gap_n_points;
    m_xgb_gap_n_bad                    = ti.gap_n_bad;
    m_xgb_gap_energy                   = ti.gap_energy;
    m_xgb_gap_num_valid_tracks         = ti.gap_num_valid_tracks;
    m_xgb_gap_flag_single_shower       = ti.gap_flag_single_shower;
    m_xgb_hol_1_n_valid_tracks         = ti.hol_1_n_valid_tracks;
    m_xgb_hol_1_min_angle              = ti.hol_1_min_angle;
    m_xgb_hol_1_energy                 = ti.hol_1_energy;
    m_xgb_hol_1_min_length             = ti.hol_1_min_length;
    m_xgb_hol_2_min_angle              = ti.hol_2_min_angle;
    m_xgb_hol_2_medium_dQ_dx           = ti.hol_2_medium_dQ_dx;
    m_xgb_hol_2_ncount                 = ti.hol_2_ncount;
    m_xgb_lol_3_angle_beam             = ti.lol_3_angle_beam;
    m_xgb_lol_3_n_valid_tracks         = ti.lol_3_n_valid_tracks;
    m_xgb_lol_3_min_angle              = ti.lol_3_min_angle;
    m_xgb_lol_3_vtx_n_segs             = ti.lol_3_vtx_n_segs;
    m_xgb_lol_3_shower_main_length     = ti.lol_3_shower_main_length;
    m_xgb_lol_3_n_out                  = ti.lol_3_n_out;
    m_xgb_lol_3_n_sum                  = ti.lol_3_n_sum;
    m_xgb_hol_1_flag_all_shower        = ti.hol_1_flag_all_shower;
    m_xgb_mgo_energy                   = ti.mgo_energy;
    m_xgb_mgo_max_energy               = ti.mgo_max_energy;
    m_xgb_mgo_total_energy             = ti.mgo_total_energy;
    m_xgb_mgo_n_showers                = ti.mgo_n_showers;
    m_xgb_mgo_max_energy_1             = ti.mgo_max_energy_1;
    m_xgb_mgo_max_energy_2             = ti.mgo_max_energy_2;
    m_xgb_mgo_total_other_energy       = ti.mgo_total_other_energy;
    m_xgb_mgo_n_total_showers          = ti.mgo_n_total_showers;
    m_xgb_mgo_total_other_energy_1     = ti.mgo_total_other_energy_1;
    m_xgb_mgt_flag_single_shower       = ti.mgt_flag_single_shower;
    m_xgb_mgt_max_energy               = ti.mgt_max_energy;
    m_xgb_mgt_total_other_energy       = ti.mgt_total_other_energy;
    m_xgb_mgt_max_energy_1             = ti.mgt_max_energy_1;
    m_xgb_mgt_e_indirect_max_energy    = ti.mgt_e_indirect_max_energy;
    m_xgb_mgt_e_direct_max_energy      = ti.mgt_e_direct_max_energy;
    m_xgb_mgt_n_direct_showers         = ti.mgt_n_direct_showers;
    m_xgb_mgt_e_direct_total_energy    = ti.mgt_e_direct_total_energy;
    m_xgb_mgt_flag_indirect_max_pio    = ti.mgt_flag_indirect_max_pio;
    m_xgb_mgt_e_indirect_total_energy  = ti.mgt_e_indirect_total_energy;
    m_xgb_mip_quality_energy           = ti.mip_quality_energy;
    m_xgb_mip_quality_overlap          = ti.mip_quality_overlap;
    m_xgb_mip_quality_n_showers        = ti.mip_quality_n_showers;
    m_xgb_mip_quality_n_tracks         = ti.mip_quality_n_tracks;
    m_xgb_mip_quality_flag_inside_pi0  = ti.mip_quality_flag_inside_pi0;
    m_xgb_mip_quality_n_pi0_showers    = ti.mip_quality_n_pi0_showers;
    m_xgb_mip_quality_shortest_length  = ti.mip_quality_shortest_length;
    m_xgb_mip_quality_acc_length       = ti.mip_quality_acc_length;
    m_xgb_mip_quality_shortest_angle   = ti.mip_quality_shortest_angle;
    m_xgb_mip_quality_flag_proton      = ti.mip_quality_flag_proton;
    m_xgb_br1_1_shower_type            = ti.br1_1_shower_type;
    m_xgb_br1_1_vtx_n_segs             = ti.br1_1_vtx_n_segs;
    m_xgb_br1_1_energy                 = ti.br1_1_energy;
    m_xgb_br1_1_n_segs                 = ti.br1_1_n_segs;
    m_xgb_br1_1_flag_sg_topology       = ti.br1_1_flag_sg_topology;
    m_xgb_br1_1_flag_sg_trajectory     = ti.br1_1_flag_sg_trajectory;
    m_xgb_br1_1_sg_length              = ti.br1_1_sg_length;
    m_xgb_br1_2_n_connected            = ti.br1_2_n_connected;
    m_xgb_br1_2_max_length             = ti.br1_2_max_length;
    m_xgb_br1_2_n_connected_1          = ti.br1_2_n_connected_1;
    m_xgb_br1_2_n_shower_segs          = ti.br1_2_n_shower_segs;
    m_xgb_br1_2_max_length_ratio       = ti.br1_2_max_length_ratio;
    m_xgb_br1_2_shower_length          = ti.br1_2_shower_length;
    m_xgb_br1_3_n_connected_p          = ti.br1_3_n_connected_p;
    m_xgb_br1_3_max_length_p           = ti.br1_3_max_length_p;
    m_xgb_br1_3_n_shower_main_segs     = ti.br1_3_n_shower_main_segs;
    m_xgb_br3_1_energy                 = ti.br3_1_energy;
    m_xgb_br3_1_n_shower_segments      = ti.br3_1_n_shower_segments;
    m_xgb_br3_1_sg_flag_trajectory     = ti.br3_1_sg_flag_trajectory;
    m_xgb_br3_1_sg_direct_length       = ti.br3_1_sg_direct_length;
    m_xgb_br3_1_sg_length              = ti.br3_1_sg_length;
    m_xgb_br3_1_total_main_length      = ti.br3_1_total_main_length;
    m_xgb_br3_1_total_length           = ti.br3_1_total_length;
    m_xgb_br3_1_iso_angle              = ti.br3_1_iso_angle;
    m_xgb_br3_1_sg_flag_topology       = ti.br3_1_sg_flag_topology;
    m_xgb_br3_2_n_ele                  = ti.br3_2_n_ele;
    m_xgb_br3_2_n_other                = ti.br3_2_n_other;
    m_xgb_br3_2_other_fid              = ti.br3_2_other_fid;
    m_xgb_br3_4_acc_length             = ti.br3_4_acc_length;
    m_xgb_br3_4_total_length           = ti.br3_4_total_length;
    m_xgb_br3_7_min_angle              = ti.br3_7_min_angle;
    m_xgb_br3_8_max_dQ_dx              = ti.br3_8_max_dQ_dx;
    m_xgb_br3_8_n_main_segs            = ti.br3_8_n_main_segs;
    m_xgb_vis_1_n_vtx_segs             = ti.vis_1_n_vtx_segs;
    m_xgb_vis_1_energy                 = ti.vis_1_energy;
    m_xgb_vis_1_num_good_tracks        = ti.vis_1_num_good_tracks;
    m_xgb_vis_1_max_angle              = ti.vis_1_max_angle;
    m_xgb_vis_1_max_shower_angle       = ti.vis_1_max_shower_angle;
    m_xgb_vis_1_tmp_length1            = ti.vis_1_tmp_length1;
    m_xgb_vis_1_tmp_length2            = ti.vis_1_tmp_length2;
    m_xgb_vis_2_n_vtx_segs             = ti.vis_2_n_vtx_segs;
    m_xgb_vis_2_min_angle              = ti.vis_2_min_angle;
    m_xgb_vis_2_min_weak_track         = ti.vis_2_min_weak_track;
    m_xgb_vis_2_angle_beam             = ti.vis_2_angle_beam;
    m_xgb_vis_2_min_angle1             = ti.vis_2_min_angle1;
    m_xgb_vis_2_iso_angle1             = ti.vis_2_iso_angle1;
    m_xgb_vis_2_min_medium_dQ_dx       = ti.vis_2_min_medium_dQ_dx;
    m_xgb_vis_2_min_length             = ti.vis_2_min_length;
    m_xgb_vis_2_sg_length              = ti.vis_2_sg_length;
    m_xgb_vis_2_max_angle              = ti.vis_2_max_angle;
    m_xgb_vis_2_max_weak_track         = ti.vis_2_max_weak_track;
    m_xgb_pio_1_mass                   = ti.pio_1_mass;
    m_xgb_pio_1_pio_type               = ti.pio_1_pio_type;
    m_xgb_pio_1_energy_1               = ti.pio_1_energy_1;
    m_xgb_pio_1_energy_2               = ti.pio_1_energy_2;
    m_xgb_pio_1_dis_1                  = ti.pio_1_dis_1;
    m_xgb_pio_1_dis_2                  = ti.pio_1_dis_2;
    m_xgb_pio_mip_id                   = ti.pio_mip_id;
    m_xgb_stem_dir_flag_single_shower  = ti.stem_dir_flag_single_shower;
    m_xgb_stem_dir_angle               = ti.stem_dir_angle;
    m_xgb_stem_dir_energy              = ti.stem_dir_energy;
    m_xgb_stem_dir_angle1              = ti.stem_dir_angle1;
    m_xgb_stem_dir_angle2              = ti.stem_dir_angle2;
    m_xgb_stem_dir_angle3              = ti.stem_dir_angle3;
    m_xgb_stem_dir_ratio               = ti.stem_dir_ratio;
    m_xgb_br2_num_valid_tracks         = ti.br2_num_valid_tracks;
    m_xgb_br2_n_shower_main_segs       = ti.br2_n_shower_main_segs;
    m_xgb_br2_max_angle                = ti.br2_max_angle;
    m_xgb_br2_sg_length                = ti.br2_sg_length;
    m_xgb_br2_flag_sg_trajectory       = ti.br2_flag_sg_trajectory;
    m_xgb_stem_len_energy              = ti.stem_len_energy;
    m_xgb_stem_len_length              = ti.stem_len_length;
    m_xgb_stem_len_flag_avoid_muon_check = ti.stem_len_flag_avoid_muon_check;
    m_xgb_stem_len_num_daughters       = ti.stem_len_num_daughters;
    m_xgb_stem_len_daughter_length     = ti.stem_len_daughter_length;
    m_xgb_brm_n_mu_segs                = ti.brm_n_mu_segs;
    m_xgb_brm_Ep                       = ti.brm_Ep;
    m_xgb_brm_acc_length               = ti.brm_acc_length;
    m_xgb_brm_shower_total_length      = ti.brm_shower_total_length;
    m_xgb_brm_connected_length         = ti.brm_connected_length;
    m_xgb_brm_n_size                   = ti.brm_n_size;
    m_xgb_brm_n_shower_main_segs       = ti.brm_n_shower_main_segs;
    m_xgb_brm_n_mu_main                = ti.brm_n_mu_main;
    m_xgb_lem_shower_main_length       = ti.lem_shower_main_length;
    m_xgb_lem_n_3seg                   = ti.lem_n_3seg;
    m_xgb_lem_e_charge                 = ti.lem_e_charge;
    m_xgb_lem_e_dQdx                   = ti.lem_e_dQdx;
    m_xgb_lem_shower_num_main_segs     = ti.lem_shower_num_main_segs;
    m_xgb_brm_acc_direct_length        = ti.brm_acc_direct_length;
    m_xgb_stw_1_energy                 = ti.stw_1_energy;
    m_xgb_stw_1_dis                    = ti.stw_1_dis;
    m_xgb_stw_1_dQ_dx                  = ti.stw_1_dQ_dx;
    m_xgb_stw_1_flag_single_shower     = ti.stw_1_flag_single_shower;
    m_xgb_stw_1_n_pi0                  = ti.stw_1_n_pi0;
    m_xgb_stw_1_num_valid_tracks       = ti.stw_1_num_valid_tracks;
    m_xgb_spt_shower_main_length       = ti.spt_shower_main_length;
    m_xgb_spt_shower_total_length      = ti.spt_shower_total_length;
    m_xgb_spt_angle_beam               = ti.spt_angle_beam;
    m_xgb_spt_angle_vertical           = ti.spt_angle_vertical;
    m_xgb_spt_max_dQ_dx                = ti.spt_max_dQ_dx;
    m_xgb_spt_angle_beam_1             = ti.spt_angle_beam_1;
    m_xgb_spt_angle_drift              = ti.spt_angle_drift;
    m_xgb_spt_angle_drift_1            = ti.spt_angle_drift_1;
    m_xgb_spt_num_valid_tracks         = ti.spt_num_valid_tracks;
    m_xgb_spt_n_vtx_segs               = ti.spt_n_vtx_segs;
    m_xgb_spt_max_length               = ti.spt_max_length;
    m_xgb_mip_energy                   = ti.mip_energy;
    m_xgb_mip_n_end_reduction          = ti.mip_n_end_reduction;
    m_xgb_mip_n_first_mip              = ti.mip_n_first_mip;
    m_xgb_mip_n_first_non_mip          = ti.mip_n_first_non_mip;
    m_xgb_mip_n_first_non_mip_1        = ti.mip_n_first_non_mip_1;
    m_xgb_mip_n_first_non_mip_2        = ti.mip_n_first_non_mip_2;
    m_xgb_mip_vec_dQ_dx_0              = ti.mip_vec_dQ_dx_0;
    m_xgb_mip_vec_dQ_dx_1              = ti.mip_vec_dQ_dx_1;
    m_xgb_mip_max_dQ_dx_sample         = ti.mip_max_dQ_dx_sample;
    m_xgb_mip_n_below_threshold        = ti.mip_n_below_threshold;
    m_xgb_mip_n_below_zero             = ti.mip_n_below_zero;
    m_xgb_mip_n_lowest                 = ti.mip_n_lowest;
    m_xgb_mip_n_highest                = ti.mip_n_highest;
    m_xgb_mip_lowest_dQ_dx             = ti.mip_lowest_dQ_dx;
    m_xgb_mip_highest_dQ_dx            = ti.mip_highest_dQ_dx;
    m_xgb_mip_medium_dQ_dx             = ti.mip_medium_dQ_dx;
    m_xgb_mip_stem_length              = ti.mip_stem_length;
    m_xgb_mip_length_main              = ti.mip_length_main;
    m_xgb_mip_length_total             = ti.mip_length_total;
    m_xgb_mip_angle_beam               = ti.mip_angle_beam;
    m_xgb_mip_iso_angle                = ti.mip_iso_angle;
    m_xgb_mip_n_vertex                 = ti.mip_n_vertex;
    m_xgb_mip_n_good_tracks            = ti.mip_n_good_tracks;
    m_xgb_mip_E_indirect_max_energy    = ti.mip_E_indirect_max_energy;
    m_xgb_mip_flag_all_above           = ti.mip_flag_all_above;
    m_xgb_mip_min_dQ_dx_5              = ti.mip_min_dQ_dx_5;
    m_xgb_mip_n_other_vertex           = ti.mip_n_other_vertex;
    m_xgb_mip_n_stem_size              = ti.mip_n_stem_size;
    m_xgb_mip_flag_stem_trajectory     = ti.mip_flag_stem_trajectory;
    m_xgb_mip_min_dis                  = ti.mip_min_dis;
    m_xgb_mip_vec_dQ_dx_2              = ti.mip_vec_dQ_dx_2;
    m_xgb_mip_vec_dQ_dx_3              = ti.mip_vec_dQ_dx_3;
    m_xgb_mip_vec_dQ_dx_4              = ti.mip_vec_dQ_dx_4;
    m_xgb_mip_vec_dQ_dx_5              = ti.mip_vec_dQ_dx_5;
    m_xgb_mip_vec_dQ_dx_6              = ti.mip_vec_dQ_dx_6;
    m_xgb_mip_vec_dQ_dx_7              = ti.mip_vec_dQ_dx_7;
    m_xgb_mip_vec_dQ_dx_8              = ti.mip_vec_dQ_dx_8;
    m_xgb_mip_vec_dQ_dx_9              = ti.mip_vec_dQ_dx_9;
    m_xgb_mip_vec_dQ_dx_10             = ti.mip_vec_dQ_dx_10;
    m_xgb_mip_vec_dQ_dx_11             = ti.mip_vec_dQ_dx_11;
    m_xgb_mip_vec_dQ_dx_12             = ti.mip_vec_dQ_dx_12;
    m_xgb_mip_vec_dQ_dx_13             = ti.mip_vec_dQ_dx_13;
    m_xgb_mip_vec_dQ_dx_14             = ti.mip_vec_dQ_dx_14;
    m_xgb_mip_vec_dQ_dx_15             = ti.mip_vec_dQ_dx_15;
    m_xgb_mip_vec_dQ_dx_16             = ti.mip_vec_dQ_dx_16;
    m_xgb_mip_vec_dQ_dx_17             = ti.mip_vec_dQ_dx_17;
    m_xgb_mip_vec_dQ_dx_18             = ti.mip_vec_dQ_dx_18;
    m_xgb_mip_vec_dQ_dx_19             = ti.mip_vec_dQ_dx_19;
    m_xgb_br3_3_score                  = ti.br3_3_score;
    m_xgb_br3_5_score                  = ti.br3_5_score;
    m_xgb_br3_6_score                  = ti.br3_6_score;
    m_xgb_pio_2_score                  = ti.pio_2_score;
    m_xgb_stw_2_score                  = ti.stw_2_score;
    m_xgb_stw_3_score                  = ti.stw_3_score;
    m_xgb_stw_4_score                  = ti.stw_4_score;
    m_xgb_sig_1_score                  = ti.sig_1_score;
    m_xgb_sig_2_score                  = ti.sig_2_score;
    m_xgb_lol_1_score                  = ti.lol_1_score;
    m_xgb_lol_2_score                  = ti.lol_2_score;
    m_xgb_tro_1_score                  = ti.tro_1_score;
    m_xgb_tro_2_score                  = ti.tro_2_score;
    m_xgb_tro_4_score                  = ti.tro_4_score;
    m_xgb_tro_5_score                  = ti.tro_5_score;
    m_xgb_br4_1_shower_main_length     = ti.br4_1_shower_main_length;
    m_xgb_br4_1_shower_total_length    = ti.br4_1_shower_total_length;
    m_xgb_br4_1_min_dis                = ti.br4_1_min_dis;
    m_xgb_br4_1_energy                 = ti.br4_1_energy;
    m_xgb_br4_1_flag_avoid_muon_check  = ti.br4_1_flag_avoid_muon_check;
    m_xgb_br4_1_n_vtx_segs             = ti.br4_1_n_vtx_segs;
    m_xgb_br4_2_ratio_45               = ti.br4_2_ratio_45;
    m_xgb_br4_2_ratio_35               = ti.br4_2_ratio_35;
    m_xgb_br4_2_ratio_25               = ti.br4_2_ratio_25;
    m_xgb_br4_2_ratio_15               = ti.br4_2_ratio_15;
    m_xgb_br4_2_ratio1_45              = ti.br4_2_ratio1_45;
    m_xgb_br4_2_ratio1_35              = ti.br4_2_ratio1_35;
    m_xgb_br4_2_ratio1_25              = ti.br4_2_ratio1_25;
    m_xgb_br4_2_ratio1_15              = ti.br4_2_ratio1_15;
    m_xgb_br4_2_iso_angle              = ti.br4_2_iso_angle;
    m_xgb_br4_2_iso_angle1             = ti.br4_2_iso_angle1;
    m_xgb_br4_2_angle                  = ti.br4_2_angle;
    m_xgb_tro_3_stem_length            = ti.tro_3_stem_length;
    m_xgb_tro_3_n_muon_segs            = ti.tro_3_n_muon_segs;
    m_xgb_br4_1_n_main_segs            = ti.br4_1_n_main_segs;

    // ---- Step 4: evaluate and transform -------------------------------------
    if (ti.br_filled == 1) {
        float val1 = m_rdr_xgboost->EvaluateMVA("MyBDT");
        // Clamp to avoid division by zero in log-odds transformation
        val1 = std::max(-0.9999f, std::min(0.9999f, val1));
        ti.nue_score = static_cast<float>(std::log10((1.0 + val1) / (1.0 - val1)));
    } else {
        ti.nue_score = -15.f;  // background-like default; matches prototype default_val = -15
    }
}
