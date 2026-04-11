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
    auto resolve = [](const std::string& p) {
        return p.empty() ? p : Persist::resolve(p);
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

void UbooneNueBDTScorer::visit(Clus::Facade::Ensemble& ensemble) const
{
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
    float val = 1e9f;

    float br3_3_v_energy;
    float br3_3_v_angle;
    float br3_3_v_dir_length;
    float br3_3_v_length;

    TMVA::Reader reader;
    reader.AddVariable("br3_3_v_energy",     &br3_3_v_energy);
    reader.AddVariable("br3_3_v_angle",      &br3_3_v_angle);
    reader.AddVariable("br3_3_v_dir_length", &br3_3_v_dir_length);
    reader.AddVariable("br3_3_v_length",     &br3_3_v_length);

    reader.BookMVA("MyBDT", m_br3_3_xml);

    for (size_t i = 0; i != ti.br3_3_v_energy.size(); ++i) {
        br3_3_v_energy     = ti.br3_3_v_energy.at(i);
        br3_3_v_angle      = ti.br3_3_v_angle.at(i);
        br3_3_v_dir_length = ti.br3_3_v_dir_length.at(i);
        br3_3_v_length     = ti.br3_3_v_length.at(i);

        float tmp_val = reader.EvaluateMVA("MyBDT");
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
    float val = 1e9f;

    float br3_5_v_dir_length;
    float br3_5_v_total_length;
    float br3_5_v_flag_avoid_muon_check;
    float br3_5_v_n_seg;
    float br3_5_v_angle;
    float br3_5_v_sg_length;
    float br3_5_v_energy;
    // br3_5_v_n_main_segs is NOT added (commented out in prototype)
    float br3_5_v_n_segs;
    float br3_5_v_shower_main_length;
    float br3_5_v_shower_total_length;

    TMVA::Reader reader;
    reader.AddVariable("br3_5_v_dir_length",            &br3_5_v_dir_length);
    reader.AddVariable("br3_5_v_total_length",          &br3_5_v_total_length);
    reader.AddVariable("br3_5_v_flag_avoid_muon_check", &br3_5_v_flag_avoid_muon_check);
    reader.AddVariable("br3_5_v_n_seg",                 &br3_5_v_n_seg);
    reader.AddVariable("br3_5_v_angle",                 &br3_5_v_angle);
    reader.AddVariable("br3_5_v_sg_length",             &br3_5_v_sg_length);
    reader.AddVariable("br3_5_v_energy",                &br3_5_v_energy);
    reader.AddVariable("br3_5_v_n_segs",                &br3_5_v_n_segs);
    reader.AddVariable("br3_5_v_shower_main_length",    &br3_5_v_shower_main_length);
    reader.AddVariable("br3_5_v_shower_total_length",   &br3_5_v_shower_total_length);

    reader.BookMVA("MyBDT", m_br3_5_xml);

    for (size_t i = 0; i != ti.br3_5_v_dir_length.size(); ++i) {
        br3_5_v_dir_length            = ti.br3_5_v_dir_length.at(i);
        br3_5_v_total_length          = ti.br3_5_v_total_length.at(i);
        br3_5_v_flag_avoid_muon_check = ti.br3_5_v_flag_avoid_muon_check.at(i);
        br3_5_v_n_seg                 = ti.br3_5_v_n_seg.at(i);
        br3_5_v_angle                 = ti.br3_5_v_angle.at(i);
        br3_5_v_sg_length             = ti.br3_5_v_sg_length.at(i);
        br3_5_v_energy                = ti.br3_5_v_energy.at(i);
        br3_5_v_n_segs                = ti.br3_5_v_n_segs.at(i);
        br3_5_v_shower_main_length    = ti.br3_5_v_shower_main_length.at(i);
        br3_5_v_shower_total_length   = ti.br3_5_v_shower_total_length.at(i);

        float tmp_val = reader.EvaluateMVA("MyBDT");
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
    float val = 1e9f;

    float br3_6_v_angle;
    float br3_6_v_angle1;
    float br3_6_v_flag_shower_trajectory;
    float br3_6_v_direct_length;
    float br3_6_v_length;
    float br3_6_v_n_other_vtx_segs;
    float br3_6_v_energy;

    TMVA::Reader reader;
    reader.AddVariable("br3_6_v_angle",                  &br3_6_v_angle);
    reader.AddVariable("br3_6_v_angle1",                 &br3_6_v_angle1);
    reader.AddVariable("br3_6_v_flag_shower_trajectory", &br3_6_v_flag_shower_trajectory);
    reader.AddVariable("br3_6_v_direct_length",          &br3_6_v_direct_length);
    reader.AddVariable("br3_6_v_length",                 &br3_6_v_length);
    reader.AddVariable("br3_6_v_n_other_vtx_segs",       &br3_6_v_n_other_vtx_segs);
    reader.AddVariable("br3_6_v_energy",                 &br3_6_v_energy);

    reader.BookMVA("MyBDT", m_br3_6_xml);

    for (size_t i = 0; i != ti.br3_6_v_angle.size(); ++i) {
        br3_6_v_angle                  = ti.br3_6_v_angle.at(i);
        br3_6_v_angle1                 = ti.br3_6_v_angle1.at(i);
        br3_6_v_flag_shower_trajectory = ti.br3_6_v_flag_shower_trajectory.at(i);
        br3_6_v_direct_length          = ti.br3_6_v_direct_length.at(i);
        br3_6_v_length                 = ti.br3_6_v_length.at(i);
        br3_6_v_n_other_vtx_segs       = ti.br3_6_v_n_other_vtx_segs.at(i);
        br3_6_v_energy                 = ti.br3_6_v_energy.at(i);

        float tmp_val = reader.EvaluateMVA("MyBDT");
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
    float val = 1e9f;

    float pio_2_v_dis2;
    float pio_2_v_angle2;
    float pio_2_v_acc_length;

    TMVA::Reader reader;
    reader.AddVariable("pio_2_v_dis2",       &pio_2_v_dis2);
    reader.AddVariable("pio_2_v_angle2",     &pio_2_v_angle2);
    reader.AddVariable("pio_2_v_acc_length", &pio_2_v_acc_length);
    reader.AddVariable("pio_mip_id",         &ti.pio_mip_id);

    reader.BookMVA("MyBDT", m_pio_2_xml);

    if (ti.pio_filled == 1 && ti.pio_flag_pio == 0) {
        for (size_t i = 0; i != ti.pio_2_v_dis2.size(); ++i) {
            pio_2_v_dis2       = ti.pio_2_v_dis2.at(i);
            pio_2_v_angle2     = ti.pio_2_v_angle2.at(i);
            pio_2_v_acc_length = ti.pio_2_v_acc_length.at(i);

            float tmp_val = reader.EvaluateMVA("MyBDT");
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
    float val = 1e9f;

    float stw_2_v_medium_dQ_dx;
    float stw_2_v_energy;
    float stw_2_v_angle;
    float stw_2_v_dir_length;
    float stw_2_v_max_dQ_dx;

    TMVA::Reader reader;
    reader.AddVariable("stw_2_v_medium_dQ_dx", &stw_2_v_medium_dQ_dx);
    reader.AddVariable("stw_2_v_energy",        &stw_2_v_energy);
    reader.AddVariable("stw_2_v_angle",         &stw_2_v_angle);
    reader.AddVariable("stw_2_v_dir_length",    &stw_2_v_dir_length);
    reader.AddVariable("stw_2_v_max_dQ_dx",     &stw_2_v_max_dQ_dx);

    reader.BookMVA("MyBDT", m_stw_2_xml);

    for (size_t i = 0; i != ti.stw_2_v_medium_dQ_dx.size(); ++i) {
        stw_2_v_medium_dQ_dx = ti.stw_2_v_medium_dQ_dx.at(i);
        stw_2_v_energy        = ti.stw_2_v_energy.at(i);
        stw_2_v_angle         = ti.stw_2_v_angle.at(i);
        stw_2_v_dir_length    = ti.stw_2_v_dir_length.at(i);
        stw_2_v_max_dQ_dx     = ti.stw_2_v_max_dQ_dx.at(i);

        float tmp_val = reader.EvaluateMVA("MyBDT");
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
    float val = 1e9f;

    float stw_3_v_angle;
    float stw_3_v_dir_length;
    float stw_3_v_energy;
    float stw_3_v_medium_dQ_dx;

    TMVA::Reader reader;
    reader.AddVariable("stw_3_v_angle",          &stw_3_v_angle);
    reader.AddVariable("stw_3_v_dir_length",     &stw_3_v_dir_length);
    reader.AddVariable("stw_3_v_energy",         &stw_3_v_energy);
    reader.AddVariable("stw_3_v_medium_dQ_dx",   &stw_3_v_medium_dQ_dx);

    reader.BookMVA("MyBDT", m_stw_3_xml);

    for (size_t i = 0; i != ti.stw_3_v_angle.size(); ++i) {
        stw_3_v_angle        = ti.stw_3_v_angle.at(i);
        stw_3_v_dir_length   = ti.stw_3_v_dir_length.at(i);
        stw_3_v_energy       = ti.stw_3_v_energy.at(i);
        stw_3_v_medium_dQ_dx = ti.stw_3_v_medium_dQ_dx.at(i);

        float tmp_val = reader.EvaluateMVA("MyBDT");
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
    float val = 1e9f;

    float stw_4_v_angle;
    float stw_4_v_dis;
    float stw_4_v_energy;

    TMVA::Reader reader;
    reader.AddVariable("stw_4_v_angle",  &stw_4_v_angle);
    reader.AddVariable("stw_4_v_dis",    &stw_4_v_dis);
    reader.AddVariable("stw_4_v_energy", &stw_4_v_energy);

    reader.BookMVA("MyBDT", m_stw_4_xml);

    for (size_t i = 0; i != ti.stw_4_v_angle.size(); ++i) {
        stw_4_v_angle  = ti.stw_4_v_angle.at(i);
        stw_4_v_dis    = ti.stw_4_v_dis.at(i);
        stw_4_v_energy = ti.stw_4_v_energy.at(i);

        float tmp_val = reader.EvaluateMVA("MyBDT");
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
    float val = 1e9f;

    float sig_1_v_angle;
    float sig_1_v_flag_single_shower;
    float sig_1_v_energy;
    float sig_1_v_energy_1;

    TMVA::Reader reader;
    reader.AddVariable("sig_1_v_angle",               &sig_1_v_angle);
    reader.AddVariable("sig_1_v_flag_single_shower",  &sig_1_v_flag_single_shower);
    reader.AddVariable("sig_1_v_energy",              &sig_1_v_energy);
    reader.AddVariable("sig_1_v_energy_1",            &sig_1_v_energy_1);

    reader.BookMVA("MyBDT", m_sig_1_xml);

    for (size_t i = 0; i != ti.sig_1_v_angle.size(); ++i) {
        sig_1_v_angle              = ti.sig_1_v_angle.at(i);
        sig_1_v_flag_single_shower = ti.sig_1_v_flag_single_shower.at(i);
        sig_1_v_energy             = ti.sig_1_v_energy.at(i);
        sig_1_v_energy_1           = ti.sig_1_v_energy_1.at(i);

        float tmp_val = reader.EvaluateMVA("MyBDT");
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
    float val = 1e9f;

    float sig_2_v_energy;
    float sig_2_v_shower_angle;
    float sig_2_v_flag_single_shower;
    float sig_2_v_medium_dQ_dx;
    float sig_2_v_start_dQ_dx;

    TMVA::Reader reader;
    reader.AddVariable("sig_2_v_energy",              &sig_2_v_energy);
    reader.AddVariable("sig_2_v_shower_angle",        &sig_2_v_shower_angle);
    reader.AddVariable("sig_2_v_flag_single_shower",  &sig_2_v_flag_single_shower);
    reader.AddVariable("sig_2_v_medium_dQ_dx",        &sig_2_v_medium_dQ_dx);
    reader.AddVariable("sig_2_v_start_dQ_dx",         &sig_2_v_start_dQ_dx);

    reader.BookMVA("MyBDT", m_sig_2_xml);

    for (size_t i = 0; i != ti.sig_2_v_energy.size(); ++i) {
        sig_2_v_energy             = ti.sig_2_v_energy.at(i);
        sig_2_v_shower_angle       = ti.sig_2_v_shower_angle.at(i);
        sig_2_v_flag_single_shower = ti.sig_2_v_flag_single_shower.at(i);
        sig_2_v_medium_dQ_dx       = ti.sig_2_v_medium_dQ_dx.at(i);
        sig_2_v_start_dQ_dx        = ti.sig_2_v_start_dQ_dx.at(i);

        float tmp_val = reader.EvaluateMVA("MyBDT");
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
    float val = 1e9f;

    float lol_1_v_energy;
    float lol_1_v_vtx_n_segs;
    float lol_1_v_nseg;
    float lol_1_v_angle;

    TMVA::Reader reader;
    reader.AddVariable("lol_1_v_energy",      &lol_1_v_energy);
    reader.AddVariable("lol_1_v_vtx_n_segs",  &lol_1_v_vtx_n_segs);
    reader.AddVariable("lol_1_v_nseg",        &lol_1_v_nseg);
    reader.AddVariable("lol_1_v_angle",       &lol_1_v_angle);

    reader.BookMVA("MyBDT", m_lol_1_xml);

    for (size_t i = 0; i != ti.lol_1_v_energy.size(); ++i) {
        lol_1_v_energy     = ti.lol_1_v_energy.at(i);
        lol_1_v_vtx_n_segs = ti.lol_1_v_vtx_n_segs.at(i);
        lol_1_v_nseg       = ti.lol_1_v_nseg.at(i);
        lol_1_v_angle      = ti.lol_1_v_angle.at(i);

        float tmp_val = reader.EvaluateMVA("MyBDT");
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
    float val = 1e9f;

    float lol_2_v_length;
    float lol_2_v_angle;
    float lol_2_v_type;
    float lol_2_v_vtx_n_segs;
    float lol_2_v_energy;
    float lol_2_v_shower_main_length;
    float lol_2_v_flag_dir_weak;

    TMVA::Reader reader;
    reader.AddVariable("lol_2_v_length",             &lol_2_v_length);
    reader.AddVariable("lol_2_v_angle",              &lol_2_v_angle);
    reader.AddVariable("lol_2_v_type",               &lol_2_v_type);
    reader.AddVariable("lol_2_v_vtx_n_segs",         &lol_2_v_vtx_n_segs);
    reader.AddVariable("lol_2_v_energy",             &lol_2_v_energy);
    reader.AddVariable("lol_2_v_shower_main_length", &lol_2_v_shower_main_length);
    reader.AddVariable("lol_2_v_flag_dir_weak",      &lol_2_v_flag_dir_weak);

    reader.BookMVA("MyBDT", m_lol_2_xml);

    for (size_t i = 0; i != ti.lol_2_v_length.size(); ++i) {
        lol_2_v_length             = ti.lol_2_v_length.at(i);
        lol_2_v_angle              = ti.lol_2_v_angle.at(i);
        lol_2_v_type               = ti.lol_2_v_type.at(i);
        lol_2_v_vtx_n_segs         = ti.lol_2_v_vtx_n_segs.at(i);
        lol_2_v_energy             = ti.lol_2_v_energy.at(i);
        lol_2_v_shower_main_length = ti.lol_2_v_shower_main_length.at(i);
        lol_2_v_flag_dir_weak      = ti.lol_2_v_flag_dir_weak.at(i);

        float tmp_val = reader.EvaluateMVA("MyBDT");
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
    float val = 1e9f;

    float tro_1_v_particle_type;
    float tro_1_v_flag_dir_weak;
    float tro_1_v_min_dis;
    float tro_1_v_sg1_length;
    float tro_1_v_shower_main_length;
    float tro_1_v_max_n_vtx_segs;
    float tro_1_v_tmp_length;
    float tro_1_v_medium_dQ_dx;
    float tro_1_v_dQ_dx_cut;
    float tro_1_v_flag_shower_topology;

    TMVA::Reader reader;
    reader.AddVariable("tro_1_v_particle_type",       &tro_1_v_particle_type);
    reader.AddVariable("tro_1_v_flag_dir_weak",       &tro_1_v_flag_dir_weak);
    reader.AddVariable("tro_1_v_min_dis",             &tro_1_v_min_dis);
    reader.AddVariable("tro_1_v_sg1_length",          &tro_1_v_sg1_length);
    reader.AddVariable("tro_1_v_shower_main_length",  &tro_1_v_shower_main_length);
    reader.AddVariable("tro_1_v_max_n_vtx_segs",      &tro_1_v_max_n_vtx_segs);
    reader.AddVariable("tro_1_v_tmp_length",          &tro_1_v_tmp_length);
    reader.AddVariable("tro_1_v_medium_dQ_dx",        &tro_1_v_medium_dQ_dx);
    reader.AddVariable("tro_1_v_dQ_dx_cut",           &tro_1_v_dQ_dx_cut);
    reader.AddVariable("tro_1_v_flag_shower_topology",&tro_1_v_flag_shower_topology);

    reader.BookMVA("MyBDT", m_tro_1_xml);

    for (size_t i = 0; i != ti.tro_1_v_particle_type.size(); ++i) {
        tro_1_v_particle_type      = ti.tro_1_v_particle_type.at(i);
        tro_1_v_flag_dir_weak      = ti.tro_1_v_flag_dir_weak.at(i);
        tro_1_v_min_dis            = ti.tro_1_v_min_dis.at(i);
        tro_1_v_sg1_length         = ti.tro_1_v_sg1_length.at(i);
        tro_1_v_shower_main_length = ti.tro_1_v_shower_main_length.at(i);
        tro_1_v_max_n_vtx_segs     = ti.tro_1_v_max_n_vtx_segs.at(i);
        tro_1_v_tmp_length         = ti.tro_1_v_tmp_length.at(i);
        tro_1_v_medium_dQ_dx       = ti.tro_1_v_medium_dQ_dx.at(i);
        tro_1_v_dQ_dx_cut          = ti.tro_1_v_dQ_dx_cut.at(i);
        tro_1_v_flag_shower_topology = ti.tro_1_v_flag_shower_topology.at(i);

        float tmp_val = reader.EvaluateMVA("MyBDT");
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
    float val = 1e9f;

    float tro_2_v_energy;
    float tro_2_v_stem_length;
    float tro_2_v_iso_angle;
    float tro_2_v_max_length;
    float tro_2_v_angle;

    TMVA::Reader reader;
    reader.AddVariable("tro_2_v_energy",      &tro_2_v_energy);
    reader.AddVariable("tro_2_v_stem_length", &tro_2_v_stem_length);
    reader.AddVariable("tro_2_v_iso_angle",   &tro_2_v_iso_angle);
    reader.AddVariable("tro_2_v_max_length",  &tro_2_v_max_length);
    reader.AddVariable("tro_2_v_angle",       &tro_2_v_angle);

    reader.BookMVA("MyBDT", m_tro_2_xml);

    for (size_t i = 0; i != ti.tro_2_v_energy.size(); ++i) {
        tro_2_v_energy      = ti.tro_2_v_energy.at(i);
        tro_2_v_stem_length = ti.tro_2_v_stem_length.at(i);
        tro_2_v_iso_angle   = ti.tro_2_v_iso_angle.at(i);
        tro_2_v_max_length  = ti.tro_2_v_max_length.at(i);
        tro_2_v_angle       = ti.tro_2_v_angle.at(i);

        float tmp_val = reader.EvaluateMVA("MyBDT");
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
    float val = 1e9f;

    float tro_4_v_dir2_mag;
    float tro_4_v_angle;
    float tro_4_v_angle1;
    float tro_4_v_angle2;
    float tro_4_v_length;
    float tro_4_v_length1;
    float tro_4_v_medium_dQ_dx;
    float tro_4_v_end_dQ_dx;
    float tro_4_v_energy;
    float tro_4_v_shower_main_length;
    float tro_4_v_flag_shower_trajectory;

    TMVA::Reader reader;
    reader.AddVariable("tro_4_v_dir2_mag",              &tro_4_v_dir2_mag);
    reader.AddVariable("tro_4_v_angle",                 &tro_4_v_angle);
    reader.AddVariable("tro_4_v_angle1",                &tro_4_v_angle1);
    reader.AddVariable("tro_4_v_angle2",                &tro_4_v_angle2);
    reader.AddVariable("tro_4_v_length",                &tro_4_v_length);
    reader.AddVariable("tro_4_v_length1",               &tro_4_v_length1);
    reader.AddVariable("tro_4_v_medium_dQ_dx",          &tro_4_v_medium_dQ_dx);
    reader.AddVariable("tro_4_v_end_dQ_dx",             &tro_4_v_end_dQ_dx);
    reader.AddVariable("tro_4_v_energy",                &tro_4_v_energy);
    reader.AddVariable("tro_4_v_shower_main_length",    &tro_4_v_shower_main_length);
    reader.AddVariable("tro_4_v_flag_shower_trajectory",&tro_4_v_flag_shower_trajectory);

    reader.BookMVA("MyBDT", m_tro_4_xml);

    for (size_t i = 0; i != ti.tro_4_v_dir2_mag.size(); ++i) {
        tro_4_v_dir2_mag              = ti.tro_4_v_dir2_mag.at(i);
        tro_4_v_angle                 = ti.tro_4_v_angle.at(i);
        tro_4_v_angle1                = ti.tro_4_v_angle1.at(i);
        tro_4_v_angle2                = ti.tro_4_v_angle2.at(i);
        tro_4_v_length                = ti.tro_4_v_length.at(i);
        tro_4_v_length1               = ti.tro_4_v_length1.at(i);
        tro_4_v_medium_dQ_dx          = ti.tro_4_v_medium_dQ_dx.at(i);
        tro_4_v_end_dQ_dx             = ti.tro_4_v_end_dQ_dx.at(i);
        tro_4_v_energy                = ti.tro_4_v_energy.at(i);
        tro_4_v_shower_main_length    = ti.tro_4_v_shower_main_length.at(i);
        tro_4_v_flag_shower_trajectory= ti.tro_4_v_flag_shower_trajectory.at(i);

        float tmp_val = reader.EvaluateMVA("MyBDT");
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
    float val = 1e9f;

    float tro_5_v_max_angle;
    float tro_5_v_min_angle;
    float tro_5_v_max_length;
    float tro_5_v_iso_angle;
    float tro_5_v_n_vtx_segs;
    float tro_5_v_min_count;
    float tro_5_v_max_count;
    float tro_5_v_energy;

    TMVA::Reader reader;
    reader.AddVariable("tro_5_v_max_angle",   &tro_5_v_max_angle);
    reader.AddVariable("tro_5_v_min_angle",   &tro_5_v_min_angle);
    reader.AddVariable("tro_5_v_max_length",  &tro_5_v_max_length);
    reader.AddVariable("tro_5_v_iso_angle",   &tro_5_v_iso_angle);
    reader.AddVariable("tro_5_v_n_vtx_segs", &tro_5_v_n_vtx_segs);
    reader.AddVariable("tro_5_v_min_count",   &tro_5_v_min_count);
    reader.AddVariable("tro_5_v_max_count",   &tro_5_v_max_count);
    reader.AddVariable("tro_5_v_energy",      &tro_5_v_energy);

    reader.BookMVA("MyBDT", m_tro_5_xml);

    for (size_t i = 0; i != ti.tro_5_v_max_angle.size(); ++i) {
        tro_5_v_max_angle   = ti.tro_5_v_max_angle.at(i);
        tro_5_v_min_angle   = ti.tro_5_v_min_angle.at(i);
        tro_5_v_max_length  = ti.tro_5_v_max_length.at(i);
        tro_5_v_iso_angle   = ti.tro_5_v_iso_angle.at(i);
        tro_5_v_n_vtx_segs  = ti.tro_5_v_n_vtx_segs.at(i);
        tro_5_v_min_count   = ti.tro_5_v_min_count.at(i);
        tro_5_v_max_count   = ti.tro_5_v_max_count.at(i);
        tro_5_v_energy      = ti.tro_5_v_energy.at(i);

        float tmp_val = reader.EvaluateMVA("MyBDT");
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
// Step 3: Build a TMVA reader with ~160 features, evaluate if br_filled==1.
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

    // ---- Step 3: build final XGBoost TMVA reader ----------------------------
    // Variable order must match the training XML exactly.
    float nue_kine_reco_Enu = static_cast<float>(ki.kine_reco_Enu);

    TMVA::Reader reader;
    reader.AddVariable("match_isFC",                   &ti.match_isFC);
    reader.AddVariable("kine_reco_Enu",                &nue_kine_reco_Enu);
    reader.AddVariable("cme_mu_energy",                &ti.cme_mu_energy);
    reader.AddVariable("cme_energy",                   &ti.cme_energy);
    reader.AddVariable("cme_mu_length",                &ti.cme_mu_length);
    reader.AddVariable("cme_length",                   &ti.cme_length);
    reader.AddVariable("cme_angle_beam",               &ti.cme_angle_beam);
    reader.AddVariable("anc_angle",                    &ti.anc_angle);
    reader.AddVariable("anc_max_angle",                &ti.anc_max_angle);
    reader.AddVariable("anc_max_length",               &ti.anc_max_length);
    reader.AddVariable("anc_acc_forward_length",       &ti.anc_acc_forward_length);
    reader.AddVariable("anc_acc_backward_length",      &ti.anc_acc_backward_length);
    reader.AddVariable("anc_acc_forward_length1",      &ti.anc_acc_forward_length1);
    reader.AddVariable("anc_shower_main_length",       &ti.anc_shower_main_length);
    reader.AddVariable("anc_shower_total_length",      &ti.anc_shower_total_length);
    reader.AddVariable("anc_flag_main_outside",        &ti.anc_flag_main_outside);
    reader.AddVariable("gap_flag_prolong_u",           &ti.gap_flag_prolong_u);
    reader.AddVariable("gap_flag_prolong_v",           &ti.gap_flag_prolong_v);
    reader.AddVariable("gap_flag_prolong_w",           &ti.gap_flag_prolong_w);
    reader.AddVariable("gap_flag_parallel",            &ti.gap_flag_parallel);
    reader.AddVariable("gap_n_points",                 &ti.gap_n_points);
    reader.AddVariable("gap_n_bad",                    &ti.gap_n_bad);
    reader.AddVariable("gap_energy",                   &ti.gap_energy);
    reader.AddVariable("gap_num_valid_tracks",         &ti.gap_num_valid_tracks);
    reader.AddVariable("gap_flag_single_shower",       &ti.gap_flag_single_shower);
    reader.AddVariable("hol_1_n_valid_tracks",         &ti.hol_1_n_valid_tracks);
    reader.AddVariable("hol_1_min_angle",              &ti.hol_1_min_angle);
    reader.AddVariable("hol_1_energy",                 &ti.hol_1_energy);
    reader.AddVariable("hol_1_min_length",             &ti.hol_1_min_length);
    reader.AddVariable("hol_2_min_angle",              &ti.hol_2_min_angle);
    reader.AddVariable("hol_2_medium_dQ_dx",           &ti.hol_2_medium_dQ_dx);
    reader.AddVariable("hol_2_ncount",                 &ti.hol_2_ncount);
    reader.AddVariable("lol_3_angle_beam",             &ti.lol_3_angle_beam);
    reader.AddVariable("lol_3_n_valid_tracks",         &ti.lol_3_n_valid_tracks);
    reader.AddVariable("lol_3_min_angle",              &ti.lol_3_min_angle);
    reader.AddVariable("lol_3_vtx_n_segs",             &ti.lol_3_vtx_n_segs);
    reader.AddVariable("lol_3_shower_main_length",     &ti.lol_3_shower_main_length);
    reader.AddVariable("lol_3_n_out",                  &ti.lol_3_n_out);
    reader.AddVariable("lol_3_n_sum",                  &ti.lol_3_n_sum);
    reader.AddVariable("hol_1_flag_all_shower",        &ti.hol_1_flag_all_shower);
    reader.AddVariable("mgo_energy",                   &ti.mgo_energy);
    reader.AddVariable("mgo_max_energy",               &ti.mgo_max_energy);
    reader.AddVariable("mgo_total_energy",             &ti.mgo_total_energy);
    reader.AddVariable("mgo_n_showers",                &ti.mgo_n_showers);
    reader.AddVariable("mgo_max_energy_1",             &ti.mgo_max_energy_1);
    reader.AddVariable("mgo_max_energy_2",             &ti.mgo_max_energy_2);
    reader.AddVariable("mgo_total_other_energy",       &ti.mgo_total_other_energy);
    reader.AddVariable("mgo_n_total_showers",          &ti.mgo_n_total_showers);
    reader.AddVariable("mgo_total_other_energy_1",     &ti.mgo_total_other_energy_1);
    reader.AddVariable("mgt_flag_single_shower",       &ti.mgt_flag_single_shower);
    reader.AddVariable("mgt_max_energy",               &ti.mgt_max_energy);
    reader.AddVariable("mgt_total_other_energy",       &ti.mgt_total_other_energy);
    reader.AddVariable("mgt_max_energy_1",             &ti.mgt_max_energy_1);
    reader.AddVariable("mgt_e_indirect_max_energy",    &ti.mgt_e_indirect_max_energy);
    reader.AddVariable("mgt_e_direct_max_energy",      &ti.mgt_e_direct_max_energy);
    reader.AddVariable("mgt_n_direct_showers",         &ti.mgt_n_direct_showers);
    reader.AddVariable("mgt_e_direct_total_energy",    &ti.mgt_e_direct_total_energy);
    reader.AddVariable("mgt_flag_indirect_max_pio",    &ti.mgt_flag_indirect_max_pio);
    reader.AddVariable("mgt_e_indirect_total_energy",  &ti.mgt_e_indirect_total_energy);
    reader.AddVariable("mip_quality_energy",           &ti.mip_quality_energy);
    reader.AddVariable("mip_quality_overlap",          &ti.mip_quality_overlap);
    reader.AddVariable("mip_quality_n_showers",        &ti.mip_quality_n_showers);
    reader.AddVariable("mip_quality_n_tracks",         &ti.mip_quality_n_tracks);
    reader.AddVariable("mip_quality_flag_inside_pi0",  &ti.mip_quality_flag_inside_pi0);
    reader.AddVariable("mip_quality_n_pi0_showers",    &ti.mip_quality_n_pi0_showers);
    reader.AddVariable("mip_quality_shortest_length",  &ti.mip_quality_shortest_length);
    reader.AddVariable("mip_quality_acc_length",       &ti.mip_quality_acc_length);
    reader.AddVariable("mip_quality_shortest_angle",   &ti.mip_quality_shortest_angle);
    reader.AddVariable("mip_quality_flag_proton",      &ti.mip_quality_flag_proton);
    reader.AddVariable("br1_1_shower_type",            &ti.br1_1_shower_type);
    reader.AddVariable("br1_1_vtx_n_segs",             &ti.br1_1_vtx_n_segs);
    reader.AddVariable("br1_1_energy",                 &ti.br1_1_energy);
    reader.AddVariable("br1_1_n_segs",                 &ti.br1_1_n_segs);
    reader.AddVariable("br1_1_flag_sg_topology",       &ti.br1_1_flag_sg_topology);
    reader.AddVariable("br1_1_flag_sg_trajectory",     &ti.br1_1_flag_sg_trajectory);
    reader.AddVariable("br1_1_sg_length",              &ti.br1_1_sg_length);
    reader.AddVariable("br1_2_n_connected",            &ti.br1_2_n_connected);
    reader.AddVariable("br1_2_max_length",             &ti.br1_2_max_length);
    reader.AddVariable("br1_2_n_connected_1",          &ti.br1_2_n_connected_1);
    reader.AddVariable("br1_2_n_shower_segs",          &ti.br1_2_n_shower_segs);
    reader.AddVariable("br1_2_max_length_ratio",       &ti.br1_2_max_length_ratio);
    reader.AddVariable("br1_2_shower_length",          &ti.br1_2_shower_length);
    reader.AddVariable("br1_3_n_connected_p",          &ti.br1_3_n_connected_p);
    reader.AddVariable("br1_3_max_length_p",           &ti.br1_3_max_length_p);
    reader.AddVariable("br1_3_n_shower_main_segs",     &ti.br1_3_n_shower_main_segs);
    reader.AddVariable("br3_1_energy",                 &ti.br3_1_energy);
    reader.AddVariable("br3_1_n_shower_segments",      &ti.br3_1_n_shower_segments);
    reader.AddVariable("br3_1_sg_flag_trajectory",     &ti.br3_1_sg_flag_trajectory);
    reader.AddVariable("br3_1_sg_direct_length",       &ti.br3_1_sg_direct_length);
    reader.AddVariable("br3_1_sg_length",              &ti.br3_1_sg_length);
    reader.AddVariable("br3_1_total_main_length",      &ti.br3_1_total_main_length);
    reader.AddVariable("br3_1_total_length",           &ti.br3_1_total_length);
    reader.AddVariable("br3_1_iso_angle",              &ti.br3_1_iso_angle);
    reader.AddVariable("br3_1_sg_flag_topology",       &ti.br3_1_sg_flag_topology);
    reader.AddVariable("br3_2_n_ele",                  &ti.br3_2_n_ele);
    reader.AddVariable("br3_2_n_other",                &ti.br3_2_n_other);
    reader.AddVariable("br3_2_other_fid",              &ti.br3_2_other_fid);
    reader.AddVariable("br3_4_acc_length",             &ti.br3_4_acc_length);
    reader.AddVariable("br3_4_total_length",           &ti.br3_4_total_length);
    reader.AddVariable("br3_7_min_angle",              &ti.br3_7_min_angle);
    reader.AddVariable("br3_8_max_dQ_dx",              &ti.br3_8_max_dQ_dx);
    reader.AddVariable("br3_8_n_main_segs",            &ti.br3_8_n_main_segs);
    reader.AddVariable("vis_1_n_vtx_segs",             &ti.vis_1_n_vtx_segs);
    reader.AddVariable("vis_1_energy",                 &ti.vis_1_energy);
    reader.AddVariable("vis_1_num_good_tracks",        &ti.vis_1_num_good_tracks);
    reader.AddVariable("vis_1_max_angle",              &ti.vis_1_max_angle);
    reader.AddVariable("vis_1_max_shower_angle",       &ti.vis_1_max_shower_angle);
    reader.AddVariable("vis_1_tmp_length1",            &ti.vis_1_tmp_length1);
    reader.AddVariable("vis_1_tmp_length2",            &ti.vis_1_tmp_length2);
    reader.AddVariable("vis_2_n_vtx_segs",             &ti.vis_2_n_vtx_segs);
    reader.AddVariable("vis_2_min_angle",              &ti.vis_2_min_angle);
    reader.AddVariable("vis_2_min_weak_track",         &ti.vis_2_min_weak_track);
    reader.AddVariable("vis_2_angle_beam",             &ti.vis_2_angle_beam);
    reader.AddVariable("vis_2_min_angle1",             &ti.vis_2_min_angle1);
    reader.AddVariable("vis_2_iso_angle1",             &ti.vis_2_iso_angle1);
    reader.AddVariable("vis_2_min_medium_dQ_dx",       &ti.vis_2_min_medium_dQ_dx);
    reader.AddVariable("vis_2_min_length",             &ti.vis_2_min_length);
    reader.AddVariable("vis_2_sg_length",              &ti.vis_2_sg_length);
    reader.AddVariable("vis_2_max_angle",              &ti.vis_2_max_angle);
    reader.AddVariable("vis_2_max_weak_track",         &ti.vis_2_max_weak_track);
    reader.AddVariable("pio_1_mass",                   &ti.pio_1_mass);
    reader.AddVariable("pio_1_pio_type",               &ti.pio_1_pio_type);
    reader.AddVariable("pio_1_energy_1",               &ti.pio_1_energy_1);
    reader.AddVariable("pio_1_energy_2",               &ti.pio_1_energy_2);
    reader.AddVariable("pio_1_dis_1",                  &ti.pio_1_dis_1);
    reader.AddVariable("pio_1_dis_2",                  &ti.pio_1_dis_2);
    reader.AddVariable("pio_mip_id",                   &ti.pio_mip_id);
    reader.AddVariable("stem_dir_flag_single_shower",  &ti.stem_dir_flag_single_shower);
    reader.AddVariable("stem_dir_angle",               &ti.stem_dir_angle);
    reader.AddVariable("stem_dir_energy",              &ti.stem_dir_energy);
    reader.AddVariable("stem_dir_angle1",              &ti.stem_dir_angle1);
    reader.AddVariable("stem_dir_angle2",              &ti.stem_dir_angle2);
    reader.AddVariable("stem_dir_angle3",              &ti.stem_dir_angle3);
    reader.AddVariable("stem_dir_ratio",               &ti.stem_dir_ratio);
    reader.AddVariable("br2_num_valid_tracks",         &ti.br2_num_valid_tracks);
    reader.AddVariable("br2_n_shower_main_segs",       &ti.br2_n_shower_main_segs);
    reader.AddVariable("br2_max_angle",                &ti.br2_max_angle);
    reader.AddVariable("br2_sg_length",                &ti.br2_sg_length);
    reader.AddVariable("br2_flag_sg_trajectory",       &ti.br2_flag_sg_trajectory);
    reader.AddVariable("stem_len_energy",              &ti.stem_len_energy);
    reader.AddVariable("stem_len_length",              &ti.stem_len_length);
    reader.AddVariable("stem_len_flag_avoid_muon_check",&ti.stem_len_flag_avoid_muon_check);
    reader.AddVariable("stem_len_num_daughters",       &ti.stem_len_num_daughters);
    reader.AddVariable("stem_len_daughter_length",     &ti.stem_len_daughter_length);
    reader.AddVariable("brm_n_mu_segs",                &ti.brm_n_mu_segs);
    reader.AddVariable("brm_Ep",                       &ti.brm_Ep);
    reader.AddVariable("brm_acc_length",               &ti.brm_acc_length);
    reader.AddVariable("brm_shower_total_length",      &ti.brm_shower_total_length);
    reader.AddVariable("brm_connected_length",         &ti.brm_connected_length);
    reader.AddVariable("brm_n_size",                   &ti.brm_n_size);
    reader.AddVariable("brm_n_shower_main_segs",       &ti.brm_n_shower_main_segs);
    reader.AddVariable("brm_n_mu_main",                &ti.brm_n_mu_main);
    reader.AddVariable("lem_shower_main_length",       &ti.lem_shower_main_length);
    reader.AddVariable("lem_n_3seg",                   &ti.lem_n_3seg);
    reader.AddVariable("lem_e_charge",                 &ti.lem_e_charge);
    reader.AddVariable("lem_e_dQdx",                   &ti.lem_e_dQdx);
    reader.AddVariable("lem_shower_num_main_segs",     &ti.lem_shower_num_main_segs);
    reader.AddVariable("brm_acc_direct_length",        &ti.brm_acc_direct_length);
    reader.AddVariable("stw_1_energy",                 &ti.stw_1_energy);
    reader.AddVariable("stw_1_dis",                    &ti.stw_1_dis);
    reader.AddVariable("stw_1_dQ_dx",                  &ti.stw_1_dQ_dx);
    reader.AddVariable("stw_1_flag_single_shower",     &ti.stw_1_flag_single_shower);
    reader.AddVariable("stw_1_n_pi0",                  &ti.stw_1_n_pi0);
    reader.AddVariable("stw_1_num_valid_tracks",       &ti.stw_1_num_valid_tracks);
    reader.AddVariable("spt_shower_main_length",       &ti.spt_shower_main_length);
    reader.AddVariable("spt_shower_total_length",      &ti.spt_shower_total_length);
    reader.AddVariable("spt_angle_beam",               &ti.spt_angle_beam);
    reader.AddVariable("spt_angle_vertical",           &ti.spt_angle_vertical);
    reader.AddVariable("spt_max_dQ_dx",                &ti.spt_max_dQ_dx);
    reader.AddVariable("spt_angle_beam_1",             &ti.spt_angle_beam_1);
    reader.AddVariable("spt_angle_drift",              &ti.spt_angle_drift);
    reader.AddVariable("spt_angle_drift_1",            &ti.spt_angle_drift_1);
    reader.AddVariable("spt_num_valid_tracks",         &ti.spt_num_valid_tracks);
    reader.AddVariable("spt_n_vtx_segs",               &ti.spt_n_vtx_segs);
    reader.AddVariable("spt_max_length",               &ti.spt_max_length);
    reader.AddVariable("mip_energy",                   &ti.mip_energy);
    reader.AddVariable("mip_n_end_reduction",          &ti.mip_n_end_reduction);
    reader.AddVariable("mip_n_first_mip",              &ti.mip_n_first_mip);
    reader.AddVariable("mip_n_first_non_mip",          &ti.mip_n_first_non_mip);
    reader.AddVariable("mip_n_first_non_mip_1",        &ti.mip_n_first_non_mip_1);
    reader.AddVariable("mip_n_first_non_mip_2",        &ti.mip_n_first_non_mip_2);
    reader.AddVariable("mip_vec_dQ_dx_0",              &ti.mip_vec_dQ_dx_0);
    reader.AddVariable("mip_vec_dQ_dx_1",              &ti.mip_vec_dQ_dx_1);
    reader.AddVariable("mip_max_dQ_dx_sample",         &ti.mip_max_dQ_dx_sample);
    reader.AddVariable("mip_n_below_threshold",        &ti.mip_n_below_threshold);
    reader.AddVariable("mip_n_below_zero",             &ti.mip_n_below_zero);
    reader.AddVariable("mip_n_lowest",                 &ti.mip_n_lowest);
    reader.AddVariable("mip_n_highest",                &ti.mip_n_highest);
    reader.AddVariable("mip_lowest_dQ_dx",             &ti.mip_lowest_dQ_dx);
    reader.AddVariable("mip_highest_dQ_dx",            &ti.mip_highest_dQ_dx);
    reader.AddVariable("mip_medium_dQ_dx",             &ti.mip_medium_dQ_dx);
    reader.AddVariable("mip_stem_length",              &ti.mip_stem_length);
    reader.AddVariable("mip_length_main",              &ti.mip_length_main);
    reader.AddVariable("mip_length_total",             &ti.mip_length_total);
    reader.AddVariable("mip_angle_beam",               &ti.mip_angle_beam);
    reader.AddVariable("mip_iso_angle",                &ti.mip_iso_angle);
    reader.AddVariable("mip_n_vertex",                 &ti.mip_n_vertex);
    reader.AddVariable("mip_n_good_tracks",            &ti.mip_n_good_tracks);
    reader.AddVariable("mip_E_indirect_max_energy",    &ti.mip_E_indirect_max_energy);
    reader.AddVariable("mip_flag_all_above",           &ti.mip_flag_all_above);
    reader.AddVariable("mip_min_dQ_dx_5",              &ti.mip_min_dQ_dx_5);
    reader.AddVariable("mip_n_other_vertex",           &ti.mip_n_other_vertex);
    reader.AddVariable("mip_n_stem_size",              &ti.mip_n_stem_size);
    reader.AddVariable("mip_flag_stem_trajectory",     &ti.mip_flag_stem_trajectory);
    reader.AddVariable("mip_min_dis",                  &ti.mip_min_dis);
    reader.AddVariable("mip_vec_dQ_dx_2",              &ti.mip_vec_dQ_dx_2);
    reader.AddVariable("mip_vec_dQ_dx_3",              &ti.mip_vec_dQ_dx_3);
    reader.AddVariable("mip_vec_dQ_dx_4",              &ti.mip_vec_dQ_dx_4);
    reader.AddVariable("mip_vec_dQ_dx_5",              &ti.mip_vec_dQ_dx_5);
    reader.AddVariable("mip_vec_dQ_dx_6",              &ti.mip_vec_dQ_dx_6);
    reader.AddVariable("mip_vec_dQ_dx_7",              &ti.mip_vec_dQ_dx_7);
    reader.AddVariable("mip_vec_dQ_dx_8",              &ti.mip_vec_dQ_dx_8);
    reader.AddVariable("mip_vec_dQ_dx_9",              &ti.mip_vec_dQ_dx_9);
    reader.AddVariable("mip_vec_dQ_dx_10",             &ti.mip_vec_dQ_dx_10);
    reader.AddVariable("mip_vec_dQ_dx_11",             &ti.mip_vec_dQ_dx_11);
    reader.AddVariable("mip_vec_dQ_dx_12",             &ti.mip_vec_dQ_dx_12);
    reader.AddVariable("mip_vec_dQ_dx_13",             &ti.mip_vec_dQ_dx_13);
    reader.AddVariable("mip_vec_dQ_dx_14",             &ti.mip_vec_dQ_dx_14);
    reader.AddVariable("mip_vec_dQ_dx_15",             &ti.mip_vec_dQ_dx_15);
    reader.AddVariable("mip_vec_dQ_dx_16",             &ti.mip_vec_dQ_dx_16);
    reader.AddVariable("mip_vec_dQ_dx_17",             &ti.mip_vec_dQ_dx_17);
    reader.AddVariable("mip_vec_dQ_dx_18",             &ti.mip_vec_dQ_dx_18);
    reader.AddVariable("mip_vec_dQ_dx_19",             &ti.mip_vec_dQ_dx_19);
    reader.AddVariable("br3_3_score",                  &ti.br3_3_score);
    reader.AddVariable("br3_5_score",                  &ti.br3_5_score);
    reader.AddVariable("br3_6_score",                  &ti.br3_6_score);
    reader.AddVariable("pio_2_score",                  &ti.pio_2_score);
    reader.AddVariable("stw_2_score",                  &ti.stw_2_score);
    reader.AddVariable("stw_3_score",                  &ti.stw_3_score);
    reader.AddVariable("stw_4_score",                  &ti.stw_4_score);
    reader.AddVariable("sig_1_score",                  &ti.sig_1_score);
    reader.AddVariable("sig_2_score",                  &ti.sig_2_score);
    reader.AddVariable("lol_1_score",                  &ti.lol_1_score);
    reader.AddVariable("lol_2_score",                  &ti.lol_2_score);
    reader.AddVariable("tro_1_score",                  &ti.tro_1_score);
    reader.AddVariable("tro_2_score",                  &ti.tro_2_score);
    reader.AddVariable("tro_4_score",                  &ti.tro_4_score);
    reader.AddVariable("tro_5_score",                  &ti.tro_5_score);
    reader.AddVariable("br4_1_shower_main_length",     &ti.br4_1_shower_main_length);
    reader.AddVariable("br4_1_shower_total_length",    &ti.br4_1_shower_total_length);
    reader.AddVariable("br4_1_min_dis",                &ti.br4_1_min_dis);
    reader.AddVariable("br4_1_energy",                 &ti.br4_1_energy);
    reader.AddVariable("br4_1_flag_avoid_muon_check",  &ti.br4_1_flag_avoid_muon_check);
    reader.AddVariable("br4_1_n_vtx_segs",             &ti.br4_1_n_vtx_segs);
    reader.AddVariable("br4_2_ratio_45",               &ti.br4_2_ratio_45);
    reader.AddVariable("br4_2_ratio_35",               &ti.br4_2_ratio_35);
    reader.AddVariable("br4_2_ratio_25",               &ti.br4_2_ratio_25);
    reader.AddVariable("br4_2_ratio_15",               &ti.br4_2_ratio_15);
    reader.AddVariable("br4_2_ratio1_45",              &ti.br4_2_ratio1_45);
    reader.AddVariable("br4_2_ratio1_35",              &ti.br4_2_ratio1_35);
    reader.AddVariable("br4_2_ratio1_25",              &ti.br4_2_ratio1_25);
    reader.AddVariable("br4_2_ratio1_15",              &ti.br4_2_ratio1_15);
    reader.AddVariable("br4_2_iso_angle",              &ti.br4_2_iso_angle);
    reader.AddVariable("br4_2_iso_angle1",             &ti.br4_2_iso_angle1);
    reader.AddVariable("br4_2_angle",                  &ti.br4_2_angle);
    reader.AddVariable("tro_3_stem_length",            &ti.tro_3_stem_length);
    reader.AddVariable("tro_3_n_muon_segs",            &ti.tro_3_n_muon_segs);
    reader.AddVariable("br4_1_n_main_segs",            &ti.br4_1_n_main_segs);

    reader.BookMVA("MyBDT", m_nue_xgboost_xml);

    // ---- Step 4: evaluate and transform -------------------------------------
    if (ti.br_filled == 1) {
        float val1 = reader.EvaluateMVA("MyBDT");
        ti.nue_score = static_cast<float>(std::log10((1.0 + val1) / (1.0 - val1)));
    } else {
        ti.nue_score = -15.f;  // background-like default; matches prototype default_val = -15
    }
}
