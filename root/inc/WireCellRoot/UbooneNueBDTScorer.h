/** IEnsembleVisitor that runs TMVA BDT scoring for the nueCC tagger.
 *
 * Must run AFTER TaggerCheckNeutrino in the visitor pipeline.
 * Reads TaggerInfo + KineInfo from TrackFitting, evaluates all sub-BDTs,
 * and writes the resulting score into TaggerInfo::nue_score.
 *
 * Only cal_bdts_xgboost() and its 30 sub-BDTs are ported here.
 * (cal_bdts() — the older TMVA combination variant — is not ported.)
 *
 * XML weight files are expected in wire-cell-data under uboone/weights/.
 */

#ifndef WIRECELLROOT_UBOONENUEBDTSCORER_H
#define WIRECELLROOT_UBOONENUEBDTSCORER_H

#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellClus/NeutrinoTaggerInfo.h"
#include "WireCellUtil/Logging.h"

#include "TMVA/Reader.h"
#include <memory>

namespace WireCell {
    namespace Root {

        class UbooneNueBDTScorer : public IConfigurable,
                                   public Clus::IEnsembleVisitor
        {
        public:
            UbooneNueBDTScorer();
            virtual ~UbooneNueBDTScorer() = default;

            virtual void configure(const WireCell::Configuration& cfg);
            virtual Configuration default_configuration() const;

            /// Read TaggerInfo/KineInfo from grouping's TrackFitting,
            /// run xgboost BDT scoring, write nue_score back.
            virtual void visit(Clus::Facade::Ensemble& ensemble) const;

        private:
            Log::logptr_t log;
            std::string m_grouping_name{"live"};

            // Paths to TMVA XML weight files (resolved at configure time).
            // Scalar-feature sub-BDTs:
            std::string m_mipid_xml;        // mipid_BDT.weights.xml
            std::string m_gap_xml;          // gap_BDT.weights.xml
            std::string m_hol_lol_xml;      // hol_lol_BDT.weights.xml
            std::string m_cme_anc_xml;      // cme_anc_BDT.weights.xml
            std::string m_mgo_mgt_xml;      // mgo_mgt_BDT.weights.xml
            std::string m_br1_xml;          // br1_BDT.weights.xml
            std::string m_br3_xml;          // br3_BDT.weights.xml
            std::string m_stemdir_br2_xml;  // stem_dir_br2_BDT.weights.xml
            std::string m_trimuon_xml;      // stl_lem_brm_BDT.weights.xml
            std::string m_br4_tro_xml;      // br4_tro_BDT.weights.xml
            std::string m_mipquality_xml;   // mipquality_BDT.weights.xml
            std::string m_pio_1_xml;        // pio_1_BDT.weights.xml
            std::string m_stw_spt_xml;      // stw_spt_BDT.weights.xml
            std::string m_vis_1_xml;        // vis_1_BDT.weights.xml
            std::string m_vis_2_xml;        // vis_2_BDT.weights.xml
            // Vector-feature sub-BDTs:
            std::string m_br3_3_xml;        // br3_3_BDT.weights.xml
            std::string m_br3_5_xml;        // br3_5_BDT.weights.xml
            std::string m_br3_6_xml;        // br3_6_BDT.weights.xml
            std::string m_pio_2_xml;        // pio_2_BDT.weights.xml
            std::string m_stw_2_xml;        // stw_2_BDT.weights.xml
            std::string m_stw_3_xml;        // stw_3_BDT.weights.xml
            std::string m_stw_4_xml;        // stw_4_BDT.weights.xml
            std::string m_sig_1_xml;        // sig_1_BDT.weights.xml
            std::string m_sig_2_xml;        // sig_2_BDT.weights.xml
            std::string m_lol_1_xml;        // lol_1_BDT.weights.xml
            std::string m_lol_2_xml;        // lol_2_BDT.weights.xml
            std::string m_tro_1_xml;        // tro_1_BDT.weights.xml
            std::string m_tro_2_xml;        // tro_2_BDT.weights.xml
            std::string m_tro_4_xml;        // tro_4_BDT.weights.xml
            std::string m_tro_5_xml;        // tro_5_BDT.weights.xml
            // Final XGBoost combiner:
            std::string m_nue_xgboost_xml;  // XGB_nue_seed2_0923.xml

            // --------------- scalar-feature sub-BDTs -------------------
            // Gate on ti.*_filled==1 (or pio_filled + pio_flag_pio).
            // Return default_val when the fill condition is not met.
            float cal_mipid_bdt      (Clus::PR::TaggerInfo& ti, float default_val) const;
            float cal_gap_bdt        (Clus::PR::TaggerInfo& ti, float default_val) const;
            float cal_hol_lol_bdt    (Clus::PR::TaggerInfo& ti, float default_val) const;
            float cal_cme_anc_bdt    (Clus::PR::TaggerInfo& ti, float default_val) const;
            float cal_mgo_mgt_bdt    (Clus::PR::TaggerInfo& ti, float default_val) const;
            float cal_br1_bdt        (Clus::PR::TaggerInfo& ti, float default_val) const;
            float cal_br3_bdt        (Clus::PR::TaggerInfo& ti, float default_val) const;
            float cal_stemdir_br2_bdt(Clus::PR::TaggerInfo& ti, float default_val) const;
            float cal_trimuon_bdt    (Clus::PR::TaggerInfo& ti, float default_val) const;
            float cal_br4_tro_bdt    (Clus::PR::TaggerInfo& ti, float default_val) const;
            float cal_mipquality_bdt (Clus::PR::TaggerInfo& ti, float default_val) const;
            float cal_pio_1_bdt      (Clus::PR::TaggerInfo& ti, float default_val) const;
            float cal_stw_spt_bdt    (Clus::PR::TaggerInfo& ti, float default_val) const;
            float cal_vis_1_bdt      (Clus::PR::TaggerInfo& ti, float default_val) const;
            float cal_vis_2_bdt      (Clus::PR::TaggerInfo& ti, float default_val) const;

            // --------------- vector-feature sub-BDTs -------------------
            // Evaluate per-element, return minimum score; return default_val if vector empty.
            float cal_br3_3_bdt(Clus::PR::TaggerInfo& ti, float default_val) const;
            float cal_br3_5_bdt(Clus::PR::TaggerInfo& ti, float default_val) const;
            float cal_br3_6_bdt(Clus::PR::TaggerInfo& ti, float default_val) const;
            float cal_pio_2_bdt(Clus::PR::TaggerInfo& ti, float default_val) const;
            float cal_stw_2_bdt(Clus::PR::TaggerInfo& ti, float default_val) const;
            float cal_stw_3_bdt(Clus::PR::TaggerInfo& ti, float default_val) const;
            float cal_stw_4_bdt(Clus::PR::TaggerInfo& ti, float default_val) const;
            float cal_sig_1_bdt(Clus::PR::TaggerInfo& ti, float default_val) const;
            float cal_sig_2_bdt(Clus::PR::TaggerInfo& ti, float default_val) const;
            float cal_lol_1_bdt(Clus::PR::TaggerInfo& ti, float default_val) const;
            float cal_lol_2_bdt(Clus::PR::TaggerInfo& ti, float default_val) const;
            float cal_tro_1_bdt(Clus::PR::TaggerInfo& ti, float default_val) const;
            float cal_tro_2_bdt(Clus::PR::TaggerInfo& ti, float default_val) const;
            float cal_tro_4_bdt(Clus::PR::TaggerInfo& ti, float default_val) const;
            float cal_tro_5_bdt(Clus::PR::TaggerInfo& ti, float default_val) const;

            /// Top-level scorer: fills all *_score fields then runs the XGBoost
            /// TMVA model; writes ti.nue_score.
            void cal_bdts_xgboost(Clus::PR::TaggerInfo& ti,
                                  const Clus::PR::KineInfo& ki) const;

            /// Initialize all TMVA readers at configure time.
            void init_readers();

            // ====== Cached TMVA readers and their float variable buffers ======
            // All mutable because visit() and the cal_*_bdt methods are const.

            // --- br3_3 vector sub-BDT ---
            mutable std::unique_ptr<TMVA::Reader> m_rdr_br3_3;
            mutable float m_br3_3_v_energy{0};
            mutable float m_br3_3_v_angle{0};
            mutable float m_br3_3_v_dir_length{0};
            mutable float m_br3_3_v_length{0};

            // --- br3_5 vector sub-BDT ---
            mutable std::unique_ptr<TMVA::Reader> m_rdr_br3_5;
            mutable float m_br3_5_v_dir_length{0};
            mutable float m_br3_5_v_total_length{0};
            mutable float m_br3_5_v_flag_avoid_muon_check{0};
            mutable float m_br3_5_v_n_seg{0};
            mutable float m_br3_5_v_angle{0};
            mutable float m_br3_5_v_sg_length{0};
            mutable float m_br3_5_v_energy{0};
            mutable float m_br3_5_v_n_segs{0};
            mutable float m_br3_5_v_shower_main_length{0};
            mutable float m_br3_5_v_shower_total_length{0};

            // --- br3_6 vector sub-BDT ---
            mutable std::unique_ptr<TMVA::Reader> m_rdr_br3_6;
            mutable float m_br3_6_v_angle{0};
            mutable float m_br3_6_v_angle1{0};
            mutable float m_br3_6_v_flag_shower_trajectory{0};
            mutable float m_br3_6_v_direct_length{0};
            mutable float m_br3_6_v_length{0};
            mutable float m_br3_6_v_n_other_vtx_segs{0};
            mutable float m_br3_6_v_energy{0};

            // --- pio_2 vector sub-BDT ---
            mutable std::unique_ptr<TMVA::Reader> m_rdr_pio_2;
            mutable float m_pio_2_v_dis2{0};
            mutable float m_pio_2_v_angle2{0};
            mutable float m_pio_2_v_acc_length{0};
            mutable float m_pio_2_v_pio_mip_id{0};

            // --- stw_2 vector sub-BDT ---
            mutable std::unique_ptr<TMVA::Reader> m_rdr_stw_2;
            mutable float m_stw_2_v_medium_dQ_dx{0};
            mutable float m_stw_2_v_energy{0};
            mutable float m_stw_2_v_angle{0};
            mutable float m_stw_2_v_dir_length{0};
            mutable float m_stw_2_v_max_dQ_dx{0};

            // --- stw_3 vector sub-BDT ---
            mutable std::unique_ptr<TMVA::Reader> m_rdr_stw_3;
            mutable float m_stw_3_v_angle{0};
            mutable float m_stw_3_v_dir_length{0};
            mutable float m_stw_3_v_energy{0};
            mutable float m_stw_3_v_medium_dQ_dx{0};

            // --- stw_4 vector sub-BDT ---
            mutable std::unique_ptr<TMVA::Reader> m_rdr_stw_4;
            mutable float m_stw_4_v_angle{0};
            mutable float m_stw_4_v_dis{0};
            mutable float m_stw_4_v_energy{0};

            // --- sig_1 vector sub-BDT ---
            mutable std::unique_ptr<TMVA::Reader> m_rdr_sig_1;
            mutable float m_sig_1_v_angle{0};
            mutable float m_sig_1_v_flag_single_shower{0};
            mutable float m_sig_1_v_energy{0};
            mutable float m_sig_1_v_energy_1{0};

            // --- sig_2 vector sub-BDT ---
            mutable std::unique_ptr<TMVA::Reader> m_rdr_sig_2;
            mutable float m_sig_2_v_energy{0};
            mutable float m_sig_2_v_shower_angle{0};
            mutable float m_sig_2_v_flag_single_shower{0};
            mutable float m_sig_2_v_medium_dQ_dx{0};
            mutable float m_sig_2_v_start_dQ_dx{0};

            // --- lol_1 vector sub-BDT ---
            mutable std::unique_ptr<TMVA::Reader> m_rdr_lol_1;
            mutable float m_lol_1_v_energy{0};
            mutable float m_lol_1_v_vtx_n_segs{0};
            mutable float m_lol_1_v_nseg{0};
            mutable float m_lol_1_v_angle{0};

            // --- lol_2 vector sub-BDT ---
            mutable std::unique_ptr<TMVA::Reader> m_rdr_lol_2;
            mutable float m_lol_2_v_length{0};
            mutable float m_lol_2_v_angle{0};
            mutable float m_lol_2_v_type{0};
            mutable float m_lol_2_v_vtx_n_segs{0};
            mutable float m_lol_2_v_energy{0};
            mutable float m_lol_2_v_shower_main_length{0};
            mutable float m_lol_2_v_flag_dir_weak{0};

            // --- tro_1 vector sub-BDT ---
            mutable std::unique_ptr<TMVA::Reader> m_rdr_tro_1;
            mutable float m_tro_1_v_particle_type{0};
            mutable float m_tro_1_v_flag_dir_weak{0};
            mutable float m_tro_1_v_min_dis{0};
            mutable float m_tro_1_v_sg1_length{0};
            mutable float m_tro_1_v_shower_main_length{0};
            mutable float m_tro_1_v_max_n_vtx_segs{0};
            mutable float m_tro_1_v_tmp_length{0};
            mutable float m_tro_1_v_medium_dQ_dx{0};
            mutable float m_tro_1_v_dQ_dx_cut{0};
            mutable float m_tro_1_v_flag_shower_topology{0};

            // --- tro_2 vector sub-BDT ---
            mutable std::unique_ptr<TMVA::Reader> m_rdr_tro_2;
            mutable float m_tro_2_v_energy{0};
            mutable float m_tro_2_v_stem_length{0};
            mutable float m_tro_2_v_iso_angle{0};
            mutable float m_tro_2_v_max_length{0};
            mutable float m_tro_2_v_angle{0};

            // --- tro_4 vector sub-BDT ---
            mutable std::unique_ptr<TMVA::Reader> m_rdr_tro_4;
            mutable float m_tro_4_v_dir2_mag{0};
            mutable float m_tro_4_v_angle{0};
            mutable float m_tro_4_v_angle1{0};
            mutable float m_tro_4_v_angle2{0};
            mutable float m_tro_4_v_length{0};
            mutable float m_tro_4_v_length1{0};
            mutable float m_tro_4_v_medium_dQ_dx{0};
            mutable float m_tro_4_v_end_dQ_dx{0};
            mutable float m_tro_4_v_energy{0};
            mutable float m_tro_4_v_shower_main_length{0};
            mutable float m_tro_4_v_flag_shower_trajectory{0};

            // --- tro_5 vector sub-BDT ---
            mutable std::unique_ptr<TMVA::Reader> m_rdr_tro_5;
            mutable float m_tro_5_v_max_angle{0};
            mutable float m_tro_5_v_min_angle{0};
            mutable float m_tro_5_v_max_length{0};
            mutable float m_tro_5_v_iso_angle{0};
            mutable float m_tro_5_v_n_vtx_segs{0};
            mutable float m_tro_5_v_min_count{0};
            mutable float m_tro_5_v_max_count{0};
            mutable float m_tro_5_v_energy{0};

            // --- Final xgboost reader and its ~160 float variable buffers ---
            mutable std::unique_ptr<TMVA::Reader> m_rdr_xgboost;
            mutable float m_xgb_match_isFC{0};
            mutable float m_xgb_kine_reco_Enu{0};
            mutable float m_xgb_cme_mu_energy{0};
            mutable float m_xgb_cme_energy{0};
            mutable float m_xgb_cme_mu_length{0};
            mutable float m_xgb_cme_length{0};
            mutable float m_xgb_cme_angle_beam{0};
            mutable float m_xgb_anc_angle{0};
            mutable float m_xgb_anc_max_angle{0};
            mutable float m_xgb_anc_max_length{0};
            mutable float m_xgb_anc_acc_forward_length{0};
            mutable float m_xgb_anc_acc_backward_length{0};
            mutable float m_xgb_anc_acc_forward_length1{0};
            mutable float m_xgb_anc_shower_main_length{0};
            mutable float m_xgb_anc_shower_total_length{0};
            mutable float m_xgb_anc_flag_main_outside{0};
            mutable float m_xgb_gap_flag_prolong_u{0};
            mutable float m_xgb_gap_flag_prolong_v{0};
            mutable float m_xgb_gap_flag_prolong_w{0};
            mutable float m_xgb_gap_flag_parallel{0};
            mutable float m_xgb_gap_n_points{0};
            mutable float m_xgb_gap_n_bad{0};
            mutable float m_xgb_gap_energy{0};
            mutable float m_xgb_gap_num_valid_tracks{0};
            mutable float m_xgb_gap_flag_single_shower{0};
            mutable float m_xgb_hol_1_n_valid_tracks{0};
            mutable float m_xgb_hol_1_min_angle{0};
            mutable float m_xgb_hol_1_energy{0};
            mutable float m_xgb_hol_1_min_length{0};
            mutable float m_xgb_hol_2_min_angle{0};
            mutable float m_xgb_hol_2_medium_dQ_dx{0};
            mutable float m_xgb_hol_2_ncount{0};
            mutable float m_xgb_lol_3_angle_beam{0};
            mutable float m_xgb_lol_3_n_valid_tracks{0};
            mutable float m_xgb_lol_3_min_angle{0};
            mutable float m_xgb_lol_3_vtx_n_segs{0};
            mutable float m_xgb_lol_3_shower_main_length{0};
            mutable float m_xgb_lol_3_n_out{0};
            mutable float m_xgb_lol_3_n_sum{0};
            mutable float m_xgb_hol_1_flag_all_shower{0};
            mutable float m_xgb_mgo_energy{0};
            mutable float m_xgb_mgo_max_energy{0};
            mutable float m_xgb_mgo_total_energy{0};
            mutable float m_xgb_mgo_n_showers{0};
            mutable float m_xgb_mgo_max_energy_1{0};
            mutable float m_xgb_mgo_max_energy_2{0};
            mutable float m_xgb_mgo_total_other_energy{0};
            mutable float m_xgb_mgo_n_total_showers{0};
            mutable float m_xgb_mgo_total_other_energy_1{0};
            mutable float m_xgb_mgt_flag_single_shower{0};
            mutable float m_xgb_mgt_max_energy{0};
            mutable float m_xgb_mgt_total_other_energy{0};
            mutable float m_xgb_mgt_max_energy_1{0};
            mutable float m_xgb_mgt_e_indirect_max_energy{0};
            mutable float m_xgb_mgt_e_direct_max_energy{0};
            mutable float m_xgb_mgt_n_direct_showers{0};
            mutable float m_xgb_mgt_e_direct_total_energy{0};
            mutable float m_xgb_mgt_flag_indirect_max_pio{0};
            mutable float m_xgb_mgt_e_indirect_total_energy{0};
            mutable float m_xgb_mip_quality_energy{0};
            mutable float m_xgb_mip_quality_overlap{0};
            mutable float m_xgb_mip_quality_n_showers{0};
            mutable float m_xgb_mip_quality_n_tracks{0};
            mutable float m_xgb_mip_quality_flag_inside_pi0{0};
            mutable float m_xgb_mip_quality_n_pi0_showers{0};
            mutable float m_xgb_mip_quality_shortest_length{0};
            mutable float m_xgb_mip_quality_acc_length{0};
            mutable float m_xgb_mip_quality_shortest_angle{0};
            mutable float m_xgb_mip_quality_flag_proton{0};
            mutable float m_xgb_br1_1_shower_type{0};
            mutable float m_xgb_br1_1_vtx_n_segs{0};
            mutable float m_xgb_br1_1_energy{0};
            mutable float m_xgb_br1_1_n_segs{0};
            mutable float m_xgb_br1_1_flag_sg_topology{0};
            mutable float m_xgb_br1_1_flag_sg_trajectory{0};
            mutable float m_xgb_br1_1_sg_length{0};
            mutable float m_xgb_br1_2_n_connected{0};
            mutable float m_xgb_br1_2_max_length{0};
            mutable float m_xgb_br1_2_n_connected_1{0};
            mutable float m_xgb_br1_2_n_shower_segs{0};
            mutable float m_xgb_br1_2_max_length_ratio{0};
            mutable float m_xgb_br1_2_shower_length{0};
            mutable float m_xgb_br1_3_n_connected_p{0};
            mutable float m_xgb_br1_3_max_length_p{0};
            mutable float m_xgb_br1_3_n_shower_main_segs{0};
            mutable float m_xgb_br3_1_energy{0};
            mutable float m_xgb_br3_1_n_shower_segments{0};
            mutable float m_xgb_br3_1_sg_flag_trajectory{0};
            mutable float m_xgb_br3_1_sg_direct_length{0};
            mutable float m_xgb_br3_1_sg_length{0};
            mutable float m_xgb_br3_1_total_main_length{0};
            mutable float m_xgb_br3_1_total_length{0};
            mutable float m_xgb_br3_1_iso_angle{0};
            mutable float m_xgb_br3_1_sg_flag_topology{0};
            mutable float m_xgb_br3_2_n_ele{0};
            mutable float m_xgb_br3_2_n_other{0};
            mutable float m_xgb_br3_2_other_fid{0};
            mutable float m_xgb_br3_4_acc_length{0};
            mutable float m_xgb_br3_4_total_length{0};
            mutable float m_xgb_br3_7_min_angle{0};
            mutable float m_xgb_br3_8_max_dQ_dx{0};
            mutable float m_xgb_br3_8_n_main_segs{0};
            mutable float m_xgb_vis_1_n_vtx_segs{0};
            mutable float m_xgb_vis_1_energy{0};
            mutable float m_xgb_vis_1_num_good_tracks{0};
            mutable float m_xgb_vis_1_max_angle{0};
            mutable float m_xgb_vis_1_max_shower_angle{0};
            mutable float m_xgb_vis_1_tmp_length1{0};
            mutable float m_xgb_vis_1_tmp_length2{0};
            mutable float m_xgb_vis_2_n_vtx_segs{0};
            mutable float m_xgb_vis_2_min_angle{0};
            mutable float m_xgb_vis_2_min_weak_track{0};
            mutable float m_xgb_vis_2_angle_beam{0};
            mutable float m_xgb_vis_2_min_angle1{0};
            mutable float m_xgb_vis_2_iso_angle1{0};
            mutable float m_xgb_vis_2_min_medium_dQ_dx{0};
            mutable float m_xgb_vis_2_min_length{0};
            mutable float m_xgb_vis_2_sg_length{0};
            mutable float m_xgb_vis_2_max_angle{0};
            mutable float m_xgb_vis_2_max_weak_track{0};
            mutable float m_xgb_pio_1_mass{0};
            mutable float m_xgb_pio_1_pio_type{0};
            mutable float m_xgb_pio_1_energy_1{0};
            mutable float m_xgb_pio_1_energy_2{0};
            mutable float m_xgb_pio_1_dis_1{0};
            mutable float m_xgb_pio_1_dis_2{0};
            mutable float m_xgb_pio_mip_id{0};
            mutable float m_xgb_stem_dir_flag_single_shower{0};
            mutable float m_xgb_stem_dir_angle{0};
            mutable float m_xgb_stem_dir_energy{0};
            mutable float m_xgb_stem_dir_angle1{0};
            mutable float m_xgb_stem_dir_angle2{0};
            mutable float m_xgb_stem_dir_angle3{0};
            mutable float m_xgb_stem_dir_ratio{0};
            mutable float m_xgb_br2_num_valid_tracks{0};
            mutable float m_xgb_br2_n_shower_main_segs{0};
            mutable float m_xgb_br2_max_angle{0};
            mutable float m_xgb_br2_sg_length{0};
            mutable float m_xgb_br2_flag_sg_trajectory{0};
            mutable float m_xgb_stem_len_energy{0};
            mutable float m_xgb_stem_len_length{0};
            mutable float m_xgb_stem_len_flag_avoid_muon_check{0};
            mutable float m_xgb_stem_len_num_daughters{0};
            mutable float m_xgb_stem_len_daughter_length{0};
            mutable float m_xgb_brm_n_mu_segs{0};
            mutable float m_xgb_brm_Ep{0};
            mutable float m_xgb_brm_acc_length{0};
            mutable float m_xgb_brm_shower_total_length{0};
            mutable float m_xgb_brm_connected_length{0};
            mutable float m_xgb_brm_n_size{0};
            mutable float m_xgb_brm_n_shower_main_segs{0};
            mutable float m_xgb_brm_n_mu_main{0};
            mutable float m_xgb_lem_shower_main_length{0};
            mutable float m_xgb_lem_n_3seg{0};
            mutable float m_xgb_lem_e_charge{0};
            mutable float m_xgb_lem_e_dQdx{0};
            mutable float m_xgb_lem_shower_num_main_segs{0};
            mutable float m_xgb_brm_acc_direct_length{0};
            mutable float m_xgb_stw_1_energy{0};
            mutable float m_xgb_stw_1_dis{0};
            mutable float m_xgb_stw_1_dQ_dx{0};
            mutable float m_xgb_stw_1_flag_single_shower{0};
            mutable float m_xgb_stw_1_n_pi0{0};
            mutable float m_xgb_stw_1_num_valid_tracks{0};
            mutable float m_xgb_spt_shower_main_length{0};
            mutable float m_xgb_spt_shower_total_length{0};
            mutable float m_xgb_spt_angle_beam{0};
            mutable float m_xgb_spt_angle_vertical{0};
            mutable float m_xgb_spt_max_dQ_dx{0};
            mutable float m_xgb_spt_angle_beam_1{0};
            mutable float m_xgb_spt_angle_drift{0};
            mutable float m_xgb_spt_angle_drift_1{0};
            mutable float m_xgb_spt_num_valid_tracks{0};
            mutable float m_xgb_spt_n_vtx_segs{0};
            mutable float m_xgb_spt_max_length{0};
            mutable float m_xgb_mip_energy{0};
            mutable float m_xgb_mip_n_end_reduction{0};
            mutable float m_xgb_mip_n_first_mip{0};
            mutable float m_xgb_mip_n_first_non_mip{0};
            mutable float m_xgb_mip_n_first_non_mip_1{0};
            mutable float m_xgb_mip_n_first_non_mip_2{0};
            mutable float m_xgb_mip_vec_dQ_dx_0{0};
            mutable float m_xgb_mip_vec_dQ_dx_1{0};
            mutable float m_xgb_mip_max_dQ_dx_sample{0};
            mutable float m_xgb_mip_n_below_threshold{0};
            mutable float m_xgb_mip_n_below_zero{0};
            mutable float m_xgb_mip_n_lowest{0};
            mutable float m_xgb_mip_n_highest{0};
            mutable float m_xgb_mip_lowest_dQ_dx{0};
            mutable float m_xgb_mip_highest_dQ_dx{0};
            mutable float m_xgb_mip_medium_dQ_dx{0};
            mutable float m_xgb_mip_stem_length{0};
            mutable float m_xgb_mip_length_main{0};
            mutable float m_xgb_mip_length_total{0};
            mutable float m_xgb_mip_angle_beam{0};
            mutable float m_xgb_mip_iso_angle{0};
            mutable float m_xgb_mip_n_vertex{0};
            mutable float m_xgb_mip_n_good_tracks{0};
            mutable float m_xgb_mip_E_indirect_max_energy{0};
            mutable float m_xgb_mip_flag_all_above{0};
            mutable float m_xgb_mip_min_dQ_dx_5{0};
            mutable float m_xgb_mip_n_other_vertex{0};
            mutable float m_xgb_mip_n_stem_size{0};
            mutable float m_xgb_mip_flag_stem_trajectory{0};
            mutable float m_xgb_mip_min_dis{0};
            mutable float m_xgb_mip_vec_dQ_dx_2{0};
            mutable float m_xgb_mip_vec_dQ_dx_3{0};
            mutable float m_xgb_mip_vec_dQ_dx_4{0};
            mutable float m_xgb_mip_vec_dQ_dx_5{0};
            mutable float m_xgb_mip_vec_dQ_dx_6{0};
            mutable float m_xgb_mip_vec_dQ_dx_7{0};
            mutable float m_xgb_mip_vec_dQ_dx_8{0};
            mutable float m_xgb_mip_vec_dQ_dx_9{0};
            mutable float m_xgb_mip_vec_dQ_dx_10{0};
            mutable float m_xgb_mip_vec_dQ_dx_11{0};
            mutable float m_xgb_mip_vec_dQ_dx_12{0};
            mutable float m_xgb_mip_vec_dQ_dx_13{0};
            mutable float m_xgb_mip_vec_dQ_dx_14{0};
            mutable float m_xgb_mip_vec_dQ_dx_15{0};
            mutable float m_xgb_mip_vec_dQ_dx_16{0};
            mutable float m_xgb_mip_vec_dQ_dx_17{0};
            mutable float m_xgb_mip_vec_dQ_dx_18{0};
            mutable float m_xgb_mip_vec_dQ_dx_19{0};
            mutable float m_xgb_br3_3_score{0};
            mutable float m_xgb_br3_5_score{0};
            mutable float m_xgb_br3_6_score{0};
            mutable float m_xgb_pio_2_score{0};
            mutable float m_xgb_stw_2_score{0};
            mutable float m_xgb_stw_3_score{0};
            mutable float m_xgb_stw_4_score{0};
            mutable float m_xgb_sig_1_score{0};
            mutable float m_xgb_sig_2_score{0};
            mutable float m_xgb_lol_1_score{0};
            mutable float m_xgb_lol_2_score{0};
            mutable float m_xgb_tro_1_score{0};
            mutable float m_xgb_tro_2_score{0};
            mutable float m_xgb_tro_4_score{0};
            mutable float m_xgb_tro_5_score{0};
            mutable float m_xgb_br4_1_shower_main_length{0};
            mutable float m_xgb_br4_1_shower_total_length{0};
            mutable float m_xgb_br4_1_min_dis{0};
            mutable float m_xgb_br4_1_energy{0};
            mutable float m_xgb_br4_1_flag_avoid_muon_check{0};
            mutable float m_xgb_br4_1_n_vtx_segs{0};
            mutable float m_xgb_br4_2_ratio_45{0};
            mutable float m_xgb_br4_2_ratio_35{0};
            mutable float m_xgb_br4_2_ratio_25{0};
            mutable float m_xgb_br4_2_ratio_15{0};
            mutable float m_xgb_br4_2_ratio1_45{0};
            mutable float m_xgb_br4_2_ratio1_35{0};
            mutable float m_xgb_br4_2_ratio1_25{0};
            mutable float m_xgb_br4_2_ratio1_15{0};
            mutable float m_xgb_br4_2_iso_angle{0};
            mutable float m_xgb_br4_2_iso_angle1{0};
            mutable float m_xgb_br4_2_angle{0};
            mutable float m_xgb_tro_3_stem_length{0};
            mutable float m_xgb_tro_3_n_muon_segs{0};
            mutable float m_xgb_br4_1_n_main_segs{0};
        };

    }  // namespace Root
}  // namespace WireCell

#endif  // WIRECELLROOT_UBOONENUEBDTSCORER_H
