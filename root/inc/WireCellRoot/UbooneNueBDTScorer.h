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

#ifndef WIRECELLROOT_UBOONENUEBTDSCORER_H
#define WIRECELLROOT_UBOONENUEBTDSCORER_H

#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellClus/NeutrinoTaggerInfo.h"
#include "WireCellUtil/Logging.h"

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
        };

    }  // namespace Root
}  // namespace WireCell

#endif  // WIRECELLROOT_UBOONENUEBTDSCORER_H
