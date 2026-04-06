/** IEnsembleVisitor that runs TMVA BDT scoring for the numu CC tagger.
 *
 * Must run AFTER TaggerCheckNeutrino in the visitor pipeline.
 * Reads TaggerInfo + KineInfo from TrackFitting, evaluates BDTs,
 * and writes the resulting scores back into TaggerInfo.
 *
 * Only cal_numu_bdts_xgboost() and its sub-BDTs are ported here.
 * (cal_numu_bdts() — the older TMVA variant — is not ported.)
 *
 * XML weight files are expected in wire-cell-data under uboone/weights/.
 */

#ifndef WIRECELLROOT_UBOONENUMUBDTSCORER_H
#define WIRECELLROOT_UBOONENUMUBDTSCORER_H

#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellClus/NeutrinoTaggerInfo.h"
#include "WireCellUtil/Logging.h"

namespace WireCell {
    namespace Root {

        class UbooneNumuBDTScorer : public IConfigurable,
                                    public Clus::IEnsembleVisitor
        {
        public:
            UbooneNumuBDTScorer();
            virtual ~UbooneNumuBDTScorer() = default;

            virtual void configure(const WireCell::Configuration& cfg);
            virtual Configuration default_configuration() const;

            /// Read TaggerInfo/KineInfo from grouping's TrackFitting,
            /// run xgboost BDT scoring, write scores back.
            virtual void visit(Clus::Facade::Ensemble& ensemble) const;

        private:
            Log::logptr_t log;

            // Paths to TMVA XML weight files (resolved via wc.resolve at configure time).
            std::string m_numu1_xml;        // numu_tagger1.weights.xml
            std::string m_numu2_xml;        // numu_tagger2.weights.xml
            std::string m_numu3_xml;        // numu_tagger3.weights.xml
            std::string m_cosmict10_xml;    // cos_tagger_10.weights.xml
            std::string m_numu_xgboost_xml; // numu_scalars_scores_0923.xml

            std::string m_grouping_name{"live"};

            // Sub-BDT scorers — each returns the BDT score for one sub-classifier.
            // TaggerInfo is passed by mutable reference so the caller can write
            // the intermediate scores (numu_1_score etc.) before the final BDT.
            float cal_cosmict_10_bdt(Clus::PR::TaggerInfo& ti, float default_val) const;
            float cal_numu_1_bdt   (Clus::PR::TaggerInfo& ti, float default_val) const;
            float cal_numu_2_bdt   (Clus::PR::TaggerInfo& ti, float default_val) const;
            float cal_numu_3_bdt   (Clus::PR::TaggerInfo& ti, float default_val) const;

            /// Top-level scorer: calls the four sub-BDTs, then runs the xgboost
            /// TMVA model. Writes all *_score fields and numu_score into ti.
            void cal_numu_bdts_xgboost(Clus::PR::TaggerInfo& ti,
                                       const Clus::PR::KineInfo& ki) const;
        };

    }  // namespace Root
}  // namespace WireCell

#endif  // WIRECELLROOT_UBOONENUMUBDTSCORER_H
