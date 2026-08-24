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

#include "WireCellRoot/TmvaGradForest.h"

#include "TMVA/Reader.h"
#include <memory>
#include <mutex>
#include <string>
#include <vector>

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

            /// Create all TMVA readers and BookMVA their XML weights.  Called
            /// LAZILY, once, from the first visit() that actually has a
            /// TrackFitting to score (sbnd_xin/docs/76 round 1) -- see
            /// UbooneNueBDTScorer.h for the rationale; same contract here
            /// (~30 MB of XML, ~0.4 s, unchanged scores).
            void ensure_readers() const;
            void init_readers_impl() const;
            mutable std::once_flag m_readers_once;

            /// fast_xgb_forest (C++ default false = TMVA::Reader, legacy).
            /// true: numu_xgboost_xml is booked by TmvaGradForest -- same
            /// score, far cheaper to load.  See UbooneNueBDTScorer.h and
            /// sbnd_xin/docs/76 round 2.  The four sub-BDTs stay on TMVA.
            bool m_fast_xgb_forest{false};
            mutable std::unique_ptr<TmvaGradForest> m_xgb_forest;
            mutable std::vector<std::string> m_xgb_names;
            mutable std::vector<const float*> m_xgb_ptrs;
            bool have_xgb() const { return m_reader_xgboost || m_xgb_forest; }
            double eval_xgb() const;

            // --- Persistent TMVA readers (created once, on first scoring visit) ---
            // Mutable because visit() and sub-BDT methods are const.
            mutable std::unique_ptr<TMVA::Reader> m_reader_cosmict10;
            mutable std::unique_ptr<TMVA::Reader> m_reader_numu1;
            mutable std::unique_ptr<TMVA::Reader> m_reader_numu2;
            mutable std::unique_ptr<TMVA::Reader> m_reader_numu3;
            mutable std::unique_ptr<TMVA::Reader> m_reader_xgboost;

            // --- Float variable buffers bound to each reader via AddVariable ---
            // cosmict_10: 5 vars
            mutable float m_cosmict10_vtx_z{0};
            mutable float m_cosmict10_flag_shower{0};
            mutable float m_cosmict10_flag_dir_weak{0};
            mutable float m_cosmict10_angle_beam{0};
            mutable float m_cosmict10_length{0};

            // numu_1: 7 vars
            mutable float m_numu1_particle_type{0};
            mutable float m_numu1_length{0};
            mutable float m_numu1_medium_dQ_dx{0};
            mutable float m_numu1_dQ_dx_cut{0};
            mutable float m_numu1_direct_length{0};
            mutable float m_numu1_n_daughter_tracks{0};
            mutable float m_numu1_n_daughter_all{0};

            // numu_2: 4 vars
            mutable float m_numu2_length{0};
            mutable float m_numu2_total_length{0};
            mutable float m_numu2_n_daughter_tracks{0};
            mutable float m_numu2_n_daughter_all{0};

            // numu_3: 7 vars
            mutable float m_numu3_particle_type{0};
            mutable float m_numu3_max_length{0};
            mutable float m_numu3_track_length{0};
            mutable float m_numu3_max_length_all{0};
            mutable float m_numu3_max_muon_length{0};
            mutable float m_numu3_n_daughter_tracks{0};
            mutable float m_numu3_n_daughter_all{0};

            // xgboost final: ~72 vars (all scalar inputs)
            // Indexed by the order they are added in init_readers_impl().
            mutable std::vector<float> m_xgb_vars;
        };

    }  // namespace Root
}  // namespace WireCell

#endif  // WIRECELLROOT_UBOONENUMUBDTSCORER_H
