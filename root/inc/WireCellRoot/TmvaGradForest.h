/** TmvaGradForest -- a compact, exact re-evaluation of a TMVA "BDT::BDT"
 *  gradient-boosted forest weight file, for the xgboost2TMVA-style models
 *  the uBooNE/SBND neutrino BDT scorers book (sbnd_xin/docs/76 round 2).
 *
 *  WHY.  UbooneNueBDTScorer's top-level combiner is a 199 MB TMVA XML
 *  (600 trees, 1.49 M nodes).  TMVA::Reader::BookMVA parses it into a full
 *  TXMLEngine DOM -- twice (GetMethodTypeFromFile, then ReadStateFromFile) --
 *  then walks every node through stringstream attribute reads: ~4.4 s of
 *  CPU and ~0.9 GB of heap per wire-cell process, i.e. per SBND production
 *  event, before a single event is processed.  This class scans the same
 *  file once (~0.3 s), keeps 24 bytes per node, and evaluates the forest
 *  with the SAME arithmetic TMVA uses, so the score is bit-identical.
 *
 *  WHAT IS MIRRORED (ROOT 6.32.02, tmva/tmva/src; every step cited):
 *   - Tools::ReadAttr(float&)  : value = atof(str)  (double, then narrowed
 *                                to Float_t)                   Tools.cxx:1779
 *   - Tools::ReadAttr(int/short&): atoi(str)                  Tools.cxx:1791
 *   - Tools::ReadAttr<Bool_t>  : std::stringstream >> value   Tools.h:338
 *   - Node::ReadXML            : children attach by pos 'l'/'r' Node.cxx
 *   - DecisionTreeNode::GoesRight (NCoef==0):
 *         result = (value[IVar] >= Cut); return cType ? result : !result;
 *                                                    DecisionTreeNode.cxx
 *   - DecisionTree::CheckEvent : descend while nType==0; regression tree
 *                                (Weights AnalysisType="1") returns the
 *                                leaf's Float_t response  DecisionTree.cxx
 *   - MethodBDT::GetGradBoostMVA: sum (double) of every tree's response;
 *                                return 2.0/(1.0+exp(-2.0*sum))-1
 *                                                       MethodBDT.cxx
 *   - Reader::EvaluateMVA      : any NaN input variable => return -999
 *                                                          Reader.cxx:514
 *
 *  WHAT IS REFUSED (load() throws WireCell::ValueError, so the caller can
 *  keep TMVA for it): any Transformations, BoostType other than "Grad",
 *  DoPreselection, Fisher-cut nodes (NCoef != 0), Weights AnalysisType
 *  other than 1 (regression leaves), a PreselectionLowBkgVar0 attribute,
 *  a variable list that does not match the caller's names in order
 *  (TMVA::MethodBase::ReadVariablesFromXML's own check).
 */
#ifndef WIRECELLROOT_TMVAGRADFOREST
#define WIRECELLROOT_TMVAGRADFOREST

#include "WireCellUtil/Logging.h"

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace WireCell::Root {

    class TmvaGradForest {
      public:
        /// Parse the weight file.  `names` are the input variables in the
        /// order the caller will pass them to evaluate() -- the same order
        /// it would AddVariable() them to a TMVA::Reader.  Throws
        /// WireCell::ValueError on any unsupported feature or mismatch.
        void load(const std::string& path, const std::vector<std::string>& names);

        /// Same, from an in-memory XML document (used by the unit test).
        void load_string(const std::string& xml, const std::vector<std::string>& names);

        /// TMVA::Reader::EvaluateMVA("...") on the forest.  `values` has
        /// nvars() entries.  Returns -999 if any of them is NaN.
        double evaluate(const float* values) const;

        size_t nvars() const { return m_nvars; }
        size_t ntrees() const { return m_roots.size(); }
        size_t nnodes() const { return m_nodes.size(); }

      private:
        void parse(std::string_view doc, const std::vector<std::string>& names);

        struct Node {
            float cut{0};       // Cut   (Float_t)
            float res{0};       // res   (Float_t leaf response)
            int32_t ivar{-1};   // IVar  (Short_t in TMVA; -1 on leaves)
            int32_t left{-1};   // index into m_nodes, -1 = none
            int32_t right{-1};
            int32_t ntype{0};   // nType: 0 = intermediate, else leaf
            bool ctype{true};   // cType
        };
        std::vector<Node> m_nodes;
        std::vector<int32_t> m_roots;  // one per tree, in file order
        size_t m_nvars{0};
    };

    /// The prototype's log-odds transform of a forest output: log10((1+v)/(1-v)).
    ///
    /// doc sbnd_xin/docs/85 sec 7.2 (owner decision 2026-08-30).  Both uBooNE
    /// BDT scorers used to clamp v to +-0.9999 first, capping |score| at
    /// 4.30103; the prototype (NeutrinoID_nue_bdts.h:300-301,
    /// NeutrinoID_numu_bdts.h:94-95) carries v in a double and does not clamp,
    /// so |score| reaches 16.25562 and the -15 "not filled" sentinel sits below
    /// the real floor.  The clamp made MicroBooNE's nue_score > 7.0 working
    /// point select zero events by construction.  It is gone; this is the
    /// transform that replaced it.
    ///
    /// The one guard kept is the degenerate one the transform needs.
    /// evaluate() above returns tanh(sum), which rounds to exactly +-1 in
    /// double once |sum| > ~18.4 (log10(2/0) = inf), and returns the -999
    /// sentinel when any input variable is NaN (log10 of a negative ratio =
    /// NaN).  Either would leave a non-finite value in a ROOT float branch and
    /// in every downstream cut, so |v| >= 1 and non-finite v are pulled to the
    /// neighbouring double -- a 1.1e-16 move to the finite extremes +-16.25562,
    /// not the +-4.30103 truncation the clamp imposed.  `log` and `who` are
    /// used only to report a hit: a silent inf is the failure mode a
    /// 3000-event batch cannot show you.
    double bdt_log_odds_score(double v, const Log::logptr_t& log, const char* who);

}  // namespace WireCell::Root

#endif
