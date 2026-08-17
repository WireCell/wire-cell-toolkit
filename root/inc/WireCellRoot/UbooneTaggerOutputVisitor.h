/** Visitor that writes tagger variables (T_tagger, T_kine) to a ROOT file.
 *
 * This runs as an IEnsembleVisitor inside the MABC pipeline,
 * AFTER UbooneNumuBDTScorer and UbooneNueBDTScorer (which fill BDT scores)
 * and AFTER UbooneMagnifyTrackingVisitor (which creates the ROOT file).
 *
 * Opens the existing ROOT output file in UPDATE mode and adds two TTrees:
 *   - T_tagger: all TaggerInfo fields (BDT input features + scores)
 *   - T_kine:   all KineInfo fields (kinematic reconstruction)
 *
 * Tree names and branch names match the prototype (WCPPID) output format
 * for direct comparison.
 */

#ifndef WIRECELLROOT_UBOONETAGGEROUTPUTVISITOR
#define WIRECELLROOT_UBOONETAGGEROUTPUTVISITOR

#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Logging.h"

class TFile;

namespace WireCell {
    namespace Root {

        class UbooneTaggerOutputVisitor : public IConfigurable, public Clus::IEnsembleVisitor {
           public:
            UbooneTaggerOutputVisitor();
            virtual ~UbooneTaggerOutputVisitor();

            virtual void configure(const WireCell::Configuration& config);
            virtual Configuration default_configuration() const;
            virtual void visit(Clus::Facade::Ensemble& ensemble) const;

           private:
            Log::logptr_t log;
            std::string m_output_filename{"tracking_proj.root"};
            std::string m_grouping_name{"live"};
            // doc sbnd_xin/docs/pr/36 §10.8 (F7 = P4): book the neutrino_type
            // verdict-bitmask branch (prototype
            // wire-cell-prod-nue-port.cxx:1486, neutrino_type/I).  Same
            // config key as TaggerCheckNeutrino's computing knob.  C++
            // default false => branch not booked => the T_tagger schema is
            // byte-identical to the pre-knob tree.
            bool m_neutrino_type_bitmask{false};
        };

    }  // namespace Root
}  // namespace WireCell

#endif
