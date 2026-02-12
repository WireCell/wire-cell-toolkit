/** Visitor that writes tracking data to a ROOT file.
 *
 * This runs as an IEnsembleVisitor inside the MABC pipeline,
 * before tensor serialization, so it can access TrackFitting directly.
 * Writes T_bad_ch, T_proj_data, and T_proj trees.
 */

#ifndef WIRECELLROOT_UBOONEMAGNIFYTRACKINGVISITOR
#define WIRECELLROOT_UBOONEMAGNIFYTRACKINGVISITOR

#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellUtil/Logging.h"

class TFile;

namespace WireCell {
    namespace Root {

        class UbooneMagnifyTrackingVisitor : public IConfigurable, public Clus::IEnsembleVisitor {
           public:
            UbooneMagnifyTrackingVisitor();
            virtual ~UbooneMagnifyTrackingVisitor();

            virtual void configure(const WireCell::Configuration& config);
            virtual Configuration default_configuration() const;
            virtual void visit(Clus::Facade::Ensemble& ensemble) const;

           private:
            Log::logptr_t log;
            std::string m_output_filename;
            std::string m_grouping_name{"live"};
            int m_runNo{0};
            int m_subRunNo{0};
            int m_eventNo{0};
            std::vector<IAnodePlane::pointer> m_anodes;
            IDetectorVolumes::pointer m_dv;

            void write_bad_channels(TFile* output_tf, Clus::Facade::Grouping& grouping) const;
            void write_proj_data(TFile* output_tf, Clus::Facade::Grouping& grouping) const;
        };
    }  // namespace Root
}  // namespace WireCell

#endif
