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

        // Structure to hold reconstruction point data
        struct WCPointTree {
            double reco_x{0};
            double reco_y{0};
            double reco_z{0};
            double reco_dQ{0};
            double reco_dx{0};
            double reco_chi2{0};
            double reco_ndf{0};
            double reco_pu{0};
            double reco_pv{0};
            double reco_pw{0};
            double reco_pt{0};
            double reco_reduced_chi2{0};
            int reco_flag_vertex{0};
            int reco_flag_track_shower{0};
            double reco_rr{0};
            int reco_mother_cluster_id{0};
            int reco_cluster_id{0};
            int reco_proto_cluster_id{0};
            int reco_particle_id{0};
        };

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
            double m_dQdx_scale{0.1};
            double m_dQdx_offset{-1000};
            bool m_flag_skip_vertex{false};

            void write_bad_channels(TFile* output_tf, Clus::Facade::Grouping& grouping) const;
            void write_proj_data(TFile* output_tf, Clus::Facade::Grouping& grouping) const;
            void write_t_rec_data(TFile* output_tf, Clus::Facade::Grouping& grouping) const;
            void write_trun(TFile* output_tf) const;
        };
    }  // namespace Root
}  // namespace WireCell

#endif
