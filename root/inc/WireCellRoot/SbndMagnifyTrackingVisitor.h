/** Visitor that writes SBND STM-stage tracking data to a ROOT file.
 *
 * Fork-by-duplication of UbooneMagnifyTrackingVisitor (which reads the
 * neutrino-PR TrackFitting via grouping.get_track_fitting() and its PRGraph).
 * This variant serves the STM tagger validation chain instead:
 *
 *  - T_rec_charge is written from the per-cluster "stm_fit"/"stm_pass" point
 *    clouds persisted by TaggerCheckSTM's save_stm_fit knob (per-pass final
 *    fits, including rejected passes; there is no PRGraph at this stage).
 *    Each (cluster, pass) becomes its own track block with
 *    ndf = cluster_id*10 + pass so the Magnify-tracking GUI shows the
 *    forward and backward fits as separate "clusters".
 *  - T_proj_data is written from the TrackFitting stored in the grouping's
 *    named "stm" slot (pred/meas 2D charge accumulated over all passes).
 *  - Two-TPC channel convention: per-plane concatenation,
 *    global = plane_base[p] + apa * nch[p] + wire, with nch[p] taken from
 *    the configured anodes.  (uBooNE's fixed 0/2400/4800 offsets are the
 *    single-TPC special case of this.)
 *
 * Runs as an IEnsembleVisitor inside the MABC pipeline, after
 * tagger_check_stm (pipeline name "stm_magnify" in the SBND config).
 */

#ifndef WIRECELLROOT_SBNDMAGNIFYTRACKINGVISITOR
#define WIRECELLROOT_SBNDMAGNIFYTRACKINGVISITOR

#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellUtil/Logging.h"

#include <array>

class TFile;

namespace WireCell {
    namespace Root {

        class SbndMagnifyTrackingVisitor : public IConfigurable, public Clus::IEnsembleVisitor {
           public:
            SbndMagnifyTrackingVisitor();
            virtual ~SbndMagnifyTrackingVisitor();

            virtual void configure(const WireCell::Configuration& config);
            virtual Configuration default_configuration() const;
            virtual void visit(Clus::Facade::Ensemble& ensemble) const;

           private:
            Log::logptr_t log;
            std::string m_output_filename;
            std::string m_grouping_name{"live"};
            std::string m_track_fitting_name{"stm"};
            int m_runNo{0};
            int m_subRunNo{0};
            int m_eventNo{0};
            // This event's triplet; see SbndPrMagnifyTrackingVisitor.h.
            mutable int m_evt_runNo{0};
            mutable int m_evt_subRunNo{0};
            mutable int m_evt_eventNo{0};
            std::vector<IAnodePlane::pointer> m_anodes;
            IDetectorVolumes::pointer m_dv;
            double m_dQdx_scale{0.1};
            double m_dQdx_offset{-1000};
            // Readout length in ticks, used only to clamp dead-channel time
            // ranges (an unbounded dead x extent converts to +-1e9 ticks).
            int m_nticks{3427};

            // Concatenated-per-plane global channel scheme, derived from the
            // anodes at visit time: nch[p] = per-TPC channel count of plane p,
            // base[p] = sum over planes p' < p of n_apa * nch[p'].
            struct ChanScheme {
                std::array<int, 3> nch{0, 0, 0};
                std::array<int, 3> base{0, 0, 0};
                int global(int plane, int apa, int wire) const {
                    return base[plane] + apa * nch[plane] + wire;
                }
            };
            ChanScheme chan_scheme() const;

            void write_bad_channels(TFile* output_tf, Clus::Facade::Grouping& grouping, const ChanScheme& cs) const;
            void write_proj_data(TFile* output_tf, Clus::Facade::Grouping& grouping, const ChanScheme& cs) const;
            void write_t_rec_data(TFile* output_tf, Clus::Facade::Grouping& grouping, const ChanScheme& cs) const;
            void write_stm_trees(TFile* output_tf, Clus::Facade::Grouping& grouping) const;
            void write_trun(TFile* output_tf) const;
        };
    }  // namespace Root
}  // namespace WireCell

#endif
