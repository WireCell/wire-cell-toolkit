/** Visitor that writes SBND neutrino-PR tracking data to a ROOT file.
 *
 * Fork-by-duplication of UbooneMagnifyTrackingVisitor (CLAUDE.md M10: the
 * uBooNE original serves the qlport gate chain and stays byte-for-byte
 * untouched; SbndMagnifyTrackingVisitor serves the STM chain and likewise
 * stays untouched).  Same data source as the uBooNE writer -- the unnamed
 * TrackFitting slot filled by TaggerCheckNeutrino and its PRGraph -- with the
 * SBND two-TPC deltas:
 *
 *  - Two-TPC channel convention (copied from SbndMagnifyTrackingVisitor):
 *    per-plane concatenation, global = plane_base[p] + apa * nch[p] + wire,
 *    with nch[p] taken from the configured anodes at visit time.  uBooNE's
 *    fixed 0/2400/4800 offsets are the single-TPC special case of this.
 *  - T_rec_charge pu/pv/pw/pt use the per-point (apa,face) recorded in
 *    PR::Fit::paf, so APA-1 tracks land on APA-1 channels instead of being
 *    overlaid on APA-0 (the uBooNE writer has no per-point APA concept).
 *  - T_bad_ch time ranges are clamped to [0, nticks] (dead regions with
 *    unbounded x extent convert to +-1e9 ticks; same fix as the STM fork).
 *
 * Tree schema (T_bad_ch, Trun, T_proj_data, T_rec_charge, empty T_proj) is
 * unchanged from the uBooNE writer so the existing Magnify-tracking convert
 * apps and GUI layout keep working.
 *
 * Runs as an IEnsembleVisitor inside the MABC pipeline, after
 * tagger_check_neutrino (pipeline name "tracking_visitor" in the SBND config).
 */

#ifndef WIRECELLROOT_SBNDPRMAGNIFYTRACKINGVISITOR
#define WIRECELLROOT_SBNDPRMAGNIFYTRACKINGVISITOR

#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellUtil/Logging.h"

#include <array>
#include <map>
#include <string>

class TFile;

namespace WireCell {
    namespace Root {

        class SbndPrMagnifyTrackingVisitor : public IConfigurable, public Clus::IEnsembleVisitor {
           public:
            SbndPrMagnifyTrackingVisitor();
            virtual ~SbndPrMagnifyTrackingVisitor();

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
            // The triplet actually stamped into THIS event's trees.  Equal to
            // the three above for a one-event-per-process job; taken from the
            // ensemble (which MultiAlgBlobClustering fills per event) when a
            // group runs in one process, where the configured numbers are one
            // constant for every event in the group.  Set at the top of
            // visit(), which is const.
            mutable int m_evt_runNo{0};
            mutable int m_evt_subRunNo{0};
            mutable int m_evt_eventNo{0};
            std::vector<IAnodePlane::pointer> m_anodes;
            IDetectorVolumes::pointer m_dv;
            // doc 87: when true, add a per-cluster T_cluster tree recording the
            // in-scope set (switch_scope's scope_filter) and the per-bundle
            // summary nusel_extract.py otherwise has to read out of the Bee zip
            // and the pctree.  DEFAULT FALSE => tracking-pr.root byte-identical.
            bool m_save_in_scope{false};
            // doc 99: when true, T_cluster's flash_id/flash_time_us/flash_pe
            // come from Cluster::get_matched_flash() (matched_flash_gid against
            // the merge-safe "opflash" PC) instead of Cluster::get_flash() (the
            // per-input row index the multi-APA merge invalidates).  DEFAULT
            // FALSE => tracking-pr.root byte-identical.
            bool m_flash_by_gid{false};
            // sbnd_xin/docs/109 group 1: write TaggerCheckNeutrino's selection
            // census (Grouping::get_nu_census) as T_bundle (one row per
            // in-beam-window flash bundle, with the reason it did or did not
            // yield a candidate) and T_flash (one row per optical flash, with
            // its TPC and flash group), plus the event counters and the
            // selection's window in Trun.  Written on every event, including
            // those with no T_tagger row.  DEFAULT FALSE => byte-identical.
            bool m_nu_provenance{false};
            // sbnd_xin/docs/109 group 2: T_cluster tgm/stm/fc read the flags
            // SBND's taggers set (Flags::TGM/STM/FC; the legacy reads are the
            // lowercase flags only the unconfigured ClusteringTaggerFlagTransfer
            // sets, so they are always 0 here), lm reads the "lm_flag" scalar
            // Q/L matching sets, beam_flash is DERIVED with the selection's own
            // test (matched_flash_gid >= 0 and cluster_t0 in the beam window;
            // -1 without a census), and matched_flash_gid / flash_tpc are
            // added.  DEFAULT FALSE => byte-identical.
            bool m_fix_cluster_flags{false};
            // sbnd_xin/docs/109 rev 3: make T_rec_charge joinable to the
            // candidate that owns its points.
            //   - cluster_id is taken from the candidate's own TrackFitting
            //     (TaggerInfo::cluster_id, the same POINTER-derived value
            //     T_tagger/T_kine carry) instead of scanning this candidate's
            //     graph for Flags::main_cluster.  That scan returns -1 when the
            //     selected activity is a demoted main -- ClusteringUnmergeBundle
            //     clears the flag and never restores it -- and returns the
            //     PRE-swap cluster when the overall-vertex search moved the main
            //     onto a companion, so neither case joins back to T_tagger.
            //   - nu_index and point_cluster_id are added: which candidate a row
            //     belongs to, and the row's own cluster (previously readable
            //     only out of the overloaded `ndf` branch).
            //   - T_proj_data is booked even when there is nothing to write, so
            //     the tree set of the file does not vary with the event.
            // DEFAULT FALSE => tracking-pr.root byte-identical.
            bool m_rec_charge_provenance{false};
            // sbnd_xin/docs/109 group 4: provenance strings written into Trun
            // as std::string branches, one per key (sorted), plus
            // wct_version.  Empty (the default) => nothing written.
            std::map<std::string, std::string> m_provenance;
            double m_dQdx_scale{0.1};
            double m_dQdx_offset{-1000};
            bool m_flag_skip_vertex{false};
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
            void write_trun(TFile* output_tf, Clus::Facade::Grouping& grouping) const;
            void write_cluster_summary(TFile* output_tf, Clus::Facade::Grouping& grouping) const;
            void write_nu_census(TFile* output_tf, Clus::Facade::Grouping& grouping) const;
        };
    }  // namespace Root
}  // namespace WireCell

#endif
