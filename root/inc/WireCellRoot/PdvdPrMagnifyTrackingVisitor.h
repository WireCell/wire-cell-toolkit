/** Visitor that writes PDVD neutrino-PR tracking data to a ROOT file.
 *
 * Fork-by-duplication of SbndPrMagnifyTrackingVisitor (doc pdvd/25 M6; the
 * SBND file is untouched), itself a fork of UbooneMagnifyTrackingVisitor.
 * Reads the unnamed TrackFitting slot + PRGraph filled by TaggerCheckNeutrino
 * (per-bundle "nu<N>" slots when nu_per_bundle is on) and writes
 * tracking-pr.root (T_rec_charge / T_proj_data / T_bad_ch / Trun, plus
 * T_cluster with save_in_scope) with the SAME schema as the SBND/uBooNE
 * writers so the converter and the Magnify-tracking GUI keep working.
 *
 * PDVD delta: the channel coordinate.  PDVD anodes are two-sided CRPs whose
 * two faces carry DIFFERENT (partly overlapping) channel sets, so the SBND
 * (plane, apa, wire) scheme -- which ignores the face -- collides.  The
 * global coordinate here is the per-plane RANK of the wire's physical channel
 * id over the whole detector (see ChanScheme), resolved per point from
 * PR::Fit::paf (apa, face).  PDVD v6 wires: nch = 3808 / 3808 / 4672.
 *
 * Runs as an IEnsembleVisitor inside the MABC pipeline, after
 * tagger_check_neutrino (pipeline name "tracking_visitor" in pdvd pr.jsonnet).
 */

#ifndef WIRECELLROOT_PDVDPRMAGNIFYTRACKINGVISITOR
#define WIRECELLROOT_PDVDPRMAGNIFYTRACKINGVISITOR

#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellUtil/Logging.h"

#include <array>
#include <cmath>
#include <map>
#include <set>
#include <tuple>
#include <vector>

class TFile;

namespace WireCell {
    namespace Root {

        class PdvdPrMagnifyTrackingVisitor : public IConfigurable, public Clus::IEnsembleVisitor {
           public:
            PdvdPrMagnifyTrackingVisitor();
            virtual ~PdvdPrMagnifyTrackingVisitor();

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
            // doc 99: T_cluster's flash columns from Cluster::get_matched_flash()
            // (matched_flash_gid vs the merge-safe "opflash" PC) instead of
            // Cluster::get_flash().  DEFAULT FALSE => byte-identical.
            //
            // PRESENT BUT UNVALIDATED HERE.  The fix was measured on SBND, whose
            // gid side is the anode ident and therefore unique across inputs.
            // PDVD's gid encoding (opflash_phys_gid / shared_flash, per-drift-side
            // flash lists) has NOT been checked against the uniqueness
            // precondition in Grouping::flash_by_gid(), so no PDVD config wires
            // this key.  Check that first, then wire it.
            bool m_flash_by_gid{false};
            double m_dQdx_scale{0.1};
            double m_dQdx_offset{-1000};
            bool m_flag_skip_vertex{false};
            // Readout length in ticks, used only to clamp dead-channel time
            // ranges (an unbounded dead x extent converts to +-1e9 ticks).
            int m_nticks{10000};

            // PDVD channel scheme (doc pdvd/25 M6).  PDVD anodes are two-sided
            // CRPs whose two faces carry DIFFERENT (partly overlapping) channel
            // sets, so the SBND (plane, apa, wire) scheme -- face-agnostic --
            // collides.  Here the global coordinate is the RANK of the wire's
            // physical channel id among all channels of that plane over the
            // whole detector: global = base[p] + rank_p(channel(apa, face, p, wire)),
            // built from the configured anodes at visit time.  Per-plane
            // contiguous, no collisions, and a strip shared by both faces of a
            // CRP maps to ONE coordinate.  PDVD v6 wires: nch = 3808/3808/4672.
            struct ChanScheme {
                std::array<int, 3> nch{0, 0, 0};
                std::array<int, 3> base{0, 0, 0};
                // (apa, face, plane) -> wire index -> per-plane compact channel index
                std::map<std::tuple<int, int, int>, std::vector<int>> wire2idx;
                int idx(int plane, int apa, int face, int wire) const {
                    auto it = wire2idx.find(std::make_tuple(apa, face, plane));
                    if (it == wire2idx.end() || it->second.empty()) return wire;  // unknown volume: pass through
                    const int n = static_cast<int>(it->second.size());
                    const int w = std::max(0, std::min(n - 1, wire));
                    return it->second[w];
                }
                int global(int plane, int apa, int face, int wire) const {
                    return base[plane] + idx(plane, apa, face, wire);
                }
                // Fractional wire position (the fitters keep it fractional).
                double globalf(int plane, int apa, int face, double wire) const {
                    const double fl = std::floor(wire);
                    return base[plane] + idx(plane, apa, face, static_cast<int>(fl)) + (wire - fl);
                }
            };
            ChanScheme chan_scheme() const;

            void write_bad_channels(TFile* output_tf, Clus::Facade::Grouping& grouping, const ChanScheme& cs) const;
            void write_proj_data(TFile* output_tf, Clus::Facade::Grouping& grouping, const ChanScheme& cs) const;
            void write_t_rec_data(TFile* output_tf, Clus::Facade::Grouping& grouping, const ChanScheme& cs) const;
            void write_trun(TFile* output_tf) const;
            void write_cluster_summary(TFile* output_tf, Clus::Facade::Grouping& grouping) const;
        };
    }  // namespace Root
}  // namespace WireCell

#endif
