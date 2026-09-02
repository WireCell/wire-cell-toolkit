/** Visitor that writes PDVD STM-stage tracking data to a ROOT file.
 *
 * Fork-by-duplication of SbndMagnifyTrackingVisitor (doc pdvd/25 M6; the SBND
 * file is untouched), itself a fork of UbooneMagnifyTrackingVisitor (which reads the
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
 *  - PDVD channel convention: per-plane RANK of the physical channel id
 *    (face-aware; see ChanScheme below).  The SBND (plane, apa, wire)
 *    scheme ignores the face and collides on PDVD's two-sided CRPs.
 *
 * Runs as an IEnsembleVisitor inside the MABC pipeline, after
 * tagger_check_stm (pipeline name "stm_magnify" in the PDVD pr.jsonnet).
 */

#ifndef WIRECELLROOT_PDVDMAGNIFYTRACKINGVISITOR
#define WIRECELLROOT_PDVDMAGNIFYTRACKINGVISITOR

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

        class PdvdMagnifyTrackingVisitor : public IConfigurable, public Clus::IEnsembleVisitor {
           public:
            PdvdMagnifyTrackingVisitor();
            virtual ~PdvdMagnifyTrackingVisitor();

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
            void write_stm_trees(TFile* output_tf, Clus::Facade::Grouping& grouping) const;
            void write_trun(TFile* output_tf) const;
        };
    }  // namespace Root
}  // namespace WireCell

#endif
