/** Dump everything the SBND PR event display needs into ONE JSON file.
 *
 * This is the "calib mode" of the PR chain, the direct analogue of QLMatching's
 * `calib_dump` that feeds sbnd_xin/ql_scan: a single self-contained per-event
 * artifact that a Bokeh viewer can read with nothing but the standard library.
 * See sbnd_xin/docs/pr/26_pr-event-display.md.
 *
 * Why a new component rather than reading what already exists:
 *
 *  - The Bee zip (mabc-pr.zip) carries track_fit / shower_track / vertices as
 *    flat x/y/z/q arrays.  It has no segment grouping (so no polylines), no
 *    Steiner layer, and no 2-D projection at all.
 *  - tracking-pr.root carries the 2-D projection in T_proj_data -- except that
 *    tree comes out with only its `cluster_id` branch whenever the PR chain
 *    also runs `tagger_output`, which reopens the file in UPDATE mode and
 *    drops the vector<vector<int>> StreamerInfo on close.  See doc pr/26.
 *
 * So the display would need three readers, one of which is broken.  This
 * writes the union instead, in the coordinate conventions doc pr/7 fixed:
 * positions in cm, `pu/pv/pw` kept FRACTIONAL per-APA wire indices, and time as
 * a slice index (`pt / nticks_per_slice`).
 *
 * Stage 2 (doc pr/26 sec 7) adds what makes an event JUDGEABLE rather than just
 * visible: a `shower_id` on every segment, `showers`, `kine` and `tagger`.  The
 * particle-flow TREE is deliberately NOT rebuilt here -- the display reads it
 * from the Bee zip's `mc.json`, which MultiAlgBlobClustering::fill_bee_pf_tree
 * already writes.  `shower_id` uses that producer's node encoding exactly
 * (`cluster_id*1000 + start_segment id`) so it is the join key from an mc.json
 * shower node to every segment of that shower.
 *
 * DEFAULT OFF.  Like every other entry in the SBND `cm_by_name` table this is
 * only instantiated when named in `pipeline_names`, so with the name absent the
 * compiled configuration -- and hence every production output -- is unchanged.
 *
 * Runs as an IEnsembleVisitor inside the MABC pipeline, after
 * tagger_check_neutrino (pipeline name "pr_display" in the SBND config).  It
 * reads only; it mutates nothing.
 */

#ifndef WIRECELLCLUS_PRDISPLAYDUMP
#define WIRECELLCLUS_PRDISPLAYDUMP

#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellAux/Logger.h"

#include <array>
#include <string>
#include <vector>

namespace WireCell::Clus {

    class PrDisplayDump : public Aux::Logger, public IConfigurable, public Clus::IEnsembleVisitor {
       public:
        PrDisplayDump();
        virtual ~PrDisplayDump();

        // IConfigurable
        virtual void configure(const WireCell::Configuration& config);
        virtual Configuration default_configuration() const;

        // IEnsembleVisitor
        virtual void visit(Facade::Ensemble& ensemble) const;

       private:
        std::string m_output_filename{"calib-pr.json"};
        std::string m_grouping_name{"live"};
        int m_runNo{0}, m_subRunNo{0}, m_eventNo{0};
        double m_dQdx_scale{0.1};
        double m_dQdx_offset{-1000};
        int m_nticks{3427};
        // Drop 2-D cells whose measured charge is at or below this, purely to
        // keep the JSON small.  0 = keep everything.
        double m_proj_charge_min{0};
        bool m_pretty{false};
        std::vector<IAnodePlane::pointer> m_anodes;
        IDetectorVolumes::pointer m_dv{nullptr};

        // Per-plane per-TPC channel counts, the same derivation
        // SbndPrMagnifyTrackingVisitor uses.  Recorded in the dump so the
        // viewer can label its wire axes without a geometry file; the dump
        // itself stores PER-APA wire indices, not global channels.
        struct ChanScheme {
            std::array<int, 3> nch{0, 0, 0};
            std::array<int, 3> base{0, 0, 0};
        };
        ChanScheme chan_scheme() const;

        Configuration dump_meta(Facade::Grouping& grouping, const ChanScheme& cs) const;
        Configuration dump_graph(Facade::Grouping& grouping) const;      // segments + vertices
        Configuration dump_showers(Facade::Grouping& grouping) const;    // PR::Shower rows
        Configuration dump_kine(Facade::Grouping& grouping) const;       // KineInfo
        Configuration dump_tagger(Facade::Grouping& grouping) const;     // TaggerInfo
        Configuration dump_track_shower(Facade::Grouping& grouping) const;
        Configuration dump_steiner(Facade::Grouping& grouping) const;
        Configuration dump_proj(Facade::Grouping& grouping) const;
        Configuration dump_dead(Facade::Grouping& grouping) const;
    };
}  // namespace WireCell::Clus

#endif
