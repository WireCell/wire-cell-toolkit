/** CascadeDeghosting: coarse-to-fine learned deghosting of an imaging cluster (wcfm doc 14).

    Replaces ProjectionDeghosting + InSliceDeghosting in a chain

      tiling -> BlobClustering (uncut blobs) -> CascadeDeghosting -> BlobGrouping -> ChargeSolving -> ...
                                   SP frame -/  (port 1)

    It is a join node: port 0 the cluster, port 1 the frame the slices were made from.  The models were trained
    on the frame's charge (0.25 x the 4-tick sum of the gauss trace per channel and slice, all ticks), which the
    slice activity does not carry (it sums only the ticks that pass the slicing threshold); with a null frame,
    or no trace of charge_tag, the node falls back to the activity (x charge_scale) and says so.

    The input cluster's blobs (uncut) are level 0.  For each level in turn:
      1. level > 0: every survivor of the previous level is cut to that level's width by BlobCutting's
         bisection (duplicated in CascadeCut.cxx), its depth counted from the uncut blob, so the final level's
         cells are exactly those of one direct BlobCutting to the final width (the models' training tier);
      2. the level graph is built (Img::Cascade::build_level: blob nodes, super-wire nodes channel // k,
         bw / bb / bb_in / ww edges) and handed to the level's ITensorForward (a TorchScript e2-projq model
         run by TorchTensorSetService), which returns a logit and a charge estimate per blob;
      3. a coarse level prunes: blobs with logit < threshold, in ascending logit, each only if every
         (super-)wire it touches keeps another survivor (the coverage guard);
         the final level keeps blobs with logit >= its threshold and, optionally, repairs connectivity with
         Img::Cascade::steiner_repair (terminal filter, Voronoi bridges, Kruskal).
    The kept final blobs are emitted as a fresh cluster graph built exactly as BlobClustering builds one
    (s-b, b-w, w-c edges and geom_clustering b-b edges), with idents renumbered from ident_base.
    Only one level graph is alive at a time.

    This node uses only the ITensorForward interface, so img does not link torch; the job must load the
    plugin that provides the forward service (WireCellPytorch).

    Configuration:
      levels            array of {width: 0 (level 0, uncut) or the cut width in wires, superwire: k (1 = real
                        wires), forward: "<type>:<name>" of an ITensorForward, threshold: logit}
      cut_max_depth, cut_min_length, cut_nudge (10, 2, 0.01) the bisection (BlobCutting's max_depth,
                        counted from the uncut blob, min_length, nudge)
      guard             (bool, true) coverage guard on the coarse levels
      max_level_nodes   (int, 0 = off) cap on a level's blob count: survivors are cut in descending logit
                        and the lowest-logit survivors are dropped once the next level would exceed the cap
      repair            (bool, true) Steiner repair of the final keep set
      repair_p_term, repair_q_floor (0.8, 1e4 e), repair_budget (<= 0: 3 (-log 0.2) - log sigmoid(final threshold))
      policy            ("uboone") geom_clustering policy for bb and for the output graph
      charge_tag        ("") the frame's trace tag holding the charge (e.g. "gauss2"); "" = all traces
      charge_scale      (0.25) the training charge unit: 0.25 x the 4-tick sum (frame or activity)
      uncer_cut         (1e11) activity entries with larger uncertainty (dummy / masked) are ignored
      ident_base        (1 << 20) first ident of the output blobs
      dump_dir          ("" = off) write one npz per level and event (graph arrays, logits, decisions)
      nthreads          (1) threads for the C++ level work that splits into independent parts -- the bisection of
                        the survivors to the next width, the node features and the cross-slice bb search -- run over contiguous
                        blocks and merged in order, so the output does not depend on it (wcfm doc 15 round 3;
                        the forward's threads are libtorch's own, OMP_NUM_THREADS)
 */
#ifndef WIRECELLIMG_CASCADEDEGHOSTING
#define WIRECELLIMG_CASCADEDEGHOSTING

#include "WireCellAux/Logger.h"
#include "WireCellIface/IBlobSet.h"
#include "WireCellIface/ICluster.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IFrame.h"
#include "WireCellIface/IFunctionNode.h"
#include "WireCellIface/IJoinNode.h"
#include "WireCellIface/ITensorForward.h"

#include <string>
#include <vector>

namespace WireCell::Img {

    /// (cluster, frame) -> cluster
    class IClusterFrameJoin : public IJoinNode<std::tuple<ICluster, IFrame>, ICluster> {
      public:
        virtual ~IClusterFrameJoin() {}
        virtual std::string signature() { return typeid(IClusterFrameJoin).name(); }
    };

    class CascadeDeghosting : public Aux::Logger, public IClusterFrameJoin, public IConfigurable {
      public:
        CascadeDeghosting();
        virtual ~CascadeDeghosting();

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        virtual bool operator()(const input_tuple_type& intup, output_pointer& out);

        struct LevelCfg {
            std::string forward_tn;
            int width{0};
            int superwire{1};
            double threshold{0};
            ITensorForward::pointer forward;
        };

      private:
        std::vector<LevelCfg> m_levels;
        int m_cut_max_depth{10};
        int m_cut_min_length{2};
        double m_cut_nudge{0.01};
        bool m_guard{true};
        int m_max_level_nodes{0};
        bool m_repair{true};
        double m_repair_p_term{0.8};
        double m_repair_q_floor{1e4};
        double m_repair_budget{0};
        std::string m_policy{"uboone"};
        std::string m_charge_tag{""};
        double m_charge_scale{0.25};
        double m_uncer_cut{1e11};
        int m_ident_base{1 << 20};
        std::string m_dump_dir{""};
        int m_nthreads{1};
        size_t m_count{0};
    };

}  // namespace WireCell::Img

#endif
