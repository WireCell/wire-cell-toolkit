/** The coarse-to-fine learned deghosting on a set of blobs (CascadeDeghosting, wcfm doc 14), as a function.

    Extracted from CascadeDeghosting::operator() (pdvd doc 130 phase 2) so that the imaging node and the clustering
    retiler (ImproveCluster_2 retile_deghost) run one code path: level by level, build the level graph
    (CascadeGraph), run the level's model (ITensorForward), prune or cut (coarse levels) or keep, repair and
    guard (the final level).  The caller builds the SliceCharge (from a cluster graph, a frame or a list of slices)
    and assembles whatever it needs from the kept cells.  The imaging node's output (the cluster graph with
    Img::geom_clustering edges) stays in img.
 */
#ifndef WIRECELLAUX_CASCADERUN
#define WIRECELLAUX_CASCADERUN

#include "WireCellAux/CascadeGraph.h"
#include "WireCellAux/CellSteiner.h"
#include "WireCellAux/Logger.h"
#include "WireCellIface/IBlob.h"
#include "WireCellIface/ITensorForward.h"

#include <functional>
#include <string>
#include <vector>

namespace WireCell::Aux::Cascade {

    /// One level of the cascade.  width: the cut width in wires (0 = the uncut input, level 0 only);
    /// superwire: channels per (super-)wire node; threshold: the logit below which a node is a candidate
    /// (coarse level) or dropped (final level); forward: the level's model.
    struct LevelSpec {
        int width{0};
        int superwire{1};
        double threshold{0.0};
        ITensorForward::pointer forward;
        std::string name;   // for messages only
    };

    struct RunParams {
        CutParams cut;                 // nudge, max_depth, min_length (CascadeDeghosting cut_*)
        bool guard{true};              // coverage guard on the coarse levels
        int max_level_nodes{0};        // cap on the nodes of the next level (0 = none)
        bool repair{true};             // Steiner repair of the final keep set
        SteinerParams steiner;         // p_term, q_floor, budget (set budget from default_repair_budget when <= 0)
        std::string policy{"uboone"};  // geom_clustering policy of the bb edges
        int nthreads{1};
        bool iso{false};               // wcfm doc 17 dense-ambiguous-slice fallback
        IsoParams iso_params;
        bool final_guard{false};       // wcfm doc 22
        std::string dump_dir;          // per-level npz dumps (cascade-<ident>-L<n>.npz), "" = none
        int ident_base{0};             // first ident of the cut cells
        // pdvd doc 130: length of a final-level edge between two cells, for the length pricing of the repair
        // (steiner.len_cost / steiner.max_bridge_len).  Empty = no lengths, the repair of wcfm doc 14.
        std::function<double(const IBlob::pointer&, const IBlob::pointer&)> edge_length;
    };

    struct RunResult {
        std::vector<IBlob::pointer> kept;   // the final-level cells kept, in node order
        size_t nin{0};                      // input blobs
        size_t nkeep_thr{0};                // final-level cells kept by the threshold alone
        SteinerResult srep;                 // the repair's counts (empty when repair is off)
    };

    /// The repair budget CascadeDeghosting derives from the final level's threshold when none is configured:
    /// three P = 0.2 cells plus one cell at the threshold.
    double default_repair_budget(double final_threshold);

    /// Run the cascade on `blobs` with the charge `sc` (which must index every blob's slice).  `ident` names the
    /// cluster in messages and dumps; `call` is the caller's call count (messages only).  Logs on `log` exactly
    /// as the imaging node did.
    RunResult run_cascade(const std::vector<IBlob::pointer>& blobs, const SliceCharge& sc,
                          const std::vector<LevelSpec>& levels, const RunParams& par, int ident,
                          Log::logptr_t log, int call = 0);

}  // namespace WireCell::Aux::Cascade

#endif
