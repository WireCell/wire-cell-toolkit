/** CascadeDeghostingFM: CascadeDeghosting (wcfm doc 14) with the stage-A head-score columns at the final level
    (wcfm docs 27-30).

    FORKED BY DUPLICATION from CascadeDeghosting (the production node, which stays byte-for-byte untouched; the
    owner's choice of 2026-10-03, doc 30).  Everything CascadeDeghosting does is done identically here; the one
    addition is a third input port and, when `head` is configured, three more per-cell input columns for the
    final level's model:

      port 0  the cluster (uncut blobs)                         } as CascadeDeghosting
      port 1  the SP frame the slices were made from            }
      port 2  the FM feature tensor set of the same anode and frame: FMFeatureExtract's output (per plane a
              coords (N, 2) i4 tensor [channel ident, slice index] and a feat (N, 128) f4 or feat_half (N, 128)
              u2 tensor).  Only the active pixels are listed.

    At the final level (k = 1, the 4-wire cells) and only there, for every cell i and plane p:
      e_p(i)   = the mean over the cell's channels c in plane p that have an FM pixel at (c, slice ident) of that
                 pixel's feature row (half-rounded, as the training graphs' b_fm), zero when none
      has_p(i) = 1 when at least one such pixel exists
      score(i) = head(e_U, e_V, e_W, has)   -- the ITensorForward `head.forward` (the TorchScript CrossHead,
                 docs 21/27: inputs eU f4[N,128], eV, eW, has f4[N,3]; output score f4[N])
      per wire node w with an FM pixel (plane, channel, slice): the log-softmax of score over the cells on w
      ls_mean(i), ls_min(i) = the mean and the minimum of those log-softmax values over the cell's edges to such
                 wires; a cell with none gets the level's column minimum over the other cells (doc 27 sec 1)
    The level's xb (N, 15) becomes (N, 18) = [the 15 charge columns, score, ls_mean, ls_min] and is handed to the
    final level's forward, whose model was trained on those 18 columns (doc 29 AB arm, scripts/d30_export.py).
    The log-softmax runs in double, as the python reference (wcfm scripts/d27_head.py head_scores).

    With `head` absent (the default) the node is CascadeDeghosting with an ignored third port.

    Configuration: all of CascadeDeghosting's, plus
      head     (absent / null = off) {forward: "<type>:<name>" of the head ITensorForward, fm_dim: 128,
                 chunk: 0}.  chunk (wcfm doc 34) = the number of cells per head forward; 0 (the default) hands the
                 whole final level to one forward, as doc 30.  The head scores each cell from its own row, so a
                 chunked forward bounds the head's working memory (2.3 GB at 4e5 cells in one batch) at the price
                 of a batch-size-dependent BLAS path: the scores may differ in the last bits (doc 34 has the gate).
 */
#ifndef WIRECELLIMG_CASCADEDEGHOSTINGFM
#define WIRECELLIMG_CASCADEDEGHOSTINGFM

#include "WireCellAux/Logger.h"
#include "WireCellIface/IBlobSet.h"
#include "WireCellIface/ICluster.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IFrame.h"
#include "WireCellIface/IFunctionNode.h"
#include "WireCellIface/IJoinNode.h"
#include "WireCellIface/ITensorForward.h"
#include "WireCellIface/ITensorSet.h"

#include <string>
#include <vector>

namespace WireCell::Img {

    /// (cluster, frame, FM feature tensor set) -> cluster
    class IClusterFrameTensorJoin : public IJoinNode<std::tuple<ICluster, IFrame, ITensorSet>, ICluster> {
      public:
        virtual ~IClusterFrameTensorJoin() {}
        virtual std::string signature() { return typeid(IClusterFrameTensorJoin).name(); }
    };

    class CascadeDeghostingFM : public Aux::Logger, public IClusterFrameTensorJoin, public IConfigurable {
      public:
        CascadeDeghostingFM();
        virtual ~CascadeDeghostingFM();

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
        // false = legacy: the slice's first tick is (slice start - frame time) / tick.  true: slice start / tick, for
        // slicers with frame-relative slice starts (MaskSlice) on frames with a non-zero time (pdvd doc 122 sec 3).
        bool m_slice_start_relative{false};
        double m_uncer_cut{1e11};
        int m_ident_base{1 << 20};
        std::string m_dump_dir{""};
        int m_nthreads{1};
        bool m_iso{false};
        bool m_final_guard{false};
        bool m_keep_slices{false};
        double m_iso_nmin{500}, m_iso_mmin{4}, m_iso_amin{0.8}, m_iso_t{-1.5}, m_iso_amb_lo{-3.0}, m_iso_amb_hi{0.5};
        // doc 30: the head-score columns of the final level (absent key => off => CascadeDeghosting's output)
        bool m_head{false};
        std::string m_head_tn{""};
        int m_fm_dim{128};
        int m_head_chunk{0};   // doc 34: cells per head forward; 0 = one forward for the level
        ITensorForward::pointer m_head_forward;
        size_t m_count{0};
    };

}  // namespace WireCell::Img

#endif
