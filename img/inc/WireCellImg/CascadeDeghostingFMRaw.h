/** CascadeDeghostingFMRaw: CascadeDeghostingFM with the head call replaced by a PCA projection of the per-view FM
    means at the final level (wcfm doc 53 part 43, the "D3 raw" arm; wcp pr_ml doc 30 sec 9).

    FORKED BY DUPLICATION from CascadeDeghostingFM (which stays byte-for-byte untouched, as CascadeDeghosting did for
    it).  Everything CascadeDeghostingFM does is done identically here -- the same three input ports (cluster, frame,
    FMFeatureExtract tensor set; the IClusterFrameTensorJoin interface is CascadeDeghostingFM's) and the same per-view
    means -- except the final level's extra columns.  At the final level (k = 1), for every cell i and plane p:
      e_p(i)   = the mean over the cell's channels c in plane p that have an FM pixel at (c, slice ident) of that
                 pixel's feature row, half-rounded (CascadeDeghostingFM's e_p), zero when none
      has_p(i) = 1 when at least one such pixel exists
      x_p(i)   = ((e_p(i) - mean) @ basis) * has_p(i)      (mean [dim], basis [dim x k] from the `pca` file)
    The level's xb (N, 15) becomes (N, 15 + 3k + 3) = [the 15 charge columns, x_0, x_1, x_2, has_0, has_1, has_2] and
    is handed to the final level's forward, whose model was trained on those columns.  No head ITensorForward is used.

    With `fmraw` absent (the default) the node is CascadeDeghosting with an ignored third port.

    Configuration: all of CascadeDeghosting's, plus
      fmraw    (absent / null = off) {pca: "<json file on WIRECELL_PATH: {mean[dim], basis[dim][k], k, dim}>",
                 fm_dim: 128}
    With dump_dir set the final level's dump also carries fmraw_cols (N, 3k+3), fm_means (N, 3, dim), fm_has (N, 3).
 */
#ifndef WIRECELLIMG_CASCADEDEGHOSTINGFMRAW
#define WIRECELLIMG_CASCADEDEGHOSTINGFMRAW

#include "WireCellAux/Logger.h"
#include "WireCellImg/CascadeDeghostingFM.h"   // IClusterFrameTensorJoin
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

    class CascadeDeghostingFMRaw : public Aux::Logger, public IClusterFrameTensorJoin, public IConfigurable {
      public:
        CascadeDeghostingFMRaw();
        virtual ~CascadeDeghostingFMRaw();

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
        // doc 53 part 43: the projected FM columns of the final level (absent key => off => CascadeDeghosting's output)
        bool m_fmraw{false};
        std::string m_pca_file{""};
        int m_fm_dim{128};
        int m_pca_k{0};
        std::vector<float> m_pca_mean;    // [fm_dim]
        std::vector<float> m_pca_basis;   // [fm_dim x k], row-major
        size_t m_count{0};
    };

}  // namespace WireCell::Img

#endif
