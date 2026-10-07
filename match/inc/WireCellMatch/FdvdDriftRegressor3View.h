/** FdvdDriftRegressor3View: per-cluster drift distance from charge diffusion, three views (FD-VD low energy).
 *
 * Counterpart of FdvdDriftRegressor (collection view only), forked by duplication: that file is untouched.
 * Selected by "views": 3 in the "drift" block of FdvdLowEQLMatching; absent or 1 = FdvdDriftRegressor, unchanged.
 *
 * The network is the three-view model of fdvd_sim doc 35 (DNN_ROI_SP diffusion_t0 FusionRegressor, mode "mask"),
 * served through an ITensorForward service with the TorchScript export of fdvd_sim/stage35/export_e4_ts.py.
 * That export computes the wire-crossing tables itself, so the toolkit passes only crops and their origins.
 *
 * Per selected cluster (the selection, the dominant CRM and the footprint rule are FdvdDriftRegressor's):
 *   - for each view U, V, W: the CRM's SP "gauss" frame rows of that view, |q|, float32, masked to the blobs'
 *     footprint in that view (the blob's own wire range in the view +- 2 wires, ticks +- 12);
 *   - a 256 wire x 1024 tick window centred on that view's masked charge centroid (round half to even);
 *   - an induction view without charge is all zeros and takes the collection crop's origin; a cluster without
 *     collection charge keeps three zero crops with origin 0 (it is still inferred, as in the one-view class);
 *   - U and V are moved to the collection crop's tick origin (shift by t0_view - t0_W, cut to the window);
 *   - x = log1p(crop) / 5, [n, 3, 256, 1024]; meta = (CRM ident, absolute channel of the first row of the U, V, W
 *     crop), [n, 4] float32;
 *   - network -> mu x 100 = drift distance to the response plane [cm], sigma = exp(0.5 clamp(logvar, -7, 7)) x 100.
 * This is the python of fdvd_sim/stage35/apply35.py (collect_one, cmd_infer) and DNN_ROI_SP
 * validation/score_val_uvw.py fusion_input.
 *
 * Result: rows (representative blob, mu, sigma); with dump_crops also the three own crops per cluster
 * [n, 3, 256, 1024] before the tick shift and before log1p, and their origins [n, 7]
 * (CRM, c0_U, t0_U, c0_V, t0_V, c0_W, t0_W; view-local rows, an empty induction view is flagged by -1e9 in t0).
 */
#ifndef WIRECELL_MATCH_FDVDDRIFTREGRESSOR3VIEW
#define WIRECELL_MATCH_FDVDDRIFTREGRESSOR3VIEW

#include "WireCellIface/ITensorForward.h"
#include "WireCellUtil/Configuration.h"
#include "WireCellUtil/Logging.h"

#include <array>
#include <string>
#include <vector>

namespace WireCell::Clus::Facade {
    class Grouping;
}

namespace WireCell::Match {

    class FdvdDriftRegressor3View {
      public:
        explicit FdvdDriftRegressor3View(const WireCell::Configuration& cfg);

        struct Result {
            std::vector<double> rows;      // [n, 3] rep, mu, sigma
            std::vector<float> crops;      // [n, 3, nch, ntk] when dump_crops
            std::vector<double> origins;   // [n, 7] when dump_crops
            size_t ncrop{0};
            size_t ncrm{0};
        };
        Result compute(const Clus::Facade::Grouping& grouping, const double* table, size_t nrow, size_t ncol,
                       int ident) const;
        bool dump_crops() const { return m_dump_crops; }
        int nch() const { return m_nch; }
        int ntk() const { return m_ntk; }

        // One footprint box of a view: wires [w0, w1), ticks [t0, t1), before dilation.
        using box_t = std::array<long, 4>;

        // The masked crop of one view.  img: [nrow, ntick] |q|.  out: [nch, ntk], zero filled here.
        // Returns false (out stays zero, c0 = t0 = 0) when the footprint holds no charge.
        static bool view_crop(const float* img, int nrow, size_t ntick, const std::vector<box_t>& boxes, int dw, int dt,
                              int nch, int ntk, float* out, long& c0, long& t0);

        // dst[:, lo:hi] = src[:, lo - d : hi - d], lo = max(0, d), hi = min(ntk, ntk + d); the rest zero.
        static void shift_ticks(const float* src, long d, int nch, int ntk, float* dst);

        // all |gauss| rows [nrow, ntick] of an SP frame archive (FrameFileSink layout), float32; ch0 = first channel
        static std::vector<float> load_rows(const std::string& path, const std::string& frame_tag, size_t& nrow,
                                            size_t& ntick, int& ch0);

      private:
        std::string m_forward{"TorchTensorSetService"};
        std::string m_frames{""};
        std::string m_frame_tag{"gauss"};
        double m_qmin_all{10e3};     // as FdvdDriftRegressor
        double m_qmin{45e3};
        double m_tick{500.0};        // ns
        std::array<int, 3> m_row0{0, 286, 572};    // first row of U, V, W in a CRM frame
        std::array<int, 3> m_nrow{286, 286, 292};  // wires per view
        int m_dw{2}, m_dt{12};       // footprint dilation
        int m_nch{256}, m_ntk{1024}; // crop
        int m_chunk{32}, m_batch{16};
        bool m_dump_crops{false};
        ITensorForward::pointer m_fwd;
        Log::logptr_t m_log;
    };

}  // namespace WireCell::Match

#endif
