/** FdvdDriftRegressor: per-cluster drift distance from charge diffusion (FD-VD low energy).
 *
 * The doc 10 drift regressor M3 (fdvd_sim/stageB/dl_drift.py; network DNN_ROI_SP diffusion_t0 DriftRegressor,
 * W view) applied inside WCT through an ITensorForward service (Pytorch::TorchTensorSetService with the
 * TorchScript export of fdvd_sim/stageB/export_drift_ts.py).  For each clustering cluster whose blob charge
 * passes the cuts it builds the "own" crop of dl_drift.crops_event (lines 194-218):
 *   - the cluster's blobs in its dominant CRM (largest positive blob charge);
 *   - that CRM's SP "gauss" frame, W rows, |q|, float32, masked to the blobs' footprint (W wires +- 2,
 *     ticks +- 12 around [floor(start / tick), + ceil(span / tick)));
 *   - a 256 W x 1024 tick window centred on the masked charge centroid (round half to even), zero padded;
 * then log1p(crop) / 5 -> network -> mu x 100 = drift distance to the response plane [cm],
 * sigma = exp(0.5 clamp(logvar, -7, 7)) x 100.  Float32 throughout, as python; python ran the network in bf16
 * on GPU (doc 16 gate M3 reports the precision difference).
 *
 * A helper owned by FdvdLowEQLMatching (its "drift" configuration block), not a graph node: the matcher
 * already holds the clustering tree and the blob table, and WCT has no ITensorSet fanout to feed a second
 * node.  The SP frames are read from files, as python does, only for the CRMs that host a selected cluster:
 * frames = a path pattern with one %d for the anode ident (e.g. ".../sp-anode%d.tar.gz", FrameFileSink
 * archives with a "<frame_tag><N>" frame).
 * Selection (dl_drift.collect_one + ovl13 Chain.drift): signed blob charge sum >= qmin_all and positive
 * sum >= qmin.
 * Result: rows (representative blob = smallest blob-table order, mu, sigma); with dump_crops also the
 * own crops [n, 256, 1024] before log1p.
 */
#ifndef WIRECELL_MATCH_FDVDDRIFTREGRESSOR
#define WIRECELL_MATCH_FDVDDRIFTREGRESSOR

#include "WireCellIface/ITensorForward.h"
#include "WireCellUtil/Configuration.h"
#include "WireCellUtil/Logging.h"

#include <string>
#include <vector>

namespace WireCell::Clus::Facade {
    class Grouping;
}

namespace WireCell::Match {

    class FdvdDriftRegressor {
      public:
        explicit FdvdDriftRegressor(const WireCell::Configuration& cfg);

        struct Result {
            std::vector<double> rows;     // [n, 3] rep, mu, sigma
            std::vector<float> crops;     // [n, nch, ntk] when dump_crops
            size_t ncrop{0};
            size_t ncrm{0};
        };
        Result compute(const Clus::Facade::Grouping& grouping, const double* table, size_t nrow, size_t ncol,
                       int ident) const;
        bool dump_crops() const { return m_dump_crops; }
        int nch() const { return m_nch; }
        int ntk() const { return m_ntk; }

        // the |gauss| W rows [nw, ntick] of an SP frame archive (FrameFileSink layout), float32
        static std::vector<float> load_w_rows(const std::string& path, const std::string& frame_tag, int w0, int nw,
                                              size_t& ntick);

      private:
        std::string m_forward{"TorchTensorSetService"};
        std::string m_frames{""};
        std::string m_frame_tag{"gauss"};
        double m_qmin_all{10e3};     // dl_drift.py QMIN (signed charge sum)
        double m_qmin{45e3};         // ovl13 Chain.drift qmin (positive charge sum)
        double m_tick{500.0};        // ns, dl_drift TICK_NS
        int m_w0{572}, m_nw{292};    // first W row of a CRM frame, number of W wires
        int m_dw{2}, m_dt{12};       // footprint dilation
        int m_nch{256}, m_ntk{1024}; // crop
        int m_chunk{64}, m_batch{16};
        bool m_dump_crops{false};
        ITensorForward::pointer m_fwd;
        Log::logptr_t m_log;
    };

}  // namespace WireCell::Match

#endif
