/** Numeric core of the FD-VD low-energy (solar) charge-light matcher.
 *
 * A port of the frozen python matcher of fdvd_sim docs 06-10 and 15 (wcp
 * WireCell/wcp-porting-validation, fdvd_sim/stageB):
 *   ql_m2m_proto.py  shift_table / pred_at / ks_dis      (light prediction)
 *   ql_purity.py     build_groups / Resp / features / decide
 *   ql10_drift.py    drift_mask (the "veto3" drift spec)
 *   ql08.py          ARMS (P5 E100, P8 E0; Qc 100 ke, KS 0.3, centroid 200 cm, strict)
 * Each function cites the python lines it ports.  The arithmetic follows
 * numpy's order wherever that is reproducible (pairwise sums, np.interp,
 * float32 storage, NumPy 2 weak-scalar float32 comparisons), so the C++
 * features agree with python to the last float32 bit except where numpy
 * uses its own SIMD exp (lifetime factor) or BLAS (centroid matmul);
 * fdvd_sim doc 16 gates Q1/Q2.
 *
 * Units follow the python: cm, us, electrons, PE.
 */
#ifndef WIRECELL_MATCH_FDVDLOWE
#define WIRECELL_MATCH_FDVDLOWE

#include "WireCellMatch/PhotonLibraryModel.h"

#include <array>
#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace WireCell::Clus::Facade {
    class Grouping;
}

namespace WireCell::Match::FdvdLowE {

    constexpr int NCH = 184;   // FD-VD 1x8x14 OpDets (= opflash columns)

    // ql_m2m_proto.py lines 52-59; ql_purity.py line 34 (PMIN_STORE), 36 (QMIN)
    struct Constants {
        double vcm{0.160563};            // drift speed, cm/us
        double xa{325.07}, xc{-325.0};   // anode (collection plane) and cathode, true x [cm]
        double xw_cm{325.07};
        double ext_a{2.0}, ext_c{1.2};   // shift-table drift-window tolerances [cm]
        double tau_us{10400.0};          // electron lifetime
        double pe_per_e{8.8885e-3 / 0.891};
        double shift_step_cm{5.0};
        double twin_us{4300.0};          // flash time window
        double qmin_cluster{50e3};       // cluster kept for Q-L: blob charge sum >= this [e]
        int pmin_store{3};               // flash groups kept with >= this many lit OpDets
        double x_resp_cm{306.9};         // drift regressor label origin (dl_drift.py X_RESP)
    };

    // numpy float64 add.reduce on contiguous data (pairwise, block 128)
    double np_sum(const double* a, size_t n);
    // the same for float32 data with float32 accumulators (numpy FLOAT_pairwise_sum)
    float np_sum_f(const float* a, size_t n);

    // FdvdBlobTable columns (FdvdBlobTable.h)
    enum BlobCol { kOrder = 0, kAnode = 1, kWmin = 9, kWmax = 10, kX = 11, kY = 12, kZ = 13, kQ = 14, kStart = 15, kSpan = 16 };

    // The clustering tree's clusters as blob-table rows: each tree blob joined by its key (anode, face, slice
    // index range, u/v/w wire ranges = Aux::fill_scalar_blob) to the table; rows sorted in python blob order.
    struct JoinedCluster {
        int index{-1};                   // position in grouping.children()
        std::vector<int> rows;           // blob-table rows, python order
    };
    std::vector<JoinedCluster> join_clusters(const Clus::Facade::Grouping& grouping, const double* table, size_t nrow,
                                             size_t ncol, size_t& n_tree_blobs, size_t& n_unjoined);
    // numpy.interp (compiled_base.c arr_interp)
    double np_interp(double x, const std::vector<double>& xp, const std::vector<double>& fp, double left, double right);

    // One cluster prepared for matching (ql_m2m_proto.collect_one lines 141-154).
    struct Cluster {
        int index{-1};                   // position in the input grouping's children()
        int rep{-1};                     // representative blob = smallest blob-table order index
        double Q{0};                     // sum of blob charge (negatives included), python blob order
        double x_app{0};                 // charge-weighted apparent x, max(q,0) weights
        int nblob{0};
        // shift table (ql_m2m_proto.shift_table)
        double t_lo{0}, dt{0}, t_hi{0};
        std::vector<float> tab;          // [nstep * NCH], predicted PE without s0
        size_t nstep() const { return tab.size() / NCH; }
    };

    // ql_m2m_proto.py shift_table (lines 98-116).  P: blob centres (apparent cm), q: max(q, 0).
    // Returns false when the drift window is empty (python returns None).
    bool shift_table(const PhotonLibraryModel& lib, const Constants& C, const std::vector<std::array<double, 3>>& P,
                     const std::vector<double>& q, Cluster& c);

    // ql_m2m_proto.py pred_at (lines 179-187): table interpolated at time t -> out[NCH] (float64)
    void pred_at(const Cluster& c, double t, double* out);

    // ql_purity.py Resp.__call__ / fire (lines 141-168), parameters from the calibration
    struct Resp {
        std::vector<double> xc, y, f;
        double tail{0};
        double operator()(double p) const;
        double fire(double p) const;
    };

    struct Calibration {
        double s0{1};
        double rwin_lo{0}, rwin_hi{0};
        Resp resp;
        std::vector<double> drift_mu, drift_dhat, drift_shat;   // ql10_drift.calib
    };
    Calibration load_calibration(const std::string& json_path);

    // Flash groups (ql_purity.build_groups with delta 0 + ql08.gr_one storage, lines 155-200)
    struct Groups {
        std::vector<double> t;           // anchor time [us]
        std::vector<float> pe;           // [ng * NCH], stored float32
        std::vector<int> nflash;         // members
    };
    // opflash rows: time [ns], NCH PE columns
    Groups build_groups(const std::vector<double>& time_ns, const std::vector<double>& pe_rows, const Constants& C);

    // ql_purity.features (lines 171-206), tshift 0
    struct Features {
        std::vector<int> k, g;                   // cluster index (in the cluster list), group index (window subset)
        std::vector<float> r, ks, dc, npp;
        std::vector<double> gt, gtot;            // per window group
        std::vector<int> gnpd, gidx;             // per window group; gidx = index into Groups
    };
    Features features(const std::vector<Cluster>& cl, const Groups& G, const Calibration& cal, const Constants& C,
                      const std::vector<std::array<double, 3>>& pd_pos);

    // ql08.ARMS row
    struct Arm {
        std::string name;
        double qc{100e3};
        int P{5};
        double E{100};
        double ks_c{0.3}, dc_c{200.0};
    };
    std::vector<Arm> default_arms();

    // ql10_drift.drift_mask (lines 77-91) for "none" (veto_k <= 0) or "vetoK"; mu NaN = no cut
    std::vector<bool> drift_mask(const Features& Fe, double veto_k, const std::vector<double>& mu,
                                 const std::vector<Cluster>& cl, const Calibration& cal, const Constants& C);

    // ql_purity.decide (lines 209-238), strict uniqueness, with the pair mask applied first (ql10_drift.dec)
    std::map<int, int> decide(const Features& Fe, const std::vector<bool>& mask, const std::vector<Cluster>& cl,
                              const Arm& arm, const Calibration& cal);

}  // namespace WireCell::Match::FdvdLowE

#endif
