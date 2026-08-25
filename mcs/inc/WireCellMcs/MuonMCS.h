/** MuonMCS -- muon momentum from Multiple Coulomb Scattering.
 *
 * ROOT-free port of the MicroBooNE standalone MCS estimator:
 *
 *   upstream: https://github.com/uboone/Multiple_Coulomb_Scattering
 *             src/mcs.{h,cxx} at 6aa0b9c576e22b954adf4472594a422929896c5d
 *   physics:  MicroBooNE collaboration, arXiv:2605.03048
 *
 * Upstream license (preserved per MIT terms):
 *
 *   MIT License
 *   Copyright (c) 2026 MicroBooNE Collaboration
 *
 *   Permission is hereby granted, free of charge, to any person obtaining a
 *   copy of this software and associated documentation files (the
 *   "Software"), to deal in the Software without restriction, including
 *   without limitation the rights to use, copy, modify, merge, publish,
 *   distribute, sublicense, and/or sell copies of the Software, and to
 *   permit persons to whom the Software is furnished to do so, subject to
 *   the following conditions: The above copyright notice and this
 *   permission notice shall be included in all copies or substantial
 *   portions of the Software.
 *
 * Port rules (doc sbnd_xin/docs/80_mcs-integration.md, "doc 80"):
 *  - Everything stays in upstream units: lengths cm, energies MeV internally
 *    (results in GeV where upstream published GeV).  Unit conversion is the
 *    caller's job (clus/ driver).
 *  - Only what ROOT forced is changed: TMatrixDEigen -> Eigen
 *    SelfAdjointEigenSolver, TGraph::Eval -> local Interp1D (linear with
 *    2-point extrapolation), TF1::GetMinimumX -> local grid-prescan + Brent
 *    reproducing the same recipe (npx=100, eps=1e-10 abs and rel,
 *    maxiter=100, up-to-10 bracket retries).
 *  - Known upstream defects are fixed behind per-bug McsOptions bools,
 *    default FIXED (doc 80 sec 6.4, owner decision 4).  All-false reproduces
 *    upstream bit-for-bit; McsCounters records every place a fix fired, so a
 *    zero counter set proves fixed == upstream on that input.
 *
 * This is a plain library, NOT a WCT plugin/component.  It must never appear
 * in a jsonnet `plugins:` list (doc 80 sec 3).
 */
#ifndef WIRECELLMCS_MUONMCS
#define WIRECELLMCS_MUONMCS

#include <functional>
#include <vector>

namespace WireCell::Mcs {

    /// One bool per upstream defect (doc 80 sec 6.4 numbering).  Default
    /// FIXED; flip all false for the bit-exact upstream behaviour the round-0
    /// fixtures record.
    struct McsOptions {
        bool fix_single_segment{true};      // #7: require >=2 fitted segments (>=1 real
                                            //     angle) else abort with -1; upstream
                                            //     minimises a constant and publishes it
        bool fix_vx_overflow{true};         // #8: |vx| >= 1 selects ivx=4; upstream's
                                            //     indicator sum falls through to ivx=0
        bool fix_nan_guards{true};          // #9+#10: probability floor before -log,
                                            //     acos clamp, drift-parallel projection
                                            //     plane fallback
        bool fix_bad_path_outputs{true};    // #12: mu_tracklen/emu_tracklen = -1 on
                                            //     bad_path; upstream publishes the
                                            //     start-end chord as a track length
        bool fix_remove_seg_by_index{true}; // #15: erase segment points by identity;
                                            //     upstream erases a contiguous block
                                            //     matched by exact float coordinates
    };

    /// Where each McsOptions fix actually fired on the last run.  All-zero
    /// means the fixed and upstream paths were arithmetically identical.
    struct McsCounters {
        int single_seg_abort{0};       // #7 fired: <2 fitted segments
        int ivx_overflow{0};           // #8 fired: a segment had |vx| >= 1
        int prob_floor{0};             // #9 fired: likelihood underflowed to 0
        int acos_clamp{0};             // #10 fired: |cos| > 1 clamped
        int degenerate_plane{0};       // #9 fired: drift-parallel prior segment
        int noncontiguous_removal{0};  // #15 fired: segment points not contiguous
        int bad_path{0};               // #12 relevant: trim failed start->end
    };

    struct McsResult {
        double mu_tracklen{-1};    // trimmed path length [cm]
        double emu_tracklen{-1};   // TOTAL energy [GeV] from the CSDA range of the
                                   // trimmed path (upstream's own PDG table --
                                   // deliberately NOT the toolkit range estimator,
                                   // doc 80 sec 8.3)
        double emu_MCS{-1};        // MCS TOTAL energy [GeV]; -1 = not computed
        double ambiguity_MCS{-1};  // max likelihood ratio of the two side minima
                                   // to the global one; 1 = most ambiguous
        int nsegs{0};              // fitted 14-cm segments
        bool bad_path{false};      // trim could not reach the end vertex
        McsCounters counters;
    };

    namespace detail {

        using VD = std::vector<double>;
        using VVD = std::vector<std::vector<double>>;

        /// Linear interpolation with 2-point extrapolation outside the range,
        /// reproducing TGraph::Eval (no spline) expression-for-expression.
        /// xs must be strictly increasing.
        double interp1d(const VD& xs, const VD& ys, double x);

        /// The upstream 20-point PDG CSDA table (mcs.cxx:437-446):
        /// KE [MeV] <-> residual range [cm].
        double ke_from_rr(double rr);
        double rr_from_ke(double ke);

        /// Grid-prescan + Brent minimisation reproducing the recipe behind
        /// TF1::GetMinimumX (ROOT's BrentMinimizer1D::Minimize): npx-point
        /// uniform scan brackets the minimum (strict <, earliest wins, bracket
        /// clamped to +-dx), Brent (golden section + parabolic steps,
        /// tol = eps*|x| + eps, convergence |x-m| <= 2 tol - (b-a)/2) refines,
        /// retried up to 10 times on non-convergence.
        struct Scan {
            double xmin0{0}, xmax0{0};    // range in
            int argmin{-1};               // grid argmin index
            double xmin1{0}, xmax1{0};    // bracket out
            double xmiddle{0};            // Brent start point
            VD grid_y;                    // the npx evaluations
        };
        struct Minimize1DResult {
            double x{0};        // minimiser location
            int nsearch{0};     // bracket+Brent rounds used (1 = converged first try)
            bool ok{false};
            Scan first_scan;    // the round-0 gate #3 quantity
        };
        Minimize1DResult minimize1d(const std::function<double(double)>& func,
                                    double xlow, double xup,
                                    int npx = 100, double eps = 1e-10, int maxiter = 100);

        /// Stage 1: Dijkstra-style trim of the input cloud to the muon path
        /// (upstream trim_trajectory, mcs.cxx:296-360).
        struct TrimResult {
            bool bad_path{false};
            VVD points_final;   // end vertex first, start vertex last (upstream order)
        };
        TrimResult trim_trajectory(const VVD& points,
                                   const VD& vtx_start, const VD& vtx_end);

        /// Stages 2+3: segmentation and inter-segment angles
        /// (upstream form_segs, mcs.cxx:365-424).
        struct SegResult {
            int nsegs_container{0};     // upstream segs.size() incl. a possible
                                        // trailing empty Track (bug #7's off-by-one)
            std::vector<int> seg_npoints;
            VD distance;                // from muon start to segment midpoint [cm]
            VD angle, angle_projB, angle_projC;  // [0] is the -1 sentinel
            VVD aAxes;                  // per-segment principal axis
            VVD COM;
            VVD track_axes;             // whole-track PCA axes (unflipped)
        };
        SegResult form_segs(const VVD& points, const VD& muon_start,
                            const VD& muon_end, double seg_length,
                            const McsOptions& opt, McsCounters& ctr);

        /// Stage 4: double-Gaussian likelihood minimised in KE
        /// (upstream estimate_energy, mcs.cxx:569-602).
        struct EnergyResult {
            double keguess{-1}, keguess_lower{-1}, keguess_higher{-1};  // [MeV]
            double l_keguess{0}, l_keguess_lower{0}, l_keguess_higher{0};
            double ambiguity{-1};
            Scan scans[3];      // first prescan of each of the three minimisations
        };
        EnergyResult estimate_energy(const VD& segs_distance, const VD& segs_angle_x,
                                     const VD& segs_angle_y, const VD& vx_comps,
                                     const McsOptions& opt, McsCounters& ctr);

        /// Single-angle likelihood terms, exposed for the round-4 pull test.
        double lnlikelihood_theta_xz(double angle, double T,
                                     const McsOptions& opt, McsCounters& ctr);
        double lnlikelihood_theta_yz(double angle, double T, double vx,
                                     const McsOptions& opt, McsCounters& ctr);
        /// {sigma1, sigma2, area_ratio} of the tuned double-Gaussian PDFs.
        VD pred_theta_xz_pars(double T);
        VD pred_theta_yz_pars(double T, int vx_index);

    }  // namespace detail

    class MuonMCS {
      public:
        explicit MuonMCS(const McsOptions& options = {})
            : m_options(options)
        {
        }

        /// Estimate the muon energy from trajectory shape.
        /// vtx_start/vtx_end: reconstructed muon endpoints {x,y,z} [cm].
        /// points: 3D cloud following the muon path [cm]; extra activity
        /// (delta rays, other prongs) is tolerated -- stage 1 trims it.
        McsResult run(const std::vector<double>& vtx_start,
                      const std::vector<double>& vtx_end,
                      const std::vector<std::vector<double>>& points) const;

        const McsOptions& options() const { return m_options; }

      private:
        McsOptions m_options;
    };

}  // namespace WireCell::Mcs

#endif
