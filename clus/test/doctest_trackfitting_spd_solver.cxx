// pr_ml doc 87b -- TrackFitting parameter spd_solver (default 0 = OFF).
//
// The four normal-equation solves of TrackFitting (fit_point, trajectory_fit,
// dQ_dx_multi_fit, dQ_dx_fit) have symmetric positive (semi-)definite matrices
// A = R^T M R (+ F^T F) but were solved by Eigen's BiCGSTAB at its defaults
// (tolerance at machine epsilon, cap 2n).  On a dense, ill-conditioned dQ/dx
// system (the slow FD-HD tuning event nue_b11099_e2, n = 25,024) BiCGSTAB ran to
// its cap and returned an answer off by up to 680 times the largest exact value.
// spd_solver > 0 switches the four solves to TrackFitting::spd_cg_ic_solve:
// ConjugateGradient with an IncompleteCholesky (AMD) preconditioner at
// spd_solver_tol (default 1e-10), the legacy solve running if it fails.
// These cases pin:
//   * the defaults are the legacy path (spd_solver 0, tolerance 1e-10), and the
//     set/get round trip;
//   * the defect on a small smear-plus-smoothness system of the same kind
//     (a 24 x 24 grid, n = 576): BiCGSTAB at the toolkit's settings misses the
//     direct (SimplicialLDLT) answer by more than 10 % of its largest value,
//     spd_cg_ic_solve converges and matches it to 1e-4;
//   * the guards: a column with no stored entries (no crash; that unknown keeps
//     its guess), a zero right-hand side, a NaN in A (returns false, x = x0, so
//     the caller runs the legacy solve), an empty or mis-sized system (false).
//   * singular 3 x 3 position systems (fit_point / trajectory_fit with a plane or
//     a coordinate unconstrained: rank 2 and rank 1, with empty columns): the
//     helper converges and gives the legacy BiCGSTAB solve's own answer, the
//     null-space components keeping their initial values -- no change of
//     behaviour where the legacy solve did not fail.

#include "WireCellUtil/doctest.h"
#include "WireCellClus/TrackFitting.h"

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace WireCell;
using namespace WireCell::Clus;
using SpMat = Eigen::SparseMatrix<double>;

namespace {
// A = R^T R + F^T F, b = R^T R x_true: R a Gaussian smear over a m x m grid of
// unknowns (radius rad, width sig), F a first-difference regulariser (weight lam).
void smear_system(int m, int rad, double sig, double lam, SpMat& A, Eigen::VectorXd& b)
{
    const int n = m * m;
    std::vector<Eigen::Triplet<double>> T;
    for (int a = 0; a < m; ++a)
        for (int c = 0; c < m; ++c)
            for (int da = -rad; da <= rad; ++da)
                for (int dc = -rad; dc <= rad; ++dc) {
                    const int a2 = a + da, c2 = c + dc;
                    if (a2 < 0 || c2 < 0 || a2 >= m || c2 >= m) continue;
                    T.emplace_back(a * m + c, a2 * m + c2, std::exp(-0.5 * (da * da + dc * dc) / (sig * sig)));
                }
    SpMat R(n, n);
    R.setFromTriplets(T.begin(), T.end());
    T.clear();
    for (int a = 0; a < m; ++a)
        for (int c = 0; c < m; ++c) {
            const int i = a * m + c;
            T.emplace_back(i, i, -lam);
            T.emplace_back(i, c + 1 < m ? i + 1 : i, lam);
        }
    SpMat F(n, n);
    F.setFromTriplets(T.begin(), T.end());
    Eigen::VectorXd xt(n);
    for (int i = 0; i < n; ++i) xt(i) = 1.0 + 0.5 * std::sin(i * 0.05) + (i % 97 == 0 ? 3.0 : 0.0);
    SpMat Rt = R.transpose(), Ft = F.transpose();
    A = Rt * R + Ft * F;
    b = Rt * (R * xt);
}
}  // namespace

TEST_CASE("pr_ml doc 87b spd_solver: defaults are the legacy path, set / get round trip")
{
    TrackFitting tf;
    CHECK(tf.get_parameter("spd_solver") == 0.0);
    CHECK(tf.get_parameter("spd_solver_tol") == 1e-10);
    TrackFitting::Parameters p;
    CHECK(p.spd_solver == 0.0);
    CHECK(p.spd_solver_tol == 1e-10);
    tf.set_parameter("spd_solver", 1.0);
    tf.set_parameter("spd_solver_tol", 1e-8);
    CHECK(tf.get_parameter("spd_solver") == 1.0);
    CHECK(tf.get_parameter("spd_solver_tol") == 1e-8);
}

TEST_CASE("pr_ml doc 87b spd_solver: BiCGSTAB at the toolkit's settings misses an ill-conditioned SPD system; CG + IC solves it")
{
    SpMat A;
    Eigen::VectorXd b;
    smear_system(24, 3, 1.5, 1e-2, A, b);
    const int n = A.rows();
    REQUIRE(n == 576);
    const Eigen::VectorXd x0 = Eigen::VectorXd::Ones(n);

    Eigen::SimplicialLDLT<SpMat> ldlt(A);
    REQUIRE(ldlt.info() == Eigen::Success);
    const Eigen::VectorXd xs = ldlt.solve(b);
    const double ms = xs.cwiseAbs().maxCoeff();
    REQUIRE(ms > 0);
    CHECK((A * xs - b).norm() / b.norm() < 1e-12);

    // the legacy call, as at the four sites: defaults, solveWithGuess
    Eigen::BiCGSTAB<SpMat> bi;
    bi.compute(A);
    const Eigen::VectorXd xb = bi.solveWithGuess(b, x0);
    CHECK(bi.iterations() == 2 * n);                          // ran to its cap
    CHECK((xb - xs).cwiseAbs().maxCoeff() / ms > 0.1);        // the defect: far from the direct answer

    Eigen::VectorXd x;
    long iters = -1;
    double err = -1;
    REQUIRE(TrackFitting::spd_cg_ic_solve(A, b, x0, 1e-10, x, iters, err));
    CHECK(iters > 0);
    CHECK(iters < 2 * n);                                      // converged before the cap
    CHECK(err <= 1e-10);
    CHECK((x - xs).cwiseAbs().maxCoeff() / ms < 1e-4);
}

TEST_CASE("pr_ml doc 87b spd_solver: guards")
{
    Eigen::VectorXd x;
    long iters = -1;
    double err = -1;

    SUBCASE("a column with no stored entries keeps its guess")
    {
        std::vector<Eigen::Triplet<double>> T{{0, 0, 4.0}, {0, 1, 1.0}, {1, 0, 1.0}, {1, 1, 3.0}};
        SpMat A(3, 3);
        A.setFromTriplets(T.begin(), T.end());   // row / column 2 empty
        Eigen::VectorXd b(3), x0(3);
        b << 1.0, 2.0, 0.0;
        x0 << 0.0, 0.0, 7.0;
        REQUIRE(TrackFitting::spd_cg_ic_solve(A, b, x0, 1e-12, x, iters, err));
        CHECK(x(2) == 7.0);
        CHECK(x(0) == doctest::Approx(1.0 / 11.0));     // [4 1; 1 3] x = [1 2]
        CHECK(x(1) == doctest::Approx(7.0 / 11.0));
    }
    SUBCASE("a zero right-hand side gives zero")
    {
        std::vector<Eigen::Triplet<double>> T{{0, 0, 2.0}, {1, 1, 2.0}};
        SpMat A(2, 2);
        A.setFromTriplets(T.begin(), T.end());
        Eigen::VectorXd b = Eigen::VectorXd::Zero(2), x0 = Eigen::VectorXd::Ones(2);
        REQUIRE(TrackFitting::spd_cg_ic_solve(A, b, x0, 1e-10, x, iters, err));
        CHECK(x.isZero(0));
    }
    SUBCASE("a NaN in A returns false with x = x0")
    {
        std::vector<Eigen::Triplet<double>> T{{0, 0, std::numeric_limits<double>::quiet_NaN()}, {1, 1, 2.0}};
        SpMat A(2, 2);
        A.setFromTriplets(T.begin(), T.end());
        Eigen::VectorXd b = Eigen::VectorXd::Ones(2), x0(2);
        x0 << 3.0, 4.0;
        CHECK_FALSE(TrackFitting::spd_cg_ic_solve(A, b, x0, 1e-10, x, iters, err));
        CHECK(x == x0);
    }
    SUBCASE("an empty or mis-sized system returns false")
    {
        SpMat A0(0, 0);
        Eigen::VectorXd e(0);
        CHECK_FALSE(TrackFitting::spd_cg_ic_solve(A0, e, e, 1e-10, x, iters, err));
        SpMat A(2, 2);
        A.setIdentity();
        Eigen::VectorXd b = Eigen::VectorXd::Ones(3), x0 = Eigen::VectorXd::Zero(2);
        CHECK_FALSE(TrackFitting::spd_cg_ic_solve(A, b, x0, 1e-10, x, iters, err));
        CHECK(x == x0);
    }
}

TEST_CASE("pr_ml doc 87b spd_solver: singular 3x3 position systems give the legacy solve's answer")
{
    struct Case { const char* name; double a[9]; };
    const Case cases[] = {
        {"rank 2: x from time rows, y + z from one plane", {1, 0, 0, 0, 1, 1, 0, 1, 1}},
        {"rank 1: x only (y and z columns empty)", {4, 0, 0, 0, 0, 0, 0, 0, 0}},
        {"rank 2: x and y (z column empty)", {2, 0, 0, 0, 3, 0, 0, 0, 0}},
    };
    for (const auto& c : cases) {
        CAPTURE(c.name);
        std::vector<Eigen::Triplet<double>> T;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                if (c.a[3 * i + j] != 0) T.emplace_back(i, j, c.a[3 * i + j]);
        SpMat A(3, 3);
        A.setFromTriplets(T.begin(), T.end());
        Eigen::VectorXd xt(3), x0(3);
        xt << 1.0, 2.0, 3.0;
        x0 << 0.5, 1.5, 1.0;
        const Eigen::VectorXd b = A * xt;

        Eigen::BiCGSTAB<SpMat> bi;   // the legacy call
        bi.compute(A);
        const Eigen::VectorXd xb = bi.solveWithGuess(b, x0);
        REQUIRE_FALSE(std::isnan(bi.error()));

        Eigen::VectorXd x;
        long iters = -1;
        double err = -1;
        REQUIRE(TrackFitting::spd_cg_ic_solve(A, b, x0, 1e-10, x, iters, err));
        CHECK((A * x - b).norm() <= 1e-12 * std::max(1.0, b.norm()));
        CHECK((x - xb).cwiseAbs().maxCoeff() < 1e-12);
    }
}
