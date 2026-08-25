// ROOT-free port of uboone/Multiple_Coulomb_Scattering src/mcs.{h,cxx}
// at 6aa0b9c (MIT, notice in the header).  See the header for the port rules.
//
// Fidelity discipline (doc 80 sec 6.2): expression FORMS are kept verbatim --
// the covariance triple loop, norm's pow-then-sqrt, the d+d^2 edge score --
// because a last-bit difference in an accumulation order moves an eigenvector
// at 1e-15, flips a point across a segment boundary, and shifts the energy by
// O(100 MeV).  Do not "clean up" arithmetic here.
//
// Upstream oddities kept deliberately (doc 80 sec 6.4, M7/M15):
//  - trim_trajectory scores endpoint proximity by sqrt(norm(diff(...))) --
//    a double square root (mcs.cxx:307-308).  Monotone, arg-min unchanged.
//  - fitPCA weights enter the covariance squared and the normalisation is
//    /track.N, wrong for non-unit weights (mcs.cxx:187-196).  Latent: nothing
//    ever passes a weight != 1.  NOT fixed -- it would break the golden.
//  - get_seg's half-open projection window and the flexible-length
//    segmentation loop are kept exactly.

#include "WireCellMcs/MuonMCS.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <tuple>

using namespace WireCell::Mcs;
using WireCell::Mcs::detail::VD;
using WireCell::Mcs::detail::VVD;

// ---------------------------------------------------------------------------
// constants (upstream mcs.h:236-240)
static const double seg_length = 14;  // cm; = X0(LAr)/rho, structural NOT tunable
                                      // (doc 80 sec 7.1)
static const double Mmu = 105.658;    // MeV

// ---------------------------------------------------------------------------
// vector helpers (upstream mcs.h:20-51) -- expression forms verbatim
namespace {

    VD diff(VD u, const VD& v)
    {
        // upstream (mcs.h:22-25)
        for (int i = 0, n = u.size(); i < n; i++) { u[i] -= v[i]; }
        return u;
    }

    double norm(const VD& v)
    {
        // upstream (mcs.h:27-31): pow-accumulate then sqrt, in this order
        double norm = 0;
        for (int i = 0, n = v.size(); i < n; i++) { norm += std::pow(v[i], 2); }
        return std::sqrt(norm);
    }

    double dot(const VD& u, const VD& v)
    {
        // upstream (mcs.h:33-37)
        double val = 0;
        for (int i = 0, n = u.size(); i < n; i++) { val += u[i] * v[i]; }
        return val;
    }

    VD cross(const VD& u, const VD& v)
    {
        // upstream (mcs.h:39-41)
        return { u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0] };
    }

    VD scale(VD v, double a)
    {
        // upstream (mcs.h:43-46)
        for (int i = 0, n = v.size(); i < n; i++) { v[i] *= a; }
        return v;
    }

}  // namespace

// ---------------------------------------------------------------------------
// interpolation: TGraph::Eval equivalent (gate #1: bitwise)
namespace WireCell::Mcs::detail {

    double interp1d(const VD& xs, const VD& ys, double x)
    {
        // Reproduces TGraph::Eval (no spline) for strictly-increasing xs:
        // exact-node x returns ys[i]; outside the range the nearest TWO points
        // extrapolate (this is why WireCell::irrterp, which CLAMPS, is unusable
        // here -- doc 80 sec 6.1); interpolation is anchored at the UPPER
        // neighbour: yn = y[up] + (x - x[up]) * (y[low] - y[up]) / (x[low] - x[up]).
        const int n = xs.size();
        if (n == 0) return 0;
        if (n == 1) return ys[0];
        int low = -1, up = -1;
        // binary search for the last xs[i] <= x
        {
            int lo = 0, hi = n - 1;
            if (x < xs[0]) { low = -1; }
            else if (x >= xs[n - 1]) { low = n - 1; }
            else {
                while (hi - lo > 1) {
                    int mid = (lo + hi) / 2;
                    if (xs[mid] <= x) lo = mid;
                    else hi = mid;
                }
                low = lo;
            }
        }
        if (low == -1) low = 0;  // below range: first two points extrapolate
        if (xs[low] == x) return ys[low];
        if (low == n - 1) low--;  // above range: last two points extrapolate
        up = low + 1;
        if (xs[low] == xs[up]) return ys[low];
        double yn = ys[up] + (x - xs[up]) * (ys[low] - ys[up]) / (xs[low] - xs[up]);
        return yn;
    }

    // the CSDA table (upstream setUKEfromRR, mcs.cxx:437-446).  uEnergy is
    // KINETIC energy in MeV despite upstream's "Total" comment (doc 80 sec 6.4
    // #13: a 10 MeV total is impossible for a 105.658 MeV muon, and run() adds
    // Mmu to it).
    namespace {
        const double table_rho = 1.396;  // LAr density [g/cm3] (mcs.h:240)
        struct CsdaTable {
            VD rr, ke;
            CsdaTable()
            {
                const int un = 20;
                const double uCSDA[un] = { .9833, 1.786, 3.321, 6.598, 10.58, 30.84, 42.50,
                                           67.32, 106.3, 172.5, 238.5, 493.4, 616.3, 855.2,
                                           1202., 1758., 2297., 4359., 5354., 7298. };
                const double uEnergy[un] = { 10.0,  14.0,  20.0,  30.0,   40.0,  80.0,  100.0,
                                             140.0, 200.0, 300.0, 400.0,  800.0, 1000., 1400.,
                                             2000., 3000., 4000., 8000., 10000., 14000. };
                rr.resize(un);
                ke.assign(uEnergy, uEnergy + un);
                for (int i = 0; i < un; i++) { rr[i] = uCSDA[i] / table_rho; }  // cm (mcs.cxx:443)
            }
        };
        const CsdaTable& csda()
        {
            static const CsdaTable t;  // upstream cleanUp() leaked and then bricked
                                       // the object (bug #11); compile-time-constant
                                       // tables need no lifecycle at all
            return t;
        }
    }  // namespace

    double ke_from_rr(double rr) { return interp1d(csda().rr, csda().ke, rr); }
    double rr_from_ke(double ke) { return interp1d(csda().ke, csda().rr, ke); }

}  // namespace WireCell::Mcs::detail

// ---------------------------------------------------------------------------
// 1-D minimisation reproducing TF1::GetMinimumX's recipe (gate #3/#6).
// Grid prescan + Brent (R.P. Brent, "Algorithms for Minimization Without
// Derivatives", ch. 5 / Numerical Recipes 10.2), written from the recipe:
// npx-point uniform scan, strict < so the EARLIEST grid minimum wins, bracket
// clamped to argmin +- dx; Brent with golden ratio constant, parabolic step
// acceptance |p| < |q r / 2| within (a,b), tol = eps |x| + eps, convergence
// |x - m| <= 2 tol - (b - a)/2; on non-convergence the bracket+Brent pair is
// retried (up to 10 times) on the shrunken interval.
namespace WireCell::Mcs::detail {

    namespace {

        Scan minim_step(const std::function<double(double)>& func,
                        double xmin, double xmax, int npx)
        {
            Scan sc;
            sc.xmin0 = xmin;
            sc.xmax0 = xmax;
            double dx = (xmax - xmin) / (npx - 1);
            double xxmin = xmin;
            double yymin = func(xxmin);
            sc.grid_y.reserve(npx);
            sc.grid_y.push_back(yymin);
            sc.argmin = 0;
            for (int i = 1; i <= npx - 1; i++) {
                double x = xmin + i * dx;
                double y = func(x);
                sc.grid_y.push_back(y);
                if (y < yymin) {
                    xxmin = x;
                    yymin = y;
                    sc.argmin = i;
                }
            }
            sc.xmin1 = std::max(xmin, xxmin - dx);
            sc.xmax1 = std::min(xmax, xxmin + dx);
            sc.xmiddle = std::min(xxmin, sc.xmax1);
            return sc;
        }

        double minim_brent(const std::function<double(double)>& func,
                           double& xmin, double& xmax, double xmiddle,
                           bool& ok, double epsabs, double epsrel, int itermax)
        {
            const double c = 3.81966011250105097e-01;  // (3 - sqrt(5))/2, golden section
            double u, v, w, x, fv, fu, fw, fx, e, p, q, r, t2, d = 0, m, tol;
            v = w = x = xmiddle;
            e = 0;
            double a = xmin;
            double b = xmax;
            fv = fw = fx = func(x);
            for (int i = 0; i < itermax; i++) {
                m = 0.5 * (a + b);
                tol = epsrel * (std::fabs(x)) + epsabs;
                t2 = 2 * tol;
                if (std::fabs(x - m) <= (t2 - 0.5 * (b - a))) {
                    ok = true;
                    return x;
                }
                if (std::fabs(e) > tol) {
                    // parabolic interpolation through (v,fv) (w,fw) (x,fx)
                    r = (x - w) * (fx - fv);
                    q = (x - v) * (fx - fw);
                    p = (x - v) * q - (x - w) * r;
                    q = 2 * (q - r);
                    if (q > 0) p = -p;
                    else q = -q;
                    r = e;
                    e = d;
                    if (std::fabs(p) >= std::fabs(0.5 * q * r) || p <= q * (a - x) || p >= q * (b - x)) {
                        // reject the parabolic step: golden-section instead
                        e = (x >= m ? a - x : b - x);
                        d = c * e;
                    }
                    else {
                        d = p / q;
                        u = x + d;
                        if (u - a < t2 || b - u < t2) d = (m - x >= 0) ? std::fabs(tol) : -std::fabs(tol);
                    }
                }
                else {
                    e = (x >= m ? a - x : b - x);
                    d = c * e;
                }
                u = (std::fabs(d) >= tol ? x + d : x + ((d >= 0) ? std::fabs(tol) : -std::fabs(tol)));
                fu = func(u);
                if (fu <= fx) {
                    if (u < x) b = x;
                    else a = x;
                    v = w; fv = fw;
                    w = x; fw = fx;
                    x = u; fx = fu;
                }
                else {
                    if (u < x) a = u;
                    else b = u;
                    if (fu <= fw || w == x) {
                        v = w; fv = fw;
                        w = u; fw = fu;
                    }
                    else if (fu <= fv || v == x || v == w) {
                        v = u; fv = fu;
                    }
                }
            }
            ok = false;
            xmin = a;
            xmax = b;
            return x;
        }

    }  // namespace

    Minimize1DResult minimize1d(const std::function<double(double)>& func,
                                double xlow, double xup, int npx, double eps, int maxiter)
    {
        Minimize1DResult mr;
        double xmin = xlow, xmax = xup;
        const int max_search = 10;
        int nsearch = 0;
        bool ok = false;
        double x = 0;
        while (!ok) {
            if (nsearch > max_search) break;
            Scan sc = minim_step(func, xmin, xmax, npx);
            if (nsearch == 0) mr.first_scan = sc;
            xmin = sc.xmin1;
            xmax = sc.xmax1;
            x = minim_brent(func, xmin, xmax, sc.xmiddle, ok, eps, eps, maxiter);
            nsearch++;
        }
        mr.x = x;
        mr.nsearch = nsearch;
        mr.ok = ok;
        return mr;
    }

}  // namespace WireCell::Mcs::detail

// ---------------------------------------------------------------------------
// stage-1 graph structures (upstream mcs.h:53-194)
namespace {

    // upstream (mcs.h:72-138).  Determinism additions per doc 80 sec 6.3:
    // every sort key gets the point id as tie-break.  Upstream is
    // order-deterministic but order-SENSITIVE; on the round-0 fixtures no tie
    // fires, so these are provably inert there.
    struct Point {
        int id;
        bool exhausted = false;
        double score;
        VD pos;
        VD edges_dist;
        Point* prior = nullptr;
        std::vector<Point*> edges;

        Point(int id_, const VD& pos_)
            : id(id_)
            , score(1e9)
            , pos(pos_)
        {
        }

        // upstream (mcs.h:87-92): d + d^2 favours short hops; exponents tuned
        // for ~0.6 cm point spacing (which SBND's fitter shares, doc 80 sec 9.3)
        double get_dist_score(const Point& p) const
        {
            double dist = norm(diff(pos, p.pos));
            return dist + std::pow(dist, 2);
        }

        // upstream (mcs.h:94-103)
        void update_path(int edge_index)
        {
            Point* edge = edges[edge_index];
            double new_score = score + edges_dist[edge_index];
            if (new_score < edge->score) {
                edge->score = new_score;
                edge->prior = this;
                edge->exhausted = false;
            }
        }

        // upstream (mcs.h:105-109)
        void update_paths()
        {
            for (int i = 0, n = edges.size(); i < n; i++) { update_path(i); }
            exhausted = true;
        }

        // upstream (mcs.h:111-119); tie-break on id (sec 6.3 item 2)
        void sort_edges()
        {
            std::vector<std::tuple<double, int, Point*>> edges_tuples;
            for (int i = 0, n = edges.size(); i < n; i++) {
                edges_tuples.push_back(std::make_tuple(edges_dist[i], edges[i]->id, edges[i]));
            }
            std::sort(edges_tuples.begin(), edges_tuples.end(),
                      [](const auto& t1, const auto& t2) {
                          if (std::get<0>(t1) != std::get<0>(t2)) return std::get<0>(t1) < std::get<0>(t2);
                          return std::get<1>(t1) < std::get<1>(t2);
                      });
            for (int i = 0, n = edges_tuples.size(); i < n; i++) {
                edges_dist[i] = std::get<0>(edges_tuples[i]);
                edges[i] = std::get<2>(edges_tuples[i]);
            }
        }

        // upstream (mcs.h:121-136).  The "only keep 4 closest" comment there
        // is stale -- nedges_max = 20 (doc 80 sec 6.4).
        void add_edge(Point* p, const VD& dir)
        {
            int nedges_max = 20;
            bool in_dir = dot(diff(p->pos, pos), dir) > 0;
            double dist = get_dist_score(*p);
            int nedges = edges_dist.size();
            if (in_dir && (nedges < nedges_max || dist < edges_dist.back())) {
                edges.push_back(p);
                edges_dist.push_back(dist);
                sort_edges();
                if (nedges >= nedges_max) {
                    edges.pop_back();
                    edges_dist.pop_back();
                }
            }
        }
    };

    // upstream (mcs.h:141-143) with an id tie-break among equal-score
    // non-exhausted points (sec 6.3 item 3)
    struct Comparator {
        bool operator()(const Point* p1, const Point* p2) const
        {
            if (!p1->exhausted && !p2->exhausted && p1->score == p2->score) {
                return p1->id < p2->id;
            }
            return !p1->exhausted && (p2->exhausted || (p1->score < p2->score));
        }
    };

    // upstream (mcs.h:145-194) with an id per point so segment removal can be
    // identity-based (bug #15 fix)
    struct Track {
        VVD points;
        VD weights;
        std::vector<int> ids;
        double total_weight = 0;
        int N = 0;

        void add_point(const VD& point, int id, double weight = 1.)
        {
            points.push_back(point);
            weights.push_back(weight);
            ids.push_back(id);
            total_weight += weight;
            N++;
        }

        // upstream remove_seg (mcs.h:164-175): match the segment's first point
        // by EXACT float coordinates, erase seg_n contiguous entries.
        // With fix_remove_seg_by_index the segment's points are erased by id
        // identity instead (bug #15: when the refit axis reorders points the
        // segment is not contiguous in track order and the contiguous erase
        // removes the wrong points; first_index == -1 is UB upstream).
        void remove_seg(const Track& seg, bool by_index, McsCounters& ctr)
        {
            if (by_index) {
                // contiguity check purely for the counter
                int first_index = -1;
                for (int i = 0; i < N; i++) {
                    if (ids[i] == seg.ids.front()) {
                        first_index = i;
                        break;
                    }
                }
                bool contiguous = (first_index >= 0 && first_index + seg.N <= N);
                if (contiguous) {
                    for (int j = 0; j < seg.N; j++) {
                        if (ids[first_index + j] != seg.ids[j]) {
                            contiguous = false;
                            break;
                        }
                    }
                }
                if (!contiguous) ctr.noncontiguous_removal++;

                std::vector<bool> drop(N, false);
                double removed_weight = 0;
                for (int sid : seg.ids) {
                    for (int i = 0; i < N; i++) {
                        if (!drop[i] && ids[i] == sid) {
                            drop[i] = true;
                            removed_weight += weights[i];
                            break;
                        }
                    }
                }
                VVD np;
                VD nw;
                std::vector<int> ni;
                for (int i = 0; i < N; i++) {
                    if (drop[i]) continue;
                    np.push_back(points[i]);
                    nw.push_back(weights[i]);
                    ni.push_back(ids[i]);
                }
                int removed = N - (int)np.size();
                points.swap(np);
                weights.swap(nw);
                ids.swap(ni);
                total_weight -= removed_weight;
                N -= removed;
                return;
            }

            // upstream-exact path (mcs.h:164-186)
            const VD& first_pos = seg.points.front();
            int first_index = -1;
            for (int i = 0; i < N; i++) {
                const VD& point_pos = points[i];
                if (point_pos[0] != first_pos[0] || point_pos[1] != first_pos[1] ||
                    point_pos[2] != first_pos[2]) {
                    continue;
                }
                first_index = i;
                break;
            }
            if (first_index < 0 || first_index + seg.N > N) {
                // upstream would index weights[-1] / erase past the end (UB).
                // No defined behaviour exists to reproduce; drop nothing.
                ctr.noncontiguous_removal++;
                return;
            }
            double removed_weight = 0;
            for (int i = 0; i < seg.N; i++) { removed_weight += weights[first_index + i]; }
            points.erase(points.begin() + first_index, points.begin() + first_index + seg.N);
            weights.erase(weights.begin() + first_index, weights.begin() + first_index + seg.N);
            ids.erase(ids.begin() + first_index, ids.begin() + first_index + seg.N);
            total_weight -= removed_weight;
            N -= seg.N;
        }
    };

    // upstream ComparePCAProjection (mcs.h:54-69): sort points by projection on
    // an axis.  Tie-break on the point id (sec 6.3 item 1 -- the
    // highest-leverage tie: get_seg slices a PREFIX of this order, so a tie
    // flip changes segment membership).  Points and ids sorted together;
    // upstream did not co-sort weights and neither do we (all weights are 1).
    void sort_points(Track& track, const VD& axis)
    {
        std::vector<int> order(track.N);
        for (int i = 0; i < track.N; i++) order[i] = i;
        std::vector<double> proj(track.N);
        for (int i = 0; i < track.N; i++) proj[i] = dot(track.points[i], axis);
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            if (proj[a] != proj[b]) return proj[a] < proj[b];
            return track.ids[a] < track.ids[b];
        });
        VVD np(track.N);
        std::vector<int> ni(track.N);
        for (int i = 0; i < track.N; i++) {
            np[i] = track.points[order[i]];
            ni[i] = track.ids[order[i]];
        }
        track.points.swap(np);
        track.ids.swap(ni);
    }

}  // namespace

// ---------------------------------------------------------------------------
namespace WireCell::Mcs::detail {

    TrimResult trim_trajectory(const VVD& points, const VD& vtx_start, const VD& vtx_end)
    {
        // upstream (mcs.cxx:296-360)
        TrimResult out;
        const int npoints_traj = points.size();

        VD muon_dir_reco = diff(vtx_end, vtx_start);
        int startpoint_index = -1;
        int endpoint_index = -1;
        double startpoint_dist = 1e9;
        double endpoint_dist = 1e9;
        for (int i = 0; i < npoints_traj; i++) {
            const VD& pos = points[i];
            // upstream (mcs.cxx:307-308): sqrt of norm (itself already a
            // sqrt).  Monotone, arg-min unchanged; kept verbatim.
            double dist_start = std::sqrt(norm(diff(pos, vtx_start)));
            double dist_end = std::sqrt(norm(diff(pos, vtx_end)));
            if (dist_start < startpoint_dist) {
                startpoint_index = i;
                startpoint_dist = dist_start;
            }
            if (dist_end < endpoint_dist) {
                endpoint_index = i;
                endpoint_dist = dist_end;
            }
        }

        // upstream (mcs.cxx:318-326); ids stay 1-based input indices.  The
        // never-pushed Point upstream allocates for i == startpoint_index (a
        // leak per event, doc 80 sec 6.4) is simply not allocated here.
        std::vector<std::unique_ptr<Point>> storage;
        storage.reserve(npoints_traj);
        std::vector<Point*> traj_points;
        traj_points.reserve(npoints_traj);
        storage.push_back(std::make_unique<Point>(startpoint_index + 1, points[startpoint_index]));
        traj_points.push_back(storage.back().get());
        traj_points[0]->score = 0;
        for (int i = 0; i < npoints_traj; i++) {
            if (i == startpoint_index) continue;
            storage.push_back(std::make_unique<Point>(i + 1, points[i]));
            traj_points.push_back(storage.back().get());
        }

        // upstream (mcs.cxx:327-328): O(N^2) full graph build
        for (int i = 0; i < npoints_traj; i++) {
            for (int j = 0; j < npoints_traj; j++) { traj_points[i]->add_edge(traj_points[j], muon_dir_reco); }
        }

        // upstream (mcs.cxx:329-342): relax-then-sort loop
        bool bad_path = false;
        while (true) {
            traj_points.front()->update_paths();
            std::sort(traj_points.begin(), traj_points.end(), Comparator());
            if (traj_points.front()->id == (endpoint_index + 1)) { break; }
            if (traj_points.front()->exhausted) {
                bad_path = true;
                break;
            }
        }

        // upstream (mcs.cxx:343-350): walk prior pointers back from the end
        VVD trajectory_points_final;
        Point* path_pointer = traj_points[0];
        while (true && !bad_path) {
            trajectory_points_final.push_back(path_pointer->pos);
            path_pointer = path_pointer->prior;
            if (path_pointer == nullptr) { break; }
        }

        // upstream (mcs.cxx:352-354): the true endpoints are re-inserted; the
        // list runs end vertex -> start vertex
        trajectory_points_final.insert(trajectory_points_final.begin(), vtx_end);
        trajectory_points_final.insert(trajectory_points_final.end(), vtx_start);

        out.bad_path = bad_path;
        out.points_final = std::move(trajectory_points_final);
        return out;
    }

    namespace {

        // upstream fitPCA (mcs.cxx:169-214).  TMatrixDEigen -> Eigen
        // SelfAdjointEigenSolver (the covariance is symmetric by construction);
        // Eigen returns eigenvalues ASCENDING where ROOT's are descending, so
        // reorder explicitly, and pin each eigenvector's sign with a stable
        // rule (largest-|component| positive) BEFORE any muon-direction flip.
        // The flip is what fixes the physical sign of the axes that matter.
        VVD fit_pca(const Track& track, VD& com, VD& evals)
        {
            com = { 0, 0, 0 };
            evals.assign(3, 0.0);
            for (int i = 0; i < track.N; i++) {
                com[0] += track.points[i][0] * track.weights[i];
                com[1] += track.points[i][1] * track.weights[i];
                com[2] += track.points[i][2] * track.weights[i];
            }
            com = scale(com, 1. / track.total_weight);

            VVD meanPoints(track.N);
            for (int i = 0; i < track.N; i++) {
                meanPoints[i].resize(3);
                meanPoints[i][0] = (track.points[i][0] - com[0]) * track.weights[i];
                meanPoints[i][1] = (track.points[i][1] - com[1]) * track.weights[i];
                meanPoints[i][2] = (track.points[i][2] - com[2]) * track.weights[i];
            }

            // upstream (mcs.cxx:192-198): the covariance triple loop, kept
            // verbatim (NOT X^T X -- accumulation order is load-bearing).
            // Weights enter squared and the normalisation is /N: latent bug
            // #14, kept (M15).
            double covData[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    for (int k = 0; k < track.N; k++) { covData[3 * i + j] += meanPoints[k][i] * meanPoints[k][j]; }
                    covData[3 * i + j] = covData[3 * i + j] / track.N;
                }
            }

            Eigen::Matrix3d cov;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) { cov(i, j) = covData[3 * i + j]; }
            }
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(cov);

            VVD axes = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
            for (int a = 0; a < 3; a++) {
                int col = 2 - a;  // ascending -> descending
                Eigen::Vector3d v = solver.eigenvectors().col(col);
                // stable sign: largest-|component| positive
                int imax = 0;
                for (int i = 1; i < 3; i++) {
                    if (std::abs(v[i]) > std::abs(v[imax])) imax = i;
                }
                if (v[imax] < 0) v = -v;
                for (int i = 0; i < 3; i++) { axes[a][i] = v[i]; }
                evals[a] = solver.eigenvalues()[col];
            }
            evals = scale(evals, 1. / norm(evals));
            return axes;
        }

        // upstream get_seg (mcs.cxx:217-229): select the points within
        // seg_length of the first point, measured by projection on axis.
        // Half-open window; slices a PREFIX of the current sort order.
        Track get_seg(const Track& track, const VD& axis, double seg_len)
        {
            Track track_seg;
            double first_proj = dot(track.points.front(), axis);
            for (int i = 0; i < track.N; i++) {
                double proj = dot(track.points[i], axis);
                if (proj < (first_proj + seg_len)) { track_seg.add_point(track.points[i], track.ids[i]); }
            }
            return track_seg;
        }

        // upstream fitSegPCA (mcs.cxx:238-292)
        bool fit_seg_pca(Track& track, const VD& aAxis, Track& seg, VVD& segFitVector,
                         VD& segCOM, double seg_len, double& currentFirstPointProj,
                         double& currentLastPointProj, const McsOptions& opt, McsCounters& ctr)
        {
            seg.points.clear();
            seg.weights.clear();
            seg.ids.clear();
            seg.total_weight = 0;
            seg.N = 0;

            Track track_buffer = get_seg(track, aAxis, seg_len * 2);

            double lastPointProj = dot(track_buffer.points.back(), aAxis);
            currentFirstPointProj = lastPointProj;
            for (int i = 0; i < int(track_buffer.points.size()); i++) {
                double proj = dot(track_buffer.points[i], aAxis);
                if ((proj > currentLastPointProj) && (proj < currentFirstPointProj)) { currentFirstPointProj = proj; }
            }

            seg = get_seg(track_buffer, aAxis, seg_len);
            if (seg.N <= 1) { return false; }

            VD eigenvalues;
            segFitVector = fit_pca(seg, segCOM, eigenvalues);
            if (dot(segFitVector[0], aAxis) < 0) {
                segFitVector[0][0] *= -1;
                segFitVector[0][1] *= -1;
                segFitVector[0][2] *= -1;
            }

            sort_points(track_buffer, segFitVector[0]);

            seg = get_seg(track_buffer, segFitVector[0], seg_len);
            if (seg.N <= 1) { return false; }

            segFitVector = fit_pca(seg, segCOM, eigenvalues);
            if (dot(segFitVector[0], aAxis) < 0) {
                segFitVector[0][0] *= -1;
                segFitVector[0][1] *= -1;
                segFitVector[0][2] *= -1;
            }

            sort_points(seg, segFitVector[0]);
            currentLastPointProj = dot(seg.points.back(), aAxis);

            sort_points(seg, aAxis);
            track.remove_seg(seg, opt.fix_remove_seg_by_index, ctr);
            return true;
        }

        // upstream get_angle (mcs.cxx:136-140), acos clamp behind bug #10
        double get_angle(const VD& v1, const VD& v2, const McsOptions& opt, McsCounters& ctr)
        {
            double l1 = norm(v1);
            double l2 = norm(v2);
            double c = dot(v1, v2) / (l1 * l2);
            if (opt.fix_nan_guards && (c > 1.0 || c < -1.0)) {
                ctr.acos_clamp++;
                c = (c > 1.0) ? 1.0 : -1.0;
            }
            return (std::acos(c));
        }

        // upstream setSegAngles (mcs.cxx:143-166).  vecy_plane = aAxis_prior x
        // x-hat is the ZERO vector for a drift-parallel prior segment
        // (bug #9's NaN); with fix_nan_guards the projection plane falls back
        // to aAxis_prior x y-hat, which is well-defined exactly there.
        // bAxis_prior/cAxis_prior are read and never used upstream
        // (mcs.cxx:147-148, dead per doc 80 sec 6.4 #17) -- dropped.
        void set_seg_angles(const VVD& priorSegFitVector, const VD& segFitVector,
                            VD& angle_vec, VD& angleProjX_vec, VD& angleProjY_vec,
                            const McsOptions& opt, McsCounters& ctr)
        {
            const VD& aAxis_prior = priorSegFitVector[0];

            VD vecy_plane = cross(aAxis_prior, { 1, 0, 0 });
            if (opt.fix_nan_guards && norm(vecy_plane) < 1e-12) {
                ctr.degenerate_plane++;
                vecy_plane = cross(aAxis_prior, { 0, 1, 0 });
            }
            VD vecx_plane = cross(aAxis_prior, vecy_plane);
            vecx_plane = scale(vecx_plane, 1. / norm(vecx_plane));
            vecy_plane = scale(vecy_plane, 1. / norm(vecy_plane));
            VD projX = diff(segFitVector, scale(vecy_plane, dot(segFitVector, vecy_plane)));
            VD projY = diff(segFitVector, scale(vecx_plane, dot(segFitVector, vecx_plane)));

            projX = scale(projX, 1. / norm(projX));

            int dirX = 1 - 2 * (dot(segFitVector, vecx_plane) < 0);
            int dirY = 1 - 2 * (dot(segFitVector, vecy_plane) < 0);

            angle_vec.push_back(get_angle(segFitVector, aAxis_prior, opt, ctr));
            angleProjX_vec.push_back(get_angle(projX, aAxis_prior, opt, ctr) * dirX);
            angleProjY_vec.push_back(get_angle(projY, aAxis_prior, opt, ctr) * dirY);
        }

    }  // namespace

    SegResult form_segs(const VVD& vec_points, const VD& muon_start, const VD& muon_end,
                        double seg_len, const McsOptions& opt, McsCounters& ctr)
    {
        // upstream (mcs.cxx:365-424)
        SegResult out;

        Track track;
        int npoints = vec_points.size();
        for (int i = 0; i < npoints; i++) { track.add_point(vec_points[i], i); }

        VD com = { 0, 0, 0 };
        VD eigenvalues = { 0, 0, 0 };
        VVD axes = fit_pca(track, com, eigenvalues);
        VD aAxis = axes[0];
        // axes[1]/axes[2] are assigned and never used upstream
        // (mcs.cxx:381-382); retained only in the returned track_axes.
        VD vec_muon = { muon_end[0] - muon_start[0], muon_end[1] - muon_start[1],
                        muon_end[2] - muon_start[2] };
        if (dot(vec_muon, aAxis) < 0) { aAxis = scale(aAxis, -1.); }
        sort_points(track, aAxis);

        Track track_remainder = get_seg(track, aAxis, 1e6);

        VVD segPCAFit_prev;  // segPCAFit_vec[iseg-1], all that is ever read back
        std::vector<int> seg_npoints;
        VVD segs_aAxis_vec;
        VVD segCOM_vec;
        VD distance_vec, angle_vec, angleProjB_vec, angleProjC_vec;

        double currentFirstPointProj = dot(track.points.front(), aAxis);
        double currentLastPointProj = dot(track.points.front(), aAxis);
        double end_proj = dot(track.points.back(), aAxis);

        int nsegs_container = 0;
        int iseg = 0;
        while ((currentLastPointProj + 0.5 * seg_len) < end_proj && track_remainder.N >= 10) {
            VVD segFitVector = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
            VD segCOM = { 0, 0, 0 };
            Track seg;
            nsegs_container++;  // upstream pushes the Track BEFORE fitting
                                // (mcs.cxx:406-408), which is how bug #7's
                                // trailing empty segment arises

            bool can_fit = fit_seg_pca(track_remainder, aAxis, seg, segFitVector, segCOM,
                                       seg_len, currentFirstPointProj, currentLastPointProj,
                                       opt, ctr);
            if (!can_fit) { break; }
            seg_npoints.push_back(seg.N);
            segs_aAxis_vec.push_back(segFitVector[0]);
            segCOM_vec.push_back(segCOM);
            distance_vec.push_back((currentFirstPointProj + currentLastPointProj) / 2 -
                                   dot(muon_start, aAxis));

            if (iseg == 0) {
                angle_vec.push_back(-1);
                angleProjB_vec.push_back(-1);
                angleProjC_vec.push_back(-1);
            }
            else {
                set_seg_angles(segPCAFit_prev, segs_aAxis_vec.back(), angle_vec,
                               angleProjB_vec, angleProjC_vec, opt, ctr);
            }
            segPCAFit_prev = segFitVector;

            iseg++;
        }

        out.nsegs_container = nsegs_container;
        out.seg_npoints = std::move(seg_npoints);
        out.distance = std::move(distance_vec);
        out.angle = std::move(angle_vec);
        out.angle_projB = std::move(angleProjB_vec);
        out.angle_projC = std::move(angleProjC_vec);
        out.aAxes = std::move(segs_aAxis_vec);
        out.COM = std::move(segCOM_vec);
        out.track_axes = std::move(axes);
        return out;
    }

    // ------------------------------------------------------------------
    // stage 4: likelihood machinery (upstream mcs.cxx:426-602)
    namespace {

        // the MicroBooNE tune (upstream configure(), mcs.cxx:25-50).  The
        // standalone mcs.cxx is the source of truth -- NOT the stale
        // commented-out variant in WireCellMCS_module.cc:388-421 (doc 80
        // sec 6.4 #18).
        const double res_sigma1_xz = 0.005776;
        const double res_sigma2_xz = 0.01821;
        const VD par_sigma1_xz = { -0.449144931, 0.793132642, -1.291292240, 0.536765147, -0.084910516, 0.146304242 };
        const VD par_sigma2_xz = { 0.562850793, 0.118940108, 0.000625509, 0, -0.000000100, 1.217639251 };
        const VD par_ratio_xz = { 0.839684805, 0.839684805, 0, 0 };
        const VD res_sigma1_yz = { 0.0449, 0.0206, 0.01403, 0.0131, 0.0114 };
        const VD res_sigma2_yz = { 0.1506, 0.06154, 0.03965, 0.04179, 0.07347 };
        const VVD par_sigma1_yz = { { -0.09, -0.084325217, 0.487240052, 0.395496655, -0.187184874, 0.166128734 },
                                    { 0.0, -0.575280374, 0.070151974, 0.187260875, -0.099717108, 0.160128002 },
                                    { -0.153367057, -0.583042532, 0.983374136, -0.712652874, 0.134743902, 0.465439107 },
                                    { -0.268993212, 0.103899779, -0.588953942, 0.282356661, -0.067930741, 0.176282668 },
                                    { -0.2, -0.724028910, 0.660065851, -0.327141529, 0.038426745, 0.094357571 } };
        const VVD par_sigma2_yz = { { 0.193250363, 3.046073125, 15.0, -7.519451962, -1.0, 0.5 },
                                    { 2.226214634, -5.407245785, 7.227026554, -2.411882646, -0.000001517, 0.224728907 },
                                    { 0.25, 1.670060693, -2.086639043, 1.122640876, 0.000000440, 0.268340031 },
                                    { 0.430965414, -0.079204912, 0.220949488, -0.000010000, -0.000010000, 0.3 },
                                    { 1.635001641, -2.414781711, 1.669221724, -0.320484259, 0.0, 0.157018001 } };
        const VVD par_ratio_yz = { { 0.519230936, 1.0, 1.663982350, 0.084218766 },
                                   { 0.7, 0.788705752, 5.0, 2.045055158 },
                                   { 0.822433461, 0.701381849, 5.0, 0.0 },
                                   { 0.897408009, 0.7, 2.182533745, 0.297834161 },
                                   { 0.9, 0.9, 0.0, 0.0 } };

        // upstream (mcs.cxx:464)
        double sigma_highland(double T) { return 13.6 * (T + Mmu) / T / (T + 2 * Mmu); }

        // upstream (mcs.cxx:430)
        double sigmoid1(double x) { return 1. / (1. + std::exp(-x)); }

        // upstream (mcs.cxx:467)
        double sigmoid_par(double x, const VD& par)
        {
            return par[0] + (par[1] - par[0]) * (1 - 1. / (1 + std::exp(-par[2] * (x / 1000. - par[3]))));
        }

        // upstream (mcs.cxx:471-476)
        double quartic_decay(double x, const VD& par)
        {
            double u = x / par.back() / 1000.;
            double val = 0;
            for (int i = 0, n = par.size(); i < n - 1; i++) { val += par[i] * std::pow(u, i); }
            return 1 + val * std::exp(-u);
        }

        // upstream (mcs.cxx:495-502)
        double double_gaussian(double angle, const VD& pars)
        {
            double sigma1 = pars[0];
            double sigma2 = pars[1];
            double ratio = pars[2];
            double gaussian1 = (1. / (std::sqrt(2 * M_PI) * sigma1)) * std::exp(-0.5 * std::pow(angle / sigma1, 2));
            double gaussian2 = (1. / (std::sqrt(2 * M_PI) * sigma2)) * std::exp(-0.5 * std::pow(angle / sigma2, 2));
            return ratio * gaussian1 + (1 - ratio) * gaussian2;
        }

        // -log with the bug-#9 underflow floor: -log(0) = +inf poisons the
        // likelihood sum; the floor keeps the term huge (~708) but finite and
        // fires ONLY where upstream produced inf.
        double neg_log(double probability, const McsOptions& opt, McsCounters& ctr)
        {
            if (opt.fix_nan_guards && !(probability > 0)) {
                ctr.prob_floor++;
                probability = std::numeric_limits<double>::min();
            }
            return -std::log(probability);
        }

    }  // namespace

    VD pred_theta_xz_pars(double T)
    {
        // upstream (mcs.cxx:479-484)
        double sigma1 = std::sqrt(std::pow(sigma_highland(T) * quartic_decay(T, par_sigma1_xz), 2) + std::pow(res_sigma1_xz, 2));
        double sigma2 = std::sqrt(std::pow(sigma_highland(T) * quartic_decay(T, par_sigma2_xz), 2) + std::pow(res_sigma2_xz, 2));
        double area_ratio = sigmoid_par(T, par_ratio_xz);
        return { sigma1, sigma2, area_ratio };
    }

    VD pred_theta_yz_pars(double T, int vx_index)
    {
        // upstream (mcs.cxx:487-492)
        double sigma1 = std::sqrt(std::pow(sigma_highland(T) * quartic_decay(T, par_sigma1_yz[vx_index]), 2) + std::pow(res_sigma1_yz[vx_index], 2));
        double sigma2 = std::sqrt(std::pow(sigma_highland(T) * quartic_decay(T, par_sigma2_yz[vx_index]), 2) + std::pow(res_sigma2_yz[vx_index], 2));
        double area_ratio = sigmoid_par(T, par_ratio_yz[vx_index]);
        return { sigma1, sigma2, area_ratio };
    }

    double lnlikelihood_theta_xz(double angle, double T, const McsOptions& opt, McsCounters& ctr)
    {
        // upstream (mcs.cxx:505-508)
        VD pars = pred_theta_xz_pars(T);
        return neg_log(double_gaussian(angle, pars), opt, ctr);
    }

    double lnlikelihood_theta_yz(double angle, double T, double vx, const McsOptions& opt, McsCounters& ctr)
    {
        // upstream (mcs.cxx:512-538)
        VD vx_edges = { 0, 0.1, 0.2, 0.35, 0.75, 1 };
        VD emu_edges = { 600, 950, 1300 };

        double vx_abs = std::abs(vx);
        // upstream (mcs.cxx:520): the five indicator terms cover [0,1) only,
        // so |vx| >= 1 (a drift-aligned segment, or 1+eps from rounding) makes
        // every term false and falls through to ivx = 0 -- the
        // most-PERPENDICULAR bin, whose res_sigma1_yz is 4x the correct one
        // (bug #8).  The fix clamps to the [0.75,1] bin, which is what the
        // slicing means physically.
        int ivx = 0 * (vx_abs >= vx_edges[0] && vx_abs < vx_edges[1]) + 1 * (vx_abs >= vx_edges[1] && vx_abs < vx_edges[2]) +
                  2 * (vx_abs >= vx_edges[2] && vx_abs < vx_edges[3]) + 3 * (vx_abs >= vx_edges[3] && vx_abs < vx_edges[4]) +
                  4 * (vx_abs >= vx_edges[4] && vx_abs < vx_edges[5]);
        if (opt.fix_vx_overflow && vx_abs >= 1) {
            ctr.ivx_overflow++;
            ivx = 4;
        }

        VD pvx = { double_gaussian(angle, pred_theta_yz_pars(T, 0)), double_gaussian(angle, pred_theta_yz_pars(T, 1)),
                   double_gaussian(angle, pred_theta_yz_pars(T, 2)), double_gaussian(angle, pred_theta_yz_pars(T, 3)),
                   double_gaussian(angle, pred_theta_yz_pars(T, 4)) };

        // upstream (mcs.cxx:526-535): sigmoid blending between vx slices at
        // high energy
        double probability = 0.0;
        double width = 50;
        double scale1 = sigmoid1((T - emu_edges[0]) / width);
        double scale2 = sigmoid1((T - emu_edges[2]) / width);
        if (ivx == 4) {
            if (T < emu_edges[1]) { probability = (1 - scale1) * pvx[ivx] + scale1 * pvx[ivx - 1]; }
            else if (T >= emu_edges[1]) { probability = (1 - scale2) * pvx[ivx - 1] + scale2 * pvx[ivx - 2]; }
        }
        else if (ivx == 3) {
            probability = (1 - scale2) * pvx[ivx] + scale2 * pvx[ivx - 1];
        }
        else if (ivx <= 2) {
            probability = pvx[ivx];
        }

        return neg_log(probability, opt, ctr);
    }

    namespace {

        // upstream lnlikelihood_track (mcs.cxx:541-566).  par layout exactly as
        // estimate_energy packs it: [nsegs, d_1..d_n, ax_1..ax_n, ay_1..ay_n,
        // vx_1..vx_n].  The loop starts at i=2, which is what makes the i=1
        // "-1" angle sentinel safe (doc 80 sec 6.4 #16).
        double lnlikelihood_track(const double* KE, const double* par,
                                  const McsOptions& opt, McsCounters& ctr)
        {
            int nsegs = par[0];
            double lnlikelihood = 0;
            for (int i = 2; i < nsegs + 1; i++) {
                double theta_xz = par[i + nsegs];
                double theta_yz = par[i + 2 * nsegs];
                double vx = par[i + 3 * nsegs];

                double distance1 = par[i - 1];
                double distance2 = par[i];
                double rrtot_guess = rr_from_ke(KE[0]);
                double rrguess1 = std::max(rrtot_guess - distance1, 1.);  // >= 1 cm
                double rrguess2 = std::max(rrtot_guess - distance2, 1.);
                double keguess1 = ke_from_rr(rrguess1);
                double keguess2 = ke_from_rr(rrguess2);
                double keguess = (keguess1 + keguess2) / 2;  // angle matched to the avg energy

                double lnl_xz = lnlikelihood_theta_xz(theta_xz, keguess, opt, ctr);
                double lnl_yz = lnlikelihood_theta_yz(theta_yz, keguess, vx, opt, ctr);
                lnlikelihood += lnl_xz + lnl_yz;
            }
            return lnlikelihood;
        }

    }  // namespace

    EnergyResult estimate_energy(const VD& segs_distance_in, const VD& segs_angle_x,
                                 const VD& segs_angle_y, const VD& vx_comps,
                                 const McsOptions& opt, McsCounters& ctr)
    {
        // upstream (mcs.cxx:569-602).  The parameter vector is packed once and
        // reused (upstream reallocated per call -- hoisted, doc 80 sec 6.4).
        EnergyResult out;
        VD par = segs_distance_in;
        par.insert(par.begin(), (double)segs_distance_in.size());
        par.insert(par.end(), segs_angle_x.begin(), segs_angle_x.end());
        par.insert(par.end(), segs_angle_y.begin(), segs_angle_y.end());
        par.insert(par.end(), vx_comps.begin(), vx_comps.end());

        const double emin = 0;
        const double emax = 4e3;  // 4 GeV max estimate

        auto lnl = [&par, &opt, &ctr](double ke) {
            double KE = ke;
            return lnlikelihood_track(&KE, &par[0], opt, ctr);
        };

        // TF1::GetMinimumX swaps in the full function range when handed an
        // empty interval (xmin >= xmax) -- reachable here when keguess sits at
        // the very bottom of the scan so keguess*0.8 < emin + 1e-3.
        auto get_minimum_x = [&](double xlow, double xup, Scan& scan) {
            if (xlow >= xup) {
                xlow = emin;
                xup = emax;
            }
            Minimize1DResult mr = minimize1d(lnl, xlow, xup);
            scan = mr.first_scan;
            return mr.x;
        };

        double keguess = get_minimum_x(emin + 1e-3, emax - 1e-3, out.scans[0]);
        double keguess_lower = get_minimum_x(emin + 1e-3, keguess * 0.8, out.scans[1]);
        double keguess_higher = get_minimum_x(std::min(keguess * 1.2, emax - 2e-3), emax - 1e-3, out.scans[2]);

        double l_keguess = std::exp(-lnl(keguess));
        double l_keguess_lower = std::exp(-lnl(keguess_lower));
        double l_keguess_higher = std::exp(-lnl(keguess_higher));

        double ambiguity_score = std::max(l_keguess_lower / l_keguess, l_keguess_higher / l_keguess);

        out.keguess = keguess;
        out.keguess_lower = keguess_lower;
        out.keguess_higher = keguess_higher;
        out.l_keguess = l_keguess;
        out.l_keguess_lower = l_keguess_lower;
        out.l_keguess_higher = l_keguess_higher;
        out.ambiguity = ambiguity_score;
        return out;
    }

}  // namespace WireCell::Mcs::detail

// ---------------------------------------------------------------------------
McsResult MuonMCS::run(const std::vector<double>& vtx_start,
                       const std::vector<double>& vtx_end,
                       const std::vector<std::vector<double>>& points) const
{
    // upstream MCS::run (mcs.cxx:65-132), minus the std::cout chatter
    using namespace WireCell::Mcs::detail;
    const McsOptions& opt = m_options;
    McsResult res;

    TrimResult trim = trim_trajectory(points, vtx_start, vtx_end);
    const VVD& tp = trim.points_final;
    int npoints_trajectory_final = tp.size();
    res.bad_path = trim.bad_path;
    if (trim.bad_path) res.counters.bad_path++;

    // upstream (mcs.cxx:82-87): range and its energy over the trimmed path
    double rr_path = 0;
    for (int i = 1; i < npoints_trajectory_final; i++) { rr_path += norm(diff(tp[i], tp[i - 1])); }
    double KE_rr_path = ke_from_rr(rr_path);
    if (opt.fix_bad_path_outputs && trim.bad_path) {
        // bug #12: on bad_path the "path" is just the two vertices, so
        // rr_path is the start-end chord -- not a track length.  Upstream
        // publishes it anyway; keep -1 instead.
        res.mu_tracklen = -1;
        res.emu_tracklen = -1;
    }
    else {
        res.mu_tracklen = rr_path;
        res.emu_tracklen = (KE_rr_path + Mmu) / 1000.;  // KE -> total E, MeV -> GeV
    }

    // upstream (mcs.cxx:89-95): unusable paths keep emu_MCS = -1
    if (trim.bad_path || npoints_trajectory_final < 20 ||
        norm(diff(tp.back(), vtx_end)) < 2 * seg_length) {
        return res;
    }

    SegResult segs = form_segs(tp, vtx_start, vtx_end, seg_length, opt, res.counters);
    res.nsegs = segs.distance.size();

    if (opt.fix_single_segment) {
        // bug #7: upstream guards on segs.size() < 2, but the container can
        // hold a trailing EMPTY Track, so one real segment (zero angles)
        // passes, the likelihood loop never executes, and Brent minimises a
        // constant -- publishing a grid artefact as emu_MCS with
        // ambiguity = 1.  Require >= 2 FITTED segments (>= 1 real angle).
        if (res.nsegs < 2) {
            if (segs.nsegs_container >= 2) res.counters.single_seg_abort++;
            return res;
        }
    }
    else {
        if (segs.nsegs_container < 2) { return res; }
    }

    // upstream (mcs.cxx:118-121)
    VD vx_components;
    for (const auto& seg_dir : segs.aAxes) {
        if (!seg_dir.empty()) { vx_components.push_back(seg_dir[0]); }
    }

    EnergyResult er = estimate_energy(segs.distance, segs.angle_projB, segs.angle_projC,
                                      vx_components, opt, res.counters);
    res.emu_MCS = (er.keguess + Mmu) / 1000.;  // KE -> total E, MeV -> GeV
    res.ambiguity_MCS = er.ambiguity;
    return res;
}
