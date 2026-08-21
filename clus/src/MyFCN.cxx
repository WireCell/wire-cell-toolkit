#include "WireCellClus/MyFCN.h"
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Logging.h"

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <cmath>

static auto s_log = WireCell::Log::logger("clus.NeutrinoPattern");

using namespace WireCell;
using namespace WireCell::Clus::PR;

// doc sbnd_xin/docs/pr/51 round 7 (robust vertex fit): pure geometry helpers
// for the per-leg dynamic direction-window substitution.  Declared in
// PRSegmentFunctions.h so they are doctestable without a cluster; implemented
// here next to the production math they replicate.
namespace WireCell::Clus::PR {

MvfitAnnulusPca mvfit_annulus_pca(const std::vector<WireCell::Point>& pts,
                                  const WireCell::Point& vtx,
                                  double rin, double rout)
{
    MvfitAnnulusPca out;
    double min_dis = 1e9;
    std::vector<size_t> kept;   // index order -- deterministic
    for (size_t i = 0; i != pts.size(); i++) {
        const double d = (pts[i] - vtx).magnitude();
        if (d < rin || d > rout) continue;
        kept.push_back(i);
        if (d < min_dis) {
            out.center = pts[i];
            min_dis = d;
        }
    }
    out.npts = (int)kept.size();
    if (kept.size() < 2) return out;

    // covariance about the near-vertex kept point, exactly as
    // MyFCN::AddSegment (NOT the window mean -- see docs/pr/28 sec. 3.1).
    Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
    for (size_t k : kept) {
        const double dx = pts[k].x() - out.center.x();
        const double dy = pts[k].y() - out.center.y();
        const double dz = pts[k].z() - out.center.z();
        cov(0, 0) += dx * dx;
        cov(0, 1) += dx * dy;
        cov(0, 2) += dx * dz;
        cov(1, 1) += dy * dy;
        cov(1, 2) += dy * dz;
        cov(2, 2) += dz * dz;
    }
    cov(1, 0) = cov(0, 1);
    cov(2, 0) = cov(0, 2);
    cov(2, 1) = cov(1, 2);
    cov /= (double)kept.size();

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(cov);
    const Eigen::Vector3d ev = es.eigenvalues();     // ascending
    const Eigen::Matrix3d evec = es.eigenvectors();
    for (int i = 0; i != 3; i++) {
        const int col = 2 - i;                       // descending
        out.vals[i] = ev(col) + std::pow(0.15 * units::cm, 2);
        const double norm = std::sqrt(evec(0, col) * evec(0, col) +
                                      evec(1, col) * evec(1, col) +
                                      evec(2, col) * evec(2, col));
        out.axes[i] = WireCell::Point(evec(0, col) / norm,
                                      evec(1, col) / norm,
                                      evec(2, col) / norm);
    }
    return out;
}

double mvfit_rout_dyn(double chord, double frac, double rmin, double rmax)
{
    return std::min(std::max(frac * chord, rmin), rmax);
}

bool mvfit_leg_disagrees(const WireCell::Point& inner_axis,
                         const MvfitAnnulusPca& outer,
                         double angle_deg, int min_pts, double min_aniso)
{
    if (outer.npts < min_pts) return false;
    if (outer.vals[1] <= 0) return false;
    if (std::sqrt(outer.vals[0] / outer.vals[1]) < min_aniso) return false;
    const double mi = inner_axis.magnitude();
    const double mo = outer.axes[0].magnitude();
    if (mi <= 0 || mo <= 0) return false;
    // folded angle: eigenvector signs are arbitrary
    double c = std::abs(inner_axis.dot(outer.axes[0])) / (mi * mo);
    if (c > 1.0) c = 1.0;
    const double angle = std::acos(c) * 180.0 / M_PI;
    return angle > angle_deg;
}

// doc sbnd_xin/docs/pr/104: junction-snap helpers (declared in
// PRSegmentFunctions.h, doctest_vertex_junction_snap.cxx).
int vjs_direction_classes(const std::vector<WireCell::Vector>& dirs, double collinear_deg)
{
    int n = 0;
    std::vector<bool> used(dirs.size(), false);
    for (size_t i = 0; i < dirs.size(); i++) {
        if (used[i]) continue;
        used[i] = true;
        n++;
        for (size_t j = i + 1; j < dirs.size(); j++) {
            if (used[j]) continue;
            const double mi = dirs[i].magnitude(), mj = dirs[j].magnitude();
            if (mi <= 0 || mj <= 0) continue;
            double c = dirs[i].dot(dirs[j]) / (mi * mj);
            if (c > 1.0) c = 1.0;
            if (c < -1.0) c = -1.0;
            if (std::acos(c) * 180.0 / M_PI > collinear_deg) used[j] = true;
        }
    }
    return n;
}

bool vjs_joint_fit(const std::vector<std::pair<WireCell::Point, WireCell::Vector>>& lines,
                   WireCell::Point& out, double& rms)
{
    Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
    Eigen::Vector3d b = Eigen::Vector3d::Zero();
    std::vector<std::pair<Eigen::Vector3d, Eigen::Matrix3d>> terms;
    for (const auto& [p, d] : lines) {
        const double m = d.magnitude();
        if (m <= 0) continue;
        Eigen::Vector3d u(d.x() / m, d.y() / m, d.z() / m);
        Eigen::Vector3d q(p.x(), p.y(), p.z());
        Eigen::Matrix3d Pm = Eigen::Matrix3d::Identity() - u * u.transpose();
        A += Pm;
        b += Pm * q;
        terms.emplace_back(q, Pm);
    }
    if (terms.size() < 2) return false;
    Eigen::FullPivLU<Eigen::Matrix3d> lu(A);
    if (!lu.isInvertible()) return false;
    Eigen::Vector3d x = lu.solve(b);
    if (!std::isfinite(x(0)) || !std::isfinite(x(1)) || !std::isfinite(x(2))) return false;
    double ss = 0;
    for (const auto& [q, Pm] : terms) {
        Eigen::Vector3d r = Pm * (x - q);
        ss += r.squaredNorm();
    }
    rms = std::sqrt(ss / terms.size());
    out = WireCell::Point(x(0), x(1), x(2));
    return true;
}

}  // namespace WireCell::Clus::PR

MyFCN::MyFCN(VertexPtr vtx, bool flag_vtx_constraint, double vtx_constraint_range, 
             double vertex_protect_dis, double vertex_protect_dis_short_track, double fit_dis)
    : vtx(vtx)
    , enforce_two_track_fit(false)
    , flag_vtx_constraint(flag_vtx_constraint)
    , vtx_constraint_range(vtx_constraint_range)
    , vertex_protect_dis(vertex_protect_dis)
    , vertex_protect_dis_short_track(vertex_protect_dis_short_track)
    , fit_dis(fit_dis)
{
    segments.clear();
    vec_points.clear();
}

MyFCN::~MyFCN()
{
}

void MyFCN::print_points()
{
    for (size_t i = 0; i != vec_points.size(); i++) {
        for (size_t j = 0; j != vec_points.at(i).size(); j++) {
            SPDLOG_LOGGER_TRACE(s_log, "print_points: {} {} {} {} {}",
                i, j,
                vec_points.at(i).at(j).x() / units::cm,
                vec_points.at(i).at(j).y() / units::cm,
                vec_points.at(i).at(j).z() / units::cm);
        }
    }
}

void MyFCN::AddSegment(SegmentPtr sg)
{
    // push in ...
    segments.push_back(sg);
    {
        std::vector<Facade::geo_point_t> pts;
        vec_points.push_back(pts);
    }

    Facade::geo_point_t center(0, 0, 0);
    double min_dis = 1e9;

    // prototype (NeutrinoID_improve_vertex.h:538): WCP::PointVector& pts =
    // sg->get_point_vec(), and ProtoSegment.h:16 defines get_point_vec() as
    // {return fit_pt_vec;} -- the *fitted* trajectory written by the
    // do_multi_tracking that always precedes fit_vertex, NOT the raw Steiner
    // path (that is get_wcpt_vec()).  The toolkit counterpart of fit_pt_vec is
    // Segment::fits().  Reading wcpts() here fed the vertex fit the unfitted
    // skeleton; see sbnd_xin/docs/pr/28 sec. 3.1.
    const auto& fits = sg->fits();
    if (fits.empty()) {
        Facade::geo_point_t a(0, 0, 0);
        vec_PCA_dirs.push_back(std::make_tuple(a, a, a));
        vec_PCA_vals.push_back(std::make_tuple(0, 0, 0));
        vec_centers.push_back(a);
        return;
    }

    std::vector<Facade::geo_point_t> pts;
    pts.reserve(fits.size());
    for (const auto& f : fits) {
        pts.push_back(f.point);
    }
    double length = 0;
    if (pts.size() > 1) {
        auto front = pts.front();
        auto back = pts.back();
        length = std::sqrt(std::pow(front.x() - back.x(), 2) + 
                          std::pow(front.y() - back.y(), 2) + 
                          std::pow(front.z() - back.z(), 2));
    }

    Facade::geo_point_t vtx_pt = vtx->fit().point;
    
    for (size_t i = 0; i != pts.size(); i++) {
        double dis_to_vertex = std::sqrt(std::pow(pts.at(i).x() - vtx_pt.x(), 2) + 
                                        std::pow(pts.at(i).y() - vtx_pt.y(), 2) + 
                                        std::pow(pts.at(i).z() - vtx_pt.z(), 2));
        if (length > 3.0 * units::cm) {
            if (dis_to_vertex < vertex_protect_dis || dis_to_vertex > fit_dis) continue;
        } else {
            if (dis_to_vertex < vertex_protect_dis_short_track || dis_to_vertex > fit_dis) continue;
        }
        
        vec_points.back().push_back(pts.at(i));
        if (dis_to_vertex < min_dis) {
            center = pts.at(i);
            min_dis = dis_to_vertex;
        }
    }

    // calculate the PCA ...
    if (vec_points.back().size() > 1) {
        int nsum = vec_points.back().size();
        
        // Eigen vectors ...
        std::vector<Facade::geo_point_t> PCA_axis(3, Facade::geo_point_t(0, 0, 0));
        double PCA_values[3] = {0, 0, 0};
        
        Eigen::Matrix3d cov_matrix = Eigen::Matrix3d::Zero();
        
        for (size_t k = 0; k != vec_points.back().size(); k++) {
            double dx = vec_points.back().at(k).x() - center.x();
            double dy = vec_points.back().at(k).y() - center.y();
            double dz = vec_points.back().at(k).z() - center.z();
            
            cov_matrix(0, 0) += dx * dx;
            cov_matrix(0, 1) += dx * dy;
            cov_matrix(0, 2) += dx * dz;
            cov_matrix(1, 1) += dy * dy;
            cov_matrix(1, 2) += dy * dz;
            cov_matrix(2, 2) += dz * dz;
        }
        
        cov_matrix(1, 0) = cov_matrix(0, 1);
        cov_matrix(2, 0) = cov_matrix(0, 2);
        cov_matrix(2, 1) = cov_matrix(1, 2);
        
        cov_matrix /= nsum;
        
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(cov_matrix);
        Eigen::Vector3d eigen_values = eigen_solver.eigenvalues();
        Eigen::Matrix3d eigen_vectors = eigen_solver.eigenvectors();
        
        // Eigen returns eigenvalues in ascending order; reverse to descending
        // (largest first) to match prototype convention (ROOT's TMatrixDEigen).
        // FitVertex zeros out row 0 (track direction = largest eigenvalue) and
        // weights rows 1,2 by sqrt(λ0/λk), so index 0 must be the largest.
        PCA_values[0] = eigen_values(2) + std::pow(0.15 * units::cm, 2);
        PCA_values[1] = eigen_values(1) + std::pow(0.15 * units::cm, 2);
        PCA_values[2] = eigen_values(0) + std::pow(0.15 * units::cm, 2);

        for (int i = 0; i != 3; i++) {
            int col = 2 - i;  // reverse column order to match descending eigenvalues
            double norm = std::sqrt(eigen_vectors(0, col) * eigen_vectors(0, col) +
                                   eigen_vectors(1, col) * eigen_vectors(1, col) +
                                   eigen_vectors(2, col) * eigen_vectors(2, col));
            PCA_axis[i] = Facade::geo_point_t(eigen_vectors(0, col) / norm,
                                     eigen_vectors(1, col) / norm,
                                     eigen_vectors(2, col) / norm);
        }
        
        vec_PCA_dirs.push_back(std::make_tuple(PCA_axis[0], PCA_axis[1], PCA_axis[2]));
        vec_PCA_vals.push_back(std::make_tuple(PCA_values[0], PCA_values[1], PCA_values[2]));
        vec_centers.push_back(center);
        
    } else {
        Facade::geo_point_t a(0, 0, 0);
        vec_PCA_dirs.push_back(std::make_tuple(a, a, a));
        vec_PCA_vals.push_back(std::make_tuple(0, 0, 0));
        if (vec_points.back().size() == 1) {
            vec_centers.push_back(vec_points.back().back());
        } else {
            vec_centers.push_back(a);
        }
    }

    // doc sbnd_xin/docs/pr/51 round 7: robust per-leg direction substitution.
    // Epilogue only -- every statement above is untouched production code, and
    // with the knob off (m_robust.on == false, the C++ default) it never runs.
    if (!m_robust.on) return;
    robust_maybe_substitute_leg(sg);
}

// doc sbnd_xin/docs/pr/51 round 7 (follow-ups 1-5 of the same doc): the
// production (1.5, 6] cm window reads back, for a long leg, mostly the
// straight 4 cm stretch the previous UpdateInfo wrote toward the current
// vertex -- the fit self-confirms wherever the vertex is.  This epilogue
// re-estimates the leg's direction over a re-seat-free dynamic annulus and
// substitutes it when it disagrees with the fit-visible one.  Design and
// gate rationale at MyFCN.h RobustParams.
void MyFCN::robust_maybe_substitute_leg(SegmentPtr sg)
{
    // only legs the production fit counts (valid inner PCA)
    if (vec_points.back().size() <= 1) return;

    // shower veto: the same three-clause recipe improve_vertex uses for its
    // ntracks census (NeutrinoVertexFinder.cxx flag_skip_two_legs block).  An
    // EM trunk's direction is not improved by a longer lever (it fans out);
    // the anisotropy gate below is the intrinsic guard when flags are unset.
    if (sg->flags_any(SegmentFlags::kShowerTrajectory) ||
        sg->flags_any(SegmentFlags::kShowerTopology) ||
        (sg->particle_info() && std::abs(sg->particle_info()->pdg()) == 11)) {
        return;
    }

    const auto& fits = sg->fits();
    if (fits.size() < 2) return;
    std::vector<Facade::geo_point_t> pts;
    pts.reserve(fits.size());
    for (const auto& f : fits) {
        pts.push_back(f.point);
    }
    const double chord = (pts.front() - pts.back()).magnitude();
    if (chord < m_robust.min_len) return;

    const Facade::geo_point_t vtx_pt = vtx->fit().point;

    // re-seat radius mirror of UpdateInfo (wcpt space, default_dis_cut =
    // 4 cm as at both UpdateInfo call sites): only a leg whose far wcpt end
    // is > 8 cm from the vertex carries a 4 cm straight re-seat to clear.
    double reseat = 0;
    const auto& seg_wcpts = sg->wcpts();
    if (!seg_wcpts.empty()) {
        const double default_dis_cut = 4.0 * units::cm;
        const double dis_front = (seg_wcpts.front().point - vtx_pt).magnitude();
        const double dis_back = (seg_wcpts.back().point - vtx_pt).magnitude();
        if (std::max(dis_front, dis_back) > 2 * default_dis_cut) {
            reseat = default_dis_cut;
        }
    }
    const double rin = std::max(vertex_protect_dis, reseat + m_robust.rin_margin);
    const double rout = mvfit_rout_dyn(chord, m_robust.rout_frac,
                                       m_robust.rout_min, m_robust.rout_max);

    const auto outer = mvfit_annulus_pca(pts, vtx_pt, rin, rout);
    const auto& inner_axis = std::get<0>(vec_PCA_dirs.back());
    if (!mvfit_leg_disagrees(inner_axis, outer, m_robust.angle,
                             m_robust.min_pts, m_robust.min_aniso)) {
        return;
    }

    // folded angle recomputed for the sentinel only
    double cang = std::abs(inner_axis.dot(outer.axes[0]));
    if (cang > 1.0) cang = 1.0;
    const double angle = std::acos(cang) * 180.0 / M_PI;

    // save the production tuple (exact restore for the charge-veto path),
    // then substitute PCA dirs/vals + center.  vec_points stays production:
    // ntracks, the fit gate's ntracks side and the npoints prior weight are
    // unchanged by construction.
    const size_t idx = vec_PCA_dirs.size() - 1;
    m_robust_backup.emplace_back(idx, vec_PCA_dirs.back(), vec_PCA_vals.back(),
                                 vec_centers.back());

    std::array<Facade::geo_point_t, 3> ax = outer.axes;
    const Facade::geo_point_t prod_ax[3] = {std::get<0>(vec_PCA_dirs.back()),
                                            std::get<1>(vec_PCA_dirs.back()),
                                            std::get<2>(vec_PCA_dirs.back())};
    for (int k = 0; k != 3; k++) {
        // hemisphere-orient toward the production counterpart: the FitVertex
        // pair-angle census is sign-sensitive, the 3x3 data term is not
        if (ax[k].dot(prod_ax[k]) < 0) {
            ax[k] = ax[k] * (-1.0);
        }
    }
    vec_PCA_dirs.back() = std::make_tuple(ax[0], ax[1], ax[2]);
    vec_PCA_vals.back() = std::make_tuple(outer.vals[0], outer.vals[1], outer.vals[2]);
    vec_centers.back() = outer.center;
    m_n_substituted++;

    SPDLOG_LOGGER_DEBUG(s_log,
        "mvfit: cluster {} leg {} substituted: chord={:.1f} cm window=({:.1f},{:.1f}] "
        "npts={} aniso={:.1f} angle={:.1f} deg",
        vtx->cluster() ? std::to_string(vtx->cluster()->get_cluster_id()) : std::string("unknown"),
        idx, chord / units::cm, rin / units::cm, rout / units::cm,
        outer.npts, std::sqrt(outer.vals[0] / outer.vals[1]), angle);
}

void MyFCN::restore_production_pca()
{
    if (!m_robust_backup.empty()) {
        SPDLOG_LOGGER_DEBUG(s_log,
            "mvfit: cluster {} charge-veto restore: {} leg(s) reverted to production PCA",
            vtx->cluster() ? std::to_string(vtx->cluster()->get_cluster_id()) : std::string("unknown"),
            m_robust_backup.size());
    }
    for (const auto& bk : m_robust_backup) {
        const size_t idx = std::get<0>(bk);
        vec_PCA_dirs.at(idx) = std::get<1>(bk);
        vec_PCA_vals.at(idx) = std::get<2>(bk);
        vec_centers.at(idx) = std::get<3>(bk);
    }
    m_robust_backup.clear();
    m_n_substituted = 0;
}

void MyFCN::update_fit_range(double tmp_vertex_protect_dis, double tmp_vertex_protect_dis_short_track, double tmp_fit_dis)
{
    vertex_protect_dis = tmp_vertex_protect_dis;
    vertex_protect_dis_short_track = tmp_vertex_protect_dis_short_track;
    fit_dis = tmp_fit_dis;

    std::vector<SegmentPtr> tmp_segments = segments;
    segments.clear();
    vec_points.clear();
    for (auto it = tmp_segments.begin(); it != tmp_segments.end(); it++) {
        AddSegment(*it);
    }
}

int MyFCN::get_fittable_tracks()
{
    int ncount = 0;
    for (size_t i = 0; i != vec_points.size(); i++) {
        if (vec_points.at(i).size() > 1) ncount++;
    }
    return ncount;
}

std::pair<SegmentPtr, int> MyFCN::get_seg_info(int i)
{
    if (i < (int)segments.size()) {
        return std::make_pair(segments.at(i), vec_points.at(i).size());
    }
    return std::make_pair(nullptr, 0);
}

std::pair<bool, WireCell::Clus::Facade::geo_point_t> MyFCN::FitVertex()
{
    Facade::geo_point_t fit_pos = vtx->fit().point;
    bool fit_flag = false;

    int ntracks = get_fittable_tracks();
    int npoints = 0;

    int n_large_angles = 0;
    for (size_t i = 0; i != vec_PCA_vals.size(); i++) {
        Eigen::Vector3d dir1(std::get<0>(vec_PCA_dirs.at(i)).x(), 
                            std::get<0>(vec_PCA_dirs.at(i)).y(), 
                            std::get<0>(vec_PCA_dirs.at(i)).z());
        for (size_t j = i + 1; j < vec_PCA_vals.size(); j++) {
            Eigen::Vector3d dir2(std::get<0>(vec_PCA_dirs.at(j)).x(), 
                                std::get<0>(vec_PCA_dirs.at(j)).y(), 
                                std::get<0>(vec_PCA_dirs.at(j)).z());
            // prototype (NeutrinoID_improve_vertex.h:689): dir1.Angle(dir2).
            // ROOT's TVector3::Angle returns 0 when either vector is null and
            // clamps the cosine to [-1,1].  A bare std::acos(dot) instead gives
            // 90 deg for a null direction -- which is what AddSegment pushes for
            // a segment with <2 points inside the annulus -- so a stub counted
            // as a large angle and could open a gate the prototype keeps shut;
            // and it gives NaN when FP rounding pushes the dot past -1, losing a
            // genuine anti-parallel pair.  See sbnd_xin/docs/pr/28 sec. 3.2.
            const double ptot2 = dir1.squaredNorm() * dir2.squaredNorm();
            double angle = 0.0;
            if (ptot2 > 0) {
                double arg = dir1.dot(dir2) / std::sqrt(ptot2);
                if (arg > 1.0) arg = 1.0;
                if (arg < -1.0) arg = -1.0;
                angle = std::acos(arg) * 180.0 / M_PI;
            }
            if (angle > 15) n_large_angles++;
        }
    }

    if ((ntracks > 2 && n_large_angles > 1) || (ntracks >= 2 && enforce_two_track_fit && n_large_angles >= 1)) {

        // doc sbnd_xin/docs/pr/51 round 7: a distorted 2-leg vertex needs
        // more corrective authority than the 0.43 cm polish prior allows --
        // relax it iff a leg was actually substituted.  >= 3-leg vertices
        // keep the production prior (already well braced).  Mutating the
        // member leaves the constraint block below verbatim; MyFCN is
        // per-call disposable so the mutation has no afterlife.
        if (m_robust.on && m_n_substituted > 0 && ntracks == 2) {
            SPDLOG_LOGGER_DEBUG(s_log,
                "mvfit: cluster {} 2-leg vertex with {} substituted leg(s): prior {:.2f} -> {:.2f} cm",
                vtx->cluster() ? std::to_string(vtx->cluster()->get_cluster_id()) : std::string("unknown"),
                m_n_substituted, vtx_constraint_range / units::cm,
                m_robust.prior_range / units::cm);
            vtx_constraint_range = m_robust.prior_range;
        }

        // start the fit ...
        Eigen::VectorXd temp_pos_3D_init(3), temp_pos_3D(3); // to be fitted
        temp_pos_3D_init(0) = fit_pos.x();
        temp_pos_3D_init(1) = fit_pos.y();
        temp_pos_3D_init(2) = fit_pos.z();

        Eigen::Vector3d b(0, 0, 0);
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(3, 3);

        for (size_t i = 0; i != vec_PCA_vals.size(); i++) {
            if (std::get<0>(vec_PCA_vals.at(i)) > 0) {
                npoints += vec_points.at(i).size();

                // fill the matrix ... first row, second column
                Eigen::MatrixXd R(3, 3);
                R(0, 0) = 0; R(0, 1) = 0; R(0, 2) = 0;
                double val1 = std::sqrt(std::get<0>(vec_PCA_vals.at(i)) / std::get<1>(vec_PCA_vals.at(i)));
                R(1, 0) = val1 * std::get<1>(vec_PCA_dirs.at(i)).x();
                R(1, 1) = val1 * std::get<1>(vec_PCA_dirs.at(i)).y();
                R(1, 2) = val1 * std::get<1>(vec_PCA_dirs.at(i)).z();
                val1 = std::sqrt(std::get<0>(vec_PCA_vals.at(i)) / std::get<2>(vec_PCA_vals.at(i)));
                R(2, 0) = val1 * std::get<2>(vec_PCA_dirs.at(i)).x();
                R(2, 1) = val1 * std::get<2>(vec_PCA_dirs.at(i)).y();
                R(2, 2) = val1 * std::get<2>(vec_PCA_dirs.at(i)).z();

                Eigen::Vector3d data(vec_centers.at(i).x(), vec_centers.at(i).y(), vec_centers.at(i).z());
                data = R * data;
                Eigen::MatrixXd RT = R.transpose();

                b += RT * data;
                A += RT * R;
            }
        }

        // add constraint ...
        if (flag_vtx_constraint) {
            Eigen::MatrixXd R = Eigen::MatrixXd::Zero(3, 3);
            R(0, 0) = 1.0 / vtx_constraint_range * std::sqrt(npoints);
            R(1, 1) = 1.0 / vtx_constraint_range * std::sqrt(npoints);
            R(2, 2) = 1.0 / vtx_constraint_range * std::sqrt(npoints);

            Eigen::Vector3d data(fit_pos.x() / vtx_constraint_range * std::sqrt(npoints), 
                                fit_pos.y() / vtx_constraint_range * std::sqrt(npoints), 
                                fit_pos.z() / vtx_constraint_range * std::sqrt(npoints));
            Eigen::MatrixXd RT = R.transpose();
            b += RT * data;
            A += RT * R;
        }

        Eigen::BiCGSTAB<Eigen::MatrixXd> solver;
        solver.compute(A);
        temp_pos_3D = solver.solveWithGuess(b, temp_pos_3D_init);

        if (!std::isnan(solver.error())) {
            fit_pos = Facade::geo_point_t(temp_pos_3D(0), temp_pos_3D(1), temp_pos_3D(2));
            fit_flag = true;
        } else {
            SPDLOG_LOGGER_TRACE(s_log, "FitVertex: Cluster: {} Fit Vertex Failed!",
                vtx->cluster() ? std::to_string(vtx->cluster()->get_cluster_id()) : std::string("unknown"));
        }
    }

    return std::make_pair(fit_flag, fit_pos);
}

bool MyFCN::UpdateInfo(Facade::geo_point_t fit_pos, Facade::Cluster& temp_cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, double default_dis_cut){
    
    // Get PC transform and cluster parameters
    auto pcts = track_fitter.get_pc_transforms();
    const auto transform = pcts->pc_transform(temp_cluster.get_scope_transform(temp_cluster.get_default_scope()));
    double cluster_t0 = temp_cluster.get_cluster_t0();
    
    // Get APA/face for the fit position
    auto test_wpid = dv->contained_by(fit_pos);
    if (test_wpid.apa() == -1 || test_wpid.face() == -1) {
        SPDLOG_LOGGER_TRACE(s_log, "UpdateInfo: Warning: fit_pos not contained in detector volume");
        return false;
    }
    int apa = test_wpid.apa();
    int face = test_wpid.face();
    
    // Get geometry parameters from TrackFitting
    const auto& wpid_offsets = track_fitter.get_wpid_offsets();
    const auto& wpid_slopes = track_fitter.get_wpid_slopes();
    
    WirePlaneId wpid(kAllLayers, face, apa);
    auto offset_it = wpid_offsets.find(wpid);
    auto slope_it = wpid_slopes.find(wpid);
    
    if (offset_it == wpid_offsets.end() || slope_it == wpid_slopes.end()) {
        SPDLOG_LOGGER_TRACE(s_log, "UpdateInfo: Warning: geometry parameters not found for APA {} Face {}", apa, face);
        return false;
    }
    
    auto offset_t = std::get<0>(offset_it->second);
    auto offset_u = std::get<1>(offset_it->second);
    auto offset_v = std::get<2>(offset_it->second);
    auto offset_w = std::get<3>(offset_it->second);
    auto slope_x = std::get<0>(slope_it->second);
    auto slope_yu = std::get<1>(slope_it->second).first;
    auto slope_zu = std::get<1>(slope_it->second).second;
    auto slope_yv = std::get<2>(slope_it->second).first;
    auto slope_zv = std::get<2>(slope_it->second).second;
    auto slope_yw = std::get<3>(slope_it->second).first;
    auto slope_zw = std::get<3>(slope_it->second).second;
    
    // Transform to raw coordinates for wire plane calculations
    auto fit_pos_raw = transform->backward(fit_pos, cluster_t0, face, apa);
    auto vtx_pos_raw = transform->backward(vtx->fit().point, cluster_t0, face, apa);
    
    // Print update information if vertex moved significantly
    if ((fit_pos - vtx->fit().point).magnitude() > 0.01 * units::cm) {
        SPDLOG_LOGGER_TRACE(s_log, "UpdateInfo: Cluster: {} Update Vertex: ({}, {}, {}, {}) <- ({}, {}, {}, {})",
            vtx->cluster() ? std::to_string(vtx->cluster()->get_cluster_id()) : std::string("unknown"),
            offset_u + (slope_yu * fit_pos_raw.y() + slope_zu * fit_pos_raw.z()),
            offset_v + (slope_yv * fit_pos_raw.y() + slope_zv * fit_pos_raw.z()),
            offset_w + (slope_yw * fit_pos_raw.y() + slope_zw * fit_pos_raw.z()),
            offset_t + slope_x * fit_pos_raw.x(),
            offset_u + (slope_yu * vtx_pos_raw.y() + slope_zu * vtx_pos_raw.z()),
            offset_v + (slope_yv * vtx_pos_raw.y() + slope_zv * vtx_pos_raw.z()),
            offset_w + (slope_yw * vtx_pos_raw.y() + slope_zw * vtx_pos_raw.z()),
            offset_t + slope_x * vtx_pos_raw.x());
    }
    
    // Get steiner point cloud from cluster.  size_major() term: Dataset::size()
    // counts ARRAYS, not points (Clustering_Util.cxx:81) -- an empty steiner
    // cloud with zero-length coordinate arrays passes the .size() test, the
    // kNN calls below then return empty, and the unguarded [0]s read past the
    // end.  Same predicate fix already applied at Clustering_Util.cxx:87 and
    // NeutrinoPatternBase.cxx:212.
    if (!temp_cluster.has_pc("steiner_pc") || temp_cluster.get_pc("steiner_pc").size() == 0
        || temp_cluster.get_pc("steiner_pc").size_major() == 0) {
        SPDLOG_LOGGER_TRACE(s_log, "UpdateInfo: Warning: steiner_pc not found in cluster");
        return false;
    }
    const auto& steiner_pc = temp_cluster.get_pc("steiner_pc");
    const auto& coords = temp_cluster.get_default_scope().coords;
    const auto& x_coords = steiner_pc.get(coords.at(0))->elements<double>();
    const auto& y_coords = steiner_pc.get(coords.at(1))->elements<double>();
    const auto& z_coords = steiner_pc.get(coords.at(2))->elements<double>();
    
    // Find closest steiner point to new fit position
    auto vtx_knn_results = temp_cluster.kd_steiner_knn(1, fit_pos, "steiner_pc");
    size_t vtx_new_idx = vtx_knn_results[0].first;
    Facade::geo_point_t vtx_new_pt(x_coords[vtx_new_idx], y_coords[vtx_new_idx], z_coords[vtx_new_idx]);
    
    // Update each segment connected to this vertex
    for (size_t i = 0; i != segments.size(); i++) {
        auto& seg_wcpts = segments.at(i)->wcpts();
        if (seg_wcpts.empty()) continue;
        
        // Determine if vertex is at front or back of segment
        bool flag_front = false;
        double dis_front = (seg_wcpts.front().point - vtx->fit().point).magnitude();
        double dis_back = (seg_wcpts.back().point - vtx->fit().point).magnitude();
        flag_front = (dis_front < dis_back);
        
        // Find closest wcpt to the PCA center, excluding points too close to vertex
        double max_dis = std::max(dis_front, dis_back);
        double dis_cut = (max_dis > 2 * default_dis_cut) ? default_dis_cut : 0;
        
        size_t min_idx = 0;
        double min_dis = 1e9;
        for (size_t j = 0; j != seg_wcpts.size(); j++) {
            double dis_to_center = (seg_wcpts.at(j).point - vec_centers.at(i)).magnitude();
            double dis_to_vtx = (seg_wcpts.at(j).point - vtx->fit().point).magnitude();
            if (dis_to_center < min_dis && dis_to_vtx > dis_cut) {
                min_idx = j;
                min_dis = dis_to_center;
            }
        }
        
        auto& min_wcp = seg_wcpts.at(min_idx);
        
        // Create new path from new vertex position to closest point using steiner graph
        std::list<WCPoint> new_list;
        new_list.push_back(WCPoint{vtx_new_pt});
        
        // Interpolate intermediate points
        double dis_step = 2.0 * units::cm;
        double total_dis = (vtx_new_pt - min_wcp.point).magnitude();
        int ncount = std::round(total_dis / dis_step);
        if (ncount < 2) ncount = 2;
        
        // double cumulative_ray_length = 0.0;
        for (int qx = 1; qx < ncount; qx++) {
            Facade::geo_point_t tmp_p(
                vtx_new_pt.x() + (min_wcp.point.x() - vtx_new_pt.x()) / ncount * qx,
                vtx_new_pt.y() + (min_wcp.point.y() - vtx_new_pt.y()) / ncount * qx,
                vtx_new_pt.z() + (min_wcp.point.z() - vtx_new_pt.z()) / ncount * qx
            );
            auto tmp_knn_results = temp_cluster.kd_steiner_knn(1, tmp_p, "steiner_pc");
            size_t tmp_idx = tmp_knn_results[0].first;
            Facade::geo_point_t tmp_pt(x_coords[tmp_idx], y_coords[tmp_idx], z_coords[tmp_idx]);
            
            // Skip if too far from target point or duplicate
            if ((tmp_pt - tmp_p).magnitude() > 0.3 * units::cm) continue;
            double dis_to_last = (tmp_pt - new_list.back().point).magnitude();
            if (dis_to_last > 0.01 * units::cm && (tmp_pt - min_wcp.point).magnitude() > 0.01 * units::cm) {
                // cumulative_ray_length += dis_to_last;
                new_list.push_back(WCPoint{tmp_pt});
            }
        }
        // cumulative_ray_length += (min_wcp.point - new_list.back().point).magnitude();
        new_list.push_back(WCPoint{min_wcp.point});
        
        // Replace segment path
        std::list<WCPoint> old_list(seg_wcpts.begin(), seg_wcpts.end());
        
        if (flag_front) {
            // Remove old path from front to min_wcp
            while (!old_list.empty() && (old_list.front().point - min_wcp.point).magnitude() > 0.01 * units::cm) {
                old_list.pop_front();
            }
            if (!old_list.empty()) old_list.pop_front();
            
            // Prepend new path
            for (auto it = new_list.rbegin(); it != new_list.rend(); ++it) {
                old_list.push_front(*it);
            }
        } else {
            // Remove old path from back to min_wcp
            while (!old_list.empty() && (old_list.back().point - min_wcp.point).magnitude() > 0.01 * units::cm) {
                old_list.pop_back();
            }
            if (!old_list.empty()) old_list.pop_back();
            
            // Append new path
            for (auto it = new_list.rbegin(); it != new_list.rend(); ++it) {
                old_list.push_back(*it);
            }
        }
        
        // Update segment wcpts
        seg_wcpts.clear();
        seg_wcpts.reserve(old_list.size());
        std::copy(old_list.begin(), old_list.end(), std::back_inserter(seg_wcpts));
        
        // Clear fit data - will need to be recalculated
        segments.at(i)->clear_fit(dv);
    }
    
    // Update vertex with new fit position and wire coordinates
    auto& vtx_fit = vtx->fit();
    vtx_fit.point = fit_pos;
    vtx_fit.pu = offset_u + (slope_yu * fit_pos_raw.y() + slope_zu * fit_pos_raw.z());
    vtx_fit.pv = offset_v + (slope_yv * fit_pos_raw.y() + slope_zv * fit_pos_raw.z());
    vtx_fit.pw = offset_w + (slope_yw * fit_pos_raw.y() + slope_zw * fit_pos_raw.z());
    vtx_fit.pt = offset_t + slope_x * fit_pos_raw.x();
    vtx_fit.paf = std::make_pair(apa, face);
    vtx_fit.flag_fix = true;
    
    // Update vertex wcpt
    vtx->wcpt().point = vtx_new_pt;

    return true;
}