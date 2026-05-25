#include "WireCellMatch/SemiAnalyticalModel.h"

#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Units.h"

// Boost 1.82's special-functions headers transitively include the deprecated
// boost/functional.hpp which uses std::unary_function/binary_function; gcc 12+
// flags them as -Werror=deprecated-declarations. Silence the warning around
// the boost includes only.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <boost/math/policies/policy.hpp>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#pragma GCC diagnostic pop

#include <algorithm>
#include <cmath>
#include <limits>

using namespace WireCell;
using namespace WireCell::Match;

namespace {

    using noLDoublePromote =
        boost::math::policies::policy<boost::math::policies::promote_double<false>>;

    constexpr double kPi = M_PI;

    // -----------------------------------------------------------------------
    // The following helpers are ported verbatim from
    // larsim/PhotonPropagation/PhotonPropagationUtils.{h,cxx}.
    // They are kept inside the anonymous namespace so they don't leak.
    double fast_acos(double xin)
    {
        constexpr double a3 = -2.08730442907856008e-02;
        constexpr double a2 =  7.68769404161671888e-02;
        constexpr double a1 = -2.12871094165645952e-01;
        constexpr double a0 =  1.57075835365209659e+00;
        const double x = std::abs(xin);
        double ret = a3;
        ret *= x; ret += a2;
        ret *= x; ret += a1;
        ret *= x; ret += a0;
        ret *= std::sqrt(1.0 - x);
        if (xin >= 0) return ret;
        return M_PI - ret;
    }

    double interpolate(const std::vector<double>& xData,
                       const std::vector<double>& yData,
                       const double x,
                       const bool extrapolate,
                       std::size_t i = 0)
    {
        if (i == 0) {
            const std::size_t n = xData.size();
            if (x >= xData[n - 2]) { i = n - 2; }
            else {
                while (x > xData[i + 1]) i++;
            }
        }
        const double xL = xData[i];
        const double xR = xData[i + 1];
        const double yL = yData[i];
        const double yR = yData[i + 1];
        if (!extrapolate) {
            if (x < xL) return yL;
            if (x > xR) return yL;
        }
        const double dydx = (yR - yL) / (xR - xL);
        return yL + dydx * (x - xL);
    }

    double interpolate2(const std::vector<double>& xDistances,
                        const std::vector<double>& rDistances,
                        const std::vector<std::vector<std::vector<double>>>& parameters,
                        const double x,
                        const double r,
                        const std::size_t k)
    {
        const std::size_t nbins_r = parameters[k].size();
        std::vector<double> interp_vals(nbins_r, 0.);

        std::size_t idx = 0;
        const std::size_t n = xDistances.size();
        if (x >= xDistances[n - 2]) { idx = n - 2; }
        else {
            while (x > xDistances[idx + 1]) idx++;
        }
        for (std::size_t i = 0; i < nbins_r; ++i) {
            interp_vals[i] = interpolate(xDistances, parameters[k][i], x, false, idx);
        }
        return interpolate(rDistances, interp_vals, r, false);
    }

    template <typename T>
    inline bool isApproximatelyZero(T a, T tol = std::numeric_limits<T>::epsilon())
    {
        return std::fabs(a) <= tol;
    }
    template <typename T>
    inline bool isApproximatelyEqual(T a, T b, T tol = std::numeric_limits<T>::epsilon())
    {
        const T diff = std::fabs(a - b);
        if (diff <= tol) return true;
        if (diff < std::fmax(std::fabs(a), std::fabs(b)) * tol) return true;
        return false;
    }
    template <typename T>
    inline bool isDefinitelyGreaterThan(T a, T b, T tol = std::numeric_limits<T>::epsilon())
    {
        const T diff = a - b;
        if (diff > tol) return true;
        if (diff > std::fmax(std::fabs(a), std::fabs(b)) * tol) return true;
        return false;
    }
    template <typename T>
    inline bool isDefinitelyLessThan(T a, T b, T tol = std::numeric_limits<T>::epsilon())
    {
        const T diff = a - b;
        if (diff < tol) return true;
        if (diff < std::fmax(std::fabs(a), std::fabs(b)) * tol) return true;
        return false;
    }

    // -----------------------------------------------------------------------
    // JSON -> nested std::vector helpers. Configuration is Json::Value.
    std::vector<double> as_vd(const Configuration& v)
    {
        std::vector<double> out;
        if (v.isArray()) {
            out.reserve(v.size());
            for (const auto& x : v) out.push_back(x.asDouble());
        }
        return out;
    }
    std::vector<std::vector<double>> as_vvd(const Configuration& v)
    {
        std::vector<std::vector<double>> out;
        if (v.isArray()) {
            out.reserve(v.size());
            for (const auto& row : v) out.push_back(as_vd(row));
        }
        return out;
    }
    std::vector<std::vector<std::vector<double>>> as_vvvd(const Configuration& v)
    {
        std::vector<std::vector<std::vector<double>>> out;
        if (v.isArray()) {
            out.reserve(v.size());
            for (const auto& m : v) out.push_back(as_vvd(m));
        }
        return out;
    }

} // namespace


SemiAnalyticalModel::SemiAnalyticalModel(const Configuration& VUVHits,
                                         const Configuration& VISHits,
                                         const Geometry& geom,
                                         const std::vector<OpticalDetector>& opdets,
                                         bool doReflectedLight)
    : m_geom(geom)
    , m_opdets(opdets)
    , m_doReflectedLight(doReflectedLight)
{
    m_isFlatPDCorr = get<bool>(VUVHits, "FlatPDCorr", false);
    m_isDomePDCorr = get<bool>(VUVHits, "DomePDCorr", false);
    m_delta_angulo_vuv = get<double>(VUVHits, "delta_angulo_vuv", 10.0);
    m_pmt_radius = get<double>(VUVHits, "PMT_radius", 10.16);
    m_maxPDDistance = get<double>(VUVHits, "MaxPDDistance", 1000.0);

    if (!m_isFlatPDCorr && !m_isDomePDCorr) {
        raise<ValueError>("SemiAnalyticalModel: neither FlatPDCorr nor DomePDCorr enabled");
    }
    if (m_isFlatPDCorr) {
        m_GHvuvpars_flat = as_vvd(VUVHits["GH_PARS_flat"]);
        m_border_corr_angulo_flat = as_vd(VUVHits["GH_border_angulo_flat"]);
        m_border_corr_flat = as_vvd(VUVHits["GH_border_flat"]);
    }
    if (m_isDomePDCorr) {
        m_GHvuvpars_dome = as_vvd(VUVHits["GH_PARS_dome"]);
        m_border_corr_angulo_dome = as_vd(VUVHits["GH_border_angulo_dome"]);
        m_border_corr_dome = as_vvd(VUVHits["GH_border_dome"]);
    }

    if (m_doReflectedLight) {
        m_delta_angulo_vis = get<double>(VISHits, "delta_angulo_vis", 10.0);
        if (m_isFlatPDCorr) {
            m_vis_distances_x_flat = as_vd(VISHits["VIS_distances_x_flat"]);
            m_vis_distances_r_flat = as_vd(VISHits["VIS_distances_r_flat"]);
            m_vispars_flat = as_vvvd(VISHits["VIS_correction_flat"]);
        }
        if (m_isDomePDCorr) {
            m_vis_distances_x_dome = as_vd(VISHits["VIS_distances_x_dome"]);
            m_vis_distances_r_dome = as_vd(VISHits["VIS_distances_r_dome"]);
            m_vispars_dome = as_vvvd(VISHits["VIS_correction_dome"]);
        }
        m_cathode_plane.h = m_geom.active_size_y;
        m_cathode_plane.w = m_geom.active_size_z;
        m_plane_depth = std::abs(m_geom.cathode_x);
    }
}

// ---------------------------------------------------------------------------
// Direct (VUV) visibilities
void SemiAnalyticalModel::detectedDirectVisibilities(std::vector<double>& vis,
                                                     const WireCell::Point& scintPoint) const
{
    vis.assign(m_opdets.size(), 0.);
    for (std::size_t i = 0; i < m_opdets.size(); ++i) {
        const auto& od = m_opdets[i];
        // same-TPC visibility (port of SBNDOpticalPath_tool): only OpDets on
        // the same X-sign as the scintillation point are visible.
        if ((scintPoint.x() < 0.) != (od.center.x() < 0.)) continue;

        const Vector relative = scintPoint - od.center;
        const double distance = relative.magnitude();
        if (distance > m_maxPDDistance) continue;
        vis[i] = VUVVisibility(scintPoint, od);
    }
}

double SemiAnalyticalModel::VUVVisibility(const WireCell::Point& scintPoint,
                                          const OpticalDetector& opDet) const
{
    const Vector relative = scintPoint - opDet.center;
    const double distance = relative.magnitude();
    // orientation 0 (anode/cathode) only
    const double cosine = std::abs(relative.x()) / distance;
    const double theta = fast_acos(cosine) * 180. / kPi;

    double solid_angle = 0.0;
    if (opDet.type == 0) {  // Arapuca / rectangle
        const Vector abs_rel(std::abs(relative.x()),
                             std::abs(relative.y()),
                             std::abs(relative.z()));
        solid_angle = Rectangle_SolidAngle(Dims{opDet.h, opDet.w}, abs_rel, opDet.orientation);
    }
    else if (opDet.type == 1) {  // Dome PMT
        solid_angle = Omega_Dome_Model(distance, theta);
    }
    else {
        raise<ValueError>("SemiAnalyticalModel::VUVVisibility: unsupported OpDet type %d", opDet.type);
    }

    const double visibility_geo = std::exp(-distance / m_geom.vuv_absorption_length) *
                                  (solid_angle / (4. * kPi));

    const std::size_t j = static_cast<std::size_t>(theta / m_delta_angulo_vuv);

    // border radius in cathode-plane coordinates (orientation 0).
    const double r = std::hypot(scintPoint.y() - m_geom.active_center_y,
                                scintPoint.z() - m_geom.active_center_z);

    double pars_ini[4] = {0., 0., 0., 0.};
    double s1 = 0., s2 = 0., s3 = 0.;
    if (opDet.type == 0 && m_isFlatPDCorr) {
        pars_ini[0] = m_GHvuvpars_flat[0][j];
        pars_ini[1] = m_GHvuvpars_flat[1][j];
        pars_ini[2] = m_GHvuvpars_flat[2][j];
        pars_ini[3] = m_GHvuvpars_flat[3][j];
        s1 = interpolate(m_border_corr_angulo_flat, m_border_corr_flat[0], theta, true);
        s2 = interpolate(m_border_corr_angulo_flat, m_border_corr_flat[1], theta, true);
        s3 = interpolate(m_border_corr_angulo_flat, m_border_corr_flat[2], theta, true);
    }
    else if (opDet.type == 1 && m_isDomePDCorr) {
        pars_ini[0] = m_GHvuvpars_dome[0][j];
        pars_ini[1] = m_GHvuvpars_dome[1][j];
        pars_ini[2] = m_GHvuvpars_dome[2][j];
        pars_ini[3] = m_GHvuvpars_dome[3][j];
        s1 = interpolate(m_border_corr_angulo_dome, m_border_corr_dome[0], theta, true);
        s2 = interpolate(m_border_corr_angulo_dome, m_border_corr_dome[1], theta, true);
        s3 = interpolate(m_border_corr_angulo_dome, m_border_corr_dome[2], theta, true);
    }
    else {
        raise<ValueError>("SemiAnalyticalModel::VUVVisibility: missing corrections for type %d", opDet.type);
    }

    pars_ini[0] += s1 * r;
    pars_ini[1] += s2 * r;
    pars_ini[2] += s3 * r;

    double GH = Gaisser_Hillas(distance, pars_ini);
    if (!std::isfinite(GH) || GH < 0 || GH > 10) GH = 0;

    return GH * visibility_geo / cosine;
}

// ---------------------------------------------------------------------------
// Reflected (VIS) visibilities
void SemiAnalyticalModel::detectedReflectedVisibilities(std::vector<double>& vis,
                                                        const WireCell::Point& scintPoint) const
{
    vis.assign(m_opdets.size(), 0.);
    if (!m_doReflectedLight) return;

    // Step 1: visibility of VUV photons on the cathode-plane hotspot.
    const double plane_depth = scintPoint.x() < 0. ? -m_plane_depth : m_plane_depth;
    const Vector scint_rel(std::abs(scintPoint.x() - plane_depth),
                           std::abs(scintPoint.y() - m_geom.active_center_y),
                           std::abs(scintPoint.z() - m_geom.active_center_z));

    const double solid_angle_cathode = Rectangle_SolidAngle(m_cathode_plane, scint_rel, 0);

    const double distance_cathode = std::abs(plane_depth - scintPoint.x());
    const double cathode_vis_geo = std::exp(-distance_cathode / m_geom.vuv_absorption_length) *
                                   (solid_angle_cathode / (4. * kPi));

    // GH correction at theta=0 (on-axis), use flat PD parameters.
    const double r = std::hypot(scintPoint.y() - m_geom.active_center_y,
                                scintPoint.z() - m_geom.active_center_z);
    double pars_ini[4] = {0., 0., 0., 0.};
    double s1 = 0., s2 = 0., s3 = 0.;
    if (!m_isFlatPDCorr) {
        raise<ValueError>("SemiAnalyticalModel: flat PD VUV correction required for VIS calculation");
    }
    pars_ini[0] = m_GHvuvpars_flat[0][0];
    pars_ini[1] = m_GHvuvpars_flat[1][0];
    pars_ini[2] = m_GHvuvpars_flat[2][0];
    pars_ini[3] = m_GHvuvpars_flat[3][0];
    s1 = interpolate(m_border_corr_angulo_flat, m_border_corr_flat[0], 0., true);
    s2 = interpolate(m_border_corr_angulo_flat, m_border_corr_flat[1], 0., true);
    s3 = interpolate(m_border_corr_angulo_flat, m_border_corr_flat[2], 0., true);
    pars_ini[0] += s1 * r;
    pars_ini[1] += s2 * r;
    pars_ini[2] += s3 * r;

    double GH = Gaisser_Hillas(distance_cathode, pars_ini);
    if (!std::isfinite(GH) || GH < 0 || GH > 10) GH = 0;
    const double cathode_visibility_rec = GH * cathode_vis_geo;

    // Step 2: per-PD visibility from the hotspot.
    const Point hotspot(plane_depth, scintPoint.y(), scintPoint.z());
    for (std::size_t i = 0; i < m_opdets.size(); ++i) {
        const auto& od = m_opdets[i];
        if ((scintPoint.x() < 0.) != (od.center.x() < 0.)) continue;
        vis[i] = VISVisibility(scintPoint, od, cathode_visibility_rec, hotspot);
    }
}

double SemiAnalyticalModel::VISVisibility(const WireCell::Point& scintPoint,
                                          const OpticalDetector& opDet,
                                          double cathode_visibility,
                                          const WireCell::Point& hotspot) const
{
    const double plane_depth = scintPoint.x() < 0. ? -m_plane_depth : m_plane_depth;

    const Vector emission_rel = hotspot - opDet.center;
    const double distance_vis = emission_rel.magnitude();
    const double cosine_vis = std::abs(emission_rel.x()) / distance_vis;
    const double theta_vis = fast_acos(cosine_vis) * 180. / kPi;

    double solid_angle_detector = 0.0;
    if (opDet.type == 0) {  // Arapuca
        const Vector abs_emission_rel(std::abs(emission_rel.x()),
                                      std::abs(emission_rel.y()),
                                      std::abs(emission_rel.z()));
        solid_angle_detector = Rectangle_SolidAngle(Dims{opDet.h, opDet.w}, abs_emission_rel, opDet.orientation);
    }
    else if (opDet.type == 1) {  // Dome PMT
        solid_angle_detector = Omega_Dome_Model(distance_vis, theta_vis);
    }
    else {
        raise<ValueError>("SemiAnalyticalModel::VISVisibility: unsupported OpDet type %d", opDet.type);
    }

    const double visibility_geo = (solid_angle_detector / (2. * kPi)) * cathode_visibility;

    const std::size_t k = static_cast<std::size_t>(theta_vis / m_delta_angulo_vis);
    const double rxy = std::hypot(scintPoint.y() - m_geom.active_center_y,
                                  scintPoint.z() - m_geom.active_center_z);
    const double d_c = std::abs(scintPoint.x() - plane_depth);

    double border_correction = 0.;
    if (opDet.type == 0 && m_isFlatPDCorr) {
        border_correction =
            interpolate2(m_vis_distances_x_flat, m_vis_distances_r_flat, m_vispars_flat, d_c, rxy, k);
    }
    else if (opDet.type == 1 && m_isDomePDCorr) {
        border_correction =
            interpolate2(m_vis_distances_x_dome, m_vis_distances_r_dome, m_vispars_dome, d_c, rxy, k);
    }
    else {
        raise<ValueError>("SemiAnalyticalModel::VISVisibility: missing corrections for type %d", opDet.type);
    }
    return border_correction * visibility_geo / cosine_vis;
}

// ---------------------------------------------------------------------------
// Solid-angle / Gaisser-Hillas helpers, ported from larsim.
double SemiAnalyticalModel::Gaisser_Hillas(double x, const double* par) const
{
    const double X_mu_0 = par[3];
    const double Norm = par[0];
    const double Diff = par[1] - X_mu_0;
    const double Term = std::pow((x - X_mu_0) / Diff, Diff / par[2]);
    const double Exp = std::exp((par[1] - x) / par[2]);
    return Norm * Term * Exp;
}

double SemiAnalyticalModel::Disk_SolidAngle(double d, double h, double b) const
{
    if (b <= 0. || d < 0. || h <= 0.) return 0.;
    const double leg2 = (b + d) * (b + d);
    const double aa = std::sqrt(h * h / (h * h + leg2));
    if (isApproximatelyZero(d)) return 2. * kPi * (1. - aa);
    const double bb = 2. * std::sqrt(b * d / (h * h + leg2));
    const double cc = 4. * b * d / leg2;

    if (isDefinitelyGreaterThan(d, b)) {
        try {
            return 2. * aa *
                   (std::sqrt(1. - cc) * boost::math::ellint_3(bb, cc, noLDoublePromote()) -
                    boost::math::ellint_1(bb, noLDoublePromote()));
        }
        catch (std::domain_error&) {
            if (isApproximatelyEqual(d, b, 1e-9)) {
                return kPi - 2. * aa * boost::math::ellint_1(bb, noLDoublePromote());
            }
            return 0.;
        }
    }
    if (isDefinitelyLessThan(d, b)) {
        try {
            return 2. * kPi -
                   2. * aa *
                       (boost::math::ellint_1(bb, noLDoublePromote()) +
                        std::sqrt(1. - cc) * boost::math::ellint_3(bb, cc, noLDoublePromote()));
        }
        catch (std::domain_error&) {
            if (isApproximatelyEqual(d, b, 1e-9)) {
                return kPi - 2. * aa * boost::math::ellint_1(bb, noLDoublePromote());
            }
            return 0.;
        }
    }
    if (isApproximatelyEqual(d, b)) {
        return kPi - 2. * aa * boost::math::ellint_1(bb, noLDoublePromote());
    }
    return 0.;
}

double SemiAnalyticalModel::Rectangle_SolidAngle(double a, double b, double d) const
{
    const double aa = a / (2. * d);
    const double bb = b / (2. * d);
    const double aux = (1. + aa * aa + bb * bb) / ((1. + aa * aa) * (1. + bb * bb));
    return 4. * fast_acos(std::sqrt(aux));
}

double SemiAnalyticalModel::Rectangle_SolidAngle(const Dims& o,
                                                 const WireCell::Vector& v,
                                                 int OpDetOrientation) const
{
    double d1, d2;
    if (OpDetOrientation == 2) {
        d1 = std::abs(v.x());
        d2 = std::abs(v.z());
    }
    else if (OpDetOrientation == 1) {
        d1 = std::abs(v.x());
        d2 = std::abs(v.y());
    }
    else {
        // anode/cathode PD, plane fixed in x [default]
        d1 = std::abs(v.y());
        d2 = std::abs(v.x());
    }
    if (isApproximatelyZero(d1) && isApproximatelyZero(v.z())) {
        return Rectangle_SolidAngle(o.h, o.w, d2);
    }
    if (isDefinitelyGreaterThan(d1, o.h * .5) &&
        isDefinitelyGreaterThan(std::abs(v.z()), o.w * .5)) {
        const double A = d1 - o.h * .5;
        const double B = std::abs(v.z()) - o.w * .5;
        return (Rectangle_SolidAngle(2. * (A + o.h), 2. * (B + o.w), d2) -
                Rectangle_SolidAngle(2. * A, 2. * (B + o.w), d2) -
                Rectangle_SolidAngle(2. * (A + o.h), 2. * B, d2) +
                Rectangle_SolidAngle(2. * A, 2. * B, d2)) *
               .25;
    }
    if (d1 <= o.h * .5 && std::abs(v.z()) <= o.w * .5) {
        const double A = -d1 + o.h * .5;
        const double B = -std::abs(v.z()) + o.w * .5;
        return (Rectangle_SolidAngle(2. * (o.h - A), 2. * (o.w - B), d2) +
                Rectangle_SolidAngle(2. * A, 2. * (o.w - B), d2) +
                Rectangle_SolidAngle(2. * (o.h - A), 2. * B, d2) +
                Rectangle_SolidAngle(2. * A, 2. * B, d2)) *
               .25;
    }
    if (isDefinitelyGreaterThan(d1, o.h * .5) && std::abs(v.z()) <= o.w * .5) {
        const double A = d1 - o.h * .5;
        const double B = -std::abs(v.z()) + o.w * .5;
        return (Rectangle_SolidAngle(2. * (A + o.h), 2. * (o.w - B), d2) -
                Rectangle_SolidAngle(2. * A, 2. * (o.w - B), d2) +
                Rectangle_SolidAngle(2. * (A + o.h), 2. * B, d2) -
                Rectangle_SolidAngle(2. * A, 2. * B, d2)) *
               .25;
    }
    if (d1 <= o.h * .5 && isDefinitelyGreaterThan(std::abs(v.z()), o.w * .5)) {
        const double A = -d1 + o.h * .5;
        const double B = std::abs(v.z()) - o.w * .5;
        return (Rectangle_SolidAngle(2. * (o.h - A), 2. * (B + o.w), d2) -
                Rectangle_SolidAngle(2. * (o.h - A), 2. * B, d2) +
                Rectangle_SolidAngle(2. * A, 2. * (B + o.w), d2) -
                Rectangle_SolidAngle(2. * A, 2. * B, d2)) *
               .25;
    }
    return 0.;
}

double SemiAnalyticalModel::Omega_Dome_Model(double distance, double theta) const
{
    constexpr double par0[9] = {0., 0., 0., 0., 0., 0.597542, 1.00872, 1.46993, 2.04221};
    constexpr double par1[9] = {0., 0., 0.19569, 0.300449, 0.555598, 0.854939, 1.39166, 2.19141, 2.57732};
    constexpr double delta_theta = 10.;
    const int j = static_cast<int>(theta / delta_theta);
    const double b = m_pmt_radius;
    const double d_break = 5. * b;

    if (distance >= d_break) {
        const double R = b - par1[j];
        const double ratio_sq = (R * R) / (distance * distance);
        return 2. * kPi * (1. - std::sqrt(1. - ratio_sq));
    }
    const double R = b - par0[j];
    const double ratio_sq = (R * R) / (distance * distance);
    return 2. * kPi * (1. - std::sqrt(1. - ratio_sq));
}
