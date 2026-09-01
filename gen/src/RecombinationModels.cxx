#include "WireCellGen/RecombinationModels.h"

#include "WireCellUtil/NamedFactory.h"

#include <cmath>

WIRECELL_FACTORY(MipRecombination, WireCell::Gen::MipRecombination, WireCell::IRecombinationModel,
                 WireCell::IConfigurable)
WIRECELL_FACTORY(BirksRecombination, WireCell::Gen::BirksRecombination, WireCell::IRecombinationModel,
                 WireCell::IConfigurable)
WIRECELL_FACTORY(BoxRecombination, WireCell::Gen::BoxRecombination, WireCell::IRecombinationModel,
                 WireCell::IConfigurable)
WIRECELL_FACTORY(PowerBoxRecombination, WireCell::Gen::PowerBoxRecombination, WireCell::IRecombinationModel,
                 WireCell::IConfigurable)

using namespace WireCell;

/*
  MIP recombination model
*/
Gen::MipRecombination::MipRecombination(double Rmip, double Wi)
  : m_rmip(Rmip)
  , m_wi(Wi)
{
}
Gen::MipRecombination::~MipRecombination() {}
double Gen::MipRecombination::operator()(double dE, double dX) { return m_rmip * dE / m_wi; }
double Gen::MipRecombination::dE(double dQ, double dX) { return dQ * m_wi / m_rmip; }
void Gen::MipRecombination::configure(const WireCell::Configuration& config)
{
    m_rmip = get(config, "Rmip", m_rmip);
    m_wi = get(config, "Wi", m_wi);
}
WireCell::Configuration Gen::MipRecombination::default_configuration() const
{
    Configuration cfg;
    cfg["Rmip"] = m_rmip;
    cfg["Wi"] = m_wi;
    return cfg;
}

/*
  Birks Recombination Model
 */
Gen::BirksRecombination::BirksRecombination(double Efield, double A3t, double k3t, double rho, double Wi)
  : m_a3t(A3t)
  , m_k3t(k3t)
  , m_efield(Efield)
  , m_rho(rho)
  , m_wi(Wi)
{
}
Gen::BirksRecombination::~BirksRecombination() {}
double Gen::BirksRecombination::operator()(double dE, double dX)
{
    // dE/dX and k3t/(Efield*rho) are both expressed in the WCT system of
    // units so their product is already dimensionless.  (An earlier version
    // multiplied dE/dX by a spurious units::cm, inflating the quenching
    // term 10x and yielding R~0.32 instead of ~0.70 for a MIP at 500 V/cm.)
    const double R = m_a3t / (1 + (dE / dX) * m_k3t / (m_efield * m_rho));
    return R * dE / m_wi;
}
double Gen::BirksRecombination::dE(double dQ, double dX)
{
    const double numerator = dQ;
    const double denominator = m_a3t/m_wi - dQ/dX * m_k3t/(m_efield*m_rho);

    return numerator / denominator;
}
void Gen::BirksRecombination::configure(const WireCell::Configuration& config)
{
    m_a3t = get(config, "A3t", m_a3t);
    m_k3t = get(config, "k3t", m_k3t);
    m_efield = get(config, "Efield", m_efield);
    m_rho = get(config, "rho", m_rho);
    m_wi = get(config, "Wi", m_wi);
}
WireCell::Configuration Gen::BirksRecombination::default_configuration() const
{
    Configuration cfg;
    cfg["A3t"] = m_a3t;
    cfg["k3t"] = m_k3t;
    cfg["Efield"] = m_efield;
    cfg["rho"] = m_rho;
    cfg["Wi"] = m_wi;
    return cfg;
}

/*
  Modified Box Model
*/
Gen::BoxRecombination::BoxRecombination(double Efield, double A, double B, double rho, double Wi)
  : m_efield(Efield)
  , m_a(A)
  , m_b(B)
  , m_rho(rho)
  , m_wi(Wi)
{
}
Gen::BoxRecombination::~BoxRecombination() {}
double Gen::BoxRecombination::operator()(double dE, double dX)
{
    // See the units note in BirksRecombination::operator() -- dE/dX times
    // B/(Efield*rho) is already dimensionless in the WCT system of units.
    const double tmp = (dE / dX) * m_b / (m_efield * m_rho);
    const double R = std::log(m_a + tmp) / tmp;
    return R * dE / m_wi;
}
double Gen::BoxRecombination::dE(double dQ, double dX)
{
    const double coeff = m_b / (m_efield * m_rho);
    const double a_exp = std::exp(dQ/dX * coeff * m_wi);
    const double numerator = (a_exp - m_a) * dX;
    const double denominator = coeff;

    return numerator / denominator;
}
void Gen::BoxRecombination::configure(const WireCell::Configuration& config)
{
    m_efield = get(config, "Efield", m_efield);
    m_a = get(config, "A", m_a);
    m_b = get(config, "B", m_b);
    m_rho = get(config, "rho", m_rho);
    m_wi = get(config, "Wi", m_wi);
}
WireCell::Configuration Gen::BoxRecombination::default_configuration() const
{
    Configuration cfg;
    cfg["Efield"] = m_efield;
    cfg["A"] = m_a;
    cfg["B"] = m_b;
    cfg["rho"] = m_rho;
    cfg["Wi"] = m_wi;
    return cfg;
}

/*
  Modified Box Model with a free power on dE/dx (SBND stopping-track fit)
*/
Gen::PowerBoxRecombination::PowerBoxRecombination(double A, double k, double p, double C, double pivot, double Wi,
                                                  double dedx_max)
  : m_a(A)
  , m_k(k)
  , m_p(p)
  , m_c(C)
  , m_pivot(pivot)
  , m_wi(Wi)
  , m_dedx_max(dedx_max)
{
}
Gen::PowerBoxRecombination::~PowerBoxRecombination() {}
double Gen::PowerBoxRecombination::forward_dqdx(double dedx) const
{
    if (dedx <= 0) return 0.0;
    const double u = m_k * std::pow(dedx / m_pivot, m_p);
    const double R = std::log(m_a + u) / u;
    return m_c * R * dedx / m_wi;  // e/cm; negative below the A<1 zero crossing (~0.75 MeV/cm)
}
double Gen::PowerBoxRecombination::operator()(double dE, double dX)
{
    const double dedx = dE / units::MeV * units::cm / dX;
    return forward_dqdx(dedx) * (dX / units::cm);
}
double Gen::PowerBoxRecombination::dE(double dQ, double dX)
{
    const double dqdx = dQ / dX * units::cm;
    if (!(dqdx > 0)) return 0.0;
    if (dqdx >= forward_dqdx(m_dedx_max)) {  // above the monotone branch: saturate
        return m_dedx_max * units::MeV / units::cm * dX;
    }
    // forward_dqdx is monotone increasing on (0, m_dedx_max], from -inf at 0+.
    // Fixed-count bisection keeps the solve deterministic.
    double lo = 0.0, hi = m_dedx_max;
    for (int i = 0; i != 60; ++i) {
        const double mid = 0.5 * (lo + hi);
        if (forward_dqdx(mid) < dqdx) { lo = mid; }
        else { hi = mid; }
    }
    return 0.5 * (lo + hi) * units::MeV / units::cm * dX;
}
void Gen::PowerBoxRecombination::configure(const WireCell::Configuration& config)
{
    m_a = get(config, "A", m_a);
    m_k = get(config, "k", m_k);
    m_p = get(config, "p", m_p);
    m_c = get(config, "C", m_c);
    m_pivot = get(config, "pivot", m_pivot);
    m_wi = get(config, "Wi", m_wi);
    m_dedx_max = get(config, "dedx_max", m_dedx_max);
}
WireCell::Configuration Gen::PowerBoxRecombination::default_configuration() const
{
    Configuration cfg;
    cfg["A"] = m_a;
    cfg["k"] = m_k;
    cfg["p"] = m_p;
    cfg["C"] = m_c;
    cfg["pivot"] = m_pivot;
    cfg["Wi"] = m_wi;
    cfg["dedx_max"] = m_dedx_max;
    return cfg;
}
