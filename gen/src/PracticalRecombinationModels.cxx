#include "WireCellGen/PracticalRecombinationModels.h"

#include "WireCellUtil/NamedFactory.h"

#include <cmath>

WIRECELL_FACTORY(PracticalBoxRecombination, WireCell::Gen::PracticalBoxRecombination,
                 WireCell::IRecombinationModel, WireCell::IConfigurable)

using namespace WireCell;

/*
  Modified Box model, parameters in practical units.

  The two method bodies are the pre-e6fb7ef3 Gen::BoxRecombination arithmetic
  transcribed VERBATIM -- same expression shapes, same operation order -- so
  that a config already written in practical units yields bit-identical
  doubles to what it yielded before that commit.  Do not "simplify" them:
  folding units::MeV (== 1) away or reassociating the products changes the
  last ulp and silently breaks the byte-identity gate (doc 88 sec 3).

  prototype (wcp-uboone-bdt/inc/WCPLEEANA/cuts.h lines 394-395):
      float beta = 0.255;
      float median_dedx = (exp((median_dqdx*43e3) * 23.6e-6*beta/1.38/0.273)
                           - alpha)/(beta/1.38/0.273);
  i.e. dQ/dx in e-/cm and dE/dx in MeV/cm, which is what dE() below computes.
*/
Gen::PracticalBoxRecombination::PracticalBoxRecombination(double Efield, double A, double B, double rho, double Wi)
  : m_efield(Efield)
  , m_a(A)
  , m_b(B)
  , m_rho(rho)
  , m_wi(Wi)
{
}
Gen::PracticalBoxRecombination::~PracticalBoxRecombination() {}
double Gen::PracticalBoxRecombination::operator()(double dE, double dX)
{
    const double tmp = (dE /units::MeV*units::cm/ dX) * m_b / (m_efield * m_rho);
    const double R = std::log(m_a + tmp) / tmp;
    return R * dE / m_wi;
}
double Gen::PracticalBoxRecombination::dE(double dQ, double dX)
{
    const double coeff = m_b / (m_efield * m_rho);
    const double a_exp = std::exp(dQ/dX*units::cm * coeff * m_wi);
    const double numerator = (a_exp - m_a) * units::MeV/units::cm *dX;
    const double denominator = coeff;

    return numerator / denominator;
}
void Gen::PracticalBoxRecombination::configure(const WireCell::Configuration& config)
{
    m_efield = get(config, "Efield", m_efield);
    m_a = get(config, "A", m_a);
    m_b = get(config, "B", m_b);
    m_rho = get(config, "rho", m_rho);
    m_wi = get(config, "Wi", m_wi);
}
WireCell::Configuration Gen::PracticalBoxRecombination::default_configuration() const
{
    Configuration cfg;
    cfg["Efield"] = m_efield;
    cfg["A"] = m_a;
    cfg["B"] = m_b;
    cfg["rho"] = m_rho;
    cfg["Wi"] = m_wi;
    return cfg;
}
