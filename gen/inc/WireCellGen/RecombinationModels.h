#ifndef WIRECELLGEN_RECOMBINATIONMODEL
#define WIRECELLGEN_RECOMBINATIONMODEL

#include "WireCellUtil/Units.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IRecombinationModel.h"

#include "WireCellUtil/Units.h"

namespace WireCell {
    namespace Gen {

        /// Model for a MIP, dQ = (Rmip/Wi)*dE
        class MipRecombination : public IRecombinationModel, public IConfigurable {
            double m_rmip, m_wi;

           public:
            MipRecombination(double Rmip = 0.7, double Wi = 23.6 * units::eV / (-1 * units::eplus));
            virtual ~MipRecombination();
            virtual double operator()(double dE, double dX = 0.0);
            virtual double dE(double dQ, double dX);
            virtual void configure(const WireCell::Configuration& config);
            virtual WireCell::Configuration default_configuration() const;
        };

        /// Birks model, R = a/(1+b*dE/dX), dQ = (R/Wi)*dE
        /// a = A3t, b=k3t/(Efield*rho) as defined in:
        /// http://lar.bnl.gov/properties/pass.html#recombination
        class BirksRecombination : public IRecombinationModel, public IConfigurable {
            double m_a3t, m_k3t, m_efield, m_rho, m_wi;

           public:
            BirksRecombination(double Efield = 500 * units::volt / units::cm, double A3t = 0.8,
                               double k3t = 0.0486 * units::gram / (units::MeV * units::cm2) *
                                            (units::kilovolt / units::cm),
                               double rho = 1.396 * units::gram / units::cm3,
                               double Wi = 23.6 * units::eV / (-1 * units::eplus));
            virtual ~BirksRecombination();
            virtual double operator()(double dE, double dX = 0.0);
            virtual double dE(double dQ, double dX);
            virtual void configure(const WireCell::Configuration& config);
            virtual WireCell::Configuration default_configuration() const;
        };

        /// Modified Box model, R = ln(a+b*dE/dX)/(b*dE/dX), dQ = (R/Wi)*dE
        /// a=A, b=B/(Efield*rho) as defined in:
        /// http://lar.bnl.gov/properties/pass.html#recombination
        class BoxRecombination : public IRecombinationModel, public IConfigurable {
            double m_efield, m_a, m_b, m_rho, m_wi;

           public:
            BoxRecombination(double Efield = 500 * units::volt / units::cm, double A = 0.930,
                             double B = 0.212 * units::gram / (units::MeV * units::cm2) * (units::kilovolt / units::cm),
                             double rho = 1.396 * units::gram / units::cm3,
                             double Wi = 23.6 * units::eV / (-1 * units::eplus));
            virtual ~BoxRecombination();
            virtual double operator()(double dE, double dX = 0.0);
            virtual double dE(double dQ, double dX);
            virtual void configure(const WireCell::Configuration& config);
            virtual WireCell::Configuration default_configuration() const;
        };

        /// Modified Box model with a free power on dE/dx:
        ///     u = k * (dEdx/pivot)^p,  R = ln(A + u)/u,  dQ = C * (R/Wi) * dE
        /// Fitted to SBND stopping-track dQ/dx vs residual range at 0.5 kV/cm
        /// (wcp-porting-validation sbnd_xin/docs/55 sec 7g; canonical parameters in
        /// nusel_display/stm_ref_dqdx.json).  The drift field is not an explicit
        /// parameter: it is absorbed into the fitted (k, p).
        ///
        /// Unlike the classes above, the defaults follow the *practical-unit*
        /// convention the jsonnet configs actually use (plain numbers: Wi in MeV,
        /// pivot in MeV/cm), so an unconfigured instance is the canonical SBND fit
        /// rather than silently broken (the unit-bearing ctor defaults above make an
        /// unconfigured BoxRecombination::dE() return 0 for any input).
        ///
        /// There is no closed-form inverse.  dE() solves the forward relation by
        /// fixed-count bisection on the monotone branch dEdx in (0, dedx_max]; the
        /// forward peaks at ~77.4 MeV/cm at the canonical parameters, so dQ/dx above
        /// the peak saturates the inverse at dedx_max.  Callers in clus/ clamp their
        /// dE/dx to [0, 50] MeV/cm, well inside the branch.
        class PowerBoxRecombination : public IRecombinationModel, public IConfigurable {
            double m_a, m_k, m_p, m_c, m_pivot, m_wi, m_dedx_max;

            double forward_dqdx(double dedx) const;  // e/cm from MeV/cm

           public:
            PowerBoxRecombination(double A = 0.93, double k = 0.282371, double p = 1.362179,
                                  double C = 0.855175,
                                  double pivot = 2.1,        // MeV/cm
                                  double Wi = 23.6e-6,       // MeV per electron-ion pair
                                  double dedx_max = 77.0);   // MeV/cm, end of the monotone branch
            virtual ~PowerBoxRecombination();
            virtual double operator()(double dE, double dX = 0.0);
            virtual double dE(double dQ, double dX);
            virtual void configure(const WireCell::Configuration& config);
            virtual WireCell::Configuration default_configuration() const;
        };

    }  // namespace Gen
}  // namespace WireCell

#endif
