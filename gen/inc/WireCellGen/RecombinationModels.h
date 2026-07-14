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

    }  // namespace Gen
}  // namespace WireCell

#endif
