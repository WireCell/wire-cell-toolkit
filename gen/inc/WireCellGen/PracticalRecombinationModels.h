#ifndef WIRECELLGEN_PRACTICALRECOMBINATIONMODELS
#define WIRECELLGEN_PRACTICALRECOMBINATIONMODELS

#include "WireCellUtil/Units.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IRecombinationModel.h"

namespace WireCell {
    namespace Gen {

        /// Modified Box model with parameters in PRACTICAL units:
        /// Efield in kV/cm, B in (kV/cm)(g/cm^2)/MeV, rho in g/cm^3, and the
        /// quenching term formed from dE/dX expressed in MeV/cm.  This is the
        /// convention of the LArSoft ModBox constants (kModBoxA/kModBoxB), of
        /// the WCP prototype (wcp-uboone-bdt/inc/WCPLEEANA/cuts.h lines
        /// 394-395), and of every recombination block in this tree's detector
        /// configs, which pass bare numbers such as B=0.255, Efield=0.273.
        ///
        /// Gen::BoxRecombination is the same physics with parameters in the
        /// WCT system of units (its own defaults carry unit factors, e.g.
        /// B = 0.212 * gram/(MeV*cm2) * (kilovolt/cm)).  Upstream e6fb7ef3
        /// made that class consistently WCT-unit; a config that supplies
        /// practical-unit numbers to it is then wrong by a factor of
        /// units::cm/units::MeV = 10 in the quenching term.  The two classes
        /// are NOT interchangeable -- pick the one matching how the
        /// parameters in your config are expressed.
        ///
        /// The bodies below are the pre-e6fb7ef3 arithmetic transcribed
        /// verbatim so that configs already written in practical units keep
        /// producing bit-identical numbers (doc 88).
        class PracticalBoxRecombination : public IRecombinationModel, public IConfigurable {
            double m_efield, m_a, m_b, m_rho, m_wi;

           public:
            /// Defaults are the canonical LArSoft ModBox point: R = 0.705 for
            /// a 2.1 MeV/cm MIP at 500 V/cm.  Wi is unchanged from the
            /// WCT-unit siblings because energy is MeV in both conventions.
            PracticalBoxRecombination(double Efield = 0.5,   // kV/cm
                                      double A = 0.930,      // dimensionless
                                      double B = 0.212,      // (kV/cm)(g/cm^2)/MeV
                                      double rho = 1.396,    // g/cm^3
                                      double Wi = 23.6 * units::eV / (-1 * units::eplus));
            virtual ~PracticalBoxRecombination();
            virtual double operator()(double dE, double dX = 0.0);
            virtual double dE(double dQ, double dX);
            virtual void configure(const WireCell::Configuration& config);
            virtual WireCell::Configuration default_configuration() const;
        };

    }  // namespace Gen
}  // namespace WireCell

#endif
