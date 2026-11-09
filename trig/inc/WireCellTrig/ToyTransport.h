/** A toy transport of MeV-scale particles in liquid argon to energy deposits.

    This stands in for Geant4 to turn the final-state particles of supernova
    neutrino interactions and radioactive decays into energy deposits.  It
    aims to get the topology right at the level that matters for TPs (how far
    MeV gammas carry energy, how long electron tracks are, how energy is
    shared among them), not precision.  The medium is infinite liquid argon:
    deposits outside the detector are dropped later, by the drifter.

    - e-, e+ :: steps of fixed length losing energy at the collision stopping
      power (a table, approximate ESTAR values), with bremsstrahlung photons
      above brems_cut emitted per step (complete-screening spectrum, rate
      from the radiation length) and Gaussian multiple scattering (Highland).
      Below electron_cut the rest is deposited at once.  A positron then
      annihilates to two back-to-back 511 keV photons.
    - gamma :: exponential steps with Compton (Klein-Nishina, exact),
      photoelectric (power-law approximation) and pair (table) attenuation;
      Compton angles from Klein-Nishina.  Below photon_cut the photon is
      absorbed.
    - p, alpha, other ions :: one deposit of all the kinetic energy over a
      short length (a Bragg-Kleeman range estimate), for recombination.
    - neutrons, neutrinos, heavier nuclei recoils :: not transported; their
      energy is reported as lost.

    Times: particles move at c and time of flight is ignored, deposits take
    the time of the primary.

    Units: WCT (MeV = 1, mm = 1).
 */

#ifndef WIRECELLTRIG_TOYTRANSPORT
#define WIRECELLTRIG_TOYTRANSPORT

#include "WireCellUtil/Point.h"
#include "WireCellUtil/Units.h"

#include <limits>
#include <random>
#include <vector>

namespace WireCell::Trig {

    struct Deposit {
        Point pos;
        double energy{0};  // deposited energy
        double length{0};  // the path length over which it was deposited
        int pdg{0};        // the particle that deposited it
    };

    class ToyTransport {
       public:
        using rng_t = std::mt19937_64;

        struct Config {
            double density{1.396 * units::g / units::cm3};
            double radiation_length{14.0 * units::cm};
            double step{0.5 * units::mm};
            double electron_cut{10 * units::keV};
            double photon_cut{10 * units::keV};
            double brems_cut{100 * units::keV};
            bool multiple_scattering{true};
        };

        /// What a transport did with the energy.
        struct Tally {
            double deposited{0};
            double lost{0};  // not transported (neutrons, ...)
        };

        ToyTransport();
        explicit ToyTransport(const Config& cfg);

        /// Transport one particle of the given PDG code and kinetic energy
        /// from pos along dir (need not be unit).  Deposits are appended.
        Tally operator()(rng_t& rng, int pdg, double kinetic, const Vector& dir, const Point& pos,
                         std::vector<Deposit>& out) const;

        // The physics, public for testing.  Per length.

        /// Electron collision stopping power, energy per length.
        double electron_dedx(double kinetic) const;
        /// Photon attenuation coefficients, per length.
        double compton_mu(double energy) const;
        double photo_mu(double energy) const;
        double pair_mu(double energy) const;
        double total_mu(double energy) const;

        const Config& config() const { return m_cfg; }

       private:
        Config m_cfg;
        double m_electron_density;  // per volume

        struct Particle {
            int pdg;
            double kinetic;
            Vector dir;
            Point pos;
        };
        void electron(rng_t& rng, const Particle& p, std::vector<Particle>& stack, std::vector<Deposit>& out,
                      Tally& tally) const;
        void photon(rng_t& rng, const Particle& p, std::vector<Particle>& stack, std::vector<Deposit>& out,
                    Tally& tally) const;
        void ion(const Particle& p, std::vector<Deposit>& out, Tally& tally) const;
    };

    /// Return a unit vector at polar angle theta (via its cosine) and
    /// azimuth phi about dir (a unit vector).
    Vector rotate(const Vector& dir, double cos_theta, double phi);

    /// An isotropic unit vector.
    Vector isotropic(ToyTransport::rng_t& rng);

}  // namespace WireCell::Trig

#endif
