/** Primary particles of low-energy interactions: radioactive decays and
    MARLEY supernova neutrino events.

    Units: WCT (MeV = 1).
 */

#ifndef WIRECELLTRIG_KINEMATICS
#define WIRECELLTRIG_KINEMATICS

#include "WireCellTrig/ToyTransport.h"
#include "WireCellUtil/Configuration.h"

#include <string>
#include <vector>

namespace WireCell::Trig {

    /// A final-state particle of an interaction.
    struct Primary {
        int pdg{0};
        double kinetic{0};
        Vector dir{0, 0, 1};
    };

    /// An interaction: its primaries and the energy that characterizes it
    /// (the neutrino energy, or the decay energy release).
    struct Interaction {
        double energy{0};
        std::vector<Primary> primaries;
    };

    /** The kinetic energy spectrum of the electron of a beta-minus decay.

        dN/dT ~ p E (Q - T)^2 F(Z, E) S(T), with the non-relativistic Fermi
        function F = 2 pi eta / (1 - exp(-2 pi eta)), eta = alpha Z E / p, for
        Z of the daughter, and the shape factor S = 1 (allowed) or
        p^2 + q^2 with q = Q - T (unique first forbidden, eg 39Ar, 85Kr).
     */
    class BetaSpectrum {
       public:
        enum class Shape { allowed, unique_first };

        BetaSpectrum(double q, int z_daughter, Shape shape = Shape::allowed, size_t npoints = 1000);

        /// The unnormalized density.
        double density(double kinetic) const;
        /// Sample a kinetic energy.
        double operator()(ToyTransport::rng_t& rng) const;
        double mean() const { return m_mean; }
        double q() const { return m_q; }

       private:
        double m_q;
        int m_z;
        Shape m_shape;
        std::vector<double> m_t, m_cdf;
        double m_mean{0};
    };

    /** A radioactive decay as a set of branches, each emitting some of a
        beta (Q, daughter Z, shape), an alpha (kinetic energy) and gammas
        (energies), isotropic and uncorrelated.

        JSON: {branches: [{probability, beta: {q, z, shape: "allowed" or
        "unique_first"}, alpha, gammas: [...]}]}, energies in WCT units.
        Probabilities are normalized.
     */
    class Decay {
       public:
        struct Branch {
            double probability{1};
            double alpha{0};
            std::vector<double> gammas;
            bool has_beta{false};
            BetaSpectrum beta{1, 1};
        };

        explicit Decay(const Configuration& cfg);

        Interaction operator()(ToyTransport::rng_t& rng) const;

        /// Mean energy release (to beta, alpha and gammas), for checks.
        double mean_energy() const;

       private:
        std::vector<Branch> m_branches;
        std::vector<double> m_cumulative;
    };

    /** Events from a MARLEY JSON output file ("format": "json").

        Each event's final particles become primaries with their kinetic
        energy; the energy is that of the incident neutrino.
     */
    class MarleyEvents {
       public:
        explicit MarleyEvents(const std::string& filename);
        size_t size() const { return m_events.size(); }
        const Interaction& operator[](size_t ind) const { return m_events[ind]; }

        /// Parse the "events" array of a MARLEY JSON document.
        static std::vector<Interaction> parse(const Configuration& doc);

       private:
        std::vector<Interaction> m_events;
    };

}  // namespace WireCell::Trig

#endif
