/** Make depo sets of low-energy interactions: supernova neutrino events
    and radioactive decays, for trigger studies.

    Each call makes one depo set for one "frame": a readout window
    [start_time, start_time + readout_time) preceded by pre_time so that
    charge drifting in from the whole volume is in steady state over the
    readout.  Every frame covers the same window: frames are independent.

    Interactions of each source are Poisson in number for the source rate
    over [start_time - pre_time, start_time + readout_time), uniform in that
    time and uniform in the sensitive volume of the anode faces.  Their
    primaries are transported by ToyTransport and each deposit becomes a depo
    of charge given by the recombination model.

    Truth is carried in the depos:

    - id :: 64 * n + s, for the n-th interaction of the job (over all
      sources) and the source index s (< 64).
    - pdg, energy :: of the particle that deposited and the energy deposited.
    - One marker depo per interaction, of zero charge and pdg 0, at the
      vertex, with the interaction energy (the neutrino energy, or the decay
      energy release).  Zero charge depos are dropped by the drifter.

    Configuration:

    - anode :: IAnodePlane type:name, whose faces' sensitive volumes are
      filled.
    - recombination :: IRecombinationModel type:name.
    - seed :: random seed (default 0).
    - nframes :: number of depo sets before end of stream (default 1).
    - start_time, readout_time, pre_time :: the window, in WCT units.
    - transport :: ToyTransport configuration: density, radiation_length,
      step, electron_cut, photon_cut, brems_cut, multiple_scattering.
    - sources :: array of objects with
      - name :: a label.
      - kind :: "marley" or "decay".
      - signal :: whether it is signal (default true for marley, false for
        decay).  Recorded for downstream truth.
      - rate_hz :: interactions per second in the volume, or
      - activity_bq_per_kg :: decays per second per kg of argon (decay only).
      - file :: the MARLEY JSON output file (marley only).  Events are used
        in order, wrapping around with a warning.
      - branches :: the decay branches (decay only), see Kinematics.h.
      - daughter :: optional {name, mean_life}: each decay of this source is
        followed at the same place, after an exponential delay, by a decay
        of the named source (eg 214Bi -> 214Po).  Give the daughter
        rate_hz 0 if it has no independent decays.  A daughter decay after
        the window is dropped.
 */

#ifndef WIRECELLTRIG_LOWEDEPOSOURCE
#define WIRECELLTRIG_LOWEDEPOSOURCE

#include "WireCellAux/Logger.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IDepoSetSource.h"
#include "WireCellIface/IRecombinationModel.h"
#include "WireCellTrig/Kinematics.h"
#include "WireCellUtil/BoundingBox.h"

#include <memory>

namespace WireCell::Trig {

    class LowEDepoSource : public Aux::Logger, public IDepoSetSource, public IConfigurable {
       public:
        LowEDepoSource();
        virtual ~LowEDepoSource();

        virtual bool operator()(IDepoSet::pointer& out);

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        /// The argon mass of the sensitive volume.
        double mass() const { return m_mass; }

        static constexpr int max_sources = 64;

       private:
        struct Source {
            std::string name;
            bool signal{false};
            double rate{0};  // per WCT time
            std::unique_ptr<MarleyEvents> marley;
            size_t next_event{0};
            std::unique_ptr<Decay> decay;
            std::string daughter_name;
            int daughter{-1};
            double daughter_mean_life{0};
        };
        std::vector<Source> m_sources;

        IAnodePlane::pointer m_anode;
        IRecombinationModel::pointer m_recomb;
        std::unique_ptr<ToyTransport> m_transport;
        std::vector<BoundingBox> m_volumes;
        std::vector<double> m_volume_cdf;
        double m_mass{0};

        ToyTransport::rng_t m_rng;
        int m_nframes{1};
        double m_start_time{0}, m_readout_time{0}, m_pre_time{0};
        int m_count{0};
        long m_ninteractions{0};
        bool m_eos_sent{false};

        Point random_position();
        Interaction next_interaction(Source& src);
        double add_interaction(size_t isrc, double time, const Point& vertex, IDepo::vector& depos);
    };

}  // namespace WireCell::Trig

#endif
