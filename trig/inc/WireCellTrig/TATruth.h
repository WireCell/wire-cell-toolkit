/** Truth tables for the TAs of a frame, from the frame's depos.

    Inputs, synchronously per frame: the depos as made by the generator
    (before drifting, see LowEDepoSource for the truth they carry) and the
    TA set found in the frame.  Output: a tensor set of three tables.

    Depos are projected to the anode without diffusion: a depo in a face's
    sensitive volume arrives at the response plane after the drift at
    drift_speed, attenuated by the lifetime, and lands on the wire of each
    plane under it.  The time of a hit on a plane is the arrival time plus
    that plane's plane_offsets entry, which should be the delay from the
    response plane to the TP start.

    A TA collects the hits within channel_tolerance of each of its TPs'
    channel and within time_tolerance of the TP's time span (each hit once
    per TA).  TAs without stored TPs use their channel and time box.

    Tables (float64, one row each, "columns" in the tensor metadata):

    - tas :: per TA: frame, plane (index), face, its times, channels and
      ADC, ninputs, then q_total, q_signal (matched charge, all and from
      signal sources), n_contrib (interactions matched), dominant_id,
      dominant_source, dominant_q, and offset_est (the charge-weighted mean
      of TP start minus arrival time, to calibrate plane_offsets).
    - interactions :: per interaction (from its marker depo): frame, id,
      source, signal, time, x, y, z, energy, edep (energy of its depos in
      the sensitive volume), q_collection (arriving charge on the collection
      plane), t_collection (its mean time there, with the plane offset),
      ndepos.
    - links :: per (TA, interaction) with matched charge: frame, ta (row of
      the TA in this frame), id, q.

    The tensor set metadata has the frame ident and coverage, nchannels and
    ident of the anode, drift_speed, period, pitch (per plane of the first
    face), plane_channels (the channels of each plane, first output only),
    mass_kg (argon in the
    sensitive volume at "density"), the source names and signal flags, the
    plane offsets and "metadata" from the configuration.

    Configuration:

    - anode :: IAnodePlane type:name.
    - drift_speed, lifetime (0 = none), period (TP sample period).
    - plane_offsets :: [u, v, w] in WCT time.
    - time_tolerance, channel_tolerance.
    - sources :: as LowEDepoSource (name, kind, signal), indexed by id % 64.
    - metadata :: anything, passed to the output.
 */

#ifndef WIRECELLTRIG_TATRUTH
#define WIRECELLTRIG_TATRUTH

#include "WireCellAux/Logger.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/ITriggerActivityTruth.h"
#include "WireCellUtil/Units.h"

#include <map>

namespace WireCell::Trig {

    class TATruth : public Aux::Logger, public ITriggerActivityTruth, public IConfigurable {
       public:
        TATruth();
        virtual ~TATruth();

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        virtual bool operator()(const input_tuple_type& intup, output_pointer& out);

        /// Column names of the tables.
        static const std::vector<std::string>& ta_columns();
        static const std::vector<std::string>& interaction_columns();
        static const std::vector<std::string>& link_columns();

       private:
        IAnodePlane::pointer m_anode;
        double m_drift_speed, m_lifetime{0}, m_period;
        std::vector<double> m_plane_offsets{0, 0, 0};
        double m_time_tolerance;
        int m_channel_tolerance{1};
        std::vector<std::string> m_source_names;
        std::vector<bool> m_source_signal;
        Configuration m_metadata;
        double m_density{1.396 * units::g / units::cm3};
        double m_mass{0};
        int m_count{0};

        struct Hit {
            double time;
            double charge;
            int id;
        };
        using hits_t = std::map<int, std::vector<Hit>>;  // by channel
    };

}  // namespace WireCell::Trig

#endif
