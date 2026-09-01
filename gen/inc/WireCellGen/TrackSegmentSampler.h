// Sample track segments into point depositions, applying an ionization model
// to turn deposited energy into a (negative) count of ionization electrons.

#ifndef WIRECELLGEN_TRACKSEGMENTSAMPLER
#define WIRECELLGEN_TRACKSEGMENTSAMPLER

#include "WireCellIface/ITrackSegmentSampler.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"

#include <functional>
#include <string>

namespace WireCell::Gen {

    /** TrackSegmentSampler converts each segment of an input
     *  ITrackSegmentSet into point depos placed uniformly along the
     *  start->stop line (energy is deposited uniformly along a segment)
     *  and emits them as one IDepoSet.
     *
     *  The ionization model turning a segment into a total electron count
     *  is selected ONCE at configure() time and held as a function object
     *  (no per-call mode dispatch).  Configuration:
     *
     *  - ionization: "quanta" (default) or "recombination".
     *
     *    "quanta" (the default) takes the electron count directly from the
     *    producer's per-segment n_electrons() and applies NO recombination
     *    model (eg the edep-sim 19.5 eV quanta split, already fluctuated -- it
     *    is not re-fluctuated here).  This keeps the emitted electrons
     *    consistent with the producer's photon/charge partition (quanta
     *    conservation).  A segment whose n_electrons() is a clearly negative
     *    "not provided" sentinel is an error in this mode; tiny sub-count
     *    negatives (floating-point noise) are clamped to zero.
     *
     *    "recombination" instead recomputes the electron count from the
     *    segment energy and dE/dX = energy/track_length via the
     *    IRecombinationModel component named by the "recombination" parameter
     *    (default "BoxRecombination"; also "BirksRecombination",
     *    "MipRecombination").  The E-field, density and W-value are configured
     *    on that component.  Use this when the producer does not supply a
     *    trustworthy n_electrons().
     *
     *  - step_size: sample spacing along the segment (default 1 mm).  A
     *    segment of straight-line length L becomes N = max(1, ceil(L/step))
     *    depos at the centers of N equal sub-intervals, each carrying 1/N of
     *    the segment's charge and energy; times interpolate likewise.
     *
     *  Emitted depo charge follows the WCT convention that ionization
     *  electrons are NEGATIVE (as the recombination models' Wi = 23.6
     *  eV/-eplus yields); depo id/pdg come from the segment and extents are
     *  zero (drifting adds diffusion).
     */
    class TrackSegmentSampler : public Aux::Logger, public ITrackSegmentSampler, public IConfigurable {
       public:
        TrackSegmentSampler();
        virtual ~TrackSegmentSampler();

        // ITrackSegmentSampler
        virtual bool operator()(const input_pointer& in, output_pointer& out);

        /// WireCell::IConfigurable interface.
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;

       private:
        // Return the SIGNED total ionization charge (electrons negative)
        // for one segment, per the configured ionization mode.
        std::function<double(const ITrackSegment&)> m_ionize;

        double m_step_size;
        std::string m_ionization;
        std::string m_recombination;
        std::size_t m_count{0};
    };

}  // namespace WireCell::Gen

#endif
