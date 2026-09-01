#ifndef WIRECELL_ITRACKSEGMENT
#define WIRECELL_ITRACKSEGMENT

#include "WireCellIface/IData.h"
#include "WireCellUtil/Point.h"

namespace WireCell {

    /** An interface to information about a track segment.
     *
     *  A track segment is a straight-line, calorimetric summary of
     *  energy deposited along the path of a charged particle track.
     *  A producer typically aggregates many tracking (eg Geant4)
     *  steps into one segment and, depending on the producer's
     *  configuration, may fold in energy deposited by secondary (eg
     *  delta-ray) tracks.  The id() and pdg() thus identify a single,
     *  possibly heuristic, "most important" associated track and not
     *  necessarily the only contributor.
     *
     *  Values are expressed in the WCT system of units.
     */
    class ITrackSegment : public IData<ITrackSegment> {
       public:
        virtual ~ITrackSegment();

        /// The location of the start of the segment.
        virtual const Point& start() const = 0;

        /// The location of the end of the segment.
        virtual const Point& stop() const = 0;

        /// The time at which the track was at start().
        virtual double start_time() const = 0;

        /// The time at which the track was at stop().
        virtual double stop_time() const = 0;

        /// The total energy deposited over the segment.
        virtual double energy() const = 0;

        /// The part of energy() lost to secondary (eg scintillation)
        /// processes and thus not producing ionization, as determined
        /// by the producer's quenching/recombination model.  Zero if
        /// the producer does not provide the partition.
        virtual double secondary() const = 0;

        /// The number of ionization electrons produced over the
        /// segment as determined by the producer.  Negative if the
        /// producer does not provide the count.
        virtual double n_electrons() const = 0;

        /// The total charged particle path length contributing to
        /// the segment.  This includes the path of any secondary
        /// tracks folded into the segment and so may exceed the
        /// start() to stop() distance.
        virtual double track_length() const = 0;

        /// The track ID of the associated track.
        virtual int id() const = 0;

        /// The PDG code of the associated track.
        virtual int pdg() const = 0;
    };

}  // namespace WireCell

#endif
