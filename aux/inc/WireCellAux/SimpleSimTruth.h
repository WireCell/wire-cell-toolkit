#ifndef WIRECELLAUX_SIMPLESIMTRUTH
#define WIRECELLAUX_SIMPLESIMTRUTH

#include "WireCellIface/ISimTruth.h"

namespace WireCell::Aux {

    // A sim-truth aggregate that simply holds shared pointers to its
    // constituent products (any of which may be null).
    class SimpleSimTruth : public WireCell::ISimTruth {
        int m_ident;
        WireCell::IPrimaryVertexSet::pointer m_primaries;
        WireCell::ITrajectorySet::pointer m_trajectories;
        WireCell::ITrackSegmentSet::pointer m_segments;
        WireCell::IPhotonHitSet::pointer m_photons;
        WireCell::Configuration m_metadata;

       public:
        SimpleSimTruth(int ident, WireCell::IPrimaryVertexSet::pointer primaries,
                       WireCell::ITrajectorySet::pointer trajectories,
                       WireCell::ITrackSegmentSet::pointer segments,
                       WireCell::IPhotonHitSet::pointer photons,
                       const WireCell::Configuration& metadata = WireCell::Configuration())
          : m_ident(ident)
          , m_primaries(primaries)
          , m_trajectories(trajectories)
          , m_segments(segments)
          , m_photons(photons)
          , m_metadata(metadata)
        {
        }
        virtual ~SimpleSimTruth();
        virtual int ident() const { return m_ident; }
        virtual WireCell::IPrimaryVertexSet::pointer primaries() const { return m_primaries; }
        virtual WireCell::ITrajectorySet::pointer trajectories() const { return m_trajectories; }
        virtual WireCell::ITrackSegmentSet::pointer segments() const { return m_segments; }
        virtual WireCell::IPhotonHitSet::pointer photons() const { return m_photons; }
        virtual WireCell::Configuration metadata() const { return m_metadata; }
    };

}  // namespace WireCell::Aux

#endif
