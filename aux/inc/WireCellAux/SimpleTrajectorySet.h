#ifndef WIRECELLAUX_SIMPLETRAJECTORYSET
#define WIRECELLAUX_SIMPLETRAJECTORYSET

#include "WireCellIface/ITrajectorySet.h"

namespace WireCell::Aux {

    class SimpleTrajectorySet : public WireCell::ITrajectorySet {
        int m_ident;
        WireCell::ITrajectory::shared_vector m_trajectories;
        WireCell::Configuration m_metadata;

       public:
        SimpleTrajectorySet(int ident, const WireCell::ITrajectory::vector& trajectories,
                            const WireCell::Configuration& metadata = WireCell::Configuration())
          : m_ident(ident)
          , m_trajectories(
                std::make_shared<WireCell::ITrajectory::vector>(trajectories.begin(), trajectories.end()))
          , m_metadata(metadata)
        {
        }
        virtual ~SimpleTrajectorySet();
        virtual int ident() const { return m_ident; }
        virtual WireCell::ITrajectory::shared_vector trajectories() const { return m_trajectories; }
        virtual WireCell::Configuration metadata() const { return m_metadata; }
    };

}  // namespace WireCell::Aux

#endif
