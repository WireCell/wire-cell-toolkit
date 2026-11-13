#ifndef WIRECELLAUX_SIMPLETRAJECTORY
#define WIRECELLAUX_SIMPLETRAJECTORY

#include "WireCellIface/ITrajectory.h"

namespace WireCell::Aux {

    // A trajectory that simply holds all the data it presents.
    class SimpleTrajectory : public WireCell::ITrajectory {
        int m_id, m_parent, m_pdg;
        WireCell::Point m_start;
        double m_start_time;
        WireCell::Vector m_momentum;
        double m_energy;
        WireCell::Configuration m_metadata;

       public:
        SimpleTrajectory(int id, int parent, int pdg, const WireCell::Point& start, double start_time,
                         const WireCell::Vector& momentum, double energy,
                         const WireCell::Configuration& metadata = WireCell::Configuration())
          : m_id(id)
          , m_parent(parent)
          , m_pdg(pdg)
          , m_start(start)
          , m_start_time(start_time)
          , m_momentum(momentum)
          , m_energy(energy)
          , m_metadata(metadata)
        {
        }
        virtual ~SimpleTrajectory();
        virtual int id() const { return m_id; }
        virtual int parent() const { return m_parent; }
        virtual int pdg() const { return m_pdg; }
        virtual WireCell::Point start() const { return m_start; }
        virtual double start_time() const { return m_start_time; }
        virtual WireCell::Vector momentum() const { return m_momentum; }
        virtual double energy() const { return m_energy; }
        virtual WireCell::Configuration metadata() const { return m_metadata; }
    };

}  // namespace WireCell::Aux

#endif
