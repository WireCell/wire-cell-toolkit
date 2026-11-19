#ifndef WIRECELLAUX_SIMPLEPHOTONHIT
#define WIRECELLAUX_SIMPLEPHOTONHIT

#include "WireCellIface/IPhotonHit.h"

namespace WireCell::Aux {

    // A photon hit that simply holds all the data it presents.
    class SimplePhotonHit : public WireCell::IPhotonHit {
        WireCell::Point m_start, m_stop;
        double m_start_time, m_stop_time;
        double m_energy;
        int m_id;
        WireCell::Configuration m_metadata;

       public:
        SimplePhotonHit(const WireCell::Point& start, const WireCell::Point& stop, double start_time,
                        double stop_time, double energy, int id,
                        const WireCell::Configuration& metadata = WireCell::Configuration())
          : m_start(start)
          , m_stop(stop)
          , m_start_time(start_time)
          , m_stop_time(stop_time)
          , m_energy(energy)
          , m_id(id)
          , m_metadata(metadata)
        {
        }
        virtual ~SimplePhotonHit();
        virtual WireCell::Point start() const { return m_start; }
        virtual WireCell::Point stop() const { return m_stop; }
        virtual double start_time() const { return m_start_time; }
        virtual double stop_time() const { return m_stop_time; }
        virtual double energy() const { return m_energy; }
        virtual int id() const { return m_id; }
        virtual WireCell::Configuration metadata() const { return m_metadata; }
    };

}  // namespace WireCell::Aux

#endif
