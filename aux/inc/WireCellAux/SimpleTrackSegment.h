#include "WireCellIface/ITrackSegment.h"

namespace WireCell::Aux {

    // A track segment that simply holds all the data it presents.
    class SimpleTrackSegment : public WireCell::ITrackSegment {
    public:
        // If track_length is negative it is taken as the start to
        // stop distance.
        SimpleTrackSegment(const WireCell::Point& start, const WireCell::Point& stop,
                           double start_time, double stop_time,
                           double energy,
                           double secondary = 0.0,
                           double n_electrons = -1.0,
                           double track_length = -1.0,
                           int id = 0, int pdg = 0);

        virtual const WireCell::Point& start() const;
        virtual const WireCell::Point& stop() const;
        virtual double start_time() const;
        virtual double stop_time() const;
        virtual double energy() const;
        virtual double secondary() const;
        virtual double n_electrons() const;
        virtual double track_length() const;
        virtual int id() const;
        virtual int pdg() const;

    private:
        // bag o' data
        WireCell::Point m_start, m_stop;
        double m_start_time, m_stop_time;
        double m_energy;
        double m_secondary;
        double m_n_electrons;
        double m_track_length;
        int m_id;
        int m_pdg;
    };

}  // namespace WireCell::Aux
