#include "WireCellAux/SimpleTrackSegment.h"

using namespace WireCell::Aux;

SimpleTrackSegment::SimpleTrackSegment(const WireCell::Point& start, const WireCell::Point& stop, double start_time,
                                       double stop_time, double energy, double secondary, double n_electrons,
                                       double track_length, int id, int pdg)
  : m_start(start)
  , m_stop(stop)
  , m_start_time(start_time)
  , m_stop_time(stop_time)
  , m_energy(energy)
  , m_secondary(secondary)
  , m_n_electrons(n_electrons)
  , m_track_length(track_length < 0 ? (stop - start).magnitude() : track_length)
  , m_id(id)
  , m_pdg(pdg)
{
}

const WireCell::Point& SimpleTrackSegment::start() const { return m_start; }
const WireCell::Point& SimpleTrackSegment::stop() const { return m_stop; }
double SimpleTrackSegment::start_time() const { return m_start_time; }
double SimpleTrackSegment::stop_time() const { return m_stop_time; }
double SimpleTrackSegment::energy() const { return m_energy; }
double SimpleTrackSegment::secondary() const { return m_secondary; }
double SimpleTrackSegment::n_electrons() const { return m_n_electrons; }
double SimpleTrackSegment::track_length() const { return m_track_length; }
int SimpleTrackSegment::id() const { return m_id; }
int SimpleTrackSegment::pdg() const { return m_pdg; }
