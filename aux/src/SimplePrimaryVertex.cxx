#include "WireCellAux/SimplePrimaryVertex.h"

using namespace WireCell::Aux;

SimplePrimaryVertex::SimplePrimaryVertex(const WireCell::Point& position, double time,
                                         const WireCell::IPrimaryParticle::vector& particles,
                                         const WireCell::Configuration& metadata)
  : m_position(position)
  , m_time(time)
  , m_particles(std::make_shared<WireCell::IPrimaryParticle::vector>(particles.begin(), particles.end()))
  , m_metadata(metadata)
{
}

SimplePrimaryVertex::~SimplePrimaryVertex() {}

WireCell::Point SimplePrimaryVertex::position() const { return m_position; }
double SimplePrimaryVertex::time() const { return m_time; }
WireCell::IPrimaryParticle::shared_vector SimplePrimaryVertex::particles() const { return m_particles; }
WireCell::Configuration SimplePrimaryVertex::metadata() const { return m_metadata; }
