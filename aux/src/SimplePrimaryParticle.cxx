#include "WireCellAux/SimplePrimaryParticle.h"

using namespace WireCell::Aux;

SimplePrimaryParticle::SimplePrimaryParticle(int pdg, const WireCell::Vector& momentum, double energy, int id,
                                             const WireCell::Configuration& metadata)
  : m_pdg(pdg)
  , m_momentum(momentum)
  , m_energy(energy)
  , m_id(id)
  , m_metadata(metadata)
{
}

SimplePrimaryParticle::~SimplePrimaryParticle() {}

int SimplePrimaryParticle::pdg() const { return m_pdg; }
WireCell::Vector SimplePrimaryParticle::momentum() const { return m_momentum; }
double SimplePrimaryParticle::energy() const { return m_energy; }
int SimplePrimaryParticle::id() const { return m_id; }
WireCell::Configuration SimplePrimaryParticle::metadata() const { return m_metadata; }
