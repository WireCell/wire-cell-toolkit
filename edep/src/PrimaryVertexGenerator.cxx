#include "WireCellEdep/PrimaryVertexGenerator.h"

#include <G4Event.hh>
#include <G4PrimaryParticle.hh>
#include <G4PrimaryVertex.hh>

using namespace WireCell;
using namespace WireCell::Edep;

PrimaryVertexGenerator::PrimaryVertexGenerator() = default;
PrimaryVertexGenerator::~PrimaryVertexGenerator() = default;

void PrimaryVertexGenerator::feed(WireCell::IPrimaryVertexSet::pointer pvs) { m_pvs = pvs; }

void PrimaryVertexGenerator::GeneratePrimaryVertex(G4Event* event)
{
    if (!m_pvs) return;

    // WCT and Geant4 share the CLHEP system of units, so no scaling (ddm-f4q.6).
    for (const auto& ivtx : *m_pvs->vertices()) {
        const Point pos = ivtx->position();
        auto* g4vtx = new G4PrimaryVertex(pos.x(), pos.y(), pos.z(), ivtx->time());

        for (const auto& ipart : *ivtx->particles()) {
            const Vector p = ipart->momentum();
            g4vtx->SetPrimary(
                new G4PrimaryParticle(ipart->pdg(), p.x(), p.y(), p.z(), ipart->energy()));
        }

        event->AddPrimaryVertex(g4vtx);
    }
}
