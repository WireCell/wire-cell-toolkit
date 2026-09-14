#ifndef WIRECELLEDEP_PRIMARYVERTEXGENERATOR
#define WIRECELLEDEP_PRIMARYVERTEXGENERATOR

#include "WireCellIface/IPrimaryVertexSet.h"

#include <G4VPrimaryGenerator.hh>

class G4Event;

namespace WireCell::Edep {

    /** A Geant4 primary generator driven by a WCT IPrimaryVertexSet.
     *
     *  It is installed as the EDepSim::TrackingService custom generator (via a
     *  GeneratorFactory), so it is CONSTRUCTED on the dedicated Geant4 thread.
     *  Before each simulate() the owner feeds the current IPrimaryVertexSet;
     *  GeneratePrimaryVertex() then converts its vertices/particles into
     *  G4PrimaryVertex/G4PrimaryParticle.
     *
     *  WCT and Geant4/edep-sim share the CLHEP system of units (mm = ns = MeV =
     *  1), so values pass through with no conversion (see decision ddm-f4q.6).
     */
    class PrimaryVertexGenerator : public G4VPrimaryGenerator {
       public:
        PrimaryVertexGenerator();
        ~PrimaryVertexGenerator() override;

        /// Set the primaries to be built by the next GeneratePrimaryVertex.
        /// Non-owning share: the set must outlive the subsequent beamOn.
        void feed(WireCell::IPrimaryVertexSet::pointer pvs);

        /// G4VPrimaryGenerator: build the G4 primaries for one event.
        void GeneratePrimaryVertex(G4Event* event) override;

       private:
        WireCell::IPrimaryVertexSet::pointer m_pvs;
    };

}  // namespace WireCell::Edep

#endif
