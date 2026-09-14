#ifndef WIRECELLEDEP_IONIZATION
#define WIRECELLEDEP_IONIZATION

#include "WireCellIface/IPrimaryTracker.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"

#include <memory>
#include <mutex>
#include <string>

namespace EDepSim {
    class TrackingService;
}

namespace WireCell::Edep {

    class PrimaryVertexGenerator;

    /** Ionization: IPrimaryVertexSet -> ITrackSegmentSet.
     *
     *  Runs edep-sim (Geant4) on the input primaries via a thread-affine
     *  EDepSim::TrackingService and summarizes the resulting TG4HitSegments as
     *  an ITrackSegmentSet (which the existing Gen::TrackSegmentSampler turns
     *  into depos).  The primaries are delivered to Geant4 by our
     *  PrimaryVertexGenerator, installed as the service's custom generator.
     *
     *  Configuration:
     *   - gdml         (string): detector GDML file.
     *   - physics_list (string): Geant4 reference physics list (default: edep-sim's).
     *   - macro        (string): inline Geant4 macro text (physics tunings).
     *   - w_quanta     (double): energy per ionization/scintillation quantum
     *                            used to split deposited energy into electrons
     *                            (default 19.5 eV, edep-sim's value).
     *
     *  Only one G4RunManager may exist per process, so use a single Ionization
     *  (and drive it serially).
     */
    class Ionization : public Aux::Logger, public IPrimaryTracker, public IConfigurable {
       public:
        Ionization();
        virtual ~Ionization();

        // IPrimaryTracker (IFunctionNode<IPrimaryVertexSet, ITrackSegmentSet>)
        virtual bool operator()(const input_pointer& in, output_pointer& out);

        // IConfigurable
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;

       private:
        void start();  // create + initialize the service, once

        std::string m_physics_list;
        std::string m_gdml;
        std::string m_macro;
        double m_w_quanta;

        std::unique_ptr<EDepSim::TrackingService> m_service;
        PrimaryVertexGenerator* m_gen{nullptr};  // borrowed; owned by the service worker
        std::once_flag m_started;
        std::mutex m_feed;  // pair feed(pvs)+simulate() atomically
        std::size_t m_count{0};
    };

}  // namespace WireCell::Edep

#endif
