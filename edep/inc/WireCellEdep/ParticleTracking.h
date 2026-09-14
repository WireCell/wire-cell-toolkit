#ifndef WIRECELLEDEP_PARTICLETRACKING
#define WIRECELLEDEP_PARTICLETRACKING

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

    /** ParticleTracking: IPrimaryVertexSet -> ISimTruth.
     *
     *  Runs edep-sim (Geant4) on the input primaries via a thread-affine
     *  EDepSim::TrackingService and returns the event's full simulation truth.
     *  Which constituents of the ISimTruth are populated is configurable; a
     *  disabled constituent is left null.
     *
     *  Configuration:
     *   - gdml         (string): detector GDML file.
     *   - physics_list (string): Geant4 reference physics list.
     *   - macro        (string): inline Geant4 macro text.
     *   - w_quanta     (double): energy per quantum for the segment
     *                            electron/photon split (default 19.5 eV).
     *   - primaries    (bool, default false): echo the input primaries.
     *   - trajectories (bool, default false): the particle trajectory tree.
     *   - segments     (bool, default true):  the ionization track segments.
     *   - photons      (bool, default false): the optical-photon hits.
     *
     *  Only one G4RunManager may exist per process, so use a single
     *  ParticleTracking and drive it serially.
     *
     *  TODO: a configurable trajectory filter to winnow the (potentially very
     *  many) trajectories as they are recorded.
     */
    class ParticleTracking : public Aux::Logger, public IPrimaryTracker, public IConfigurable {
       public:
        ParticleTracking();
        virtual ~ParticleTracking();

        // IPrimaryTracker (IFunctionNode<IPrimaryVertexSet, ISimTruth>)
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
        bool m_do_primaries{false};
        bool m_do_trajectories{false};
        bool m_do_segments{true};
        bool m_do_photons{false};

        std::unique_ptr<EDepSim::TrackingService> m_service;
        PrimaryVertexGenerator* m_gen{nullptr};  // borrowed; owned by the service worker
        std::once_flag m_started;
        std::mutex m_feed;
        std::size_t m_count{0};
    };

}  // namespace WireCell::Edep

#endif
