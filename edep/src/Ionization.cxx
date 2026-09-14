#include "WireCellEdep/Ionization.h"
#include "WireCellEdep/PrimaryVertexGenerator.h"

#include "WireCellAux/SimpleTrackSegment.h"
#include "WireCellAux/SimpleTrackSegmentSet.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/Units.h"

#include "EDepSimTrackingService.hh"
#include "TG4Event.h"

#include <sstream>
#include <unordered_map>
#include <vector>

WIRECELL_FACTORY(Ionization, WireCell::Edep::Ionization, WireCell::INamed, WireCell::IPrimaryTracker,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Edep;

namespace {
    // Break inline macro text into one command per line; the TrackingService
    // applies each and skips blank/'#' lines.
    std::vector<std::string> split_lines(const std::string& text)
    {
        std::vector<std::string> lines;
        std::istringstream iss(text);
        std::string line;
        while (std::getline(iss, line)) {
            lines.push_back(line);
        }
        return lines;
    }
}  // namespace

Ionization::Ionization()
  : Aux::Logger("Ionization", "edep")
  , m_w_quanta(19.5 * units::eV)
{
}

Ionization::~Ionization() {}

WireCell::Configuration Ionization::default_configuration() const
{
    Configuration cfg;
    cfg["gdml"] = m_gdml;
    cfg["physics_list"] = m_physics_list;
    cfg["macro"] = m_macro;
    cfg["w_quanta"] = m_w_quanta;
    return cfg;
}

void Ionization::configure(const WireCell::Configuration& cfg)
{
    m_gdml = get<std::string>(cfg, "gdml", m_gdml);
    m_physics_list = get<std::string>(cfg, "physics_list", m_physics_list);
    m_macro = get<std::string>(cfg, "macro", m_macro);
    m_w_quanta = get<double>(cfg, "w_quanta", m_w_quanta);
}

void Ionization::start()
{
    std::call_once(m_started, [this] {
        m_service = std::make_unique<EDepSim::TrackingService>();

        // The factory runs ON the Geant4 thread (affinity); it publishes the
        // generator so operator() can feed it.  initialize() blocks until the
        // factory has run, so m_gen is set on return.
        auto holder = std::make_shared<PrimaryVertexGenerator*>(nullptr);
        m_service->initialize(
            m_physics_list, m_gdml,
            EDepSim::GeneratorFactory([holder] {
                auto* g = new PrimaryVertexGenerator();
                *holder = g;
                return g;
            }),
            split_lines(m_macro));
        m_gen = *holder;
        log->debug("initialized edep-sim tracking service (gdml={})", m_gdml);
    });
}

bool Ionization::operator()(const input_pointer& in, output_pointer& out)
{
    out = nullptr;
    if (!in) {  // end of stream
        log->debug("EOS at call {}", m_count);
        return true;
    }

    start();

    // Track one event on the dedicated Geant4 thread; the id also seeds the RNG.
    std::shared_ptr<TG4Event> event;
    {
        std::lock_guard<std::mutex> lk(m_feed);
        m_gen->feed(in);
        event = m_service->simulate(static_cast<unsigned long>(in->ident()));
    }

    // TrackId -> PDG, so a hit segment can be labeled with its primary's PDG.
    std::unordered_map<int, int> trackid_pdg;
    for (const auto& traj : event->Trajectories) {
        trackid_pdg[traj.GetTrackId()] = traj.GetPDGCode();
    }

    // Convert every TG4HitSegment (across all sensitive detectors) to an
    // ITrackSegment.  Field mapping mirrors the existing Arrow chain.  WCT and
    // Geant4 share CLHEP units, so no conversion (ddm-f4q.6).
    ITrackSegment::vector segments;
    for (const auto& [sdname, hits] : event->SegmentDetectors) {
        (void) sdname;
        for (const auto& seg : hits) {
            const double edep = seg.GetEnergyDeposit();
            const double sdep = seg.GetSecondaryDeposit();
            // edep-sim quanta model: N_q = E/W; N_ph = N_q * Esec/E; N_e = N_q - N_ph.
            const double n_q = edep / m_w_quanta;
            const double n_ph = (edep > 0.0) ? n_q * sdep / edep : 0.0;
            const double n_e = n_q - n_ph;

            const TLorentzVector& p0 = seg.GetStart();
            const TLorentzVector& p1 = seg.GetStop();

            int pdg = 0;
            auto it = trackid_pdg.find(seg.GetPrimaryId());
            if (it != trackid_pdg.end()) {
                pdg = it->second;
            }

            segments.push_back(std::make_shared<Aux::SimpleTrackSegment>(
                Point(p0.X(), p0.Y(), p0.Z()), Point(p1.X(), p1.Y(), p1.Z()), p0.T(), p1.T(), edep,
                sdep, n_e, seg.GetTrackLength(), seg.GetPrimaryId(), pdg));
        }
    }

    out = std::make_shared<Aux::SimpleTrackSegmentSet>(in->ident(), segments);
    log->debug("event {}: {} segments from {} detector(s)", in->ident(), segments.size(),
               event->SegmentDetectors.size());
    ++m_count;
    return true;
}
