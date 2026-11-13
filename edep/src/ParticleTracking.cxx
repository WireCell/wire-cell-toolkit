#include "WireCellEdep/ParticleTracking.h"
#include "WireCellEdep/PrimaryVertexGenerator.h"

#include "WireCellAux/SimpleTrackSegment.h"
#include "WireCellAux/SimpleTrackSegmentSet.h"
#include "WireCellAux/SimpleTrajectory.h"
#include "WireCellAux/SimpleTrajectorySet.h"
#include "WireCellAux/SimplePhotonHit.h"
#include "WireCellAux/SimplePhotonHitSet.h"
#include "WireCellAux/SimpleSimTruth.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/Units.h"

#include "EDepSimTrackingService.hh"
#include "TG4Event.h"

#include <sstream>
#include <unordered_map>
#include <vector>

WIRECELL_FACTORY(ParticleTracking, WireCell::Edep::ParticleTracking, WireCell::INamed,
                 WireCell::IPrimaryTracker, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Edep;

namespace {
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
    inline Point vect(const TLorentzVector& v) { return Point(v.X(), v.Y(), v.Z()); }
}  // namespace

ParticleTracking::ParticleTracking()
  : Aux::Logger("ParticleTracking", "edep")
  , m_w_quanta(19.5 * units::eV)
{
}

ParticleTracking::~ParticleTracking() {}

WireCell::Configuration ParticleTracking::default_configuration() const
{
    Configuration cfg;
    cfg["gdml"] = m_gdml;
    cfg["physics_list"] = m_physics_list;
    cfg["macro"] = m_macro;
    cfg["w_quanta"] = m_w_quanta;
    cfg["primaries"] = m_do_primaries;
    cfg["trajectories"] = m_do_trajectories;
    cfg["segments"] = m_do_segments;
    cfg["photons"] = m_do_photons;
    return cfg;
}

void ParticleTracking::configure(const WireCell::Configuration& cfg)
{
    m_gdml = get<std::string>(cfg, "gdml", m_gdml);
    m_physics_list = get<std::string>(cfg, "physics_list", m_physics_list);
    m_macro = get<std::string>(cfg, "macro", m_macro);
    m_w_quanta = get<double>(cfg, "w_quanta", m_w_quanta);
    m_do_primaries = get<bool>(cfg, "primaries", m_do_primaries);
    m_do_trajectories = get<bool>(cfg, "trajectories", m_do_trajectories);
    m_do_segments = get<bool>(cfg, "segments", m_do_segments);
    m_do_photons = get<bool>(cfg, "photons", m_do_photons);
}

void ParticleTracking::start()
{
    std::call_once(m_started, [this] {
        m_service = std::make_unique<EDepSim::TrackingService>();
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

bool ParticleTracking::operator()(const input_pointer& in, output_pointer& out)
{
    out = nullptr;
    if (!in) {  // end of stream
        log->debug("EOS at call {}", m_count);
        return true;
    }

    start();

    std::shared_ptr<TG4Event> event;
    {
        std::lock_guard<std::mutex> lk(m_feed);
        m_gen->feed(in);
        event = m_service->simulate(static_cast<unsigned long>(in->ident()));
    }

    const int ident = in->ident();

    // TrackId -> PDG (used to label segments and as the trajectory PDG).  Also
    // TrackId -> trajectory-start position (first point), for segment-free need.
    std::unordered_map<int, int> trackid_pdg;
    for (const auto& traj : event->Trajectories) {
        trackid_pdg[traj.GetTrackId()] = traj.GetPDGCode();
    }

    // --- segments ---
    ITrackSegmentSet::pointer segments;
    if (m_do_segments) {
        ITrackSegment::vector segs;
        for (const auto& [sdname, hits] : event->SegmentDetectors) {
            (void) sdname;
            for (const auto& seg : hits) {
                const double edep = seg.GetEnergyDeposit();
                const double sdep = seg.GetSecondaryDeposit();
                const double n_q = edep / m_w_quanta;
                const double n_ph = (edep > 0.0) ? n_q * sdep / edep : 0.0;
                const double n_e = n_q - n_ph;
                const TLorentzVector& p0 = seg.GetStart();
                const TLorentzVector& p1 = seg.GetStop();
                int pdg = 0;
                auto it = trackid_pdg.find(seg.GetPrimaryId());
                if (it != trackid_pdg.end()) pdg = it->second;
                segs.push_back(std::make_shared<Aux::SimpleTrackSegment>(
                    vect(p0), vect(p1), p0.T(), p1.T(), edep, sdep, n_e, seg.GetTrackLength(),
                    seg.GetPrimaryId(), pdg));
            }
        }
        segments = std::make_shared<Aux::SimpleTrackSegmentSet>(ident, segs);
    }

    // --- trajectory tree ---
    ITrajectorySet::pointer trajectories;
    if (m_do_trajectories) {
        ITrajectory::vector trajs;
        for (const auto& traj : event->Trajectories) {
            const TLorentzVector& mom = traj.GetInitialMomentum();
            Point start(0, 0, 0);
            double start_time = 0;
            if (!traj.Points.empty()) {
                const TLorentzVector& p = traj.Points.front().GetPosition();
                start = vect(p);
                start_time = p.T();
            }
            Configuration md;
            md["name"] = std::string(traj.GetName());
            trajs.push_back(std::make_shared<Aux::SimpleTrajectory>(
                traj.GetTrackId(), traj.GetParentId(), traj.GetPDGCode(), start, start_time,
                Vector(mom.X(), mom.Y(), mom.Z()), mom.E(), md));
        }
        trajectories = std::make_shared<Aux::SimpleTrajectorySet>(ident, trajs);
    }

    // --- photon hits ---
    IPhotonHitSet::pointer photons;
    if (m_do_photons) {
        IPhotonHit::vector hits;
        for (const auto& [sdname, phits] : event->PhotonDetectors) {
            for (const auto& hit : phits) {
                const TLorentzVector& p0 = hit.GetStart();
                const TLorentzVector& p1 = hit.GetStop();
                Configuration md;
                md["sensitive_detector"] = sdname;
                md["process"] = hit.GetProcess();
                md["wavelength"] = hit.GetWavelength();
                hits.push_back(std::make_shared<Aux::SimplePhotonHit>(
                    vect(p0), vect(p1), p0.T(), p1.T(), hit.GetEnergyDeposit(), hit.GetPrimaryId(),
                    md));
            }
        }
        photons = std::make_shared<Aux::SimplePhotonHitSet>(ident, hits);
    }

    // --- primaries (echo the input) ---
    IPrimaryVertexSet::pointer primaries;
    if (m_do_primaries) {
        primaries = in;
    }

    out = std::make_shared<Aux::SimpleSimTruth>(ident, primaries, trajectories, segments, photons);
    log->debug("event {}: truth [primaries={} trajectories={} segments={} photons={}]", ident,
               (bool) primaries, (bool) trajectories, (bool) segments, (bool) photons);
    ++m_count;
    return true;
}
