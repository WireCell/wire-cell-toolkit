#include "WireCellGen/ParticleGun.h"

#include "WireCellAux/ParticleInfo.h"
#include "WireCellAux/SimplePrimaryParticle.h"
#include "WireCellAux/SimplePrimaryVertex.h"
#include "WireCellAux/SimplePrimaryVertexSet.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Exceptions.h"

#include <cmath>

WIRECELL_FACTORY(ParticleGun, WireCell::Gen::ParticleGun, WireCell::INamed,
                 WireCell::IPrimaryVertexSetSource, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Gen;

namespace {
    constexpr double PI = 3.14159265358979323846;
}

// Spec parsing is done inline in configure() via a lambda so it can flip the
// object's m_needs_rng flag.

ParticleGun::ParticleGun()
  : Aux::Logger("ParticleGun", "gen")
{
    m_mass = Aux::ParticleInfo::pdg_to_mass(m_pdg);
}

ParticleGun::~ParticleGun() {}

WireCell::Configuration ParticleGun::default_configuration() const
{
    Configuration cfg;
    cfg["particle"] = "muon";
    cfg["energy"] = 1000.0 * units::MeV;
    cfg["energy_type"] = "kinetic";
    cfg["position"][0] = 0.0;
    cfg["position"][1] = 0.0;
    cfg["position"][2] = 0.0;
    cfg["direction"]["theta"] = 0.0;
    cfg["direction"]["phi"] = 0.0;
    cfg["direction"]["max_angle"] = 0.0;
    cfg["count"] = 1;
    cfg["event0"] = 1;
    cfg["rng"] = m_rng_tn;
    return cfg;
}

void ParticleGun::configure(const WireCell::Configuration& cfg)
{
    m_rng_tn = get<std::string>(cfg, "rng", m_rng_tn);
    m_needs_rng = false;

    auto parse_spec = [this](const Configuration& v) -> Spec {
        Spec s;
        if (v.isNumeric()) {
            s.lo = s.hi = v.asDouble();
        }
        else if (v.isArray() && v.size() == 2) {
            s.range = true;
            s.lo = v[0].asDouble();
            s.hi = v[1].asDouble();
        }
        else if (v.isObject() && v.isMember("min") && v.isMember("max")) {
            s.range = true;
            s.lo = v["min"].asDouble();
            s.hi = v["max"].asDouble();
        }
        if (s.range) m_needs_rng = true;
        return s;
    };

    // particle: PDG int or canonical name
    auto part = cfg["particle"];
    if (part.isIntegral()) {
        m_pdg = part.asInt();
    }
    else if (part.isString()) {
        m_pdg = Aux::ParticleInfo::name_to_pdg(part.asString());
    }
    if (m_pdg == 0) {
        THROW(ValueError() << errmsg{"ParticleGun: unknown particle"});
    }
    m_mass = Aux::ParticleInfo::pdg_to_mass(m_pdg);
    m_name = Aux::ParticleInfo::pdg_to_name(m_pdg);

    m_energy_type = get<std::string>(cfg, "energy_type", m_energy_type);
    m_energy = parse_spec(cfg["energy"]);

    auto pos = cfg["position"];
    for (int i = 0; i < 3; ++i) {
        m_pos[i] = parse_spec(pos.isArray() ? pos[i] : Configuration());
    }

    auto dir = cfg["direction"];
    if (dir.isMember("dir")) {
        m_dir0 = convert<Point>(dir["dir"]).norm();
    }
    else {
        const double th = get<double>(dir, "theta", 0.0);
        const double ph = get<double>(dir, "phi", 0.0);
        m_dir0 = Vector(std::sin(th) * std::cos(ph), std::sin(th) * std::sin(ph), std::cos(th));
    }
    m_max_angle = get<double>(dir, "max_angle", 0.0);
    if (m_max_angle > 0.0) m_needs_rng = true;

    m_number = get<int>(cfg, "count", m_number);
    m_event0 = get<int>(cfg, "event0", m_event0);

    if (m_needs_rng) {
        m_rng = Factory::find_tn<IRandom>(m_rng_tn);
        if (!m_rng) {
            THROW(KeyError() << errmsg{"ParticleGun: failed to get IRandom " + m_rng_tn});
        }
    }

    log->debug("particle={} ({}) mass={} MeV count={}", m_name, m_pdg, m_mass / units::MeV, m_number);
}

double ParticleGun::sample(const Spec& s) const
{
    if (!s.range) return s.lo;
    return m_rng->uniform(s.lo, s.hi);
}

WireCell::Vector ParticleGun::direction() const
{
    if (m_max_angle <= 0.0 || !m_rng) {
        return m_dir0;
    }
    // Uniform in solid angle within the cone: cos-theta uniform in [cos(max),1].
    const double ct = m_rng->uniform(std::cos(m_max_angle), 1.0);
    const double st = std::sqrt(std::max(0.0, 1.0 - ct * ct));
    const double ph = m_rng->uniform(0.0, 2.0 * PI);
    const Vector local(st * std::cos(ph), st * std::sin(ph), ct);

    // Rotate the local +z frame onto the central direction.
    const Vector w = m_dir0;
    const Vector a = (std::abs(w.x()) < 0.9) ? Vector(1, 0, 0) : Vector(0, 1, 0);
    const Vector u = a.cross(w).norm();
    const Vector v = w.cross(u);
    return u * local.x() + v * local.y() + w * local.z();
}

bool ParticleGun::operator()(output_pointer& out)
{
    if (m_eos) {
        return false;
    }
    if (m_emitted >= m_number) {
        out = nullptr;
        m_eos = true;
        return true;
    }

    const double eval = sample(m_energy);
    double total = 0.0, pmag = 0.0;
    if (m_energy_type == "momentum") {
        pmag = eval;
        total = std::sqrt(pmag * pmag + m_mass * m_mass);
    }
    else if (m_energy_type == "total") {
        total = eval;
        pmag = std::sqrt(std::max(0.0, total * total - m_mass * m_mass));
    }
    else {  // kinetic (default)
        total = eval + m_mass;
        pmag = std::sqrt(std::max(0.0, total * total - m_mass * m_mass));
    }

    const Vector mom = direction() * pmag;
    const Point pos(sample(m_pos[0]), sample(m_pos[1]), sample(m_pos[2]));

    const int ident = m_event0 + m_emitted;
    auto part = std::make_shared<Aux::SimplePrimaryParticle>(m_pdg, mom, total);
    IPrimaryParticle::vector parts{part};
    auto vtx = std::make_shared<Aux::SimplePrimaryVertex>(pos, 0.0, parts);
    IPrimaryVertex::vector vtxs{vtx};
    out = std::make_shared<Aux::SimplePrimaryVertexSet>(ident, vtxs);

    log->debug("event {}: {} |p|={} MeV", ident, m_name, pmag / units::MeV);
    ++m_emitted;
    return true;
}
