#include "WireCellTrig/LowEDepoSource.h"

#include "WireCellAux/SimpleDepo.h"
#include "WireCellAux/SimpleDepoSet.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"

#include <algorithm>
#include <cmath>
#include <limits>

WIRECELL_FACTORY(LowEDepoSource, WireCell::Trig::LowEDepoSource, WireCell::INamed, WireCell::IDepoSetSource,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Trig;

LowEDepoSource::LowEDepoSource()
  : Aux::Logger("LowEDepoSource", "trig")
  , m_readout_time(3 * units::ms)
  , m_pre_time(2.5 * units::ms)
{
}

LowEDepoSource::~LowEDepoSource() {}

WireCell::Configuration LowEDepoSource::default_configuration() const
{
    Configuration cfg;
    cfg["anode"] = "";
    cfg["recombination"] = "";
    cfg["seed"] = 0;
    cfg["nframes"] = m_nframes;
    cfg["start_time"] = m_start_time;
    cfg["readout_time"] = m_readout_time;
    cfg["pre_time"] = m_pre_time;
    ToyTransport::Config tc;
    cfg["transport"]["density"] = tc.density;
    cfg["transport"]["radiation_length"] = tc.radiation_length;
    cfg["transport"]["step"] = tc.step;
    cfg["transport"]["electron_cut"] = tc.electron_cut;
    cfg["transport"]["photon_cut"] = tc.photon_cut;
    cfg["transport"]["brems_cut"] = tc.brems_cut;
    cfg["transport"]["multiple_scattering"] = tc.multiple_scattering;
    cfg["sources"] = Json::arrayValue;
    return cfg;
}

void LowEDepoSource::configure(const WireCell::Configuration& cfg)
{
    m_anode = Factory::find_tn<IAnodePlane>(get<std::string>(cfg, "anode", ""));
    m_recomb = Factory::find_tn<IRecombinationModel>(get<std::string>(cfg, "recombination", ""));
    m_rng.seed(get<int>(cfg, "seed", 0));
    m_nframes = get(cfg, "nframes", m_nframes);
    m_start_time = get(cfg, "start_time", m_start_time);
    m_readout_time = get(cfg, "readout_time", m_readout_time);
    m_pre_time = get(cfg, "pre_time", m_pre_time);

    ToyTransport::Config tc;
    const auto& jtc = cfg["transport"];
    tc.density = get(jtc, "density", tc.density);
    tc.radiation_length = get(jtc, "radiation_length", tc.radiation_length);
    tc.step = get(jtc, "step", tc.step);
    tc.electron_cut = get(jtc, "electron_cut", tc.electron_cut);
    tc.photon_cut = get(jtc, "photon_cut", tc.photon_cut);
    tc.brems_cut = get(jtc, "brems_cut", tc.brems_cut);
    tc.multiple_scattering = get(jtc, "multiple_scattering", tc.multiple_scattering);
    m_transport = std::make_unique<ToyTransport>(tc);

    m_volumes.clear();
    m_volume_cdf.clear();
    double volume = 0;
    for (const auto& face : m_anode->faces()) {
        if (!face) continue;
        const auto bb = face->sensitive();
        if (bb.empty()) continue;
        const auto& ray = bb.bounds();
        const Vector diag = ray.second - ray.first;
        const double vol = std::abs(diag.x() * diag.y() * diag.z());
        if (vol <= 0) continue;
        m_volumes.push_back(bb);
        volume += vol;
        m_volume_cdf.push_back(volume);
    }
    if (m_volumes.empty()) {
        THROW(ValueError() << errmsg{"LowEDepoSource: the anode has no sensitive volume"});
    }
    for (auto& c : m_volume_cdf) c /= volume;
    m_mass = volume * tc.density;

    m_sources.clear();
    for (const auto& js : cfg["sources"]) {
        Source src;
        src.name = get<std::string>(js, "name", "");
        const std::string kind = get<std::string>(js, "kind", "");
        if (kind == "marley") {
            src.marley = std::make_unique<MarleyEvents>(get<std::string>(js, "file", ""));
            src.signal = get(js, "signal", true);
        }
        else if (kind == "decay") {
            src.decay = std::make_unique<Decay>(js);
            src.signal = get(js, "signal", false);
        }
        else {
            THROW(ValueError() << errmsg{"LowEDepoSource: unknown source kind \"" + kind + "\" for " + src.name});
        }
        if (js.isMember("rate_hz")) {
            src.rate = js["rate_hz"].asDouble() / units::s;
        }
        else if (js.isMember("activity_bq_per_kg") && src.decay) {
            src.rate = js["activity_bq_per_kg"].asDouble() / units::s * (m_mass / units::kg);
        }
        else {
            THROW(ValueError() << errmsg{"LowEDepoSource: source " + src.name + " needs rate_hz or activity_bq_per_kg"});
        }
        if (js["daughter"].isObject()) {
            src.daughter_name = get<std::string>(js["daughter"], "name", "");
            src.daughter_mean_life = get(js["daughter"], "mean_life", 0.0);
        }
        log->debug("source {}: {} {}, rate {} Hz", m_sources.size(), src.name, kind, src.rate * units::s);
        m_sources.push_back(std::move(src));
    }
    if (m_sources.size() > (size_t) max_sources) {
        THROW(ValueError() << errmsg{"LowEDepoSource: too many sources"});
    }
    for (auto& src : m_sources) {
        if (src.daughter_name.empty()) continue;
        for (size_t ind = 0; ind < m_sources.size(); ++ind) {
            if (m_sources[ind].name == src.daughter_name) src.daughter = ind;
        }
        if (src.daughter < 0 || !m_sources[src.daughter].decay) {
            THROW(ValueError() << errmsg{"LowEDepoSource: no decay source named " + src.daughter_name});
        }
        if (src.daughter_mean_life <= 0) {
            THROW(ValueError() << errmsg{"LowEDepoSource: daughter mean_life must be positive for " + src.name});
        }
    }
    log->debug("sensitive volume {} m3, mass {} t, {} sources, {} frames", volume / units::m3,
               m_mass / (1000 * units::kg), m_sources.size(), m_nframes);
    m_count = 0;
    m_eos_sent = false;
}

Point LowEDepoSource::random_position()
{
    std::uniform_real_distribution<double> uni(0, 1);
    const double u = uni(m_rng);
    size_t ind = std::upper_bound(m_volume_cdf.begin(), m_volume_cdf.end(), u) - m_volume_cdf.begin();
    ind = std::min(ind, m_volumes.size() - 1);
    const auto& ray = m_volumes[ind].bounds();
    const Vector diag = ray.second - ray.first;
    return ray.first + Vector(uni(m_rng) * diag.x(), uni(m_rng) * diag.y(), uni(m_rng) * diag.z());
}

Interaction LowEDepoSource::next_interaction(Source& src)
{
    if (src.decay) return (*src.decay)(m_rng);
    if (src.next_event >= src.marley->size()) {
        log->warn("source {}: reusing MARLEY events after {}", src.name, src.marley->size());
        src.next_event = 0;
    }
    auto inter = (*src.marley)[src.next_event++];
    // MARLEY gives a direction per its configuration; randomize the frame so
    // reuse does not repeat topologies.
    const Vector axis = isotropic(m_rng);
    const double angle = std::uniform_real_distribution<double>(0, 2 * M_PI)(m_rng);
    for (auto& prim : inter.primaries) {
        // Rodrigues rotation of the direction about axis by angle.
        const Vector& d = prim.dir;
        prim.dir = (d * std::cos(angle) + axis.cross(d) * std::sin(angle) +
                    axis * (axis.dot(d) * (1 - std::cos(angle))))
                       .norm();
    }
    return inter;
}

double LowEDepoSource::add_interaction(size_t isrc, double time, const Point& vertex, IDepo::vector& depos)
{
    if (m_ninteractions >= std::numeric_limits<int>::max() / max_sources) {
        THROW(ValueError() << errmsg{"LowEDepoSource: too many interactions for the depo id"});
    }
    const int id = (int) (m_ninteractions++ * max_sources + isrc);
    const auto inter = next_interaction(m_sources[isrc]);
    depos.push_back(std::make_shared<Aux::SimpleDepo>(time, vertex, 0.0, nullptr, 0, 0, id, 0, inter.energy));
    double qtot = 0;
    std::vector<Deposit> deposits;
    for (const auto& prim : inter.primaries) {
        deposits.clear();
        (*m_transport)(m_rng, prim.pdg, prim.kinetic, prim.dir, vertex, deposits);
        for (const auto& dep : deposits) {
            if (dep.energy <= 0) continue;
            const double charge = (*m_recomb)(dep.energy, dep.length);
            qtot += charge;
            depos.push_back(
                std::make_shared<Aux::SimpleDepo>(time, dep.pos, charge, nullptr, 0, 0, id, dep.pdg, dep.energy));
        }
    }
    return qtot;
}

bool LowEDepoSource::operator()(IDepoSet::pointer& out)
{
    out = nullptr;
    if (m_eos_sent) return false;
    if (m_count >= m_nframes) {
        m_eos_sent = true;
        log->debug("EOS after {} frames, {} interactions", m_count, m_ninteractions);
        return true;
    }

    const double t0 = m_start_time - m_pre_time;
    const double duration = m_pre_time + m_readout_time;
    std::uniform_real_distribution<double> uni(0, 1);

    IDepo::vector depos;
    std::vector<size_t> ninter(m_sources.size(), 0);
    double qtot = 0;
    for (size_t isrc = 0; isrc < m_sources.size(); ++isrc) {
        auto& src = m_sources[isrc];
        const size_t n = std::poisson_distribution<long>(src.rate * duration)(m_rng);
        ninter[isrc] = n;
        for (size_t ind = 0; ind < n; ++ind) {
            const double time = t0 + uni(m_rng) * duration;
            const Point vertex = random_position();
            qtot += add_interaction(isrc, time, vertex, depos);
            // a chain of daughters at the same place
            size_t parent = isrc;
            double tchain = time;
            for (int depth = 0; depth < max_sources && m_sources[parent].daughter >= 0; ++depth) {
                const auto& par = m_sources[parent];
                tchain += std::exponential_distribution<double>(1 / par.daughter_mean_life)(m_rng);
                if (tchain >= t0 + duration) break;
                parent = par.daughter;
                ++ninter[parent];
                qtot += add_interaction(parent, tchain, vertex, depos);
            }
        }
    }
    std::sort(depos.begin(), depos.end(), [](const auto& a, const auto& b) { return a->time() < b->time(); });

    std::string counts;
    for (size_t isrc = 0; isrc < m_sources.size(); ++isrc) {
        counts += " " + m_sources[isrc].name + "=" + std::to_string(ninter[isrc]);
    }
    log->debug("frame {}: {} depos, charge {} e, interactions:{}", m_count, depos.size(), qtot, counts);
    out = std::make_shared<Aux::SimpleDepoSet>(m_count++, depos);
    return true;
}
