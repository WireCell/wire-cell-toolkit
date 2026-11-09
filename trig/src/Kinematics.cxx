#include "WireCellTrig/Kinematics.h"

#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Persist.h"

#include <algorithm>
#include <cmath>

using namespace WireCell;
using namespace WireCell::Trig;

namespace {
    const double electron_mass = 0.51099895 * units::MeV;
    const double fine_structure = 1 / 137.035999;
    const int pdg_electron = 11, pdg_gamma = 22, pdg_alpha = 1000020040;

    double uniform(ToyTransport::rng_t& rng)
    {
        return std::uniform_real_distribution<double>(0, 1)(rng);
    }
}  // namespace

BetaSpectrum::BetaSpectrum(double q, int z_daughter, Shape shape, size_t npoints)
  : m_q(q)
  , m_z(z_daughter)
  , m_shape(shape)
{
    if (q <= 0 || npoints < 2) {
        THROW(ValueError() << errmsg{"BetaSpectrum: need a positive Q and at least two points"});
    }
    m_t.resize(npoints);
    m_cdf.resize(npoints);
    double sum = 0, tsum = 0;
    for (size_t ind = 0; ind < npoints; ++ind) {
        m_t[ind] = q * ind / (npoints - 1);
        if (ind) {
            const double tmid = 0.5 * (m_t[ind] + m_t[ind - 1]);
            const double w = density(tmid);
            sum += w;
            tsum += w * tmid;
        }
        m_cdf[ind] = sum;
    }
    for (auto& c : m_cdf) c /= sum;
    m_mean = tsum / sum;
}

double BetaSpectrum::density(double kinetic) const
{
    if (kinetic <= 0 || kinetic >= m_q) return 0;
    const double e = kinetic + electron_mass;
    const double p = std::sqrt(kinetic * (kinetic + 2 * electron_mass));
    const double q = m_q - kinetic;
    const double eta = fine_structure * m_z * e / p;
    const double fermi = 2 * M_PI * eta / (1 - std::exp(-2 * M_PI * eta));
    const double shape = m_shape == Shape::unique_first ? (p * p + q * q) : 1.0;
    return p * e * q * q * fermi * shape;
}

double BetaSpectrum::operator()(ToyTransport::rng_t& rng) const
{
    const double u = uniform(rng);
    auto it = std::upper_bound(m_cdf.begin(), m_cdf.end(), u);
    if (it == m_cdf.begin()) return m_t.front();
    if (it == m_cdf.end()) return m_t.back();
    const size_t hi = it - m_cdf.begin();
    const double c0 = m_cdf[hi - 1], c1 = m_cdf[hi];
    const double f = c1 > c0 ? (u - c0) / (c1 - c0) : 0;
    return m_t[hi - 1] + f * (m_t[hi] - m_t[hi - 1]);
}

Decay::Decay(const Configuration& cfg)
{
    double total = 0;
    for (const auto& jb : cfg["branches"]) {
        Branch br;
        br.probability = get(jb, "probability", 1.0);
        br.alpha = get(jb, "alpha", 0.0);
        for (const auto& g : jb["gammas"]) br.gammas.push_back(g.asDouble());
        const auto& jbeta = jb["beta"];
        if (jbeta.isObject()) {
            const std::string shape = get<std::string>(jbeta, "shape", "allowed");
            if (shape != "allowed" && shape != "unique_first") {
                THROW(ValueError() << errmsg{"Decay: unknown beta shape: " + shape});
            }
            br.has_beta = true;
            br.beta = BetaSpectrum(jbeta["q"].asDouble(), jbeta["z"].asInt(),
                                   shape == "allowed" ? BetaSpectrum::Shape::allowed
                                                      : BetaSpectrum::Shape::unique_first);
        }
        if (br.probability < 0) {
            THROW(ValueError() << errmsg{"Decay: negative branch probability"});
        }
        total += br.probability;
        m_branches.push_back(br);
        m_cumulative.push_back(total);
    }
    if (m_branches.empty() || total <= 0) {
        THROW(ValueError() << errmsg{"Decay: no branches"});
    }
    for (auto& c : m_cumulative) c /= total;
}

Interaction Decay::operator()(ToyTransport::rng_t& rng) const
{
    const double u = uniform(rng);
    size_t ind = std::upper_bound(m_cumulative.begin(), m_cumulative.end(), u) - m_cumulative.begin();
    ind = std::min(ind, m_branches.size() - 1);
    const auto& br = m_branches[ind];

    Interaction inter;
    if (br.has_beta) {
        const double t = br.beta(rng);
        inter.primaries.push_back({pdg_electron, t, isotropic(rng)});
        inter.energy += t;
    }
    if (br.alpha > 0) {
        inter.primaries.push_back({pdg_alpha, br.alpha, isotropic(rng)});
        inter.energy += br.alpha;
    }
    for (double g : br.gammas) {
        inter.primaries.push_back({pdg_gamma, g, isotropic(rng)});
        inter.energy += g;
    }
    return inter;
}

double Decay::mean_energy() const
{
    double mean = 0, last = 0;
    for (size_t ind = 0; ind < m_branches.size(); ++ind) {
        const auto& br = m_branches[ind];
        const double prob = m_cumulative[ind] - last;
        last = m_cumulative[ind];
        double e = br.alpha;
        for (double g : br.gammas) e += g;
        if (br.has_beta) e += br.beta.mean();
        mean += prob * e;
    }
    return mean;
}

std::vector<Interaction> MarleyEvents::parse(const Configuration& doc)
{
    std::vector<Interaction> events;
    for (const auto& jev : doc["events"]) {
        Interaction inter;
        for (const auto& jp : jev["initial_particles"]) {
            const int pdg = jp["pdg"].asInt();
            if (std::abs(pdg) == 12 || std::abs(pdg) == 14 || std::abs(pdg) == 16) {
                inter.energy = jp["E"].asDouble() * units::MeV;
            }
        }
        for (const auto& jp : jev["final_particles"]) {
            Primary prim;
            prim.pdg = jp["pdg"].asInt();
            prim.kinetic = (jp["E"].asDouble() - jp["mass"].asDouble()) * units::MeV;
            const Vector mom(jp["px"].asDouble(), jp["py"].asDouble(), jp["pz"].asDouble());
            prim.dir = mom.magnitude() > 0 ? mom.norm() : Vector(0, 0, 1);
            inter.primaries.push_back(prim);
        }
        events.push_back(inter);
    }
    return events;
}

MarleyEvents::MarleyEvents(const std::string& filename)
{
    const auto doc = Persist::load(filename);
    m_events = parse(doc);
    if (m_events.empty()) {
        THROW(ValueError() << errmsg{"MarleyEvents: no events in " + filename});
    }
}
