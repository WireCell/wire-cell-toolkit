#include "WireCellTrig/ToyTransport.h"

#include <algorithm>
#include <cmath>

using namespace WireCell;
using namespace WireCell::Trig;

namespace {
    const double electron_mass = 0.51099895 * units::MeV;
    const double classical_electron_radius = 2.8179403e-13 * units::cm;
    const double avogadro = 6.02214076e23;
    const double argon_z_over_a = 18.0 / 39.948;  // per g/mol

    // Electron collision stopping power in argon, MeV cm2/g, approximating
    // ESTAR.  Interpolated in log-log.
    const std::vector<std::pair<double, double>> electron_stopping{
        {0.01, 15.9}, {0.02, 9.3},  {0.05, 4.70}, {0.1, 2.94}, {0.2, 2.00}, {0.5, 1.47}, {1.0, 1.37},
        {2.0, 1.39},  {5.0, 1.46},  {10.0, 1.52}, {20.0, 1.58}, {50.0, 1.66}, {100.0, 1.71},
    };

    // Pair production attenuation in argon, cm2/g, approximating XCOM.
    // Interpolated linearly in log(E).
    const std::vector<std::pair<double, double>> pair_attenuation{
        {1.022, 0.0},   {1.5, 0.0002},  {2.0, 0.0010},  {3.0, 0.0025},  {5.0, 0.0050},
        {10.0, 0.0100}, {20.0, 0.0150}, {50.0, 0.0200}, {100.0, 0.0225},
    };

    // Photoelectric attenuation in argon above the K edge, cm2/g at 100 keV.
    const double photo_attenuation_100kev = 0.075;

    const int pdg_electron = 11, pdg_gamma = 22, pdg_neutron = 2112, pdg_proton = 2212;
    const int pdg_alpha = 1000020040;

    double uniform(ToyTransport::rng_t& rng)
    {
        return std::uniform_real_distribution<double>(0, 1)(rng);
    }
}  // namespace

Vector WireCell::Trig::rotate(const Vector& dir, double cos_theta, double phi)
{
    const Vector other = std::abs(dir.x()) < 0.9 ? Vector(1, 0, 0) : Vector(0, 1, 0);
    const Vector u = dir.cross(other).norm();
    const Vector v = dir.cross(u);
    const double sin_theta = std::sqrt(std::max(0.0, 1 - cos_theta * cos_theta));
    return (cos_theta * dir + sin_theta * std::cos(phi) * u + sin_theta * std::sin(phi) * v).norm();
}

Vector WireCell::Trig::isotropic(ToyTransport::rng_t& rng)
{
    const double cos_theta = 2 * uniform(rng) - 1;
    const double phi = 2 * M_PI * uniform(rng);
    const double sin_theta = std::sqrt(1 - cos_theta * cos_theta);
    return Vector(sin_theta * std::cos(phi), sin_theta * std::sin(phi), cos_theta);
}

ToyTransport::ToyTransport()
  : ToyTransport(Config{})
{
}

ToyTransport::ToyTransport(const Config& cfg)
  : m_cfg(cfg)
{
    const double rho_gcm3 = m_cfg.density / (units::g / units::cm3);
    m_electron_density = rho_gcm3 * avogadro * argon_z_over_a / units::cm3;
}

double ToyTransport::electron_dedx(double kinetic) const
{
    const double t = kinetic / units::MeV;
    const auto& tab = electron_stopping;
    double s = 0;
    if (t <= tab.front().first) {
        s = tab.front().second * std::pow(tab.front().first / t, 0.8);
    }
    else if (t >= tab.back().first) {
        s = tab.back().second;
    }
    else {
        auto hi = std::upper_bound(tab.begin(), tab.end(), t,
                                   [](double val, const auto& one) { return val < one.first; });
        auto lo = hi - 1;
        const double f = std::log(t / lo->first) / std::log(hi->first / lo->first);
        s = std::exp(std::log(lo->second) + f * std::log(hi->second / lo->second));
    }
    return s * units::MeV * units::cm2 / units::g * m_cfg.density;
}

double ToyTransport::compton_mu(double energy) const
{
    const double k = energy / electron_mass;
    const double l = std::log(1 + 2 * k);
    const double sigma = 2 * M_PI * classical_electron_radius * classical_electron_radius *
                         ((1 + k) / (k * k) * (2 * (1 + k) / (1 + 2 * k) - l / k) + l / (2 * k) -
                          (1 + 3 * k) / ((1 + 2 * k) * (1 + 2 * k)));
    return m_electron_density * sigma;
}

double ToyTransport::photo_mu(double energy) const
{
    const double ratio = 0.1 * units::MeV / energy;
    return photo_attenuation_100kev * ratio * ratio * ratio * units::cm2 / units::g * m_cfg.density;
}

double ToyTransport::pair_mu(double energy) const
{
    const double e = energy / units::MeV;
    const auto& tab = pair_attenuation;
    if (e <= tab.front().first) return 0;
    double a = tab.back().second;
    if (e < tab.back().first) {
        auto hi = std::upper_bound(tab.begin(), tab.end(), e,
                                   [](double val, const auto& one) { return val < one.first; });
        auto lo = hi - 1;
        const double f = std::log(e / lo->first) / std::log(hi->first / lo->first);
        a = lo->second + f * (hi->second - lo->second);
    }
    return a * units::cm2 / units::g * m_cfg.density;
}

double ToyTransport::total_mu(double energy) const
{
    return compton_mu(energy) + photo_mu(energy) + pair_mu(energy);
}

ToyTransport::Tally ToyTransport::operator()(rng_t& rng, int pdg, double kinetic, const Vector& dir,
                                             const Point& pos, std::vector<Deposit>& out) const
{
    Tally tally;
    std::vector<Particle> stack{{pdg, kinetic, dir.norm(), pos}};
    while (!stack.empty()) {
        const Particle p = stack.back();
        stack.pop_back();
        if (p.kinetic <= 0) continue;
        const int apdg = std::abs(p.pdg);
        if (apdg == pdg_electron) {
            electron(rng, p, stack, out, tally);
        }
        else if (p.pdg == pdg_gamma) {
            photon(rng, p, stack, out, tally);
        }
        else if (p.pdg == pdg_proton || (p.pdg > 1000000000 && p.pdg <= pdg_alpha)) {
            ion(p, out, tally);  // p, d, t, He3, alpha
        }
        else {
            tally.lost += p.kinetic;  // neutrons, neutrinos, heavier recoils
        }
    }
    return tally;
}

void ToyTransport::electron(rng_t& rng, const Particle& p, std::vector<Particle>& stack, std::vector<Deposit>& out,
                            Tally& tally) const
{
    const double X0 = m_cfg.radiation_length;
    const double kc = m_cfg.brems_cut;
    double T = p.kinetic;
    Vector dir = p.dir;
    Point pos = p.pos;
    std::normal_distribution<double> gaus(0, 1);

    while (T > m_cfg.electron_cut) {
        // stopping power at the middle of the step
        const double dedx0 = electron_dedx(T);
        const double dedx = electron_dedx(std::max(T - 0.5 * dedx0 * m_cfg.step, 0.5 * T));
        double step = m_cfg.step;
        double dE = dedx * step;
        if (dE >= T) {
            dE = T;
            step = T / dedx;
        }

        // Bremsstrahlung photons above kc, complete screening spectrum
        // (4/3 - 4/3 y + y^2)/k with y = k/T.
        if (T - dE > kc) {
            const double Tb = T - dE;
            const double r = kc / Tb;
            const double per_x0 = 4.0 / 3.0 * std::log(1 / r) - 4.0 / 3.0 * (1 - r) + 0.5 * (1 - r * r);
            const int nphot = std::poisson_distribution<int>(step / X0 * per_x0)(rng);
            for (int ip = 0; ip < nphot; ++ip) {
                double k = 0;
                while (true) {
                    k = kc * std::pow(Tb / kc, uniform(rng));
                    const double y = k / Tb;
                    if (uniform(rng) * 4.0 / 3.0 < 4.0 / 3.0 - 4.0 / 3.0 * y + y * y) break;
                }
                k = std::min(k, T - dE);
                if (k <= 0) break;
                T -= k;
                stack.push_back({pdg_gamma, k, dir, pos});
            }
        }

        out.push_back({pos + dir * (0.5 * step), dE, step, p.pdg});
        tally.deposited += dE;
        pos = pos + dir * step;
        T -= dE;

        if (m_cfg.multiple_scattering && T > m_cfg.electron_cut) {
            const double pc = std::sqrt(T * (T + 2 * electron_mass));
            const double beta = pc / (T + electron_mass);
            const double ratio = step / X0;
            const double corr = std::max(0.1, 1 + 0.038 * std::log(ratio / (beta * beta)));
            const double theta0 = 13.6 * units::MeV / (beta * pc) * std::sqrt(ratio) * corr;
            double cos_theta = 0, phi = 2 * M_PI * uniform(rng);
            if (theta0 > 1.0) {
                cos_theta = 2 * uniform(rng) - 1;
            }
            else {
                const double tx = theta0 * gaus(rng), ty = theta0 * gaus(rng);
                cos_theta = std::cos(std::min(M_PI, std::hypot(tx, ty)));
                phi = std::atan2(ty, tx);
            }
            dir = rotate(dir, cos_theta, phi);
        }
    }
    if (T > 0) {
        out.push_back({pos, T, std::max(T / electron_dedx(T), 1e-3 * units::mm), p.pdg});
        tally.deposited += T;
    }
    if (p.pdg == -pdg_electron) {  // annihilation at rest
        const Vector d = isotropic(rng);
        stack.push_back({pdg_gamma, electron_mass, d, pos});
        stack.push_back({pdg_gamma, electron_mass, -1.0 * d, pos});
    }
}

void ToyTransport::photon(rng_t& rng, const Particle& p, std::vector<Particle>& stack, std::vector<Deposit>& out,
                          Tally& tally) const
{
    double k = p.kinetic;
    Vector dir = p.dir;
    Point pos = p.pos;
    while (true) {
        if (k <= m_cfg.photon_cut) {
            out.push_back({pos, k, 1e-3 * units::mm, pdg_gamma});
            tally.deposited += k;
            return;
        }
        const double mu_c = compton_mu(k), mu_p = photo_mu(k), mu_e = pair_mu(k);
        const double mu = mu_c + mu_p + mu_e;
        pos = pos + dir * (-std::log(1 - uniform(rng)) / mu);
        const double which = uniform(rng) * mu;

        if (which < mu_p) {  // photoelectric: all to an electron
            stack.push_back({pdg_electron, k, dir, pos});
            return;
        }
        if (which < mu_p + mu_e) {  // pair: share what is above threshold
            const double avail = k - 2 * electron_mass;
            const double share = uniform(rng);
            const double cos_theta = std::cos(std::min(M_PI, electron_mass / k));
            stack.push_back({pdg_electron, share * avail, rotate(dir, cos_theta, 2 * M_PI * uniform(rng)), pos});
            stack.push_back({-pdg_electron, (1 - share) * avail, rotate(dir, cos_theta, 2 * M_PI * uniform(rng)), pos});
            return;
        }

        // Compton, sampled from Klein-Nishina (as Geant4 G4KleinNishinaCompton).
        const double kappa = k / electron_mass;
        const double eps0 = 1 / (1 + 2 * kappa);
        const double eps0sq = eps0 * eps0;
        const double alpha1 = -std::log(eps0);
        const double alpha2 = 0.5 * (1 - eps0sq);
        double eps = 1, one_cost = 0;
        while (true) {
            if (alpha1 > (alpha1 + alpha2) * uniform(rng)) {
                eps = std::exp(-alpha1 * uniform(rng));
            }
            else {
                eps = std::sqrt(eps0sq + (1 - eps0sq) * uniform(rng));
            }
            one_cost = (1 - eps) / (eps * kappa);
            const double sint2 = one_cost * (2 - one_cost);
            const double greject = 1 - eps * sint2 / (1 + eps * eps);
            if (greject >= uniform(rng)) break;
        }
        const double kprime = eps * k;
        const Vector newdir = rotate(dir, 1 - one_cost, 2 * M_PI * uniform(rng));
        const Vector pe = k * dir - kprime * newdir;  // electron momentum
        const double te = k - kprime;
        if (te > 0 && pe.magnitude() > 0) {
            stack.push_back({pdg_electron, te, pe.norm(), pos});
        }
        k = kprime;
        dir = newdir;
    }
}

void ToyTransport::ion(const Particle& p, std::vector<Deposit>& out, Tally& tally) const
{
    // Bragg-Kleeman range R = a T^1.75 for protons in LAr, about 0.35 mm at
    // 5 MeV.  For heavier ions, scale by mass over charge squared at the
    // same energy per nucleon.
    const double a_proton = 0.35 * units::mm / std::pow(5.0, 1.75);
    double mass_number = 1, charge = 1;
    if (p.pdg > 1000000000) {
        mass_number = (p.pdg / 10) % 1000;
        charge = (p.pdg / 10000) % 1000;
    }
    const double t_per_nucleon = p.kinetic / units::MeV / mass_number;
    const double range = a_proton * std::pow(t_per_nucleon, 1.75) * mass_number / (charge * charge);
    out.push_back({p.pos, p.kinetic, std::max(range, 1e-3 * units::mm), p.pdg});
    tally.deposited += p.kinetic;
}
