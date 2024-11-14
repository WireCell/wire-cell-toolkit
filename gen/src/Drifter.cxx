#include "WireCellGen/Drifter.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/String.h"

#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IAnodeFace.h"
#include "WireCellIface/IWirePlane.h"
#include "WireCellAux/SimpleDepo.h"

#include <boost/range.hpp>

#include <sstream>

WIRECELL_FACTORY(Drifter, WireCell::Gen::Drifter,
                 WireCell::INamed,
                 WireCell::IDrifter, WireCell::IConfigurable)

using namespace std;
using namespace WireCell;

bool Gen::Drifter::DepoTimeCompare::operator()(const IDepo::pointer& lhs, const IDepo::pointer& rhs) const
{
    if (lhs->time() == rhs->time()) {
        if (lhs->pos().x() == lhs->pos().x()) {
            return lhs.get() < rhs.get();  // break tie by pointer
        }
        return lhs->pos().x() < lhs->pos().x();
    }
    return lhs->time() < rhs->time();
}

// Xregion helper

Gen::Drifter::Xregion::Xregion(Configuration cfg)
    : anode(make_coordbounds(cfg["anode"]))
    , response(make_coordbounds(cfg["response"]))
    , cathode(make_coordbounds(cfg["cathode"]))
    , bulk(*response, *cathode)
    , near(*anode, *response)
{
}

Gen::Drifter::Drifter()
  : Aux::Logger("Drifter", "gen")
  , m_rng(nullptr)
  , m_rng_tn("Random")
  , m_DL(7.2 * units::centimeter2 / units::second)   // from arXiv:1508.07059v2
  , m_DT(12.0 * units::centimeter2 / units::second)  // ditto
  , m_lifetime(8 * units::ms)                        // read off RHS of figure 6 in MICROBOONE-NOTE-1003-PUB
  , m_fluctuate(true)
  , m_speed(1.6 * units::mm / units::us)
  , m_toffset(0.0)
  , n_dropped(0)
  , n_drifted(0)
{
}

Gen::Drifter::~Drifter() {}

WireCell::Configuration Gen::Drifter::default_configuration() const
{
    Configuration cfg;
    cfg["DL"] = m_DL;
    cfg["DT"] = m_DT;
    cfg["lifetime"] = m_lifetime;
    cfg["fluctuate"] = m_fluctuate;
    cfg["drift_speed"] = m_speed;
    cfg["time_offset"] = m_toffset;

    // see comments in .h file
    cfg["xregions"] = Json::arrayValue;
    return cfg;
}

void Gen::Drifter::configure(const WireCell::Configuration& cfg)
{
    reset();

    m_rng_tn = get(cfg, "rng", m_rng_tn);
    m_rng = Factory::find_tn<IRandom>(m_rng_tn);

    m_DL = get<double>(cfg, "DL", m_DL);
    m_DT = get<double>(cfg, "DT", m_DT);
    m_lifetime = get<double>(cfg, "lifetime", m_lifetime);
    m_fluctuate = get<bool>(cfg, "fluctuate", m_fluctuate);
    m_speed = get<double>(cfg, "drift_speed", m_speed);
    m_toffset = get<double>(cfg, "time_offset", m_toffset);

    auto jxregions = cfg["xregions"];
    if (jxregions.empty()) {
        log->critical("no xregions given so I can do nothing");
        THROW(ValueError() << errmsg{"no xregions given"});
    }
    for (auto jone : jxregions) {
        if (!jone.isNull()) {
          m_xregions.emplace_back(Xregion(jone));
        }
    }
    log->debug("time offset: {} ms, drift speed: {} mm/us", m_toffset / units::ms,
             m_speed / (units::mm / units::us));
}

void Gen::Drifter::reset() { m_xregions.clear(); }

bool Gen::Drifter::insert(const input_pointer& depo)
{
    // electrical charge to drift.  Electrons should be negative
    const double Qi = depo->charge();
    if (Qi == 0.0) {
        // Yes, some silly depo sources ask us to drift nothing....
        return false;
    }

    // Find which X region to add, or reject.  Maybe there is a faster
    // way to do this than a loop.  For example, the regions could be
    // examined in order to find some regular binning of the X axis
    // and then at worse only explicitly check extent partial bins
    // near to their edges.

    double respx = 0, direction = 0.0;

    auto xrit = std::find_if(m_xregions.begin(), m_xregions.end(), Gen::Drifter::IsInsideResp(depo));
    if (xrit != m_xregions.end()) {
        // Back up in space and time.  This is a best effort fudge.  See:
        // https://github.com/WireCell/wire-cell-gen/issues/22

        respx = xrit->response->location();
        direction = -1.0;
    }
    else {
        xrit = std::find_if(m_xregions.begin(), m_xregions.end(), Gen::Drifter::IsInsideBulk(depo));
        if (xrit != m_xregions.end()) {  // in bulk
            respx = xrit->response->location();
            direction = 1.0;
        }
    }
    if (xrit == m_xregions.end()) {
        return false;  // outside both regions
    }

    Point pos = depo->pos();
    const double dt = std::abs((respx - pos.x()) / m_speed);
    pos.x(respx);

    double dL = depo->extent_long();
    double dT = depo->extent_tran();

    double Qf = Qi;
    if (direction > 0) {
        // final number of electrons after drift if no fluctuation.
        const double absorbprob = 1 - exp(-dt / m_lifetime);

        double dQ = Qi * absorbprob;

        // How many electrons remain, with fluctuation.
        // Note: fano/recomb fluctuation should be done before the depo was first made.
        if (m_fluctuate) {
            double sign = 1.0;
            if (Qi < 0) {
                sign = -1.0;
            }
            dQ = sign * m_rng->binomial((int) std::abs(Qi), absorbprob);
        }
        Qf = Qi - dQ;

        dL = sqrt(2.0 * m_DL * dt + dL * dL);
        dT = sqrt(2.0 * m_DT * dt + dT * dT);
    }

    auto newdepo = make_shared<Aux::SimpleDepo>(depo->time() + direction * dt + m_toffset,
                                                pos, Qf, depo, dL, dT, depo->id());
    xrit->depos.insert(newdepo);
    return true;
}

bool by_time(const IDepo::pointer& lhs, const IDepo::pointer& rhs) { return lhs->time() < rhs->time(); }

// save all cached depos to the output queue sorted in time order
void Gen::Drifter::flush(output_queue& outq)
{
    for (auto& xr : m_xregions) {
        outq.insert(outq.end(), xr.depos.begin(), xr.depos.end());
        xr.depos.clear();
    }
    std::sort(outq.begin(), outq.end(), by_time);
    outq.push_back(nullptr);
}

void Gen::Drifter::flush_ripe(output_queue& outq, double now)
{
    // It might be faster to use a sorted set which would avoid an
    // exhaustive iteration of each depos stash.  Or not.
    for (auto& xr : m_xregions) {
        if (xr.depos.empty()) {
            continue;
        }

        Xregion::ordered_depos_t::const_iterator depo_beg = xr.depos.begin(), depo_end = xr.depos.end();
        Xregion::ordered_depos_t::iterator depoit = depo_beg;
        while (depoit != depo_end) {
            if ((*depoit)->time() < now) {
                ++depoit;
                continue;
            }
            break;
        }
        if (depoit == depo_beg) {
            continue;
        }

        outq.insert(outq.end(), depo_beg, depoit);
        xr.depos.erase(depo_beg, depoit);
    }
    if (outq.empty()) {
        return;
    }
    std::sort(outq.begin(), outq.end(), by_time);
}

// always returns true because by hook or crook we consume the input.
bool Gen::Drifter::operator()(const input_pointer& depo, output_queue& outq)
{
    if (m_speed <= 0.0) {       // should be in configure()?
        log->error("illegal drift speed: {}", m_speed);
        return false;
    }

    if (!depo) {  // no more inputs expected, EOS flush

        flush(outq);

        if (n_dropped) {
            log->debug("see EOS, "
                       "( dropped:{} + drifted:{} ) / total:{} depos "
                       "from stream, outside of all {} drift xregions",
                       n_dropped, n_drifted, n_dropped + n_drifted,
                       m_xregions.size());
        }
        n_drifted = n_dropped = 0;
        return true;
    }

    bool ok = insert(depo);
    if (!ok) {
        ++n_dropped;  // depo is outside our regions but
        return true;  // nothing changed, so just bail
    }
    ++n_drifted;

    // At this point, time has advanced and maybe some drifted repos
    // are ripe for removal.

    flush_ripe(outq, depo->time());

    return true;
}
