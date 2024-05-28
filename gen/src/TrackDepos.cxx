#include "WireCellGen/TrackDepos.h"
#include "WireCellAux/SimpleDepo.h"
#include "WireCellUtil/NamedFactory.h"

#include "WireCellUtil/Point.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Testing.h"
#include "WireCellUtil/Persist.h"

#include <sstream>

WIRECELL_FACTORY(TrackDepos, WireCell::Gen::TrackDepos, WireCell::IDepoSource, WireCell::IConfigurable)

using namespace std;
using namespace WireCell;

Gen::TrackDepos::TrackDepos(double stepsize, double clight)
    : Aux::Logger("TrackDepos", "gen")
    , m_stepsize(stepsize)
    , m_clight(clight)
    , m_count(0)
{
}

Gen::TrackDepos::~TrackDepos() {}

Configuration Gen::TrackDepos::default_configuration() const
{
    Configuration cfg;
    cfg["step_size"] = 1.0 * units::mm;
    cfg["clight"] = 1.0;  // fraction of speed of light the track goes
    cfg["tracks"] = Json::arrayValue;
    cfg["group_time"] = -1;  // if positive then chunk the
                             // stream of output depos into
                             // groups delineated by EOS
                             // markers and such that each
                             // group spans a time period no
                             // more than group_time.
    return cfg;
}

void Gen::TrackDepos::configure(const Configuration& cfg)
{
    m_stepsize = get<double>(cfg, "step_size", m_stepsize);
    Assert(m_stepsize > 0);
    m_clight = get<double>(cfg, "clight", m_clight);
    for (auto track : cfg["tracks"]) {
        double time = get<double>(track, "time", 0.0);
        double charge = get<double>(track, "charge", -1.0);
        Ray ray = get<Ray>(track, "ray");
        add_track(time, ray, charge);
    }
    double gt = get<double>(cfg, "group_time", -1);
    if (m_depos.empty() or gt <= 0.0) {
        m_depos.push_back(nullptr);
        return;
    }
    std::deque<WireCell::IDepo::pointer> grouped;
    double now = m_depos.front()->time();
    double end = now + gt;
    for (auto depo : m_depos) {
        if (depo->time() < end) {
            grouped.push_back(depo);
            continue;
        }
        grouped.push_back(nullptr);
        now = depo->time();
        end = now + gt;
        grouped.push_back(depo);
    }
    grouped.push_back(nullptr);
    m_depos = grouped;
}

void Gen::TrackDepos::add_track(double time, const WireCell::Ray& ray, double charge)
{
    log->debug("add_track({} us, ({} -> {})cm, {})", time / units::us, ray.first / units::cm, ray.second / units::cm,
             charge);
    m_tracks.push_back(track_t(time, ray, charge));

    const WireCell::Vector dir = WireCell::ray_unit(ray);
    const double length = WireCell::ray_length(ray);
    double step = 0;
    int count = 0;

    double charge_per_depo = units::eplus;  // charge of one positron
    if (charge > 0) {
        charge_per_depo = -charge / (length / m_stepsize);
    }
    else if (charge <= 0) {
        charge_per_depo = charge;
    }

    while (step < length) {
        const double now = time + step / (m_clight * units::clight);
        const WireCell::Point here = ray.first + dir * step;
        Aux::SimpleDepo* sdepo = new Aux::SimpleDepo(now, here, charge_per_depo);
        m_depos.push_back(WireCell::IDepo::pointer(sdepo));
        step += m_stepsize;
        ++count;
    }

    // earliest first
    std::sort(m_depos.begin(), m_depos.end(), ascending_time);
    log->debug("depos: {} over {}mm", m_depos.size(), length / units::mm);
}

bool Gen::TrackDepos::operator()(output_pointer& out)
{
    if (m_depos.empty()) {
        return false;
    }
    out = m_depos.front();
    m_depos.pop_front();

    if (!out) {  // chirp
        log->debug("EOS at call {}", m_count);
    }

    ++m_count;
    return true;
}

WireCell::IDepo::vector Gen::TrackDepos::depos() { return WireCell::IDepo::vector(m_depos.begin(), m_depos.end()); }
