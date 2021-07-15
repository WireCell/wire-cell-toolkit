#include "WireCellGen/TrackDepos.h"
#include "WireCellIface/SimpleDepo.h"
#include "WireCellUtil/NamedFactory.h"

#include "WireCellUtil/Point.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Testing.h"
#include "WireCellUtil/Persist.h"

#include "WireCellUtil/nljs2jcpp.hpp" // remove when ditch JsonCPP
#include "WireCellGen/Cfg/TrackDepos/Nljs.hpp"
using nljs_t = WireCellGen::Cfg::TrackDepos::data_t;

#include <sstream>

WIRECELL_FACTORY(TrackDepos, WireCell::Gen::TrackDepos, WireCell::IDepoSource, WireCell::IConfigurable)

using namespace std;
using namespace WireCell;

Gen::TrackDepos::TrackDepos(double stepsize, double clight)
    : m_cfg{stepsize, clight}
    , m_count(0)
    , l(Log::logger("sim"))
{

}

Gen::TrackDepos::~TrackDepos()
{
}


Configuration Gen::TrackDepos::default_configuration() const
{
    nljs_t nljs = m_cfg;
    return nljs.get<Json::Value>();
}

// wart: translate from schema/codegen rep of a ray to hand-written
template <typename RAY>
WireCell::Ray ray2ray(const RAY& ray)
{
    return WireCell::Ray(
        WireCell::Point(ray.tail.x, ray.tail.y, ray.tail.z),
        WireCell::Point(ray.head.x, ray.head.y, ray.head.z));
}

void Gen::TrackDepos::configure(const Configuration& cfg)
{
    nljs_t nljs = cfg;
    m_cfg = nljs.get<config_t>();

    for (auto& track : m_cfg.tracks) {
        add_track(track.time, ray2ray(track.ray), track.charge);
    }

    apply_grouping();
}

void Gen::TrackDepos::apply_grouping()
{
    if (m_depos.empty() or m_cfg.group_time <= 0.0) {
        m_depos.push_back(nullptr);
        return;
    }
    std::deque<WireCell::IDepo::pointer> grouped;
    double now = m_depos.front()->time();
    double end = now + m_cfg.group_time;
    for (auto depo : m_depos) {
        if (depo->time() < end) {
            grouped.push_back(depo);
            continue;
        }
        grouped.push_back(nullptr);
        now = depo->time();
        end = now + m_cfg.group_time;
        grouped.push_back(depo);
    }
    grouped.push_back(nullptr);
    m_depos = grouped;
}

static std::string dump(IDepo::pointer d)
{
    std::stringstream ss;
    ss << "q=" << d->charge() / units::eplus << "eles, t=" << d->time() / units::us << "us, r=" << d->pos() / units::mm
       << "mm";
    return ss.str();
}

void Gen::TrackDepos::add_track(double time, WireCell::Ray ray, double charge)
{
    // store this only to facilitate unit tests
    m_tracks.emplace_back(time, ray, charge);

    l->debug("add_track({} us, ({} -> {})cm, {})", time / units::us, ray.first / units::cm, ray.second / units::cm,
             charge);

    const WireCell::Vector dir = WireCell::ray_unit(ray);
    const double length = WireCell::ray_length(ray);
    double step = 0;
    int count = 0;

    double charge_per_depo = units::eplus;  // charge of one positron
    if (charge > 0) {
        charge_per_depo = -charge / (length / m_cfg.step_size);
    }
    else if (charge <= 0) {
        charge_per_depo = charge;
    }

    while (step < length) {
        const double now = time + step / (m_cfg.clight * units::clight);
        const WireCell::Point here = ray.first + dir * step;
        SimpleDepo* sdepo = new SimpleDepo(now, here, charge_per_depo);
        m_depos.push_back(WireCell::IDepo::pointer(sdepo));
        step += m_cfg.step_size;
        ++count;
    }

    // earliest first
    std::sort(m_depos.begin(), m_depos.end(), ascending_time);

    l->debug("depos: {} over {}mm", m_depos.size(), length / units::mm);
}

bool Gen::TrackDepos::operator()(output_pointer& out)
{
    if (m_depos.empty()) {
        return false;
    }
    out = m_depos.front();
    m_depos.pop_front();

    if (!out) {  // chirp
        l->debug("EOS at call {}", m_count);
    }

    ++m_count;
    return true;
}

WireCell::IDepo::vector Gen::TrackDepos::depos() {
    return WireCell::IDepo::vector(m_depos.begin(), m_depos.end());
}


