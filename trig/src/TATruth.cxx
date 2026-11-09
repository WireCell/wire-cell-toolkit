#include "WireCellTrig/TATruth.h"
#include "WireCellTrig/LowEDepoSource.h"

#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <unordered_map>

WIRECELL_FACTORY(TATruth, WireCell::Trig::TATruth, WireCell::INamed, WireCell::ITriggerActivityTruth,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Trig;

const std::vector<std::string>& TATruth::ta_columns()
{
    static const std::vector<std::string> cols{
        "frame",        "plane",       "face",           "time_start", "time_end",   "time_peak",
        "channel_start", "channel_end", "channel_peak",  "adc_integral", "adc_peak", "ninputs",
        "q_total",      "q_signal",    "n_contrib",      "dominant_id", "dominant_source", "dominant_q",
        "offset_est"};
    return cols;
}

const std::vector<std::string>& TATruth::interaction_columns()
{
    static const std::vector<std::string> cols{"frame", "id", "source", "signal", "time", "x", "y",
                                               "z", "energy", "edep", "q_collection", "t_collection",
                                               "ndepos"};
    return cols;
}

const std::vector<std::string>& TATruth::link_columns()
{
    static const std::vector<std::string> cols{"frame", "ta", "id", "q"};
    return cols;
}

TATruth::TATruth()
  : Aux::Logger("TATruth", "trig")
  , m_drift_speed(1.6 * units::mm / units::us)
  , m_period(500 * units::ns)
  , m_time_tolerance(5 * units::us)
{
}

TATruth::~TATruth() {}

WireCell::Configuration TATruth::default_configuration() const
{
    Configuration cfg;
    cfg["anode"] = "";
    cfg["drift_speed"] = m_drift_speed;
    cfg["lifetime"] = m_lifetime;
    cfg["period"] = m_period;
    for (double off : m_plane_offsets) cfg["plane_offsets"].append(off);
    cfg["time_tolerance"] = m_time_tolerance;
    cfg["channel_tolerance"] = m_channel_tolerance;
    cfg["sources"] = Json::arrayValue;
    cfg["metadata"] = Json::objectValue;
    cfg["density"] = m_density;
    return cfg;
}

void TATruth::configure(const WireCell::Configuration& cfg)
{
    m_anode = Factory::find_tn<IAnodePlane>(get<std::string>(cfg, "anode", ""));
    m_drift_speed = get(cfg, "drift_speed", m_drift_speed);
    m_lifetime = get(cfg, "lifetime", m_lifetime);
    m_period = get(cfg, "period", m_period);
    m_time_tolerance = get(cfg, "time_tolerance", m_time_tolerance);
    m_channel_tolerance = get(cfg, "channel_tolerance", m_channel_tolerance);
    if (cfg["plane_offsets"].isArray()) {
        if (cfg["plane_offsets"].size() != 3) {
            THROW(ValueError() << errmsg{"TATruth: plane_offsets needs three values"});
        }
        for (int ind = 0; ind < 3; ++ind) m_plane_offsets[ind] = cfg["plane_offsets"][ind].asDouble();
    }
    m_source_names.clear();
    m_source_signal.clear();
    for (const auto& js : cfg["sources"]) {
        m_source_names.push_back(get<std::string>(js, "name", ""));
        const bool def = get<std::string>(js, "kind", "") == "marley";
        m_source_signal.push_back(get(js, "signal", def));
    }
    m_metadata = cfg["metadata"];
    m_density = get(cfg, "density", m_density);
    m_mass = 0;
    for (const auto& face : m_anode->faces()) {
        if (!face || face->sensitive().empty()) continue;
        const auto& ray = face->sensitive().bounds();
        const Vector diag = ray.second - ray.first;
        m_mass += std::abs(diag.x() * diag.y() * diag.z()) * m_density;
    }
}

namespace {
    ITensor::pointer make_table(const std::string& name, const std::vector<std::string>& cols,
                                const std::vector<double>& data)
    {
        Configuration md;
        md["name"] = name;
        for (const auto& c : cols) md["columns"].append(c);
        const size_t ncols = cols.size();
        const size_t nrows = data.size() / ncols;
        return std::make_shared<Aux::SimpleTensor>(ITensor::shape_t{nrows, ncols},
                                                   data.empty() ? (const double*) nullptr : data.data(), md);
    }
}  // namespace

bool TATruth::operator()(const input_tuple_type& intup, output_pointer& out)
{
    out = nullptr;
    const auto& depos = std::get<0>(intup);
    const auto& taset = std::get<1>(intup);
    if (!depos || !taset) {
        if (depos || taset) {
            log->warn("end of stream on one input only at call {}", m_count);
        }
        return true;
    }

    // Interactions from markers, and hits from depos.
    struct Inter {
        double time{0}, energy{0}, edep{0}, qcol{0}, qtcol{0};
        Point pos;
        size_t ndepos{0};
    };
    std::map<int, Inter> inters;
    hits_t hits;
    const auto faces = m_anode->faces();
    for (const auto& depo : *depos->depos()) {
        const int id = depo->id();
        if (depo->charge() == 0) {
            auto& in = inters[id];
            in.time = depo->time();
            in.pos = depo->pos();
            in.energy = depo->energy();
            continue;
        }
        const Point& pos = depo->pos();
        for (const auto& face : faces) {
            if (!face || !face->sensitive().inside(pos)) continue;
            const auto planes = face->planes();
            const double xresp = planes.front()->pimpos()->origin().x();
            const double dt = face->dirx() * (pos.x() - xresp) / m_drift_speed;
            double q = std::abs(depo->charge());
            if (m_lifetime > 0 && dt > 0) q *= std::exp(-dt / m_lifetime);
            const double tarr = depo->time() + dt;
            auto& in = inters[id];
            in.edep += depo->energy();
            ++in.ndepos;
            for (const auto& plane : planes) {
                const auto* pimpos = plane->pimpos();
                const int wind = pimpos->region_binning().bin(pimpos->distance(pos, 2));
                const auto& wires = plane->wires();
                if (wind < 0 || wind >= (int) wires.size()) continue;
                const int pind = plane->planeid().index();
                if (pind < 0 || pind > 2) continue;
                hits[wires[wind]->channel()].push_back({tarr + m_plane_offsets[pind], q, id});
                if (pind == 2) {
                    in.qcol += q;
                    in.qtcol += q * (tarr + m_plane_offsets[pind]);
                }
            }
            break;
        }
    }
    for (auto& [ch, chits] : hits) {
        std::sort(chits.begin(), chits.end(), [](const Hit& a, const Hit& b) { return a.time < b.time; });
    }

    const int frame = taset->ident();
    const size_t nsources = m_source_signal.size();
    auto is_signal = [&](int id) {
        const size_t src = id % LowEDepoSource::max_sources;
        return src < nsources && m_source_signal[src];
    };

    std::vector<double> ta_rows, link_rows;
    size_t ita = 0;
    for (const auto& ta : taset->activities()) {
        std::unordered_map<int, double> qby;
        std::set<std::pair<int, size_t>> used;
        double off_sum = 0, off_w = 0;

        auto match = [&](int channel, double tbeg, double tend) {
            const auto wpid = m_anode->resolve(channel);
            const int pind = wpid.index();
            const double offset = (pind >= 0 && pind <= 2) ? m_plane_offsets[pind] : 0;
            for (int dc = -m_channel_tolerance; dc <= m_channel_tolerance; ++dc) {
                auto it = hits.find(channel + dc);
                if (it == hits.end()) continue;
                const auto& chits = it->second;
                auto lo = std::lower_bound(chits.begin(), chits.end(), tbeg - m_time_tolerance,
                                           [](const Hit& h, double t) { return h.time < t; });
                for (; lo != chits.end() && lo->time <= tend + m_time_tolerance; ++lo) {
                    const size_t ind = lo - chits.begin();
                    if (!used.insert({channel + dc, ind}).second) continue;
                    qby[lo->id] += lo->charge;
                    off_sum += lo->charge * (tbeg - lo->time + offset);
                    off_w += lo->charge;
                }
            }
        };
        if (!ta.inputs.empty()) {
            for (const auto& tp : ta.inputs) {
                match(tp.channel, tp.time_start, tp.time_start + tp.samples_over_threshold * m_period);
            }
        }
        else {
            for (int ch = ta.channel_start; ch <= ta.channel_end; ++ch) match(ch, ta.time_start, ta.time_end);
        }

        double qtot = 0, qsig = 0, qdom = 0;
        int dom = -1;
        for (const auto& [id, q] : qby) {
            qtot += q;
            if (is_signal(id)) qsig += q;
            if (q > qdom) {
                qdom = q;
                dom = id;
            }
            link_rows.insert(link_rows.end(), {(double) frame, (double) ita, (double) id, q});
        }
        ta_rows.insert(ta_rows.end(),
                       {(double) frame, (double) ta.plane.index(), (double) ta.plane.face(), ta.time_start,
                        ta.time_end, ta.time_peak, (double) ta.channel_start, (double) ta.channel_end,
                        (double) ta.channel_peak, (double) ta.adc_integral, (double) ta.adc_peak,
                        (double) ta.ninputs, qtot, qsig, (double) qby.size(), (double) dom,
                        dom < 0 ? -1.0 : (double) (dom % LowEDepoSource::max_sources), qdom,
                        off_w > 0 ? off_sum / off_w : 0.0});
        ++ita;
    }

    std::vector<double> inter_rows;
    for (const auto& [id, in] : inters) {
        const int src = id % LowEDepoSource::max_sources;
        inter_rows.insert(inter_rows.end(),
                          {(double) frame, (double) id, (double) src, is_signal(id) ? 1.0 : 0.0, in.time,
                           in.pos.x(), in.pos.y(), in.pos.z(), in.energy, in.edep, in.qcol,
                           in.qcol > 0 ? in.qtcol / in.qcol : 0.0, (double) in.ndepos});
    }

    Configuration md;
    md["frame"] = frame;
    md["time_begin"] = taset->time_begin();
    md["time_end"] = taset->time_end();
    md["nchannels"] = (int) m_anode->channels().size();
    md["anode"] = m_anode->ident();
    md["drift_speed"] = m_drift_speed;
    md["mass_kg"] = m_mass / units::kg;
    for (const auto& plane : faces.front()->planes()) md["pitch"].append(plane->pimpos()->pitch());
    if (m_count == 0) {  // the channels of each plane, once
        std::vector<std::set<int>> pchans(3);
        for (const auto& face : faces) {
            if (!face) continue;
            for (const auto& plane : face->planes()) {
                const int pind = plane->planeid().index();
                if (pind < 0 || pind > 2) continue;
                for (const auto& ich : plane->channels()) pchans[pind].insert(ich->ident());
            }
        }
        for (int pind = 0; pind < 3; ++pind) {
            md["plane_channels"][pind] = Json::arrayValue;
            for (int ch : pchans[pind]) md["plane_channels"][pind].append(ch);
        }
    }
    md["period"] = m_period;
    for (size_t ind = 0; ind < nsources; ++ind) {
        Configuration js;
        js["name"] = m_source_names[ind];
        js["signal"] = (bool) m_source_signal[ind];
        md["sources"].append(js);
    }
    for (double off : m_plane_offsets) md["plane_offsets"].append(off);
    md["metadata"] = m_metadata;

    auto tensors = std::make_shared<ITensor::vector>();
    tensors->push_back(make_table("tas", ta_columns(), ta_rows));
    tensors->push_back(make_table("interactions", interaction_columns(), inter_rows));
    tensors->push_back(make_table("links", link_columns(), link_rows));
    log->debug("frame {}: {} TAs, {} interactions, {} links", frame, ita, inters.size(),
               link_rows.size() / link_columns().size());
    out = std::make_shared<Aux::SimpleTensorSet>(frame, md, tensors);
    ++m_count;
    return true;
}
