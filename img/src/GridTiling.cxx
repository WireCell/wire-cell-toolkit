#include "WireCellImg/GridTiling.h"

#include "WireCellUtil/RayTiling.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellAux/SimpleBlob.h"

WIRECELL_FACTORY(GridTiling, WireCell::Img::GridTiling,
                 WireCell::INamed,
                 WireCell::ITiling, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::RayGrid;

Img::GridTiling::GridTiling()
    : Aux::Logger("GridTiling", "img")
{
}

Img::GridTiling::~GridTiling() {}

void Img::GridTiling::configure(const WireCell::Configuration& cfg)
{
    m_anode = Factory::find_tn<IAnodePlane>(cfg["anode"].asString());
    m_face = m_anode->face(cfg["face"].asInt());
    m_nudge = get(cfg, "nudge", m_nudge);
    if (!m_face) {
        log->error("no such face: {}", cfg["face"]);
        THROW(ValueError() << errmsg{"APA lacks face"});
    }
    m_threshold = get(cfg, "threshold", m_threshold);
    log->debug("configured with anode={} face={}",
               m_anode->ident(), m_face->ident());
}

WireCell::Configuration Img::GridTiling::default_configuration() const
{
    Configuration cfg;
    cfg["anode"] = "";  // user must set
    cfg["face"] = 0;    // the face number to focus on
    // activity must be bigger than this much to consider.
    cfg["threshold"] = m_threshold;
    cfg["nudge"] = m_nudge;
    return cfg;
}

IBlobSet::pointer Img::GridTiling::make_empty(const input_pointer& slice)
{
    return std::make_shared<Aux::SimpleBlobSet>(slice->ident(), slice);
}

bool Img::GridTiling::operator()(const input_pointer& slice, output_pointer& out)
{
    out = nullptr;
    if (!slice) {
        m_blobs_seen = 0;
        log->debug("EOS");
        return true;  // eos
    }

    const auto anodeid = m_anode->ident();
    const auto faceid = m_face->ident();

    auto chvs = slice->activity();

    if (chvs.empty()) {
        log->debug("anode={} face={} slice={}, time={} ms no activity",
                   anodeid, faceid, slice->ident(),
                   slice->frame()->time()/units::ms);
        out = make_empty(slice);
        return true;
    }

    const int sbs_ident = slice->ident();
    Aux::SimpleBlobSet* sbs = new Aux::SimpleBlobSet(sbs_ident, slice);
    out = IBlobSet::pointer(sbs);

    const int nbounds_layers = 2;
    const int nlayers = m_face->nplanes() + nbounds_layers;
    std::vector<std::vector<Activity::value_t> > measures(nlayers);
    measures[0].push_back(1);  // assume first two layers in RayGrid::Coordinates
    measures[1].push_back(1);  // are for horiz/vert bounds

    const int nactivities = slice->activity().size();
    int total_activity = 0;
    if (nactivities < m_face->nplanes()) {
        SPDLOG_LOGGER_TRACE(log, "anode={} face={} slice={} too few activities n={} / nplanes={}", anodeid, faceid, slice->ident(), nactivities, m_face->nplanes());
        return true;
    }

    for (const auto& chv : slice->activity()) {
        for (const auto& wire : chv.first->wires()) {
            auto wpid = wire->planeid();
            if (wpid.face() != faceid) {
                // l->trace("anode={} face={} slice={} chan={} wip={} q={} skip {}", anodeid, faceid, slice->ident(),
                // chv.first->ident(), wire->index(), chv.second, wpid);
                continue;
            }
            // l->trace("anode={} face={} slice={} chan={} wip={} q={} keep {}", anodeid, faceid, slice->ident(),
            // chv.first->ident(), wire->index(), chv.second, wpid);
            const int pit_ind = wire->index();
            const int layer = nbounds_layers + wire->planeid().index();
            auto& m = measures[layer];
            if (pit_ind < 0) {
                continue;
            }
            if ((int) m.size() <= pit_ind) {
                m.resize(pit_ind + 1, 0.0);
            }
            m[pit_ind] += 1.0;
            ++total_activity;
        }
    }

    if (!total_activity) {
        SPDLOG_LOGGER_TRACE(log, "anode={} face={} slice={} no total activity", anodeid, faceid, slice->ident());
        return true;
    }
    size_t nactive_layers = 0;
    for (size_t ind = 0; ind < measures.size(); ++ind) {
        const auto& blah = measures[ind];
        if (blah.empty()) {
            SPDLOG_LOGGER_TRACE(log, "anode={} face={} slice={} empty active layer ind={} out of {}", anodeid, faceid, slice->ident(), ind, measures.size());
            continue;
        }
        ++nactive_layers;
    }
    if (nactive_layers != measures.size()) {
        SPDLOG_LOGGER_TRACE(log, "anode={} face={} slice={} missing active layers {}  != {}, {} activities", anodeid, faceid, slice->ident(), nactive_layers, measures.size(), nactivities);
        return true;
    }

    activities_t activities;
    for (int layer = 0; layer < nlayers; ++layer) {
        auto& m = measures[layer];
        Activity activity(layer, {m.begin(), m.end()}, 0, m_threshold);
        //SPDLOG_LOGGER_TRACE(log, "L{} A={}", layer, activity.as_string());
        activities.push_back(activity);
    }

    // SPDLOG_LOGGER_TRACE(log, "anode={} face={} slice={} making blobs",
    //                     anodeid, faceid, slice->ident());
    auto blobs = make_blobs(m_face->raygrid(), activities, m_nudge);

    const float blob_value = 0.0;  // tiling doesn't consider particular charge
    for (const auto& blob : blobs) {
        Aux::SimpleBlob* sb = new Aux::SimpleBlob(m_blobs_seen++, blob_value,
                                                  0.0, blob, slice, m_face);
        sbs->m_blobs.push_back(IBlob::pointer(sb));
    }
    SPDLOG_LOGGER_TRACE(log, "anode={} face={} slice={} produced {} blobs",
                        anodeid, faceid, slice->ident(),
                        sbs->m_blobs.size());

    return true;
}
