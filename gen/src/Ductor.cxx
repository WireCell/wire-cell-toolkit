#include "WireCellGen/Ductor.h"
#include "WireCellGen/BinnedDiffusion.h"
#include "WireCellGen/ImpactZipper.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellIface/SimpleTrace.h"
#include "WireCellIface/SimpleFrame.h"

#include <string>

#include "WireCellUtil/nljs2jcpp.hpp" // remove when ditch JsonCPP
#include "WireCellGen/Cfg/Ductor/Nljs.hpp"
using nljs_t = WireCellGen::Cfg::Ductor::data_t;

WIRECELL_FACTORY(Ductor, WireCell::Gen::Ductor, WireCell::IDuctor, WireCell::IConfigurable)

using namespace std;
using namespace WireCell;

Gen::Ductor::Ductor()
    : l(Log::logger("sim"))
{
}

WireCell::Configuration Gen::Ductor::default_configuration() const
{
    nljs_t nljs = m_cfg;
    return nljs.get<Json::Value>();
}

void Gen::Ductor::configure(const WireCell::Configuration& cfg)
{
    nljs_t nljs = cfg;
    m_cfg = nljs.get<config_t>();
    
    m_anode = Factory::find_tn<IAnodePlane>(m_cfg.anode);

    m_mode = "continuous";
    if (m_cfg.fixed) {
        m_mode = "fixed";
    }
    else if (!m_cfg.continuous) {
        m_mode = "discontinuous";
    }

    m_rng = nullptr;
    if (m_cfg.fluctuate) {
        m_rng = Factory::find_tn<IRandom>(m_cfg.rng);
    }

    if (m_cfg.pirs.empty()) {
        l->critical("must configure with some plane impace response components");
        THROW(ValueError() << errmsg{"Gen::Ductor: must configure with some plane impact response components"});
    }
    m_pirs.clear();
    for (const auto& tn : m_cfg.pirs) {
        auto pir = Factory::find_tn<IPlaneImpactResponse>(tn);
        m_pirs.push_back(pir);
    }

    m_frame_count = m_cfg.first_frame_number;
    m_start_time = m_cfg.start_time;

    l->debug("AnodePlane: {}, mode: {}, fluctuate: {}, time start: {} ms, readout time: {} ms, frame start: {}",
             m_cfg.anode, m_mode, (m_cfg.fluctuate ? "on" : "off"), m_start_time / units::ms, m_cfg.readout_time / units::ms,
             m_frame_count);
}

ITrace::vector Gen::Ductor::process_face(IAnodeFace::pointer face, const IDepo::vector& face_depos)
{
    ITrace::vector traces;

    int iplane = -1;
    for (auto plane : face->planes()) {
        ++iplane;

        const Pimpos* pimpos = plane->pimpos();

        Binning tbins(m_cfg.readout_time / m_cfg.tick, m_start_time, m_start_time + m_cfg.readout_time);

        Gen::BinnedDiffusion bindiff(*pimpos, tbins, m_cfg.nsigma, m_rng);
        for (auto depo : face_depos) {
            bindiff.add(depo, depo->extent_long() / m_cfg.drift_speed, depo->extent_tran());
        }

        auto& wires = plane->wires();

        auto pir = m_pirs.at(iplane);
        Gen::ImpactZipper zipper(pir, bindiff);

        const int nwires = pimpos->region_binning().nbins();
        for (int iwire = 0; iwire < nwires; ++iwire) {
            auto wave = zipper.waveform(iwire);

            auto mm = Waveform::edge(wave);
            if (mm.first == (int) wave.size()) {  // all zero
                continue;
            }

            int chid = wires[iwire]->channel();
            int tbin = mm.first;

            ITrace::ChargeSequence charge(wave.begin() + mm.first, wave.begin() + mm.second);
            auto trace = make_shared<SimpleTrace>(chid, tbin, charge);
            traces.push_back(trace);
        }
    }
    return traces;
}

void Gen::Ductor::process(output_queue& frames)
{
    ITrace::vector traces;

    for (auto face : m_anode->faces()) {
        // Select the depos which are in this face's sensitive volume
        IDepo::vector face_depos, dropped_depos;
        auto bb = face->sensitive();
        if (bb.empty()) {
            l->debug("anode: {} face: {} is marked insensitive, skipping", m_anode->ident(), face->ident());
            continue;
        }

        for (auto depo : m_depos) {
            if (bb.inside(depo->pos())) {
                face_depos.push_back(depo);
            }
            else {
                dropped_depos.push_back(depo);
            }
        }

        if (face_depos.size()) {
            auto ray = bb.bounds();
            l->debug(
                "anode: {}, face: {}, processing {} depos spanning "
                "t:[{},{}]ms, bb:[{}-->{}]cm",
                m_anode->ident(), face->ident(), face_depos.size(), face_depos.front()->time() / units::ms,
                face_depos.back()->time() / units::ms, ray.first / units::cm, ray.second / units::cm);
        }
        if (dropped_depos.size()) {
            auto ray = bb.bounds();
            l->debug(
                "anode: {}, face: {}, dropped {} depos spanning "
                "t:[{},{}]ms, outside bb:[{}-->{}]cm",
                m_anode->ident(), face->ident(), dropped_depos.size(), dropped_depos.front()->time() / units::ms,
                dropped_depos.back()->time() / units::ms, ray.first / units::cm, ray.second / units::cm);
        }

        auto newtraces = process_face(face, face_depos);
        traces.insert(traces.end(), newtraces.begin(), newtraces.end());
    }

    auto frame = make_shared<SimpleFrame>(m_frame_count, m_start_time, traces, m_cfg.tick);
    IFrame::trace_list_t indices(traces.size());
    for (size_t ind = 0; ind < traces.size(); ++ind) {
        indices[ind] = ind;
    }
    std::string trace_tag("ductor");
    frame->tag_traces(string("ductor") + std::to_string(m_anode->ident()), indices);
    frame->tag_frame("ductor");
    frames.push_back(frame);
    l->debug("made frame: {} with {} traces @ {}ms", m_frame_count, traces.size(), m_start_time / units::ms);

    // fixme: what about frame overflow here?  If the depos extend
    // beyond the readout where does their info go?  2nd order,
    // diffusion and finite field response can cause depos near the
    // end of the readout to have some portion of their waveforms
    // lost?
    m_depos.clear();

    if (m_mode == "continuous") {
        m_start_time += m_cfg.readout_time;
    }

    ++m_frame_count;
}

// Return true if ready to start processing and capture start time if
// in continuous mode.
bool Gen::Ductor::start_processing(const input_pointer& depo)
{
    if (!depo) {
        return true;
    }

    if (m_mode == "fixed") {
        // fixed mode waits until EOS
        return false;
    }

    if (m_mode == "discontinuous") {
        // discontinuous mode sets start time on first depo.
        if (m_depos.empty()) {
            m_start_time = depo->time();
            return false;
        }
    }

    // continuous and discontinuous modes follow Just Enough
    // Processing(TM) strategy.

    // Note: we use this depo time even if it may not actually be
    // inside our sensitive volume.
    bool ok = depo->time() > m_start_time + m_cfg.readout_time;
    return ok;
}

bool Gen::Ductor::operator()(const input_pointer& depo, output_queue& frames)
{
    if (start_processing(depo)) {
        process(frames);
    }

    if (depo) {
        m_depos.push_back(depo);
    }
    else {
        frames.push_back(nullptr);
    }

    return true;
}
