#include "WireCellGen/DepoSplat.h"
#include "WireCellGen/GaussianDiffusion.h"

#include "WireCellUtil/Binning.h"
#include "WireCellUtil/NamedFactory.h"

// from ductor
#include "WireCellGen/BinnedDiffusion.h"

#include "WireCellUtil/Units.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellAux/SimpleTrace.h"
#include "WireCellAux/SimpleFrame.h"
#include <string>

// #include <iostream>


WIRECELL_FACTORY(DepoSplat, WireCell::Gen::DepoSplat, WireCell::IDuctor, WireCell::IConfigurable)

using namespace std;
using namespace WireCell;
using WireCell::Aux::SimpleTrace;
using WireCell::Aux::SimpleFrame;

Gen::DepoSplat::~DepoSplat() {}


Gen::DepoSplat::DepoSplat()
  : m_anode_tn("AnodePlane")
  , m_start_time(0.0 * units::ns)
  , m_readout_time(5.0 * units::ms)
  , m_tick(0.5 * units::us)
  , m_drift_speed(1.0 * units::mm / units::us)
  , m_nsigma(3.0)
  , m_mode("continuous")
  , m_frame_count(0)
  , l(Log::logger("sim"))
{
}

WireCell::Configuration Gen::DepoSplat::default_configuration() const
{
    Configuration cfg;

    /// How many Gaussian sigma due to diffusion to keep before truncating.
    put(cfg, "nsigma", m_nsigma);

    /// The initial time for this ductor
    put(cfg, "start_time", m_start_time);

    /// The time span for each readout.
    put(cfg, "readout_time", m_readout_time);

    /// The sample period
    put(cfg, "tick", m_tick);

    /// If false then determine start time of each readout based on the
    /// input depos.  This option is useful when running WCT sim on a
    /// source of depos which have already been "chunked" in time.  If
    /// true then this Ductor will continuously simulate all time in
    /// "readout_time" frames leading to empty frames in the case of
    /// some readout time with no depos.
    put(cfg, "continuous", true);

    /// Fixed mode simply reads out the same time window all the time.
    /// It implies discontinuous (continuous == false).
    put(cfg, "fixed", false);

    /// The nominal speed of drifting electrons
    put(cfg, "drift_speed", m_drift_speed);

    /// Allow for a custom starting frame number
    put(cfg, "first_frame_number", m_frame_count);

    /// Name of component providing the anode plane.
    put(cfg, "anode", m_anode_tn);

    // Tag to apply to output frame if non-empty.
    cfg["frame_tag"] = "";

    return cfg;
}

void Gen::DepoSplat::configure(const WireCell::Configuration& cfg)
{
    m_anode_tn = get<string>(cfg, "anode", m_anode_tn);
    m_anode = Factory::find_tn<IAnodePlane>(m_anode_tn);

    m_nsigma = get<double>(cfg, "nsigma", m_nsigma);
    bool continuous = get<bool>(cfg, "continuous", true);
    bool fixed = get<bool>(cfg, "fixed", false);

    m_mode = "continuous";
    if (fixed) {
        m_mode = "fixed";
    }
    else if (!continuous) {
        m_mode = "discontinuous";
    }

    m_rng = nullptr;
    if (get<bool>(cfg, "fluctuate", false)) {
        auto tn = get<std::string>(cfg, "rng", "Random");
        m_rng = Factory::find_tn<IRandom>(tn);
    }

    m_readout_time = get<double>(cfg, "readout_time", m_readout_time);
    m_tick = get<double>(cfg, "tick", m_tick);
    m_start_time = get<double>(cfg, "start_time", m_start_time);
    m_drift_speed = get<double>(cfg, "drift_speed", m_drift_speed);
    m_frame_count = get<int>(cfg, "first_frame_number", m_frame_count);

    m_frame_tag = get<std::string>(cfg, "frame_tag", "");

    l->debug("DepoSplat: tagging {}, AnodePlane: {}, mode: {}, time start: {} ms, readout time: {} ms, frame start: {}, fluctuate: {}",
             m_frame_tag,
             m_anode_tn, m_mode, m_start_time / units::ms, m_readout_time / units::ms,
             m_frame_count, m_rng != nullptr);
}


void Gen::DepoSplat::process(output_queue& frames)
{
    ITrace::vector traces;

    for (auto face : m_anode->faces()) {
        // Select the depos which are in this face's sensitive volume
        IDepo::vector face_depos, dropped_depos;
        auto bb = face->sensitive();
        if (bb.empty()) {
            l->debug("DepoSplat: anode: {} face: {} is "
                     "marked insensitive, skipping",
                     m_anode->ident(), face->ident());
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
            l->debug("DepoSplat: "
                "anode: {}, face: {}, processing {} depos spanning "
                "t:[{},{}]ms, bb:[{}-->{}]cm",
                m_anode->ident(), face->ident(), face_depos.size(), face_depos.front()->time() / units::ms,
                face_depos.back()->time() / units::ms, ray.first / units::cm, ray.second / units::cm);
        }
        if (dropped_depos.size()) {
            auto ray = bb.bounds();
            l->debug("DepoSplat: "
                "anode: {}, face: {}, dropped {} depos spanning "
                "t:[{},{}]ms, outside bb:[{}-->{}]cm",
                m_anode->ident(), face->ident(), dropped_depos.size(), dropped_depos.front()->time() / units::ms,
                dropped_depos.back()->time() / units::ms, ray.first / units::cm, ray.second / units::cm);
        }

        auto newtraces = process_face(face, face_depos);
        traces.insert(traces.end(), newtraces.begin(), newtraces.end());
    }

    auto frame = make_shared<SimpleFrame>(m_frame_count, m_start_time, traces, m_tick);
    IFrame::trace_list_t indices(traces.size());
    for (size_t ind = 0; ind < traces.size(); ++ind) {
        indices[ind] = ind;
    }
    frame->tag_frame(m_frame_tag);
    frames.push_back(frame);
    l->debug("DepoSplat: made frame: {} with {} traces @ {}ms", m_frame_count, traces.size(), m_start_time / units::ms);

    // fixme: what about frame overflow here?  If the depos extend
    // beyond the readout where does their info go?  2nd order,
    // diffusion and finite field response can cause depos near the
    // end of the readout to have some portion of their waveforms
    // lost?
    m_depos.clear();

    if (m_mode == "continuous") {
        m_start_time += m_readout_time;
    }

    ++m_frame_count;
}

// Return true if ready to start processing and capture start time if
// in continuous mode.
bool Gen::DepoSplat::start_processing(const input_pointer& depo)
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
    bool ok = depo->time() > m_start_time + m_readout_time;
    return ok;
}

bool Gen::DepoSplat::operator()(const input_pointer& depo, output_queue& frames)
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


ITrace::vector Gen::DepoSplat::process_face(IAnodeFace::pointer face, const IDepo::vector& depos)

{
    // why????
    const int time_offset = 2;  // # of ticks
    // const double difusion_scaler = 6.;
    const double charge_scaler = 1.;  // 18.;

    // channel-charge map
    std::unordered_map<int, std::vector<float> > chch;

    // tick-edged bins
    Binning tbins(m_readout_time / m_tick, m_start_time, m_start_time + m_readout_time);

    int iplane = -1;
    for (auto plane : face->planes()) {
        ++iplane;

        const Pimpos* pimpos = plane->pimpos();

        // wire-centered pitch bins
        const Binning& wbins = pimpos->region_binning();

        auto& wires = plane->wires();

        // std::cerr << "splat: plane " << plane->planeid() << " "
        //           << "wbins:" << wbins << " "
        //           << "tbins:" << tbins << " "
        //           << "#wires:" << wires.size() << " "
        //           << "#depos:" << depos.size() << "\n";

        int t_dropped{0}, p_dropped{0};

        int idepo = 0;
        for (auto depo : depos) {
            // const double tsig = depo->extent_long() * difusion_scaler;
            // const double psig = depo->extent_tran() * difusion_scaler;

            double sigma_L = depo->extent_long();
            double sigma_T = depo->extent_tran();

            // l->info("dirft: sigma_L: {} sigma_T: {}", sigma_L, sigma_T);

            if (true) {
                int nrebin = 1;
                double time_slice_width = nrebin * m_drift_speed * m_tick;                         // units::mm
                double add_sigma_L = 1.428249 * time_slice_width / nrebin / (m_tick / units::us);  // units::mm
                sigma_L = sqrt(pow(depo->extent_long(), 2) + pow(add_sigma_L, 2));  // / time_slice_width;
            }

            // What are these magic numbers???
            if (true) {
                double add_sigma_T = wbins.binsize();
                if (iplane == 0)
                    add_sigma_T *= (0.402993 * 0.3);
                else if (iplane == 1)
                    add_sigma_T *= (0.402993 * 0.5);
                else if (iplane == 2)
                    add_sigma_T *= (0.188060 * 0.2);
                sigma_T = sqrt(pow(depo->extent_tran(), 2) + pow(add_sigma_T, 2));  // / wbins.binsize();
            }

            // l->info("final: sigma_L: {} sigma_T: {}", sigma_L, sigma_T);

            const double tsig = sigma_L / m_drift_speed;
            const double psig = sigma_T;

            const double pwid = m_nsigma * psig;
            const double pcen = pimpos->distance(depo->pos());

            const double twid = m_nsigma * tsig;
            const double tcen = depo->time();

            const int pbeg = std::max(wbins.bin(pcen - pwid), 0);
            const int pend = std::min(wbins.bin(pcen + pwid) + 1, (int) wires.size());
            const int tbeg = std::max(tbins.bin(tcen - twid), 0);                      // fixme what limits
            const int tend = std::min(tbins.bin(tcen + twid) + 1, tbins.nbins() - 1);  //  to enforce here?

            if (tbeg > tend) {
                ++t_dropped;
                // l->debug("DepoSplat: dropping depo times (us) "
                //          "tbeg:{} > tend:{}, "
                //          "tcen:{}, twid:{}",
                //          tbeg/units::us, tend/units::us,
                //          tcen/units::us, twid/units::us);
                continue;
            }

            if (pbeg > pend) {
                ++p_dropped;
                continue;
            }

            Gen::GausDesc time_desc(tcen, tsig);
            Gen::GausDesc pitch_desc(pcen, psig);

            auto gd = std::make_shared<Gen::GaussianDiffusion>(depo, time_desc, pitch_desc);
            gd->set_sampling(tbins, wbins, m_nsigma, m_rng, 1);
            const auto patch = gd->patch();

            // std::stringstream ss;
            // ss << "splat: depo=" << depo->pos()/units::mm << "mm "
            //           << "@" << depo->time()*m_drift_speed<< " mm "
            //           << "p=(" << pcen << "+-" << pwid << "), t=(" << tcen << "+-" << twid << ") "
            //           << "pi=[" << pbeg << " " << pend << "], ti=[" << tbeg << " " << tend << "]";
            // l->info(ss.str());
            // l->info("tsig: {} m_drift_speed: {} ", tsig, m_drift_speed);
            // l->info("psig: {}", psig);
            // l->info("tbins: ({},{}) nbin: {}", tbins.min(), tbins.max(), tbins.nbins() );
            // l->info("wbins: ({},{}) nbin: {}", wbins.min(), wbins.max(), wbins.nbins() );
            // l->info("tbin range: ({},{})", tbeg, tend );
            // l->info("wbin range: ({},{})", pbeg, pend );
            // l->info("gd->offset: {}, {}", gd->toffset_bin(), gd->poffset_bin());
            // l->info("patch bins: {}, {}", patch.cols(), patch.rows());
            // l->info("\n");

            for (int ip = pbeg; ip < pend; ++ip) {
                auto irow = ip - gd->poffset_bin();
                if (irow < 0 or irow >= patch.rows()) continue;
                auto iwire = wires[ip];
                auto& charge = chch[iwire->channel()];
                if ((int) charge.size() < tend) {
                    charge.resize(tend, 0.0);
                }
                for (int it = tbeg; it < tend; ++it) {
                    auto icol = it - gd->toffset_bin() + time_offset;
                    if (icol < 0 or icol >= patch.cols()) continue;
                    charge[it] += std::abs(patch(irow, icol) * charge_scaler);
                }
            }
            ++idepo;
        } // over depos
        l->debug("DepoSplat: plane {} "
                 "dropped {} (time) and {} (pitch) from {} total",
                 iplane, t_dropped, p_dropped, depos.size());
         (void)idepo;
    }

    // make output traces
    ITrace::vector traces;
    for (auto& chchit : chch) {
        const int chid = chchit.first;
        auto& chv = chchit.second;
        auto trace = std::make_shared<SimpleTrace>(chid, 0, chv);
        traces.push_back(trace);
    }
   
    return traces;
}
