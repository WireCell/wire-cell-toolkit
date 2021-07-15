#include "WireCellGen/DepoSplat.h"
#include "WireCellGen/GaussianDiffusion.h"
#include "WireCellIface/SimpleTrace.h"

#include "WireCellUtil/Binning.h"
#include "WireCellUtil/NamedFactory.h"

#include <iostream>

WIRECELL_FACTORY(DepoSplat, WireCell::Gen::DepoSplat, WireCell::IDuctor, WireCell::IConfigurable)

using namespace WireCell;

Gen::DepoSplat::DepoSplat()
  : Ductor()
  , l(Log::logger("sim"))
{
}

Gen::DepoSplat::~DepoSplat() {}

ITrace::vector Gen::DepoSplat::process_face(IAnodeFace::pointer face, const IDepo::vector& depos)

{
    const int time_offset = 2;  // # of ticks
    // const double difusion_scaler = 6.;
    const double charge_scaler = 1.;  // 18.;

    // channel-charge map
    std::unordered_map<int, std::vector<float> > chch;

    // tick-edged bins
    Binning tbins(m_cfg.readout_time / m_cfg.tick, m_start_time, m_start_time + m_cfg.readout_time);

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

        int idepo = 0;
        for (auto depo : depos) {
            // const double tsig = depo->extent_long() * difusion_scaler;
            // const double psig = depo->extent_tran() * difusion_scaler;

            double sigma_L = depo->extent_long();
            double sigma_T = depo->extent_tran();

            // l->info("dirft: sigma_L: {} sigma_T: {}", sigma_L, sigma_T);

            if (true) {
                int nrebin = 1;
                double time_slice_width = nrebin * m_cfg.drift_speed * m_cfg.tick; // units::mm
                double add_sigma_L = 1.428249 * time_slice_width / nrebin / (m_cfg.tick / units::us); // units::mm
                sigma_L = sqrt(pow(depo->extent_long(), 2) + pow(add_sigma_L, 2));  // / time_slice_width;
            }

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

            const double tsig = sigma_L / m_cfg.drift_speed;
            const double psig = sigma_T;

            const double pwid = m_cfg.nsigma * psig;
            const double pcen = pimpos->distance(depo->pos());

            const double twid = m_cfg.nsigma * tsig;
            const double tcen = depo->time();

            const int pbeg = std::max(wbins.bin(pcen - pwid), 0);
            const int pend = std::min(wbins.bin(pcen + pwid) + 1, (int) wires.size());
            const int tbeg = std::max(tbins.bin(tcen - twid), 0);                      // fixme what limits
            const int tend = std::min(tbins.bin(tcen + twid) + 1, tbins.nbins() - 1);  //  to enforce here?

            if (tbeg > tend) continue;

            if (pbeg > pend) continue;

            Gen::GausDesc time_desc(tcen, tsig);
            Gen::GausDesc pitch_desc(pcen, psig);

            auto gd = std::make_shared<Gen::GaussianDiffusion>(depo, time_desc, pitch_desc);
            gd->set_sampling(tbins, wbins, m_cfg.nsigma, 0, 1);
            const auto patch = gd->patch();


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
        }
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
