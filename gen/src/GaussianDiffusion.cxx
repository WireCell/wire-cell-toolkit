#include "WireCellGen/GaussianDiffusion.h"

#include <iostream>  // debugging

using namespace WireCell;
using namespace std;

std::vector<double> Gen::GausDesc::sample(double start, double step, int nsamples) const
{
    std::vector<double> ret;

    if (!sigma) {
        ret.resize(1, 0);
        ret[0] = 1;
        if (nsamples != 1) {
            cerr << "NOT one sample for true point source: " << nsamples << "\n";
        }
    }
    else {
        ret.resize(nsamples, 0.0);
        for (int ind = 0; ind < nsamples; ++ind) {
            const double rel = (start + ind * step - center) / sigma;
            ret[ind] = exp(-0.5 * rel * rel);
        }
    }

    return ret;
}

std::vector<double> Gen::GausDesc::binint(double start, double step, int nbins) const
{
    std::vector<double> bins;

    if (!sigma) {
        bins.resize(1, 0);
        bins[0] = 1;
        if (nbins != 1) {
            cerr << "NOT one bin for true point source: " << nbins << "\n";
        }
    }
    else {
        bins.resize(nbins, 0.0);
        std::vector<double> erfs(nbins + 1, 0.0);
        const double sqrt2 = sqrt(2.0);
        for (int ind = 0; ind <= nbins; ++ind) {
            double x = (start + step * ind - center) / (sqrt2 * sigma);
            erfs[ind] = 0.5 * std::erf(x);
        }

        // double tot = 0.0;
        for (int ibin = 0; ibin < nbins; ++ibin) {
            const double val = erfs[ibin + 1] - erfs[ibin];
            // tot += val;
            bins[ibin] = val;
        }
    }
    return bins;
}

// integral Normal distribution with weighting function
// a linear weighting for the charge in each pbin
// Integral of charge spectrum <pvec> done by GausDesc::binint (do not do it again using Erf())
std::vector<double> Gen::GausDesc::weight(double start, double step, int nbins, std::vector<double> pvec) const
{
    std::vector<double> wt;
    if (!sigma) {
        wt.resize(1, 0);
        wt[0] = (start + step - center) / step;
    }
    else {
        wt.resize(nbins, 0.0);
        const double pi = 4.0 * atan(1);
        double x2 = start;
        double x1 = 0;
        double gaus2 = exp(-0.5 * (start - center) / sigma * (start - center) / sigma);
        double gaus1 = 0;
        for (int ind = 0; ind < nbins; ind++) {
            x1 = x2;
            x2 = x1 + step;
            double rel = (x2 - center) / sigma;
            gaus1 = gaus2;
            gaus2 = exp(-0.5 * rel * rel);

            // weighting
            wt[ind] =
                -1.0 * sigma / (x1 - x2) * (gaus2 - gaus1) / sqrt(2.0 * pi) / pvec[ind] + (center - x2) / (x1 - x2);
            /* std::cerr<<"Gaus: "<<"1 and 2 "<<gaus1<<", "<<gaus2<<std::endl; */
            /* std::cerr<<"center, x1, x2: "<<center<<", "<<x1<<", "<<x2<<std::endl; */
            /* std::cerr<<"Total charge: "<<pvec[ind]<<std::endl; */
            /* std::cerr<<"weight: "<<ind<<" "<<wt[ind]<<std::endl; */
        }
    }
    return wt;
}

// std::pair<int,int> Gen::GausDesc::subsample_range(int nsamples, double xmin, double xmax, double nsigma) const
// {
//     const double sample_size = (xmax-xmin)/(nsamples-1);
//     // find closest sample indices
//     int imin = int(round((center - nsigma*sigma - xmin)/sample_size));
//     int imax = int(round((center + nsigma*sigma - xmin)/sample_size));

//     return std::make_pair(std::max(imin, 0), std::min(imax+1, nsamples));
// }

/// GaussianDiffusion

Gen::GaussianDiffusion::GaussianDiffusion(const IDepo::pointer& depo, const GausDesc& time_desc,
                                          const GausDesc& pitch_desc)
  : m_deposition(depo)
  , m_time_desc(time_desc)
  , m_pitch_desc(pitch_desc)
  , m_toffset_bin(-1)
  , m_poffset_bin(-1)
{
}

void Gen::GaussianDiffusion::set_sampling(const Binning& tbin,  // overall time tick binning
                                          const Binning& pbin,  // overall impact position binning
                                          double nsigma, IRandom::pointer fluctuate, unsigned int weightstrat)
{
    if (m_patch.size() > 0) {
        return;
    }

    /// Sample time dimension
    auto tval_range = m_time_desc.sigma_range(nsigma);
    auto tbin_range = tbin.sample_bin_range(tval_range.first, tval_range.second);
    const size_t ntss = tbin_range.second - tbin_range.first;
    m_toffset_bin = tbin_range.first;
    // auto tvec =  m_time_desc.sample(tbin.center(m_toffset_bin), tbin.binsize(), ntss);
    auto tvec = m_time_desc.binint(tbin.edge(m_toffset_bin), tbin.binsize(), ntss);

    if (!ntss) {
        cerr << "Gen::GaussianDiffusion: no time bins for [" << tval_range.first / units::us << ","
             << tval_range.second / units::us << "] us\n";
        return;
    }

    /// Sample pitch dimension.
    auto pval_range = m_pitch_desc.sigma_range(nsigma);
    auto pbin_range = pbin.sample_bin_range(pval_range.first, pval_range.second);
    const size_t npss = pbin_range.second - pbin_range.first;
    m_poffset_bin = pbin_range.first;
    // auto pvec = m_pitch_desc.sample(pbin.center(m_poffset_bin), pbin.binsize(), npss);
    auto pvec = m_pitch_desc.binint(pbin.edge(m_poffset_bin), pbin.binsize(), npss);

    if (!npss) {
        cerr << "No impact bins [" << pval_range.first / units::mm << "," << pval_range.second / units::mm << "] mm\n";
        return;
    }

    // weightstrat = 1;

    // make charge weights for later interpolation.
    /// fixme: for hanyu.
    if (weightstrat == 2) {
        auto wvec = m_pitch_desc.weight(pbin.edge(m_poffset_bin), pbin.binsize(), npss, pvec);
        m_qweights = wvec;
    }
    if (weightstrat == 1) {
        m_qweights.resize(npss, 0.5);
    }

    // start making the time vs impact patch of charge.
    patch_t ret = patch_t::Zero(npss, ntss);
    double raw_sum = 0.0;

    // Convolve the two independent Gaussians
    for (size_t ip = 0; ip < npss; ++ip) {
        for (size_t it = 0; it < ntss; ++it) {
            const double val = pvec[ip] * tvec[it];
            raw_sum += val;
            ret(ip, it) = (float) val;
        }
    }

    // Depo charge should be in units of "e" so negative, but
    // explicitly track sign in case positive charge is given.
    const double depo_charge = m_deposition->charge();
    const double charge_sign = depo_charge < 0 ? -1 : 1;

    // normalize to total charge
    ret *= depo_charge / raw_sum;

    double fluc_sum = 0;
    if (fluctuate) {
        // double unfluc_sum = 0;

        for (size_t ip = 0; ip < npss; ++ip) {
            for (size_t it = 0; it < ntss; ++it) {
                const double oldval = ret(ip, it);
                // unfluc_sum += oldval;
                // should be a multinomial distribution, n_i follows binomial distribution
                // but n_i, n_j has covariance -n_tot * p_i * p_j
                // normalize later to approximate this multinomial distribution (how precise?)
                // how precise? better than poisson and 10000 total charge corresponds to a <1% level error.
                double number = 0.0;
                const double relval = oldval / depo_charge;
                if (relval >= 1.0) {
                    number = (int) std::abs(depo_charge);
                }
                else {
                    number = fluctuate->binomial((int) std::abs(depo_charge), relval);
                }
                // the charge should be negative -- ionization electrons
                number *= charge_sign;
                fluc_sum += number;
                ret(ip, it) = number;
            }
        }
        if (fluc_sum == 0) {
            // cerr << "No net charge after fluctuation. Total unfluctuated = "
            //      << unfluc_sum
            //      << " Qdepo = " << m_deposition->charge()
            //      << endl;
            return;
        }
        else {
            ret *= m_deposition->charge() / fluc_sum;
        }
    }

    {  // debugging
       // double retsum=0.0;
       // for (size_t ip = 0; ip < npss; ++ip) {
       //    for (size_t it = 0; it < ntss; ++it) {
       //        retsum += ret(ip,it);
       //    }
       //}
       // cerr << "GaussianDiffusion: Q in electrons: depo=" << m_deposition->charge()/units::eplus
       //      << " rawsum=" << raw_sum/units::eplus << " flucsum=" << fluc_sum/units::eplus
       //      << " returned=" << retsum/units::eplus << endl;
    }

    m_patch = ret;
}

void Gen::GaussianDiffusion::clear_sampling()
{
    m_patch.resize(0, 0);
    m_qweights.clear();
    m_qweights.shrink_to_fit();
}

// patch = nimpacts rows X nticks columns
// patch(row,col)
const Gen::GaussianDiffusion::patch_t& Gen::GaussianDiffusion::patch() const { return m_patch; }

const std::vector<double> Gen::GaussianDiffusion::weights() const { return m_qweights; }
