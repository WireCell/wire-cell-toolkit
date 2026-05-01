/***************************************************************************/
/*             Noise filter (NF) module for ProtoDUNE VD                   */
/***************************************************************************/
// Features:                       
//  baseline correction      
//  CNR removal   
//   Add ShieldCoupling removal for TDE u plane, idea from Lardon https://github.com/dune-lardon/lardon/blob/main/noise_filter.py#L181
//
// 10/19/2025, created by X.Ning (xning@bnl.gov) 


#include "WireCellSigProc/ProtoduneVD.h"
#include "WireCellSigProc/Derivations.h"

#include "WireCellAux/DftTools.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/Units.h"

#include "WireCellUtil/Persist.h"

#include <cmath>
#include <complex>
#include <iostream>
#include <set>
#include <fstream>

WIRECELL_FACTORY(PDVDOneChannelNoise, WireCell::SigProc::PDVD::OneChannelNoise, WireCell::IChannelFilter,
                 WireCell::IConfigurable)
WIRECELL_FACTORY(PDVDCoherentNoiseSub, WireCell::SigProc::PDVD::CoherentNoiseSub, WireCell::IChannelFilter,
                 WireCell::IConfigurable)
WIRECELL_FACTORY(PDVDFEMBNoiseSub, WireCell::SigProc::PDVD::FEMBNoiseSub, WireCell::IChannelFilter,
                 WireCell::IConfigurable)
WIRECELL_FACTORY(PDVDShieldCouplingSub, WireCell::SigProc::PDVD::ShieldCouplingSub, WireCell::IChannelFilter,
                 WireCell::IConfigurable)

using namespace WireCell::SigProc;
using WireCell::Aux::DftTools::fwd_r2c;
using WireCell::Aux::DftTools::inv_c2r;


bool PDVD::Subtract_WScaling(WireCell::IChannelFilter::channel_signals_t& chansig,
                                   const WireCell::Waveform::realseq_t& medians,
                                   const WireCell::Waveform::compseq_t& respec, int res_offset,
                                   std::vector<std::vector<int> >& rois,
                                   const IDFT::pointer& dft,
                                   const WireCell::Waveform::realseq_t& time_filter_wf,
                                   const WireCell::Waveform::realseq_t& lf_filter_wf,
                                   float decon_limit1, float roi_min_max_ratio,
                                   float rms_threshold,
                                   WireCell::SigProc::CoherentNoiseDump* dump)
{
    double ave_coef = 0;
    double_t ave_coef1 = 0;
    std::map<int, double> coef_all;

    const int nbin = medians.size();

    const size_t nrois = rois.size();
    if (dump) {
        const size_t nch = chansig.size();
        dump->roi_max_per_ch.assign(nch * nrois, 0.f);
        dump->roi_min_per_ch.assign(nch * nrois, 0.f);
        dump->roi_accepted_per_ch.assign(nch * nrois, 0);
        dump->scaling_coef.assign(nch, 0.f);
    }

    for (auto it : chansig) {
        int ch = it.first;
        WireCell::IChannelFilter::signal_t& signal = it.second;

        double sum2 = 0;
        double sum3 = 0;
        double coef = 0;

        std::pair<double, double> temp = Derivations::CalcRMS(signal);
        // std::pair<double,double> temp = WireCell::Waveform::mean_rms(signal);

        // if ( abs(ch-6117)<5)
        //     std::cout << ch << " " << temp.first << " " << temp.second << " "  << std::endl;

        for (int j = 0; j != nbin; j++) {
            if (fabs(signal.at(j)) < 4 * temp.second) {
                sum2 += signal.at(j) * medians.at(j);
                sum3 += medians.at(j) * medians.at(j);
            }
        }
        if (sum3 > 0) {
            coef = sum2 / sum3;
        }
        // protect against the extreme cases
        // if (coef < 0.6) coef = 0.6;
        // if (coef > 1.5) coef = 1.5;

        coef_all[ch] = coef;
        if (coef != 0) {  // FIXME: expect some fluctuation?
            ave_coef += coef;
            ave_coef1++;
        }
    }

    if (ave_coef1 > 0) {
        ave_coef = ave_coef / ave_coef1;
    }
    // const int achannel = chansig.begin()->first;
    // std::string filename = "scale_output"+ std::to_string(achannel)+"_"+std::to_string(ave_coef)+".txt";
    // std::ofstream outfile(filename);
    
    
    int ch_idx = 0;
    for (auto it : chansig) {
        int ch = it.first;
        WireCell::IChannelFilter::signal_t& signal = it.second;
        float scaling;
        // if (ave_coef != 0) {
        if (ave_coef > 0) {
            scaling = coef_all[ch] / ave_coef;
            // add some protections ...
            if (scaling < 0) scaling = 0;
            //	    if (scaling < 0.5 && scaling > 0.3) scaling = 0.5;
            if (scaling > 1.5) scaling = 1.5;
        }
        else {
            scaling = 0;
        }

        if (dump) dump->scaling_coef[ch_idx] = scaling;
        // if ( abs(ch-6117)<5)
        // if(achannel<3072*2){

            // std::cout << ch << " " << scaling << " "  << std::endl;
            // outfile << ch << " " << scaling <<  std::endl;
        // scaling = 1.0;
        // }


        if (respec.size() > 0 && (respec.at(0).real() != 1 || respec.at(0).imag() != 0) && res_offset != 0) {
            int nbin = signal.size();
            WireCell::Waveform::realseq_t signal_roi(nbin, 0);
            if (rms_threshold) {
                signal_roi = signal;
            }
            else {
                for (auto roi : rois) {
                    const int bin0 = std::max(roi.front() - 1, 0);
                    const int binf = std::min(roi.back() + 1, nbin - 1);
                    const double m0 = signal[bin0];
                    const double mf = signal[binf];
                    const double roi_run = binf - bin0;
                    const double roi_rise = mf - m0;
                    for (auto bin : roi) {
                        const double m = m0 + (bin - bin0) / roi_run * roi_rise;
                        signal_roi.at(bin) = signal.at(bin) - m;
                    }
                }
            }

            // do the deconvolution with a very loose low-frequency filter
            WireCell::Waveform::compseq_t signal_roi_freq = fwd_r2c(dft, signal_roi);
            WireCell::Waveform::shrink(signal_roi_freq, respec);
            for (size_t i = 0; i != signal_roi_freq.size(); i++) {
                signal_roi_freq.at(i) *= time_filter_wf.at(i) * lf_filter_wf.at(i);
            }
            WireCell::Waveform::realseq_t signal_roi_decon = inv_c2r(dft, signal_roi_freq);

            if (rms_threshold) {
                std::pair<double, double> temp = Derivations::CalcRMS(signal_roi_decon);
                double mean = temp.first;
                double rms = temp.second;
                // if (ch==580) {
                //     std::cout << "[Jujube] dbg_info_ch" << ch << " mean    " << mean << std::endl;
                //     std::cout << "[Jujube] dbg_info_ch" << ch << " rms     " << rms << std::endl;
                // }
                // std::cout << "[Jujube] dfg_rms_ch" << ch << "\t" << rms << std::endl;
                for (size_t i = 0; i != signal_roi_freq.size(); i++) {
                    signal_roi_decon.at(i) -= mean;
                }
                double rms_limit = rms_threshold * rms;
                if (rms_limit < 0.02 && rms_limit > 0) rms_limit = 0.02;
                if (rms_limit > decon_limit1) rms_limit = decon_limit1;
                decon_limit1 = rms_limit;
            }

            std::map<int, bool> flag_replace;
            for (auto roi : rois) {
                flag_replace[roi.front()] = false;
            }

            // judge if any ROI is good ...
            for (size_t r = 0; r < rois.size(); ++r) {
                const auto& roi = rois[r];
                const int bin0 = std::max(roi.front() - 1, 0);
                const int binf = std::min(roi.back() + 1, nbin - 1);
                double max_val = 0;
                double min_val = 0;
                for (int i = bin0; i <= binf; i++) {
                    int time_bin = i - res_offset;
                    if (time_bin < 0) time_bin += nbin;
                    if (time_bin >= nbin) time_bin -= nbin;
                    if (i == bin0) {
                        max_val = signal_roi_decon.at(time_bin);
                        min_val = signal_roi_decon.at(time_bin);
                    }
                    else {
                        if (signal_roi_decon.at(time_bin) > max_val) max_val = signal_roi_decon.at(time_bin);
                        if (signal_roi_decon.at(time_bin) < min_val) min_val = signal_roi_decon.at(time_bin);
                    }
                }

                const bool accepted = (max_val > decon_limit1 && fabs(min_val) < max_val * roi_min_max_ratio);
                if (accepted) flag_replace[roi.front()] = true;

                if (dump) {
                    const size_t idx = (size_t)ch_idx * nrois + r;
                    dump->roi_max_per_ch[idx] = max_val;
                    dump->roi_min_per_ch[idx] = min_val;
                    dump->roi_accepted_per_ch[idx] = accepted ? 1 : 0;
                }
            }

            //    for (auto roi: rois){
            // flag_replace[roi.front()] = true;
            //    }

            WireCell::Waveform::realseq_t temp_medians = medians;

            for (auto roi : rois) {
                // original code used the bins just outside the ROI
                const int bin0 = std::max(roi.front() - 1, 0);
                const int binf = std::min(roi.back() + 1, nbin - 1);
                if (flag_replace[roi.front()]) {
                    const double m0 = temp_medians[bin0];
                    const double mf = temp_medians[binf];
                    const double roi_run = binf - bin0;
                    const double roi_rise = mf - m0;
                    for (auto bin : roi) {
                        const double m = m0 + (bin - bin0) / roi_run * roi_rise;
                        temp_medians.at(bin) = m;
                    }
                }
            }



            // collection plane, directly subtracti ...
            for (int i = 0; i != nbin; i++) {
                if (fabs(signal.at(i)) > 0.001) {
                    signal.at(i) = signal.at(i) - temp_medians.at(i) * scaling;
                }
            }
        }
        else {
            // collection plane, directly subtracti ...
            for (int i = 0; i != nbin; i++) {
                if (fabs(signal.at(i)) > 0.001) {
                    signal.at(i) = signal.at(i) - medians.at(i) * scaling;
                }
            }
        }
        chansig[ch] = signal;
        ++ch_idx;
    }

    if (dump) {
        dump->ave_coef = ave_coef;
    }

    return true;
}

std::vector<std::vector<int> > PDVD::SignalProtection(WireCell::Waveform::realseq_t& medians,
                                                            const WireCell::Waveform::compseq_t& respec,
                                                            const IDFT::pointer& dft,
                                                            int res_offset,
                                                            int pad_f, int pad_b,
                                                            const WireCell::Waveform::realseq_t& time_filter_wf,
                                                            const WireCell::Waveform::realseq_t& lf_filter_wf,
                                                            float upper_decon_limit,
                                                            float upper_adc_limit,
                                                            float protection_factor, float min_adc_limit,
                                                            WireCell::SigProc::CoherentNoiseDump* dump)
{
    // WireCell::Waveform::realseq_t temp1;
    // for (int i=0;i!=medians.size();i++){
    //   if (fabs(medians.at(i) - mean) < 4.5*rms)
    //     temp1.push_back(medians.at(i));
    // }
    // temp = WireCell::Waveform::mean_rms(temp1);
    // mean = temp.first;
    // rms = temp.second;

    // std::cout << temp.first << " " << temp.second << std::endl;
    const int nbin = medians.size();

    //    const int protection_factor = 5.0;
    // move to input ...
    // const float upper_decon_limit = 0.05;
    // const float upper_adc_limit = 15;
    // const float min_adc_limit = 50;

    std::vector<bool> signalsBool;
    signalsBool.resize(nbin, false);

    if (dump) {
        dump->signal_bool_raw.assign(nbin, 0);
    }

    // calculate the RMS
    std::pair<double, double> temp = Derivations::CalcRMS(medians);
    double mean = temp.first;
    double rms = temp.second;

    float limit;
    if (protection_factor * rms > upper_adc_limit) {
        limit = protection_factor * rms;
    }
    else {
        limit = upper_adc_limit;
    }
    if (min_adc_limit < limit) {
        limit = min_adc_limit;
    }

    if (dump) {
        dump->adc_threshold_chosen = limit;
        dump->mean_adc = mean;
        dump->rms_adc = rms;
    }

    for (int j = 0; j != nbin; j++) {
        float content = medians.at(j);
        if (fabs(content - mean) > limit) {
            signalsBool.at(j) = true;
            if (dump) dump->signal_bool_raw.at(j) = 1;
            // add the front and back padding
            for (int k = 0; k != pad_b; k++) {
                int bin = j + k + 1;
                if (bin > nbin - 1) bin = nbin - 1;
                signalsBool.at(bin) = true;
            }
            for (int k = 0; k != pad_f; k++) {
                int bin = j - k - 1;
                if (bin < 0) {
                    bin = 0;
                }
                signalsBool.at(bin) = true;
            }
        }
    }

    // the deconvolution protection code ...
    WireCell::Waveform::realseq_t medians_decon;
    bool decon_stage_ran = false;
    if (respec.size() > 0 && (respec.at(0).real() != 1 || respec.at(0).imag() != 0) && res_offset != 0) {
        decon_stage_ran = true;

        WireCell::Waveform::compseq_t medians_freq = fwd_r2c(dft, medians);
        WireCell::Waveform::shrink(medians_freq, respec);

        for (size_t i = 0; i != medians_freq.size(); i++) {
            medians_freq.at(i) *= time_filter_wf.at(i) * lf_filter_wf.at(i);
        }
        medians_decon = inv_c2r(dft, medians_freq);

        temp = Derivations::CalcRMS(medians_decon);
        mean = temp.first;
        rms = temp.second;

        if (protection_factor * rms > upper_decon_limit) {
            limit = protection_factor * rms;
        }
        else {
            limit = upper_decon_limit;
        }

        if (dump) {
            dump->decon_threshold_chosen = limit;
            dump->mean_decon = mean;
            dump->rms_decon = rms;
            dump->decon_stage_ran = 1;
            dump->medians_decon_aligned.resize(nbin);
            for (int i = 0; i < nbin; ++i) {
                int dt_bin = i - res_offset;
                if (dt_bin < 0) dt_bin += nbin;
                if (dt_bin >= nbin) dt_bin -= nbin;
                dump->medians_decon_aligned[i] = medians_decon[dt_bin];
            }
        }

        for (int j = 0; j != nbin; j++) {
            float content = medians_decon.at(j);
            if ((content - mean) > limit) {
                int time_bin = j + res_offset;
                if (time_bin >= nbin) time_bin -= nbin;
                signalsBool.at(time_bin) = true;
                if (dump) dump->signal_bool_raw.at(time_bin) = 1;
                // add the front and back padding
                for (int k = 0; k != pad_b; k++) {
                    int bin = time_bin + k + 1;
                    if (bin > nbin - 1) bin = nbin - 1;
                    signalsBool.at(bin) = true;
                }
                for (int k = 0; k != pad_f; k++) {
                    int bin = time_bin - k - 1;
                    if (bin < 0) {
                        bin = 0;
                    }
                    signalsBool.at(bin) = true;
                }
            }
        }
    }

    // {
    // partition waveform indices into consecutive regions with
    // signalsBool true.
    std::vector<std::vector<int> > rois;
    bool inside = false;
    for (int ind = 0; ind < nbin; ++ind) {
        if (inside) {
            if (signalsBool[ind]) {  // still inside
                rois.back().push_back(ind);
            }
            else {
                inside = false;
            }
        }
        else {                       // outside the Rio
            if (signalsBool[ind]) {  // just entered ROI
                std::vector<int> roi;
                roi.push_back(ind);
                rois.push_back(roi);
                inside = true;
            }
        }
    }

    std::map<int, bool> flag_replace;
    for (auto roi : rois) {
        flag_replace[roi.front()] = true;
    }

    if (respec.size() > 0 && (respec.at(0).real() != 1 || respec.at(0).imag() != 0) && res_offset != 0) {
        for (auto roi : rois) {
            flag_replace[roi.front()] = false;
        }
    }

    // Replace medians for above regions with interpolation on values
    // just outside each region.
    for (auto roi : rois) {
        // original code used the bins just outside the ROI
        const int bin0 = std::max(roi.front() - 1, 0);
        const int binf = std::min(roi.back() + 1, nbin - 1);
        if (flag_replace[roi.front()]) {
            const double m0 = medians[bin0];
            const double mf = medians[binf];
            const double roi_run = binf - bin0;
            const double roi_rise = mf - m0;
            for (auto bin : roi) {
                const double m = m0 + (bin - bin0) / roi_run * roi_rise;
                medians.at(bin) = m;
            }
        }
    }

    if (dump) {
        dump->signal_bool.assign(nbin, 0);
        for (int i = 0; i < nbin; ++i) dump->signal_bool[i] = signalsBool[i] ? 1 : 0;
        dump->roi_starts.reserve(rois.size());
        dump->roi_ends.reserve(rois.size());
        for (auto& roi : rois) {
            dump->roi_starts.push_back(roi.front());
            dump->roi_ends.push_back(roi.back());
        }
        // Per-ROI median-decon scalars: same accept formula as Subtract_WScaling
        // but applied to the median rather than per-channel signal_roi_decon.
        const float dl1 = dump->decon_limit1;
        const float rmm = dump->roi_min_max_ratio;
        dump->roi_max_median.resize(rois.size(), 0.f);
        dump->roi_min_median.resize(rois.size(), 0.f);
        dump->roi_ratio_median.resize(rois.size(), 0.f);
        dump->roi_accepted_median.resize(rois.size(), 0);
        if (decon_stage_ran) {
            for (size_t r = 0; r < rois.size(); ++r) {
                const auto& roi = rois[r];
                const int bin0 = std::max(roi.front() - 1, 0);
                const int binf = std::min(roi.back() + 1, nbin - 1);
                double mx = 0, mn = 0;
                bool first = true;
                for (int i = bin0; i <= binf; ++i) {
                    int dt_bin = i - res_offset;
                    if (dt_bin < 0) dt_bin += nbin;
                    if (dt_bin >= nbin) dt_bin -= nbin;
                    const double v = medians_decon.at(dt_bin);
                    if (first) { mx = mn = v; first = false; }
                    else { if (v > mx) mx = v; if (v < mn) mn = v; }
                }
                dump->roi_max_median[r] = mx;
                dump->roi_min_median[r] = mn;
                dump->roi_ratio_median[r] = (mx > 0) ? std::fabs(mn) / mx : 1e9f;
                dump->roi_accepted_median[r] = (mx > dl1 && std::fabs(mn) < mx * rmm) ? 1 : 0;
            }
        }
    }

    return rois;
}


bool PDVD::SignalFilter(WireCell::Waveform::realseq_t& sig)
{
    const double sigFactor = 4.0;
    const int padBins = 8;

    float rmsVal = PDVD::CalcRMSWithFlags(sig);
    float sigThreshold = sigFactor * rmsVal;

    float ADCval;
    std::vector<bool> signalRegions;
    int numBins = sig.size();

    for (int i = 0; i < numBins; i++) {
        ADCval = sig.at(i);
        if (((ADCval > sigThreshold) || (ADCval < -1.0 * sigThreshold)) && (ADCval < 16384.0 /*4096.0*/)) {
            signalRegions.push_back(true);
        }
        else {
            signalRegions.push_back(false);
        }
    }

    for (int i = 0; i < numBins; i++) {
        if (signalRegions[i] == true) {
            int bin1 = i - padBins;
            if (bin1 < 0) {
                bin1 = 0;
            }
            int bin2 = i + padBins;
            if (bin2 > numBins) {
                bin2 = numBins;
            }

            for (int j = bin1; j < bin2; j++) {
                ADCval = sig.at(j);
                if (ADCval < 16384.0 /*4096.0*/) {
                    sig.at(j) = sig.at(j) + 200000.0/*20000.0*/;
                }
            }
        }
    }
    return true;
}

bool PDVD::RawAdapativeBaselineAlg(WireCell::Waveform::realseq_t& sig)
{
    const int windowSize = 512/*20*/;
    const int numBins = sig.size();
    int minWindowBins = windowSize / 2;

    std::vector<double> baselineVec(numBins, 0.0);
    std::vector<bool> isFilledVec(numBins, false);

    int numFlaggedBins = 0;

    for (int j = 0; j < numBins; j++) {
        if (sig.at(j) == 100000.0/*10000.0*/) {
            numFlaggedBins++;
        }
    }
    if (numFlaggedBins == numBins) {
        return true;  // Eventually replace this with flag check
    }

    double baselineVal = 0.0;
    int windowBins = 0;
    // int index;
    double ADCval = 0.0;
    for (int j = 0; j <= windowSize / 2; j++) {
        ADCval = sig.at(j);
        if (ADCval < 16384.0/*4096.0*/) {
            baselineVal += ADCval;
            windowBins++;
        }
    }

    if (windowBins == 0) {
        baselineVec[0] = 0.0;
    }
    else {
        baselineVec[0] = baselineVal / ((double) windowBins);
    }

    if (windowBins < minWindowBins) {
        isFilledVec[0] = false;
    }
    else {
        isFilledVec[0] = true;
    }

    for (int j = 1; j < numBins; j++) {
        int oldIndex = j - windowSize / 2 - 1;
        int newIndex = j + windowSize / 2;

        if (oldIndex >= 0) {
            ADCval = sig.at(oldIndex);
            if (ADCval < 16384.0/*4096.0*/) {
                baselineVal -= sig.at(oldIndex);
                windowBins--;
            }
        }
        if (newIndex < numBins) {
            ADCval = sig.at(newIndex);
            if (ADCval < 16384.0 /*4096*/) {
                baselineVal += sig.at(newIndex);
                windowBins++;
            }
        }

        if (windowBins == 0) {
            baselineVec[j] = 0.0;
        }
        else {
            baselineVec[j] = baselineVal / windowBins;
        }

        if (windowBins < minWindowBins) {
            isFilledVec[j] = false;
        }
        else {
            isFilledVec[j] = true;
        }
    }

    for (int j = 0; j < numBins; j++) {
        bool downFlag = false;
        bool upFlag = false;

        ADCval = sig.at(j);
        if (ADCval != 100000.0/*10000.0*/) {
            if (isFilledVec[j] == false) {
                int downIndex = j;
                while ((isFilledVec[downIndex] == false) && (downIndex > 0) && (sig.at(downIndex) != 100000.0/*10000.0*/)) {
                    downIndex--;
                }

                if (isFilledVec[downIndex] == false) {
                    downFlag = true;
                }

                int upIndex = j;
                while ((isFilledVec[upIndex] == false) && (upIndex < numBins - 1) && (sig.at(upIndex) != 100000.0/*10000.0*/)) {
                    upIndex++;
                }

                if (isFilledVec[upIndex] == false) {
                    upFlag = true;
                }

                if ((downFlag == false) && (upFlag == false)) {
                    baselineVec[j] = (upIndex != downIndex)
                        ? ((j - downIndex) * baselineVec[downIndex] + (upIndex - j) * baselineVec[upIndex]) /
                              ((double) upIndex - downIndex)
                        : baselineVec[downIndex];
                }
                else if ((downFlag == true) && (upFlag == false)) {
                    baselineVec[j] = baselineVec[upIndex];
                }
                else if ((downFlag == false) && (upFlag == true)) {
                    baselineVec[j] = baselineVec[downIndex];
                }
                else {
                    baselineVec[j] = 0.0;
                }
            }

            sig.at(j) = ADCval - baselineVec[j];
        }
    }

    return true;
}

bool PDVD::RemoveFilterFlags(WireCell::Waveform::realseq_t& sig)
{
    int numBins = sig.size();
    for (int i = 0; i < numBins; i++) {
        double ADCval = sig.at(i);
        if (ADCval > 16384.0/*4096.0*/) {
            if (ADCval > 100000.0/*10000.0*/) {
                sig.at(i) = ADCval - 200000.0/*20000.0*/;
            }
            else {
                sig.at(i) = 0.0;
            }
        }
    }

    return true;
}


float PDVD::CalcRMSWithFlags(const WireCell::Waveform::realseq_t& sig)
{
    float theRMS = 0.0;

    WireCell::Waveform::realseq_t temp;
    temp.reserve(sig.size());
    for (size_t i = 0; i != sig.size(); i++) {
        if (sig.at(i) < 16384.0/*4096*/) temp.push_back(sig.at(i));
    }
    float par[3];
    if (temp.size() > 0) {
        par[0] = WireCell::Waveform::percentile_binned(temp, 0.5 - 0.34);
        par[1] = WireCell::Waveform::percentile_binned(temp, 0.5);
        par[2] = WireCell::Waveform::percentile_binned(temp, 0.5 + 0.34);

        theRMS = sqrt((pow(par[2] - par[1], 2) + pow(par[1] - par[0], 2)) / 2.);
    }

    return theRMS;
}

bool PDVD::NoisyFilterAlg(WireCell::Waveform::realseq_t& sig, float min_rms, float max_rms)
{
    const double rmsVal = PDVD::CalcRMSWithFlags(sig);

    if (rmsVal > max_rms || rmsVal < min_rms) {
        int numBins = sig.size();
        for (int i = 0; i < numBins; i++) {
            sig.at(i) = 100000.0/*10000.0*/;
        }

        return true;
    }

    return false;
}

float PDVD::get_rms_and_rois(const WireCell::Waveform::realseq_t& signal, std::vector<std::vector<int> >& rois)
{

    std::pair<double, double> temp = Derivations::CalcRMS(signal);

    std::vector<int> roi;
    size_t last_bin_roi=0;
    int flag_continue=0;
    int start=1;
    int IS_signal=0;
    for (size_t j = 0; j != signal.size(); j++)
    {
        if (signal.at(j) - temp.first < -3.5 * temp.second)
        {
            IS_signal=1;

            if (start == 1) {
                last_bin_roi = j;
                start=0;
                roi.push_back(j);
            }
            else {

                if ((last_bin_roi + 1) == j) {
                    flag_continue = 1;
                    last_bin_roi = j;
                    //cout<<" 1 "<<last_bin_roi<<" j "<<j<<endl;
                }
                else {
                    flag_continue = 0;
                    last_bin_roi = j;
                    rois.push_back(roi);
                    roi.clear();
                    roi.resize(0);

                    roi.push_back(j); //start of next roi
                }

                if (flag_continue == 1) { roi.push_back(j); }
            }

        }
    }
    if (IS_signal==1) {
        rois.push_back(roi);
        roi.clear();
        roi.resize(0);
    }

    return temp.second;
}

bool PDVD::Is_FEMB_noise(const WireCell::IChannelFilter::channel_signals_t& chansig, int& beg, int& end, float min_width)
{
    // project all channels to 1D signal
    int nsignals = chansig.begin()->second.size();
    WireCell::Waveform::realseq_t signal(nsignals);
    for (const auto& cs: chansig) {
        std::transform(signal.begin(), signal.end(), cs.second.begin(), signal.begin(), std::plus<float>() );
    }

    std::vector<std::vector<int>> rois;
    /*double rms =*/ PDVD::get_rms_and_rois(signal, rois);
    for(auto roi_tmp : rois){
        double width = roi_tmp.size();
        if( width > min_width ){ // found the noise
            beg = std::max(roi_tmp[0]-20, 0); // FIXME: make it configurable? 
            end = std::min(roi_tmp.back()+20, nsignals-1);
            return true;
        }
    }

    return false;
}


/*
 * Classes
 */

PDVD::OneChannelNoise::OneChannelNoise(const std::string& anode, const std::string& noisedb)
  : m_anode_tn(anode)
  , m_noisedb_tn(noisedb)
  , m_check_partial()  // fixme, here too.
  , m_log(Log::logger("sigproc"))
{
}
PDVD::OneChannelNoise::~OneChannelNoise() {}

void PDVD::OneChannelNoise::configure(const WireCell::Configuration& cfg)
{
    m_anode_tn = get(cfg, "anode", m_anode_tn);
    m_anode = Factory::find_tn<IAnodePlane>(m_anode_tn);
    m_noisedb_tn = get(cfg, "noisedb", m_noisedb_tn);
    m_noisedb = Factory::find_tn<IChannelNoiseDatabase>(m_noisedb_tn);

    std::string dft_tn = get<std::string>(cfg, "dft", "FftwDFT");
    m_dft = Factory::find_tn<IDFT>(dft_tn);

    m_adaptive_baseline = get<bool>(cfg, "adaptive_baseline", m_adaptive_baseline);
}
WireCell::Configuration PDVD::OneChannelNoise::default_configuration() const
{
    Configuration cfg;
    cfg["anode"] = m_anode_tn;
    cfg["noisedb"] = m_noisedb_tn;
    cfg["dft"] = "FftwDFT";     // type-name for the DFT to use
    cfg["adaptive_baseline"] = false;
    return cfg;
}

WireCell::Waveform::ChannelMaskMap PDVD::OneChannelNoise::apply(int ch, signal_t& signal) const
{
    WireCell::Waveform::ChannelMaskMap ret;

    // do we need a nominal baseline correction?
    // float baseline = m_noisedb->nominal_baseline(ch);

    // correct rc undershoot
    auto spectrum = fwd_r2c(m_dft, signal);
    bool is_partial = m_adaptive_baseline ? m_check_partial(spectrum) : false;
    // bool is_partial = m_check_partial(spectrum);  // Xin's "IS_RC()"
    // if(ch>3072*2) is_partial=false;

    // if (!is_partial) {
    //     auto const& spec = m_noisedb->rcrc(ch);  // set rc_layers in chndb-xxx.jsonnet
    //     WireCell::Waveform::shrink(spectrum, spec);
    // }


    // Per-channel frequency mask from chndb freqmasks field.
    // Empty spectrum (default) = no-op; channels not listed in any freqmasks
    // channel_info entry are unaffected.
    {
        auto const& spec = m_noisedb->noise(ch);
        if (!spec.empty()) {
            if (spec.size() == spectrum.size()) {
                WireCell::Waveform::scale(spectrum, spec);
            }
            else {
                m_log->warn("PDVD OneChannelNoise ch={}: freqmask size {} != FFT size {}, "
                            "mask SKIPPED — check params.nf.nsamples vs actual frame nticks",
                            ch, spec.size(), spectrum.size());
            }
        }
    }

    // remove the DC component
    spectrum.front() = 0;
    signal = inv_c2r(m_dft, spectrum);

    // Now calculate the baseline ...
    std::pair<double, double> temp = WireCell::Waveform::mean_rms(signal);
    auto temp_signal = signal;
    for (size_t i = 0; i != temp_signal.size(); i++) {
        if (fabs(temp_signal.at(i) - temp.first) > 6 * temp.second) {
            temp_signal.at(i) = temp.first;
        }
    }
    float baseline = WireCell::Waveform::median_binned(temp_signal);
    // correct baseline
    WireCell::Waveform::increase(signal, baseline * (-1));

    // Now do the adaptive baseline
    if (is_partial) {
    //    std::cout << "[PDVD] is_partical channel: " << ch << std::endl;

        auto wpid = m_anode->resolve(ch);
        const int iplane = wpid.index();
        // add something
        WireCell::Waveform::BinRange temp_bin_range;
        temp_bin_range.first = 0;
        temp_bin_range.second = signal.size();

        if (iplane != 2) {  // not collection
            ret["lf_noisy"][ch].push_back(temp_bin_range);
            // std::cout << "Partial " << ch << std::endl;
        }
        PDVD::SignalFilter(signal);
        PDVD::RawAdapativeBaselineAlg(signal);
        PDVD::RemoveFilterFlags(signal);
    }

    const float min_rms = m_noisedb->min_rms_cut(ch);
    const float max_rms = m_noisedb->max_rms_cut(ch);
    PDVD::SignalFilter(signal);
    const double rms_val = PDVD::CalcRMSWithFlags(signal);
    bool is_noisy = PDVD::NoisyFilterAlg(signal, min_rms, max_rms);
    PDVD::RemoveFilterFlags(signal);
    if (is_noisy) {
        WireCell::Waveform::BinRange temp_bin_range;
        temp_bin_range.first = 0;
        temp_bin_range.second = signal.size();
        ret["noisy"][ch].push_back(temp_bin_range);
    }
    return ret;
}

WireCell::Waveform::ChannelMaskMap PDVD::OneChannelNoise::apply(channel_signals_t& chansig) const
{
    return WireCell::Waveform::ChannelMaskMap();
}

PDVD::CoherentNoiseSub::CoherentNoiseSub(const std::string& anode, const std::string& noisedb,
                                               float rms_threshold)
  : m_anode_tn(anode)
  , m_noisedb_tn(noisedb)
  , m_rms_threshold(rms_threshold)
{
}
PDVD::CoherentNoiseSub::~CoherentNoiseSub() {}

WireCell::Waveform::ChannelMaskMap PDVD::CoherentNoiseSub::apply(channel_signals_t& chansig) const
{
    // std::cout << "grouped channel size: " << " " << chansig.size() << std::endl;
    // find the median among all
    WireCell::Waveform::realseq_t medians = Derivations::CalcMedian(chansig);



    // std::cout << medians.size() << " " << medians.at(100) << " " << medians.at(101) << std::endl;

    // For Xin: here is how you can get the response spectrum for this group.
    const int achannel = chansig.begin()->first;

    // std::cerr << "CoherentNoiseSub: ch=" << achannel << " response offset:" << m_noisedb->response_offset(achannel)
    // << std::endl;

    const Waveform::compseq_t& respec = m_noisedb->response(achannel);
    const int res_offset = m_noisedb->response_offset(achannel);
    const int pad_f = m_noisedb->pad_window_front(achannel);
    const int pad_b = m_noisedb->pad_window_back(achannel);

    // need to move these to data base, consult with Brett ...
    // also need to be time dependent ...
    const float decon_limit = m_noisedb->coherent_nf_decon_limit(achannel);  // 0.02;
    const float adc_limit = m_noisedb->coherent_nf_adc_limit(achannel);                  // 15;
    const float decon_limit1 = m_noisedb->coherent_nf_decon_limit1(achannel);            // 0.08; // loose filter
    const float roi_min_max_ratio = m_noisedb->coherent_nf_roi_min_max_ratio(achannel);  // 0.8 default

    const float protection_factor = m_noisedb->coherent_nf_protection_factor(achannel);
    const float min_adc_limit = m_noisedb->coherent_nf_min_adc_limit(achannel);

    const int plane = m_anode->resolve(achannel).index();
    const int nfbins = respec.size();
    auto time_filter_wf = Factory::find<IFilterWaveform>("HfFilter", m_time_filters.at(plane))->filter_waveform(nfbins);
    auto lf_tighter_wf = Factory::find<IFilterWaveform>("LfFilter", m_lf_tighter_filter)->filter_waveform(nfbins);
    auto lf_loose_wf = Factory::find<IFilterWaveform>("LfFilter", m_lf_loose_filter)->filter_waveform(nfbins);

    // Opt-in coherent-NF dump. Off-state cost: one .empty() check; the
    // DumpRecord is never constructed and SP/SW receive a nullptr.
    SigProc::CoherentNoiseDump dump_rec;
    SigProc::CoherentNoiseDump* dump_ptr = nullptr;
    if (!m_dump_path.empty()
        && (m_dump_groups.empty() || m_dump_groups.count(achannel))) {
        dump_ptr = &dump_rec;
        dump_rec.apa = m_anode->ident();
        dump_rec.gid = achannel;
        dump_rec.plane = plane;
        dump_rec.nbin = medians.size();
        dump_rec.res_offset = res_offset;
        dump_rec.channels.reserve(chansig.size());
        for (auto& it : chansig) dump_rec.channels.push_back(it.first);
        dump_rec.decon_limit = decon_limit;
        dump_rec.decon_limit1 = decon_limit1;
        dump_rec.roi_min_max_ratio = roi_min_max_ratio;
        dump_rec.min_adc_limit = min_adc_limit;
        dump_rec.upper_adc_limit = adc_limit;
        dump_rec.upper_decon_limit = decon_limit;
        dump_rec.protection_factor = protection_factor;
        dump_rec.pad_front = pad_f;
        dump_rec.pad_back = pad_b;
        dump_rec.time_filter_name = m_time_filters.at(plane);
        dump_rec.lf_tighter_filter_name = m_lf_tighter_filter;
        dump_rec.lf_loose_filter_name = m_lf_loose_filter;
        dump_rec.median = medians;
    }

    // do the signal protection and adaptive baseline
    std::vector<std::vector<int> > rois =
        PDVD::SignalProtection(medians, respec, m_dft,
                                     res_offset, pad_f, pad_b, time_filter_wf, lf_tighter_wf,
                                     decon_limit, adc_limit,
                                     protection_factor, min_adc_limit,
                                     dump_ptr);

    // if (achannel == 3840){
    // 	std::cout << "Xin1: " << rois.size() << std::endl;
    	// for (size_t i=0;i!=rois.size();i++){
    	//     std::cout << "ROI in SignalProtection: " << rois.at(i).front() << " " << rois.at(i).back() << std::endl;
    	// }
    // }

    // std::cerr <<"\tSigprotection done: " << chansig.size() << " " << medians.size() << " " << medians.at(100) << " "
    // << medians.at(101) << std::endl;

    // if(achannel==9438||achannel==9116 || achannel==7160 || achannel==4182 ){
    // if(true){
    // std::cout<<"print out median: "<<achannel<<std::endl;
    // std::string filename = "medians_output"+ std::to_string(achannel)+".txt";
    // std::ofstream outfile(filename);
    
    // for (size_t i = 0; i < medians.size(); ++i) {
    //     outfile << i << " " << medians[i]<<std::endl;
    // }
    
    // outfile.close();
    // int k=0;
    // for (auto it : chansig) {
    // if(k<10){
    //     std::string filename2 = "medians_output"+ std::to_string(achannel)+"_"+std::to_string(k)+".txt";
    //     std::ofstream outfile2(filename2);
    
    //     for (size_t i = 0; i < medians.size(); ++i) {
    //         outfile2 << i << " " << it.second.at(i)<<std::endl;
    //     }
    
    //     outfile2.close();
    //     k++;
    // }

    // }

    // }

    // calculate the scaling coefficient and subtract
    PDVD::Subtract_WScaling(chansig, medians, respec, res_offset, rois,
                                  m_dft, time_filter_wf, lf_loose_wf,
                                  decon_limit1, roi_min_max_ratio,
                                  m_rms_threshold,
                                  dump_ptr);

    if (dump_ptr) {
        SigProc::CoherentNoiseDumpWriter::write(m_dump_path, *dump_ptr);
    }

    // WireCell::IChannelFilter::signal_t& signal = chansig.begin()->second;
    // for (size_t i=0;i!=signal.size();i++){
    // 	signal.at(i) = medians.at(i);
    // }

    // std::cerr <<"\tSubtrace_WScaling done" << std::endl;

    // for (auto it: chansig){
    // 	std::cout << "Xin3 " << it.first << std::endl;
    // 	break;
    // }

    return WireCell::Waveform::ChannelMaskMap();  // not implemented
}
WireCell::Waveform::ChannelMaskMap PDVD::CoherentNoiseSub::apply(int channel, signal_t& sig) const
{
    return WireCell::Waveform::ChannelMaskMap();  // not implemented
}

void PDVD::CoherentNoiseSub::configure(const WireCell::Configuration& cfg)
{
    m_anode_tn = get(cfg, "anode", m_anode_tn);
    m_anode = Factory::find_tn<IAnodePlane>(m_anode_tn);
    m_noisedb_tn = get(cfg, "noisedb", m_noisedb_tn);
    m_noisedb = Factory::find_tn<IChannelNoiseDatabase>(m_noisedb_tn);

    std::string dft_tn = get<std::string>(cfg, "dft", "FftwDFT");
    m_dft = Factory::find_tn<IDFT>(dft_tn);

    m_rms_threshold = get<float>(cfg, "rms_threshold", m_rms_threshold);

    m_time_filters.clear();
    for (auto jtf : cfg["time_filters"]) {
        m_time_filters.push_back(jtf.asString());
    }
    m_lf_tighter_filter = get<std::string>(cfg, "lf_tighter_filter", m_lf_tighter_filter);
    m_lf_loose_filter = get<std::string>(cfg, "lf_loose_filter", m_lf_loose_filter);

    m_dump_path = get<std::string>(cfg, "debug_dump_path", m_dump_path);
    m_dump_groups.clear();
    for (auto jg : cfg["debug_dump_groups"]) {
        m_dump_groups.insert(jg.asInt());
    }
}
WireCell::Configuration PDVD::CoherentNoiseSub::default_configuration() const
{
    Configuration cfg;
    cfg["anode"] = m_anode_tn;
    cfg["noisedb"] = m_noisedb_tn;
    cfg["dft"] = "FftwDFT";     // type-name for the DFT to use

    cfg["rms_threshold"] = m_rms_threshold;
    cfg["time_filters"] = Json::arrayValue;
    for (const auto& f : m_time_filters) {
        cfg["time_filters"].append(f);
    }
    cfg["lf_tighter_filter"] = m_lf_tighter_filter;
    cfg["lf_loose_filter"] = m_lf_loose_filter;

    cfg["debug_dump_path"] = m_dump_path;
    cfg["debug_dump_groups"] = Json::arrayValue;

    return cfg;
}

PDVD::FEMBNoiseSub::FEMBNoiseSub(const std::string& anode, float width)
  : m_anode_tn(anode)
  , m_width(width)
{
}
PDVD::FEMBNoiseSub::~FEMBNoiseSub() {}

WireCell::Waveform::ChannelMaskMap PDVD::FEMBNoiseSub::apply(channel_signals_t& chansig) const
{
    WireCell::Waveform::ChannelMaskMap ret;

    // WireCell::Waveform::realseq_t medians = Derivations::CalcMedian(chansig);
    // const int achannel = chansig.begin()->first;
    // std::cout << "[wgu] PDVD::FEMBNoiseSub::apply first channel: " << achannel << std::endl;

    // determine if FEMB negative pulse
    WireCell::Waveform::BinRange fembnoise_bins;
    bool is_femb_noise = Is_FEMB_noise(chansig, fembnoise_bins.first, fembnoise_bins.second, m_width);
    if (is_femb_noise) {
        for (auto const& cs : chansig) {
            ret["femb_noise"][cs.first].push_back(fembnoise_bins);
            // std::cout << "[wgu] FEMB Noise channel= " << cs.first << " , time bins: "
            // << fembnoise_bins.first << " " << fembnoise_bins.second << std::endl;
        }
    }

    return ret;
}
WireCell::Waveform::ChannelMaskMap PDVD::FEMBNoiseSub::apply(int channel, signal_t& sig) const
{
    return WireCell::Waveform::ChannelMaskMap();  // not implemented
}

void PDVD::FEMBNoiseSub::configure(const WireCell::Configuration& cfg)
{
    m_anode_tn = get(cfg, "anode", m_anode_tn);
    m_anode = Factory::find_tn<IAnodePlane>(m_anode_tn);

    m_width = get<float>(cfg, "width", m_width);
}
WireCell::Configuration PDVD::FEMBNoiseSub::default_configuration() const
{
    Configuration cfg;
    cfg["anode"] = m_anode_tn;
    cfg["width"] = m_width;

    return cfg;
}
bool PDVD::Signal_mask_top_u(WireCell::Waveform::realseq_t& sig)
{
    // specific for TDE u plane shield plane noise filter
    // protect according to the positive part of the signal, cause the noise is pure negative.

    const double sigFactor = 4.0;
    const int padBins = 70;

    float rmsVal = PDVD::CalcRMSWithFlags(sig);
    float sigThreshold = sigFactor * rmsVal;

    float ADCval;
    std::vector<bool> signalRegions;
    int numBins = sig.size();

    for (int i = 0; i < numBins; i++) {
        ADCval = sig.at(i);
        if (((ADCval > sigThreshold) ) && (ADCval < 16384.0 /*4096.0*/)) {
            signalRegions.push_back(true);
        }
        else {
            signalRegions.push_back(false);
        }
    }

    for (int i = 0; i < numBins; i++) {
        if (signalRegions[i] == true) {
            int bin1 = i - padBins;
            if (bin1 < 0) {
                bin1 = 0;
            }
            int bin2 = i + padBins;
            if (bin2 > numBins) {
                bin2 = numBins;
            }

            for (int j = bin1; j < bin2; j++) {
                ADCval = sig.at(j);
                if (ADCval < 16384.0 /*4096.0*/) {
                    sig.at(j) = sig.at(j) + 200000.0/*20000.0*/;
                    // sig.at(j) = 0;
                }
            }
        }
    }
    return true;
}
WireCell::Waveform::realseq_t PDVD::CalcMedian_shieldCoupling_u(const WireCell::IChannelFilter::channel_signals_t& chansig){

    float max_rms = 0;
    float count_max_rms = 0;
    const int nchannel = chansig.size();
    const int nbins = (chansig.begin()->second).size();
    // float content[nchannel][nbins];
    std::vector<float> content(nchannel * nbins, 0.0);  // 2D array [channel*nbins + bin]

    int start_ch = 0;
    for (auto it : chansig) {
        // const int ch = it.first;
        WireCell::IChannelFilter::signal_t& signal = it.second;
        WireCell::IChannelFilter::signal_t filtered_signal;
        for (const auto& value : signal) {
            if (value <= 100000) {
                    filtered_signal.push_back(value);
            }
        } 
        std::pair<double, double> temp = WireCell::Waveform::mean_rms(filtered_signal);

        if (temp.second > 0) {
            max_rms += temp.second;
            count_max_rms++;
        }

        for (int i = 0; i != nbins; i++) {
            content[start_ch * nbins + i] = signal.at(i);
        }
        start_ch++;
    }

    if (count_max_rms > 0) {
        max_rms /= count_max_rms;
    }

    WireCell::Waveform::realseq_t medians(nbins);
    for (int ibin = 0; ibin != nbins; ibin++) {
        WireCell::Waveform::realseq_t temp;
        for (int ich = 0; ich != nchannel; ich++) {
            const float cont = content.at(ich * nbins + ibin);
            if (cont < 5 * max_rms && fabs(cont) > 0.001) {
                temp.push_back(cont);
                // std::cout<< cont<<"\t";
            }

        }
        if (temp.size() > 0) {
            medians.at(ibin) = WireCell::Waveform::median_binned(temp);
            // if(medians.at(ibin)<-18){
            //     std::cout<<"median array size: "<<temp.size()<<std::endl;
            //     std::cout<<"max_rms = "<<max_rms<<std::endl; 
            //     std::cout<<"nchannel = "<<nchannel<<std::endl; 
            //     std::cout<<"median = "<<medians.at(ibin)<<std::endl; 
            //     for (int ich = 0; ich != temp.size(); ich++) {
            //         std::cout<<temp[ich]<<"\t";
            //     }
            //     std::cout<<std::endl;
            // }
        }
        else {
            medians.at(ibin) = 0;
        }
    }

    return medians;
}

/*
 * Shield Coupling Removal Class
 */

PDVD::ShieldCouplingSub::ShieldCouplingSub(const std::string& anode, const std::string& noisedb,
                                               float rms_threshold)
  : m_anode_tn(anode)
  , m_noisedb_tn(noisedb)
  , m_rms_threshold(rms_threshold)
{
}

PDVD::ShieldCouplingSub::~ShieldCouplingSub() {}

void PDVD::ShieldCouplingSub::configure(const WireCell::Configuration& cfg)
{
    m_anode_tn = get(cfg, "anode", m_anode_tn);
    m_anode = Factory::find_tn<IAnodePlane>(m_anode_tn);
    m_noisedb_tn = get(cfg, "noisedb", m_noisedb_tn);
    m_noisedb = Factory::find_tn<IChannelNoiseDatabase>(m_noisedb_tn);

    auto strip_file = get<std::string>(cfg, "strip_length", "");
    if (!strip_file.empty()) {
        auto jdata = Persist::load(strip_file);  // Auto-handles .bz2
        for (const auto& entry : jdata) {
            int ch = entry["channel"].asInt();
            float len = entry["length"].asFloat();
            m_strip_lengths[ch] = len;
        }
    }
    std::string dft_tn = get<std::string>(cfg, "dft", "FftwDFT");
    m_dft = Factory::find_tn<IDFT>(dft_tn);

    m_rms_threshold = get<float>(cfg, "rms_threshold", m_rms_threshold);
    // m_capa_weight = get<bool>(cfg, "capa_weight", m_capa_weight);
    // m_calibrated = get<bool>(cfg, "calibrated", m_calibrated);
    // m_group_size = get<int>(cfg, "group_size", m_group_size);
    // m_min_channels = get<int>(cfg, "min_channels", m_min_channels);
}

WireCell::Configuration PDVD::ShieldCouplingSub::default_configuration() const
{
    Configuration cfg;
    cfg["anode"] = m_anode_tn;
    cfg["noisedb"] = m_noisedb_tn;
    // cfg["capa_weight"] = m_capa_weight;
    // cfg["calibrated"] = m_calibrated;
    // cfg["group_size"] = m_group_size;
    // cfg["min_channels"] = m_min_channels;
    return cfg;
}

WireCell::Waveform::ChannelMaskMap PDVD::ShieldCouplingSub::apply(
    channel_signals_t& chansig) const
{
    WireCell::Waveform::ChannelMaskMap ret;
    
    if (chansig.empty()) {
        return ret;
    }
    // std::cout << "grouped channel size: " << " " << chansig.size() << std::endl;
    
    const int nbins = (chansig.begin()->second).size();

    std::map<int, float> scale_factors;
    for (auto& cs : chansig) {
        const int ch = cs.first;
        auto& signal = cs.second;
        
        float strip_length = 1.0;
        auto it = m_strip_lengths.find(ch);
        if (it != m_strip_lengths.end() && it->second > 0) {
            strip_length = it->second;
        }else{
            std::cerr<<"error!! strip length for ch "<<ch<<" not found"<< std::endl;
        }
        
        scale_factors[ch] = strip_length;
        
        // Scale DOWN: divide by strip length (like calib/capa in Lardon)
        for (int ibin = 0; ibin < nbins; ibin++) {
            signal.at(ibin) /= strip_length;
        }
        
        Signal_mask_top_u(signal);
    }


    WireCell::Waveform::realseq_t medians = PDVD::CalcMedian_shieldCoupling_u(chansig);

    // const int achannel = chansig.begin()->first;
   
    
    // if(false){
    //     std::cout<<"print out median: "<<achannel<<std::endl;
    // std::string filename = "medians_output_tu"+ std::to_string(achannel)+".txt";
    // std::ofstream outfile(filename);
    
    // for (size_t i = 0; i < medians.size(); ++i) {
    //     outfile << i << " " << medians[i]<<std::endl;
    // }
    
    // outfile.close();
    // int k=0;
    // for (auto it : chansig) {
    // if(it.first % 40==0){
    //     std::string filename2 = "medians_output_tu"+ std::to_string(achannel)+"_"+std::to_string(k)+".txt";
    //     std::ofstream outfile2(filename2);
    
    //     for (size_t i = 0; i < medians.size(); ++i) {
    //         outfile2 << i << " " << it.second.at(i)<<std::endl;
    //     }
    
    //     outfile2.close();
    //     k++;
    // }

    // }

    // }


    for (auto it : chansig) {
        int ch = it.first;
        WireCell::IChannelFilter::signal_t& signal = it.second;
        PDVD::RemoveFilterFlags(signal); 

        int nbin = signal.size();
        float scaling=1;
                for (int i = 0; i != nbin; i++) {
                if (fabs(signal.at(i)) > 0.001) {
                    signal.at(i) = signal.at(i) - medians.at(i) * scaling;
                }
            }

         chansig[ch] = signal;
    }

    for (auto& cs : chansig) {
        const int ch = cs.first;
        auto& signal = cs.second;

        float strip_length = scale_factors[ch];
        if (strip_length > 0) {
            for (int ibin = 0; ibin < nbins; ibin++) {
                signal.at(ibin) *= strip_length;
            }
        }
    }

  return ret;
}

WireCell::Waveform::ChannelMaskMap PDVD::ShieldCouplingSub::apply(
    int channel, signal_t& sig) const
{
    return WireCell::Waveform::ChannelMaskMap();
}



// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
