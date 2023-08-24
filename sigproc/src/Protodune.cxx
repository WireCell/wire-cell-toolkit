
/***************************************************************************/
/*             Noise filter (NF) module for protoDUNE-SP                   */
/***************************************************************************/
/* features:                                                               */
/* - frequency-domain resampling (sticky code and FEMB time clock)         */
/* - undershoot correction                                                 */
/* - partial undershoot correction                                         */
/* - 50 kHz collection plane noise filter                                  */
/* - ledge waveform identification (multiple ledges)                       */

#include "WireCellSigProc/Microboone.h"
#include "WireCellSigProc/Protodune.h"
#include "WireCellSigProc/Derivations.h"

#include "WireCellAux/DftTools.h"

#include "WireCellUtil/NamedFactory.h"

#include <cmath>
#include <complex>
#include <iostream>
#include <set>

WIRECELL_FACTORY(pdStickyCodeMitig, WireCell::SigProc::Protodune::StickyCodeMitig, WireCell::IChannelFilter,
                 WireCell::IConfigurable)
WIRECELL_FACTORY(pdOneChannelNoise, WireCell::SigProc::Protodune::OneChannelNoise, WireCell::IChannelFilter,
                 WireCell::IConfigurable)
WIRECELL_FACTORY(pdRelGainCalib, WireCell::SigProc::Protodune::RelGainCalib, WireCell::IChannelFilter,
                 WireCell::IConfigurable)

using namespace WireCell::SigProc;
using WireCell::Aux::DftTools::fwd_r2c;
using WireCell::Aux::DftTools::inv_c2r;


int LedgeIdentify1(WireCell::Waveform::realseq_t& signal, double baseline, int LedgeStart[3], int LedgeEnd[3])
{
    // find a maximum of 3 ledges in one waveform
    // the number of ledges is returned
    // the start and end of ledge is stored in the array
    int UNIT = 10;    // 5;    // rebin unit
    int CONTIN = 10;  // 20; // length of the continuous region
    int JUMPS = 2;    // 4;   // how many bins can accidental jump
    // We recommend UNIT*CONTIN = 100, UNIT*JUMPS = 20
    std::vector<int> averaged;  // store the rebinned waveform
    int up = signal.size() / UNIT;
    int nticks = signal.size();
    // rebin
    for (int i = 0; i < up; i++) {
        int temp = 0;
        for (int j = 0; j < UNIT; j++) {
            temp += signal.at(i * UNIT + j);
        }
        averaged.push_back(temp);
    }
    // relax the selection cuts if there is a large signal
    auto imax = std::min_element(signal.begin(), signal.end());
    double max_value = *imax;
    if (max_value - baseline > 1000) {
        CONTIN /= 1.25;
        JUMPS *= 1.2;
    }
    // start judging
    int NumberOfLedge = 0;                                    // to be returned
    int StartOfLastLedgeCandidate = 0;                        // store the last ledge candidate
    for (int LE = 0; LE < 3; LE++) {                          // find three ledges
        if (LE > 0 && StartOfLastLedgeCandidate == 0) break;  // no ledge candidate in the last round search
        if (StartOfLastLedgeCandidate > nticks - 200) break;  // the last candidate is too late in the window
        if (NumberOfLedge > 0 && (LedgeEnd[NumberOfLedge - 1] + 200) > nticks)
            break;  // the end last ledge reaches the end of readout window
        int ledge = 0;
        int decreaseD = 0, tolerence = 0;
        int start = 0;        // end = 0; // temporary start and end
        int StartWindow = 0;  // where to start
        if (StartOfLastLedgeCandidate == 0) {
            StartWindow = 0;
        }
        else {
            if (NumberOfLedge == 0)  // no ledge found. start from 200 ticks after the last candidate
                StartWindow = (StartOfLastLedgeCandidate + 200) / UNIT;
            else {
                if (StartOfLastLedgeCandidate > LedgeEnd[NumberOfLedge - 1])  // the last candidate is not a real ledge
                    // StartWindow = (StartOfLastLedgeCandidate+200)/UNIT;
                    StartWindow = (StartOfLastLedgeCandidate + 50) / UNIT;
                else  // the last candidate is a real ledge
                    StartWindow = (LedgeEnd[NumberOfLedge - 1] + 30) / UNIT;
            }
        }
        for (int i = StartWindow + 1; i < up - 1; i++) {
            if (averaged.at(i) < averaged.at(i - 1)) {
                if (decreaseD == 0) start = i;
                decreaseD += 1;
            }
            else {
                if (averaged.at(i + 1) < averaged.at(i - 1) && tolerence < JUMPS &&
                    decreaseD > 0) {  // we can ignore several jumps in the decreasing region
                    decreaseD += 2;
                    tolerence++;
                    i = i + 1;
                }
                else {
                    if (decreaseD > CONTIN) {
                        ledge = 1;
                        StartOfLastLedgeCandidate = start * UNIT;
                        break;
                    }
                    else {
                        decreaseD = 0;
                        tolerence = 0;
                        start = 0;
                        // end = 0;
                    }
                }
            }
        }
        // find the sharp start edge
        // if(ledge == 1&&LedgeStart>30){
        //   int edge = 0;
        //   int i = LedgeStart/UNIT-1;
        //   if(averaged.at(i)>averaged.at(i-1)&&averaged.at(i-1)>averaged.at(i-2)){ // find a edge
        //           edge = 1;
        //   }
        //   if(edge == 0) ledge = 0; // if no edge, this is not ledge
        //   if((averaged.at(i)-averaged.at(i-2)<10*UNIT)&&(averaged.at(i)-averaged.at(i-3)<10*UNIT)) // slope cut
        //           ledge = 0;
        //   if(averaged.at(LedgeStart/UNIT)-baseline*UNIT>300*UNIT) ledge = 0; // ledge is close to the baseline
        // }
        // determine the end of ledge
        // case 1: find a jump of 5 ADC in the rebinned waveform
        // case 2: a continuous 20 ticks has an average close to baseline, and a RMS larger than 3 ADC
        // case 3: reaching the tick 6000 but neither 1 nor 2 occurs
        int tempLedgeEnd = 0;
        if (ledge == 1) {
            for (int i = StartOfLastLedgeCandidate / UNIT; i < up - 1; i++) {  // case 1
                if (averaged.at(i + 1) - averaged.at(i) > 5 * UNIT) {
                    tempLedgeEnd = i * UNIT + 5;
                    if (tempLedgeEnd > nticks) tempLedgeEnd = nticks - 1;
                    break;
                }
            }

            if (tempLedgeEnd == 0) {  // not find a jump, case 2
                WireCell::Waveform::realseq_t tempA(20);
                for (int i = StartOfLastLedgeCandidate + 80; i < nticks - 20; i += 20) {
                    for (int j = i; j < 20 + i; j++) {
                        tempA.at(j - i) = signal.at(j);
                    }
                    auto stat = WireCell::Waveform::mean_rms(tempA);
                    if (stat.first - baseline < 2 && stat.second > 3) {
                        tempLedgeEnd = i;
                        break;
                    }
                }
            }
            if (tempLedgeEnd == 0) tempLedgeEnd = nticks - 1;
        }
        // test the decay time
        // if(ledge == 1&&StartOfLastLedgeCandidate>20){
        if (ledge == 1) {
            double height = signal.at(StartOfLastLedgeCandidate + 1) - signal.at(tempLedgeEnd);
            if (height < 0) ledge = 0;                              // not a ledge
            if ((tempLedgeEnd - StartOfLastLedgeCandidate) > 80) {  // test the decay if the ledge length is longenough
                double height50 = 0;
                height50 = signal.at(StartOfLastLedgeCandidate + 51);
                double height50Pre = signal.at(StartOfLastLedgeCandidate + 1) -
                                     height * (1 - exp(-50 / 100.));  // minimum 100 ticks decay time
                // if the measured is much smaller than the predicted, this is not ledge
                if (height50 - height50Pre < -12 /*-8*/) ledge = 0;
                // cout << "height50-height50Pre: " << height50-height50Pre << endl;
            }
        }

        if (ledge == 1) {  // ledge is close to the baseline
            if (averaged.at(tempLedgeEnd / UNIT) - baseline * UNIT > 5. * UNIT) ledge = 0;

            if (averaged.at(StartOfLastLedgeCandidate / UNIT) - baseline * UNIT > 100 * UNIT) ledge = 0;
        }

        if (ledge == 1 && StartOfLastLedgeCandidate > 1000) {  // ledge always follows a pulse
            // std::cerr << "[wgu] ch: " << m_ch << " StartOfLastLedgeCandidate: " << StartOfLastLedgeCandidate <<
            // std::endl;
            float hmax = 0;
            for (int i = StartOfLastLedgeCandidate - 1000; i < StartOfLastLedgeCandidate - 100; i++) {
                if (signal.at(i) > hmax) hmax = signal.at(i);
            }
            if (hmax - baseline < 200) ledge = 0;  // no large pre-signal
            // std::cerr << "[wgu] hmax: " << hmax << std::endl;
        }

        if (ledge == 1) {
            LedgeStart[NumberOfLedge] = std::max(StartOfLastLedgeCandidate - 20, 0);
            LedgeEnd[NumberOfLedge] = std::min(tempLedgeEnd + 10, (int) signal.size());  // FIXME: ends at a threshold
            NumberOfLedge++;
        }
    }  // LE

    return NumberOfLedge;
}

// adapted from WCP
// FIXME: some hardcoded 6000 ticks
bool LedgeIdentify(WireCell::Waveform::realseq_t& signal /*TH1F* h2*/, double baseline, int& LedgeStart, int& LedgeEnd)
{
    int ledge = false;
    int UNIT = 10;                  // 5;    // rebin unit
    int CONTIN = 10;                // 20; // length of the continuous region
    int JUMPS = 2;                  // 4;   // how many bins can accidental jump
    std::vector<int> averaged;      // store the rebinned waveform
    int up = signal.size() / UNIT;  // h2->GetNbinsX()/UNIT;
    int nticks = signal.size();     // h2->GetNbinsX();
    // rebin
    for (int i = 0; i < up; i++) {
        int temp = 0;
        for (int j = 0; j < UNIT; j++) {
            temp += signal.at(i * UNIT + j);  // h2->GetBinContent(i*UNIT+1+j);
        }
        averaged.push_back(temp);
    }
    // refine the selection cuts if there is a large signal
    auto imax = std::min_element(signal.begin(), signal.end());
    double max_value = *imax;

    // if(h2->GetMaximum()-baseline>1000) { CONTIN = 16; JUMPS = 5; }
    if (max_value - baseline > 1000) {
        CONTIN = 16;
        JUMPS = 5;
    }
    // start judging
    int decreaseD = 0, tolerence = 0;
    int start = 0;
    // int end = 0; // never used?
    for (int i = 1; i < up - 1; i++) {
        if (averaged.at(i) < averaged.at(i - 1)) {
            if (decreaseD == 0) start = i;
            decreaseD += 1;
        }
        else {
            if (averaged.at(i + 1) < averaged.at(i - 1) && tolerence < JUMPS &&
                decreaseD > 0) {  // we can ignore several jumps in the decreasing region
                decreaseD += 2;
                tolerence++;
                i = i + 1;
            }
            else {
                if (decreaseD > CONTIN) {
                    ledge = true;
                    LedgeStart = start * UNIT;
                    break;
                }
                else {
                    decreaseD = 0;
                    tolerence = 0;
                    start = 0;
                    // end = 0;// end never used?
                }
            }
        }
    }
    // find the sharp start edge
    if (ledge && LedgeStart > 30) {
        if (averaged.at(LedgeStart / UNIT) - baseline * UNIT > 150 * UNIT)
            ledge = false;  // ledge is close to the baseline
    }
    // test the decay time
    if (ledge && LedgeStart > 20) {
        double height = 0;
        if (LedgeStart < 5750) {  // calculate the height of edge
            double tempHeight = signal.at(LedgeStart + 200) + signal.at(LedgeStart + 220) +
                                signal.at(LedgeStart + 180) + signal.at(LedgeStart + 240);
            height = signal.at(LedgeStart) - tempHeight / 4;
            height /= 0.7;
        }
        // else height =  h2 ->GetBinContent(LedgeStart+1) -  h2 ->GetBinContent(6000);
        else
            height = signal.at(LedgeStart) - signal.back();
        if (height < 0) height = 80;             // norminal value
        if (height > 30 && LedgeStart < 5900) {  // test the decay with a relatively large height
            double height50 = 0, height100 = 0;

            height50 = signal.at(LedgeStart + 50);
            height100 = signal.at(LedgeStart + 100);
            double height50Pre =
                signal.at(LedgeStart) - height * (1 - std::exp(-50 / 100.));  // minimum 100 ticks decay time
            double height100Pre =
                signal.at(LedgeStart) - height * (1 - std::exp(-100. / 100));  // minimum 100 ticks decay time
            // if the measured is much smaller than the predicted, this is not ledge
            if (height50 - height50Pre < -8) ledge = false;
            if (height100 - height100Pre < -8) ledge = false;
        }
    }

    // determine the end of ledge
    // case 1: find a jump of 10 ADC in the rebinned waveform
    // case 2: a continuous 20 ticks has an average close to baseline, and a RMS larger than 3 ADC
    // case 3: reaching the tick 6000 but neither 1 nor 2 occurs
    if (ledge) {
        LedgeEnd = 0;
        for (int i = LedgeStart / UNIT; i < up - 1; i++) {  // case 1
            if (averaged.at(i + 1) - averaged.at(i) > 50) {
                LedgeEnd = i * UNIT + 5;
                break;
            }
        }
        if (LedgeEnd == 0) {  // not find a jump, case 2
            // double tempA[20];
            WireCell::Waveform::realseq_t tempA(20);
            for (int i = LedgeStart + 80; i < nticks - 20; i += 20) {
                for (int j = i; j < 20 + i; j++) {
                    // tempA[j-i] = h2->GetBinContent(j+1);
                    tempA.at(j - i) = signal.at(j);
                }
                auto wfinfo = WireCell::Waveform::mean_rms(tempA);
                // if(TMath::Mean(20,tempA)-baseline<2&&TMath::RMS(20,tempA)>3){
                if (wfinfo.first - baseline < 2 && wfinfo.second > 3) {
                    LedgeEnd = i;
                    break;
                }
            }
        }
        if (LedgeEnd == 0) LedgeEnd = 6000;
    }
    return ledge;
}

bool Protodune::LinearInterpSticky(WireCell::Waveform::realseq_t& signal, WireCell::Waveform::BinRangeList& rng_list,
                                   float stky_sig_like_val, float stky_sig_like_rms)
{
    WireCell::Waveform::BinRangeList trivial_range_list;  // some ranges will be removed
                                                          // from sticky ranges (rng_list)

    const int nsiglen = signal.size();
    std::pair<double, double> temp = WireCell::Waveform::mean_rms(signal);
    for (auto const& rng : rng_list) {
        int start = rng.first - 1;
        int end = rng.second + 1;
        if (start >= 0 && end < nsiglen) {
            std::vector<float> digits;
            for (int i = start + 1; i < end; i++) {
                digits.push_back(signal.at(i));
            }

            auto min = std::max_element(digits.begin(), digits.end());
            auto max = std::min_element(digits.begin(), digits.end());
            double max_value = *max;
            double min_value = *min;

            double start_content = signal.at(start);
            double end_content = signal.at(end);
            bool isSignalLike = false;
            if (rng.second == rng.first)
                isSignalLike = true;  // single sticky, do FFT interp later
            else {
                // peak > stky_sig_like_val [unit in adc]
                // two ends > stky_sig_like_rms [unit in rms]
                if (max_value - temp.first > stky_sig_like_val) {
                    if ((start_content - temp.first > stky_sig_like_rms * temp.second) &&
                        (end_content - temp.first > stky_sig_like_rms * temp.second)) {
                        isSignalLike = true;
                    }
                }
                else if (temp.first - min_value > stky_sig_like_val) {
                    if ((temp.first - start_content > stky_sig_like_rms * temp.second) &&
                        (temp.first - end_content > stky_sig_like_rms * temp.second)) {
                        isSignalLike = true;
                    }
                }
            }

            if (!isSignalLike) {  // only apply linear interpolation on non-signal-like region
                for (int i = start + 1; i < end; i++) {
                    double content = start_content + (end_content - start_content) / (end - start) * (i - start);
                    signal.at(i) = content;
                }
                trivial_range_list.push_back(rng);
            }
        }
        else if (start < 0 && end < nsiglen) {  // sticky at the first tick
            for (int i = 0; i < end; i++) {
                signal.at(i) = signal.at(end);
            }
        }
        else if (start >= 0 && end == nsiglen) {  // sticky at the last tick
            for (int i = start + 1; i < end; i++) {
                signal.at(i) = signal.at(start);
            }
        }
    }

    for (auto& elem : trivial_range_list) {
        // use the remove/erase idiom
        rng_list.erase(std::remove(rng_list.begin(), rng_list.end(), elem), rng_list.end());
    }
    return true;
}

bool Protodune::FftInterpSticky(const IDFT::pointer& dft,
                                WireCell::Waveform::realseq_t& signal, WireCell::Waveform::BinRangeList& rng_list)
{
    ;
    const int nsiglen = signal.size();
    // group into two subsamples ("even" & "odd")
    int nsublen = nsiglen / 2;
    int nsublen2 = nsiglen - nsublen;
    WireCell::Waveform::realseq_t signal_even(nsublen);  // 0-th bin, 2, 4, etc.
    WireCell::Waveform::realseq_t signal_odd(nsublen2);
    for (int j = 0; j < nsiglen; j++) {
        if (j % 2 == 0)
            signal_even.at(j / 2) = signal.at(j);
        else
            signal_odd.at((j - 1) / 2) = signal.at(j);
    }

    // dft resampling for "even", see example in test_zero_padding.cxx
    auto tran_even = fwd_r2c(dft, signal_even);
    tran_even.resize(nsublen * 2);
    if (nsublen % 2 == 0) {
        std::rotate(tran_even.begin() + nsublen / 2, tran_even.begin() + nsublen, tran_even.end());
    }
    else {
        std::rotate(tran_even.begin() + (nsublen + 1) / 2, tran_even.begin() + nsublen, tran_even.end());
    }
    // inverse FFT
    auto signal_even_fc = inv_c2r(dft, tran_even);
    float scale = tran_even.size() / nsublen;
    WireCell::Waveform::scale(signal_even_fc, scale);

    // similar for "odd"
    auto tran_odd = fwd_r2c(dft, signal_odd);
    tran_odd.resize(nsublen2 * 2);
    if (nsublen2 % 2 == 0) {
        std::rotate(tran_odd.begin() + nsublen2 / 2, tran_odd.begin() + nsublen2, tran_odd.end());
    }
    else {
        std::rotate(tran_odd.begin() + (nsublen2 + 1) / 2, tran_odd.begin() + nsublen2, tran_odd.end());
    }
    auto signal_odd_fc = inv_c2r(dft, tran_odd);
    float scale2 = tran_odd.size() / nsublen2;
    WireCell::Waveform::scale(signal_odd_fc, scale2);

    // replace the linear interpolation with dft interpolation
    for (size_t j = 0; j < rng_list.size(); j++) {
        int start = rng_list.at(j).first - 1;
        int end = rng_list.at(j).second + 1;
        if (start >= 0 && end <= nsiglen) {
            for (int i = start + 1; i < end; i++) {
                if (i % 2 == 0) {  // predict "even" with "odd"
                    signal.at(i) = signal_odd_fc.at(i - 1);
                }
                else {
                    signal.at(i) = signal_even_fc.at(i);
                }
            }
        }
    }

    return true;
}

bool Protodune::FftShiftSticky(const IDFT::pointer& dft,
                               WireCell::Waveform::realseq_t& signal, double toffset,
                               std::vector<std::pair<int, int> >& st_ranges)
{
    const int nsiglen = signal.size();
    // group into two subsamples ("even" & "odd")
    int nsublen = nsiglen / 2;
    int nsublen2 = nsiglen - nsublen;
    WireCell::Waveform::realseq_t signal_even(nsublen);  // 0-th bin, 2, 4, etc.
    WireCell::Waveform::realseq_t signal_odd(nsublen2);
    for (int j = 0; j < nsiglen; j++) {
        if (j % 2 == 0)
            signal_even.at(j / 2) = signal.at(j);
        else
            signal_odd.at((j - 1) / 2) = signal.at(j);
    }

    // dft shift for "even"
    auto tran_even = fwd_r2c(dft, signal_even);
    double f0 = 1. / nsublen;
    const double PI = std::atan(1.0) * 4;
    for (size_t i = 0; i < tran_even.size(); i++) {
        double fi = i * f0;
        double omega = 2 * PI * fi;
        auto z = tran_even.at(i);
        std::complex<float> z1(0, omega * toffset);
        // std::complex<double> z2 = z * std::exp(z1);
        tran_even.at(i) = z * std::exp(z1);
    }
    // inverse FFT
    auto signal_even_fc = inv_c2r(dft, tran_even);

    // similar to "odd"
    auto tran_odd = fwd_r2c(dft, signal_odd);
    f0 = 1. / nsublen2;
    for (size_t i = 0; i < tran_odd.size(); i++) {
        double fi = i * f0;
        double omega = 2 * PI * fi;
        auto z = tran_odd.at(i);
        std::complex<float> z1(0, omega * toffset);
        // std::complex z2 = z * std::exp(z1);
        tran_odd.at(i) = z * std::exp(z1);
    }
    //
    auto signal_odd_fc = inv_c2r(dft, tran_odd);
    
    // float scale = 1./tran_odd.size();
    // WireCell::Waveform::scale(signal_odd_fc, 1./nsublen2);

    // replace the linear interpolation with dft interpolation
    for (size_t j = 0; j < st_ranges.size(); j++) {
        int start = st_ranges.at(j).first - 1;
        int end = st_ranges.at(j).second + 1;
        if (start >= 0 && end <= nsiglen) {
            for (int i = start + 1; i < end; i++) {
                if (i % 2 == 0) {  // predict "even" with "odd"
                    int ind = (i - 1) / 2. - toffset;
                    if (ind >= 0 && ind < nsublen2) signal.at(i) = signal_odd_fc.at(ind);
                    // std::cerr << "lc:fc= " << signal_lc.at(i) << " " << signel_even_fc.at() << std::endl;
                }

                else {
                    int ind = i / 2. - toffset;
                    if (ind >= 0 && ind < nsublen) signal.at(i) = signal_even_fc.at(ind);
                }
            }
        }
    }

    return true;
}

bool Protodune::FftScaling(const IDFT::pointer& dft,
                           WireCell::Waveform::realseq_t& signal, int nsamples)
{
    const int nsiglen = signal.size();
    auto tran = fwd_r2c(dft, signal);
    tran.resize(nsamples);
    if (nsiglen % 2 == 0) {  // ref test_zero_padding.cxx
        std::rotate(tran.begin() + nsiglen / 2, tran.begin() + nsiglen, tran.end());
    }
    else {
        std::rotate(tran.begin() + (nsiglen + 1) / 2, tran.begin() + nsiglen, tran.end());
    }
    // inverse FFT
    auto signal_fc = inv_c2r(dft, tran);

    WireCell::Waveform::scale(signal_fc, nsamples / nsiglen);
    signal = signal_fc;

    return true;
}

/*
 * Classes
 */

Protodune::StickyCodeMitig::StickyCodeMitig(const std::string& anode, const std::string& noisedb,
                                            float stky_sig_like_val, float stky_sig_like_rms, int stky_max_len)
  : m_anode_tn(anode)
  , m_noisedb_tn(noisedb)
  , m_stky_sig_like_val(stky_sig_like_val)
  , m_stky_sig_like_rms(stky_sig_like_rms)
  , m_stky_max_len(stky_max_len)

{
}
Protodune::StickyCodeMitig::~StickyCodeMitig() {}

void Protodune::StickyCodeMitig::configure(const WireCell::Configuration& cfg)
{
    m_anode_tn = get(cfg, "anode", m_anode_tn);
    m_anode = Factory::find_tn<IAnodePlane>(m_anode_tn);
    if (!m_anode) {
        THROW(KeyError() << errmsg{"failed to get IAnodePlane: " + m_anode_tn});
    }

    m_noisedb_tn = get(cfg, "noisedb", m_noisedb_tn);
    m_noisedb = Factory::find_tn<IChannelNoiseDatabase>(m_noisedb_tn);

    std::string dft_tn = get<std::string>(cfg, "dft", "FftwDFT");
    m_dft = Factory::find_tn<IDFT>(dft_tn);

    m_extra_stky.clear();
    auto jext = cfg["extra_stky"];
    if (!jext.isNull()) {
        for (auto jone : jext) {
            auto jchans = jone["channels"];
            for (auto jchan : jchans) {
                int channel = jchan.asInt();
                // std::cerr << "[wgu] ch# " << channel << " has " << jone["bits"].size() << " extra stky bits:" <<
                // std::endl;
                for (auto bit : jone["bits"]) {
                    // std::cerr << "[wgu] " << bit.asInt() << std::endl;
                    m_extra_stky[channel].push_back((short int) bit.asInt());
                }
            }
        }
    }

    if (cfg.isMember("stky_sig_like_val")) {
        m_stky_sig_like_val = get(cfg, "stky_sig_like_val", m_stky_sig_like_val);
    }

    if (cfg.isMember("stky_sig_like_rms")) {
        m_stky_sig_like_rms = get(cfg, "stky_sig_like_rms", m_stky_sig_like_rms);
    }

    if (cfg.isMember("stky_max_len")) {
        m_stky_max_len = get(cfg, "stky_max_len", m_stky_max_len);
    }
}
WireCell::Configuration Protodune::StickyCodeMitig::default_configuration() const
{
    Configuration cfg;
    cfg["anode"] = m_anode_tn;
    cfg["noisedb"] = m_noisedb_tn;
    cfg["stky_sig_like_val"] = m_stky_sig_like_val;
    cfg["stky_sig_like_rms"] = m_stky_sig_like_rms;
    cfg["stky_max_len"] = m_stky_max_len;
    cfg["dft"] = "FftwDFT";     // type-name for the DFT to use
    return cfg;
}

WireCell::Waveform::ChannelMaskMap Protodune::StickyCodeMitig::apply(int ch, signal_t& signal) const
{
    WireCell::Waveform::ChannelMaskMap ret;

    const int nsiglen = signal.size();
    // tag sticky bins
    // WireCell::Waveform::BinRange sticky_rng;
    WireCell::Waveform::BinRangeList sticky_rng_list;
    std::vector<short int> extra;
    if (m_extra_stky.find(ch) != m_extra_stky.end()) {
        extra = m_extra_stky.at(ch);
    }

    for (int i = 0; i < nsiglen; i++) {
        int val = signal.at(i);
        int mod = val % 64;
        auto it = std::find(extra.begin(), extra.end(), mod);
        // bool is_stky = (mod==0 || mod==1 || mod==63 || it!=extra.end());
        if (it != extra.end()) {
            if (sticky_rng_list.empty()) {
                sticky_rng_list.push_back(std::make_pair(i, i));
            }
            else if ((sticky_rng_list.back().second + 1) == i) {
                sticky_rng_list.back().second = i;
            }
            else {
                sticky_rng_list.push_back(std::make_pair(i, i));
            }
        }
    }

    // int ent_stkylen = 0;
    for (auto rng : sticky_rng_list) {
        int stkylen = rng.second - rng.first;
        if (stkylen > m_stky_max_len) {
            ret["sticky"][ch].push_back(rng);
            // ent_stkylen += stkylen;
        }
    }
    // std::cerr << "[wgu] ch: " << ch << " long_stkylen: " << long_stkylen << std::endl;

    LinearInterpSticky(signal, sticky_rng_list, m_stky_sig_like_val, m_stky_sig_like_rms);
    FftInterpSticky(m_dft, signal, sticky_rng_list);

    // Now calculate the baseline ...
    std::pair<double, double> temp = WireCell::Waveform::mean_rms(signal);
    auto temp_signal = signal;
    for (size_t i = 0; i != temp_signal.size(); i++) {
        if (fabs(temp_signal.at(i) - temp.first) > 6 * temp.second) {
            temp_signal.at(i) = temp.first;
        }
    }
    float baseline = WireCell::Waveform::median_binned(temp_signal);

    // Ledge
    auto wpid = m_anode->resolve(ch);
    const int iplane = wpid.index();
    if (iplane == 2) {
        int ledge_start[3], ledge_end[3];
        int nledge = LedgeIdentify1(signal, baseline, ledge_start, ledge_end);
        for (int LE = 0; LE < nledge; LE++) {
            // check overlap of sticky in a ledge
            float overlap = 0;
            for (auto rng : sticky_rng_list) {
                int redge = std::min(rng.second, ledge_end[LE]);
                int ledge = std::max(rng.first, ledge_start[LE]);
                if (redge >= ledge) overlap += (redge - ledge + 1);
            }
            if (overlap < 0.5 * (ledge_end[LE] - ledge_start[LE])) {
                WireCell::Waveform::BinRange ledge_bins;
                ledge_bins.first = ledge_start[LE];
                ledge_bins.second = ledge_end[LE];
                ret["ledge"][ch].push_back(ledge_bins);  // FIXME: maybe collection plane only?
                // std::cerr << "[wgu] ledge found, ch "<< ch << " , bins= [" << ledge_bins.first << " , " <<
                // ledge_bins.second << " ], sticky overlap: " << overlap/(ledge_end[LE]-ledge_start[LE]) << std::endl;
            }
        }
    }

    return ret;
}

WireCell::Waveform::ChannelMaskMap Protodune::StickyCodeMitig::apply(channel_signals_t& chansig) const
{
    return WireCell::Waveform::ChannelMaskMap();
}


Protodune::OneChannelNoise::OneChannelNoise(const std::string& anode, const std::string& noisedb)
  : m_anode_tn(anode)
  , m_noisedb_tn(noisedb)
  , m_check_partial()  // fixme, here too.
  , m_resmp()
{
}
Protodune::OneChannelNoise::~OneChannelNoise() {}

void Protodune::OneChannelNoise::configure(const WireCell::Configuration& cfg)
{
    m_anode_tn = get(cfg, "anode", m_anode_tn);
    m_anode = Factory::find_tn<IAnodePlane>(m_anode_tn);
    m_noisedb_tn = get(cfg, "noisedb", m_noisedb_tn);
    m_noisedb = Factory::find_tn<IChannelNoiseDatabase>(m_noisedb_tn);

    std::string dft_tn = get<std::string>(cfg, "dft", "FftwDFT");
    m_dft = Factory::find_tn<IDFT>(dft_tn);

    m_resmp.clear();
    auto jext = cfg["resmp"];
    if (!jext.isNull()) {
        for (auto jone : jext) {
            int smpin = jone["sample_from"].asInt();
            for (auto jch : jone["channels"]) {
                int channel = jch.asInt();
                m_resmp[channel] = smpin;
            }
        }
    }
}
WireCell::Configuration Protodune::OneChannelNoise::default_configuration() const
{
    Configuration cfg;
    cfg["anode"] = m_anode_tn;
    cfg["noisedb"] = m_noisedb_tn;
    cfg["dft"] = "FftwDFT";     // type-name for the DFT to use
    return cfg;
}

WireCell::Waveform::ChannelMaskMap Protodune::OneChannelNoise::apply(int ch, signal_t& signal) const
{
    WireCell::Waveform::ChannelMaskMap ret;

    // do we need a nominal baseline correction?
    // float baseline = m_noisedb->nominal_baseline(ch);

    // correct the FEMB 302 clock issue
    auto it = m_resmp.find(ch);
    if (it != m_resmp.end()) {
        int smpin = m_resmp.at(ch);
        int smpout = signal.size();
        signal.resize(smpin);
        FftScaling(m_dft, signal, smpout);
        // std::cerr << "[wgu] ch: " << ch << " smpin: " << smpin << " smpout: " << smpout << std::endl;
    }

    // correct rc undershoot
    auto spectrum = fwd_r2c(m_dft, signal);
    bool is_partial = m_check_partial(spectrum);  // Xin's "IS_RC()"

    if (!is_partial) {
        auto const& spec = m_noisedb->rcrc(ch);  // rc_layers set to 1 in channel noise db
        WireCell::Waveform::shrink(spectrum, spec);
    }

    // remove the "50kHz" noise in some collection channels
    // FIXME: do we need a channel list input?
    auto wpid = m_anode->resolve(ch);
    const int iplane = wpid.index();
    if (iplane == 2) {
        auto mag = WireCell::Waveform::magnitude(spectrum);
        mag.at(0) = 0;
        Microboone::RawAdapativeBaselineAlg(mag);  // subtract "linear" background in spectrum

        auto const& spec = m_noisedb->noise(ch);

        // spec -> freqBins;
        std::vector<std::pair<int, int> > freqBins;
        for (int i = 0; i < (int) spec.size(); i++) {
            float flag_re = spec.at(i).real();
            // if(i<3000) std::cout << "[wgu] spec # " << i << " : " << flag_re << std::endl;
            if (flag_re < 0.5) {  // found a tagged bin
                if (freqBins.empty()) {
                    freqBins.push_back(std::make_pair(i, i));
                }
                else if (i == freqBins.back().second + 1) {
                    freqBins.back().second = i;
                }
                else {
                    freqBins.push_back(std::make_pair(i, i));
                }
            }
        }

        // std::cout << "[wgu] freqBins.size(): " << freqBins.size() << std::endl;
        // for(auto br: freqBins) {
        //     std::cout << "[wgu] hibin: " << br.second << " lobin: " << br.first << std::endl;
        // }

        int n_harmonic = 0;
        for (int iter = 0; iter < 5; iter++) {
            // std::cout << "[wgu] iter= " << iter << std::endl;

            for (auto bin : freqBins) {
                int istart = bin.first;
                int iend = bin.second + 1;
                int nslice = iend - istart;
                // std::cout << "hibin: " << iend << " lobin: " << istart << std::endl;

                // std::cerr << istart << " " << iend << std::endl;
                WireCell::Waveform::realseq_t mag_slice(nslice);  // slice of magnitude spectrum
                std::copy(mag.begin() + istart, mag.begin() + iend, mag_slice.begin());
                std::pair<double, double> stat = WireCell::Waveform::mean_rms(mag_slice);
                // std::cerr << "[wgu] mean/rms: " << stat.first << " " << stat.second << std::endl;
                double cut = stat.first + 2.7 * stat.second;
                if (istart > 1050) {  // if(i>17){
                    cut = stat.first + 3 * stat.second;
                }

                for (int j = istart; j < iend; j++) {
                    float content = mag.at(j);
                    if (content > cut) {
                        spectrum.at(j).real(0);
                        spectrum.at(j).imag(0);
                        spectrum.at(mag.size() + 1 - j).real(0);
                        spectrum.at(mag.size() + 1 - j).imag(0);
                        mag.at(j) = 0;
                        mag.at(mag.size() + 1 - j) = 0;
                        n_harmonic++;
                    }
                }
            }
        }

        if (n_harmonic > 5) {
            WireCell::Waveform::BinRange temp_harmonic_bins;
            temp_harmonic_bins.first = 0;
            temp_harmonic_bins.second = signal.size();
            ret["harmonic"][ch].push_back(temp_harmonic_bins);
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

    // Now do the adaptive baseline for the bad RC channels
    if (is_partial) {
        // add something
        WireCell::Waveform::BinRange temp_chirped_bins;
        temp_chirped_bins.first = 0;
        temp_chirped_bins.second = signal.size();

        if (iplane != 2) {  // not collection
            ret["lf_noisy"][ch].push_back(temp_chirped_bins);
            // std::cout << "Partial " << ch << std::endl;
        }
        Microboone::SignalFilter(signal);
        Microboone::RawAdapativeBaselineAlg(signal);
        Microboone::RemoveFilterFlags(signal);
    }

    const float min_rms = m_noisedb->min_rms_cut(ch);
    const float max_rms = m_noisedb->max_rms_cut(ch);
    //
    // temp = Derivations::CalcRMS(signal);
    // if(ch==14431) std::cerr << "[wgu] channel: " << ch << " rms: " << temp.second << std::endl;
    // if(temp.second<min_rms || temp.second>max_rms){
    //     WireCell::Waveform::BinRange temp_chirped_bins;
    //     temp_chirped_bins.first = 0;
    //     temp_chirped_bins.second = signal.size();
    //     ret["noisy"][ch].push_back(temp_chirped_bins);
    // }
    // alternative RMS tagging
    Microboone::SignalFilter(signal);
    bool is_noisy = Microboone::NoisyFilterAlg(signal, min_rms, max_rms);
    Microboone::RemoveFilterFlags(signal);
    if (is_noisy) {
        WireCell::Waveform::BinRange temp_chirped_bins;
        temp_chirped_bins.first = 0;
        temp_chirped_bins.second = signal.size();
        ret["noisy"][ch].push_back(temp_chirped_bins);
    }

    return ret;
}

WireCell::Waveform::ChannelMaskMap Protodune::OneChannelNoise::apply(channel_signals_t& chansig) const
{
    return WireCell::Waveform::ChannelMaskMap();
}

Protodune::RelGainCalib::RelGainCalib(const std::string& anode, const std::string& noisedb, float gain_def,
                                      float gain_min_cut, float gain_max_cut)
  : m_anode_tn(anode)
  , m_noisedb_tn(noisedb)
  , m_gain_def(gain_def)
  , m_gain_min_cut(gain_min_cut)
  , m_gain_max_cut(gain_max_cut)
{
}
Protodune::RelGainCalib::~RelGainCalib() {}

void Protodune::RelGainCalib::configure(const WireCell::Configuration& cfg)
{
    m_anode_tn = get(cfg, "anode", m_anode_tn);
    m_anode = Factory::find_tn<IAnodePlane>(m_anode_tn);
    if (!m_anode) {
        THROW(KeyError() << errmsg{"failed to get IAnodePlane: " + m_anode_tn});
    }

    m_noisedb_tn = get(cfg, "noisedb", m_noisedb_tn);
    m_noisedb = Factory::find_tn<IChannelNoiseDatabase>(m_noisedb_tn);

    if (cfg.isMember("gain_def")) {
        m_gain_def = get(cfg, "gain_def", m_gain_def);
    }

    if (cfg.isMember("gain_min_cut")) {
        m_gain_min_cut = get(cfg, "gain_min_cut", m_gain_min_cut);
    }

    if (cfg.isMember("gain_max_cut")) {
        m_gain_max_cut = get(cfg, "gain_max_cut", m_gain_max_cut);
    }

    m_rel_gain.clear();
    for (auto jch : cfg["rel_gain"]) {
        m_rel_gain.push_back(jch.asDouble());
    }
}
WireCell::Configuration Protodune::RelGainCalib::default_configuration() const
{
    Configuration cfg;
    cfg["anode"] = m_anode_tn;
    cfg["noisedb"] = m_noisedb_tn;
    cfg["gain_def"] = m_gain_def;
    cfg["gain_min_cut"] = m_gain_min_cut;
    cfg["gain_max_cut"] = m_gain_max_cut;
    return cfg;
}

WireCell::Waveform::ChannelMaskMap Protodune::RelGainCalib::apply(int ch, signal_t& signal) const
{
    WireCell::Waveform::ChannelMaskMap ret;

    size_t nsiglen = signal.size();
    // FIXME: already subtracted in OneChannelNoise::apply(). Need to calculate again?
    // std::pair<double,double> temp = WireCell::Waveform::mean_rms(signal);
    // auto temp_signal = signal;
    // for (size_t i=0;i!=temp_signal.size();i++){
    // if (fabs(temp_signal.at(i)-temp.first)>6*temp.second){
    //     temp_signal.at(i) = temp.first;
    // }
    // }
    // float baseline = WireCell::Waveform::median_binned(temp_signal);
    float baseline = 0.;

    float gain_fac = m_rel_gain.at(ch);

    bool isDefVal = false;
    if (gain_fac > m_gain_max_cut || gain_fac < m_gain_min_cut) {
        // std::cerr << "[wgu] gain_max_cut: " << m_gain_max_cut << " gain_min_cut: " << m_gain_min_cut << std::endl;
        // std::cerr << "[wgu] ch: " << ch << " relative gain: " << gain_fac << std::endl;
        gain_fac = m_gain_def;
        isDefVal = true;
    }

    if (!isDefVal) {
        for (size_t ind = 0; ind < nsiglen; ind++) {
            float sigout = gain_fac * (signal.at(ind) - baseline);
            sigout += baseline;
            signal.at(ind) = sigout;
        }
    }

    return ret;
}

WireCell::Waveform::ChannelMaskMap Protodune::RelGainCalib::apply(channel_signals_t& chansig) const
{
    return WireCell::Waveform::ChannelMaskMap();
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
