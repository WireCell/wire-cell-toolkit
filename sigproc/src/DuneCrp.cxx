/***************************************************************************/
/*             Noise filter (NF) module for DUNE CRP                       */
/***************************************************************************/
// Features:                       
//   undershoot correction         
//   partial undershoot correction 
//
// 07/19/2023, created by W.Gu (wgu@bnl.gov) 


#include "WireCellSigProc/DuneCrp.h"
#include "WireCellSigProc/Derivations.h"

#include "WireCellAux/DftTools.h"

#include "WireCellUtil/NamedFactory.h"

#include <cmath>
#include <complex>
#include <iostream>
#include <set>

WIRECELL_FACTORY(duneCrpOneChannelNoise, WireCell::SigProc::DuneCrp::OneChannelNoise, WireCell::IChannelFilter,
                 WireCell::IConfigurable)

using namespace WireCell::SigProc;
using WireCell::Aux::DftTools::fwd_r2c;
using WireCell::Aux::DftTools::inv_c2r;

/*
 * Functions (adapted from Microboone)
 */

bool DuneCrp::SignalFilter(WireCell::Waveform::realseq_t& sig)
{
    const double sigFactor = 4.0;
    const int padBins = 8;

    float rmsVal = DuneCrp::CalcRMSWithFlags(sig);
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

bool DuneCrp::RawAdapativeBaselineAlg(WireCell::Waveform::realseq_t& sig)
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

bool DuneCrp::RemoveFilterFlags(WireCell::Waveform::realseq_t& sig)
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


float DuneCrp::CalcRMSWithFlags(const WireCell::Waveform::realseq_t& sig)
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

bool DuneCrp::NoisyFilterAlg(WireCell::Waveform::realseq_t& sig, float min_rms, float max_rms)
{
    const double rmsVal = DuneCrp::CalcRMSWithFlags(sig);

    if (rmsVal > max_rms || rmsVal < min_rms) {
        int numBins = sig.size();
        for (int i = 0; i < numBins; i++) {
            sig.at(i) = 100000.0/*10000.0*/;
        }

        return true;
    }

    return false;
}


/*
 * Classes
 */

DuneCrp::OneChannelNoise::OneChannelNoise(const std::string& anode, const std::string& noisedb)
  : m_anode_tn(anode)
  , m_noisedb_tn(noisedb)
  , m_check_partial()  // fixme, here too.
{
}
DuneCrp::OneChannelNoise::~OneChannelNoise() {}

void DuneCrp::OneChannelNoise::configure(const WireCell::Configuration& cfg)
{
    m_anode_tn = get(cfg, "anode", m_anode_tn);
    m_anode = Factory::find_tn<IAnodePlane>(m_anode_tn);
    m_noisedb_tn = get(cfg, "noisedb", m_noisedb_tn);
    m_noisedb = Factory::find_tn<IChannelNoiseDatabase>(m_noisedb_tn);

    std::string dft_tn = get<std::string>(cfg, "dft", "FftwDFT");
    m_dft = Factory::find_tn<IDFT>(dft_tn);
}
WireCell::Configuration DuneCrp::OneChannelNoise::default_configuration() const
{
    Configuration cfg;
    cfg["anode"] = m_anode_tn;
    cfg["noisedb"] = m_noisedb_tn;
    cfg["dft"] = "FftwDFT";     // type-name for the DFT to use
    return cfg;
}

WireCell::Waveform::ChannelMaskMap DuneCrp::OneChannelNoise::apply(int ch, signal_t& signal) const
{
    WireCell::Waveform::ChannelMaskMap ret;

    // do we need a nominal baseline correction?
    // float baseline = m_noisedb->nominal_baseline(ch);

    // correct rc undershoot
    auto spectrum = fwd_r2c(m_dft, signal);
    bool is_partial = m_check_partial(spectrum);  // Xin's "IS_RC()"

    if (!is_partial) {
        auto const& spec = m_noisedb->rcrc(ch);  // rc_layers set to 1 in channel noise db
        WireCell::Waveform::shrink(spectrum, spec);
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
       // std::cout << "[duneCrp] is_partical channel: " << ch << std::endl;

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
        DuneCrp::SignalFilter(signal);
        DuneCrp::RawAdapativeBaselineAlg(signal);
        DuneCrp::RemoveFilterFlags(signal);
    }

    // const float min_rms = m_noisedb->min_rms_cut(ch);
    // const float max_rms = m_noisedb->max_rms_cut(ch);
    // // alternative RMS tagging
    // DuneCrp::SignalFilter(signal);
    // bool is_noisy = DuneCrp::NoisyFilterAlg(signal, min_rms, max_rms);
    // DuneCrp::RemoveFilterFlags(signal);
    // if (is_noisy) {
    //     WireCell::Waveform::BinRange temp_bin_range;
    //     temp_bin_range.first = 0;
    //     temp_bin_range.second = signal.size();
    //     ret["noisy"][ch].push_back(temp_bin_range);
    // }

    return ret;
}

WireCell::Waveform::ChannelMaskMap DuneCrp::OneChannelNoise::apply(channel_signals_t& chansig) const
{
    return WireCell::Waveform::ChannelMaskMap();
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
