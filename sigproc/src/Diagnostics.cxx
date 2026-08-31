#include "WireCellSigProc/Diagnostics.h"

#include <cmath>
#include <complex>
//#include <iostream>

using namespace WireCell::SigProc;

Diagnostics::Partial::Partial(int nfreqs, double maxpower)
  : nfreqs(nfreqs)
  , maxpower(maxpower)
{
}

bool Diagnostics::Partial::operator()(const WireCell::Waveform::compseq_t& spec) const
{
    const double mag0 = std::abs(spec[0 + 1]);
    double sum = mag0;
    for (int ind = 1; ind <= nfreqs && ind + 1 < (int) spec.size(); ++ind) {
        const double magi = std::abs(spec[ind + 1]);
        if (mag0 <= magi) {
            return false;
        }
        sum += magi;
    }
    // if (sum/(nfreqs+1) > maxpower){
    //   std::cout << mag0 << " " << std::abs(spec[2]) << " " << std::abs(spec[3]) << " " << std::abs(spec[4]) << " " <<
    //   std::abs(spec[5]) << std::endl;
    // }

    return sum / (nfreqs + 1) > maxpower;
}

Diagnostics::Chirp::Chirp(int windowSize, double chirpMinRMS, double maxNormalNeighborFrac)
  : windowSize(windowSize)
  , chirpMinRMS(chirpMinRMS)
  , maxNormalNeighborFrac(maxNormalNeighborFrac)
{
}

bool Diagnostics::Chirp::operator()(const WireCell::Waveform::realseq_t& sig, int& beg, int& end) const
{
    beg = end = -1;

    // this alg is from Mike Mooney

    // magic numbers moved into constructor
    // const int windowSize = 20;
    // const double chirpMinRMS = 0.9;
    // const double maxNormalNeighborFrac = 0.20;

    int counter = 0;
    double runningAmpMean = 0.0;
    double runningAmpRMS = 0.0;
    int numLowRMS = 0;
    int firstLowRMSBin = -1;  // computed below assuming
    int lastLowRMSBin = -1;   //  1-based arrays
    bool lowRMSFlag = false;
    double RMSfirst = 0.0;
    double RMSsecond = 0.0;
    double RMSthird = 0.0;
    int numNormalNeighbors = 0;

    const int numBins = sig.size();

    for (int ibin = 0; ibin < numBins; ibin++) {
        double ADCval = sig[ibin];

        runningAmpMean += ADCval;
        runningAmpRMS += ADCval * ADCval;

        counter++;
        if (counter == windowSize) {
            runningAmpMean /= (double) windowSize;
            runningAmpRMS /= (double) windowSize;
            runningAmpRMS = std::sqrt(std::max(0.0, runningAmpRMS - runningAmpMean * runningAmpMean));

            RMSfirst = RMSsecond;
            RMSsecond = RMSthird;
            RMSthird = runningAmpRMS;

            if (runningAmpRMS < chirpMinRMS) {
                numLowRMS++;
            }

            if (ibin >= 3 * windowSize - 1) {
                if ((RMSsecond < chirpMinRMS) && ((RMSfirst > chirpMinRMS) || (RMSthird > chirpMinRMS))) {
                    numNormalNeighbors++;
                }

                if (lowRMSFlag == false) {
                    if ((RMSsecond < chirpMinRMS) && (RMSthird < chirpMinRMS)) {
                        lowRMSFlag = true;
                        firstLowRMSBin = ibin - 2 * windowSize + 1;
                        lastLowRMSBin = ibin - windowSize + 1;
                    }

                    if ((ibin == 3 * windowSize - 1) && (RMSfirst < chirpMinRMS) && (RMSsecond < chirpMinRMS)) {
                        lowRMSFlag = true;
                        firstLowRMSBin = ibin - 3 * windowSize + 1;
                        lastLowRMSBin = ibin - 2 * windowSize + 1;
                    }
                }
                else {
                    if ((RMSsecond < chirpMinRMS) && (RMSthird < chirpMinRMS)) {
                        lastLowRMSBin = ibin - windowSize + 1;
                    }
                }
            }

            counter = 0;
            runningAmpMean = 0.0;
            runningAmpRMS = 0.0;
        }
    }

    double chirpFrac = ((double) numLowRMS) / (((double) numBins) / ((double) windowSize));
    double normalNeighborFrac = ((double) numNormalNeighbors) / ((double) numLowRMS);
    if ((numLowRMS > 4) &&
        ((normalNeighborFrac < maxNormalNeighborFrac) ||
         ((numLowRMS < 2.0 / maxNormalNeighborFrac) && (lastLowRMSBin - firstLowRMSBin == numLowRMS * windowSize)))) {
        firstLowRMSBin = std::max(1, firstLowRMSBin - windowSize);
        lastLowRMSBin = std::min(numBins, lastLowRMSBin + 2 * windowSize);

        if ((numBins - lastLowRMSBin) < windowSize) {
            lastLowRMSBin = numBins;
        }

        if (chirpFrac > 0.99) {
            firstLowRMSBin = 1;
            lastLowRMSBin = numBins;
        }

        beg = firstLowRMSBin - 1;
        end = lastLowRMSBin;
        return true;
    }
    return false;
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
