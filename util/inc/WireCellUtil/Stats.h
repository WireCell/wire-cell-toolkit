#ifndef WIRECELLUTIL_STATS
#define WIRECELLUTIL_STATS

#include "WireCellUtil/Binning.h"

#include <ostream>
#include <vector>

namespace WireCell {

    // Keep track of count, sum and sum of squares with optional
    // simple histogramming.
    struct Stats {
        bool sample{false};
        int S0{0};
        double S1{0}, S2{0};
        std::vector<double> samples; 

        double mean() const {
            if (S0 == 0) return 0.0;
            return S1/S0;
        }

        double rms() const {
            if (S0 <= 1) return -1.0;
            const double d = S2 - S1*S1/S0;
            // round off errors can give small negative
            if (d < 0) return 0.0;
            return sqrt(d/(S0-1));
        }

        void operator()(double val) {
            S0 += 1;
            S1 += val;
            S2 += val*val;
            if (sample) samples.push_back(val);
        }

        // Return histogram of samples.  If binning min/max are both
        // zero, set and use sample min/max.  If over_under_flow is
        // true, the histogram will have nbins+2 with the two extreme
        // bins holding underflow and overflow counts.  In either case
        // the passed "bins" object is modified.
        std::vector<double> histogram(Binning& bins, 
                                      bool over_under_flow = false) const
        {
            std::vector<double> hist;
            if (not sample) return hist;

            if (bins.min() == 0 and bins.max() == 0) {
                auto mm = std::minmax_element(samples.begin(), samples.end());
                bins = Binning(bins.nbins(), *mm.first, *mm.second);
            }

            if (over_under_flow) {
                bins = Binning(bins.nbins()+2,
                               bins.min()-bins.binsize(),
                               bins.max()+bins.binsize());
            }

            const int nbins = bins.nbins();

            hist.resize(nbins, 0.0);
            for (const double s : samples) {
                int ibin = bins.bin(s);
                if (ibin < 0) {
                    if (over_under_flow) {
                        ibin = 0;
                    }
                    else continue;
                }
                if (ibin >= nbins) {
                    if (over_under_flow) {
                        ibin = nbins-1;
                    }
                    else continue;
                }
                ++hist[ibin];
            }
            return hist;
        }

    };

    std::ostream& operator<<(std::ostream& o, const Stats& s) {
        o << s.mean() << " +/- " << s.rms()
          // << std::setprecision(std::numeric_limits<double>::max_digits10 - 1) 
          << " [" << s.S0 << " " << s.S1 << " " << s.S2 << "]";
        return o;
    }
}

#endif
