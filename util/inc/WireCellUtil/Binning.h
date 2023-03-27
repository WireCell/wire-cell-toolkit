#ifndef WIRECELLUTIL_BINNING_H
#define WIRECELLUTIL_BINNING_H

#include <map>  // for std::pair
#include <cmath>
#include <iostream>  // for ostream

namespace WireCell {

    /** A binning is a uniform discretization of a linear space.

        This class largely provides methods that give semantic labels
        to various calculations related to a binned region.
    */
    class Binning {
        int m_nbins;
        double m_minval, m_maxval, m_binsize;

       public:
        /** Create a binning
            \param nbins gives the number of uniform partitions of space between bounds.
            \param minval gives the lower bound of the linear space (low edge of bin 0)
            \param maxval gives the upper bound of the linear space (high edge of bin nbins-1)
        */
        Binning(int nbins, double minval, double maxval)
          : m_nbins(0)
          , m_minval(0)
          , m_maxval(0)
          , m_binsize(0)
        {
            set(nbins, minval, maxval);
        }
        Binning()
          : m_nbins(0)
          , m_minval(0)
          , m_maxval(0)
          , m_binsize(0)
        {
        }

        // Post constructor setting
        void set(int nbins, double minval, double maxval)
        {
            m_nbins = nbins;
            m_minval = minval;
            m_maxval = maxval;
            m_binsize = ((maxval - minval) / nbins);
        }

        // Access given number of bins.
        int nbins() const { return m_nbins; }

        // Access given minimum range of binning.
        double min() const { return m_minval; }

        // Access given maximum range of binning.
        double max() const { return m_maxval; }

        /// Return the max-min
        double span() const { return m_maxval - m_minval; }

        // Binning as a range.
        std::pair<double, double> range() const { return std::make_pair(m_minval, m_maxval); }

        // Return half open range of bin indices or alternatively
        // fully closed range of edge indices.
        std::pair<int, int> irange() const { return std::make_pair(0, m_nbins); }

        // Access resulting bin size..
        double binsize() const { return m_binsize; }

        /// Return the bin containing value.  If val is in range,
        /// return value is [0,nbins-1] but no range checking is
        /// performed.
        int bin(double val) const { return int((val - m_minval) / m_binsize); }

        /// Return the center value of given bin.  Range checking is
        /// not done.
        double center(int ind) const { return m_minval + (ind + 0.5) * m_binsize; }

        /// Return the edge, nominally in [0,nbins] closest to the
        /// given value.  Range checking is not done so returned edge
        /// may be outside of range.
        int edge_index(double val) const { return int(round((val - m_minval) / m_binsize)); }

        /// Return the position of the given bin edge.  Range checking
        /// is not done.
        double edge(int ind) const { return m_minval + ind * m_binsize; }

        /// Return true value is in range.  Range is considered half
        /// open.  Ig, edge(nbins) is not inside range.
        bool inside(double val) const { return m_minval <= val && val < m_maxval; }

        /// Return true if bin is in bounds.
        bool inbounds(int bin) const { return 0 <= bin && bin < m_nbins; }

        /// Return half-open bin range which covers the range of
        /// values.  Bounds are forced to return values in [0,nbins].
        std::pair<int, int> sample_bin_range(double minval, double maxval) const
        {
            return std::make_pair(std::max(bin(minval), 0), std::min(bin(maxval) + 1, m_nbins));
        }
    };

    inline std::ostream& operator<<(std::ostream& os, const WireCell::Binning& bins)
    {
        os << bins.nbins() << "@[" << bins.min() << "," << bins.max() << "]";
        return os;
    }

    // Return a subset of the binning which contains the bounds
    inline
    Binning subset(const Binning& bins, double xmin, double xmax) {
        if (xmin > xmax) std::swap(xmin, xmax);
        const int lo = std::max(0,            bins.bin(xmin));
        const int hi = std::min(bins.nbins(), bins.bin(xmax));
        const int n = hi-lo;
        if (n <= 0) {
            return Binning();
        }
        return Binning(hi-lo, bins.edge(lo), bins.edge(hi));
    }

    // P(X<=L) for X ~ N(mean,sigma)
    inline
    double gcumulative(double L, double mean=0, double sigma=1) {
        const double scale = sqrt(2)*sigma;
        return 0.5+0.5*std::erf((L - mean)/scale);
    }
    // P(L1 <= X <= L2) for X ~ N(mean,sigma)
    inline
    double gbounds(double L1, double L2, double mean=0, double sigma=1) {
        if (L1 == L2) return 0;
        if (L1 > L2) std::swap(L1, L2);
        const double scale = sqrt(2)*sigma;
        return 0.5*( std::erf((L2 - mean)/scale) - std::erf((L1 - mean)/scale));
    }
        
    // Add the absolutely normlized, bin-integrated Gaussian
    // distribution to iterated elements.  Return sum of bin integrals
    // which will not < 1.0 given finite binning span.
    template<typename OutputIt>
    double gaussian(OutputIt out, const Binning& bins, double mean=0, double sigma=1)
    {
        const int nedges = bins.nbins() + 1;
        double last_cumu = 0;
        double total = 0;
        for (int edge = 0 ; edge < nedges; ++edge) {
            const double this_cumu = gcumulative(bins.edge(edge), mean, sigma);
            if (edge) {
                const double bin_cumu = (this_cumu - last_cumu);
                total += bin_cumu;
                *out = bin_cumu;
                ++out;
            }
            last_cumu = this_cumu;
        }
        return total;
    }

}  // namespace WireCell

#endif /* WIRECELLUTIL_BINNING_H */
