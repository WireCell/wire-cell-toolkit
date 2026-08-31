#ifndef WIRECELLGEN_BINNEDDIFFUSION
#define WIRECELLGEN_BINNEDDIFFUSION

#include "WireCellUtil/Pimpos.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/Units.h"

#include "WireCellIface/IDepo.h"
#include "WireCellIface/IDFT.h"

#include "WireCellGen/ImpactData.h"

#include <deque>

namespace WireCell {
    namespace Gen {

        /**  A BinnedDiffusion maintains an association between impact
         * positions along the pitch direction of a wire plane and
         * the diffused depositions that drift to them.
         *
         * It covers a fixed and discretely sampled time and pitch
         * domain.
         */
        class BinnedDiffusion {
           public:
            /** Create a BinnedDiffusion.

               Arguments are:

               - pimpos :: a Pimpos instance defining the wire and impact binning.

               - tbins :: a Binning instance defining the time sampling binning.

               - nsigma :: number of sigma the 2D (transverse X
                 longitudinal) Gaussian extends.

               - fluctuate :: set to an IRandom if charge-preserving
                 Poisson fluctuations are to be applied.

               - calcstrat :: set a calculation strategy that gives
                 how the microscopic distribution of charge between
                 two impacts will be interpolated toward either edge.

             */

            /// Useful to client code to mark a calculation strategy.
            enum ImpactDataCalculationStrategy { constant = 1, linear = 2 };

            BinnedDiffusion(const Pimpos& pimpos, const IDFT::pointer& dft,
                            const Binning& tbins, double nsigma = 3.0,
                            IRandom::pointer fluctuate = nullptr, ImpactDataCalculationStrategy calcstrat = linear);

            const Pimpos& pimpos() const { return m_pimpos; }
            const Binning& tbins() const { return m_tbins; }

            /// Add a deposition and its associated diffusion sigmas.
            /// Return false if no activity falls within the domain.
            bool add(IDepo::pointer deposition, double sigma_time, double sigma_pitch);

            /// Unconditionally associate an already built
            /// GaussianDiffusion to one impact.
            void add(std::shared_ptr<GaussianDiffusion> gd, int impact_index);

            /// Drop any stored ImpactData within the half open
            /// impact index range.
            void erase(int begin_impact_index, int end_impact_index);

            /// Return the data in the given impact bin.  Note, this
            /// bin represents drifted charge between two impact
            /// positions.  Take care when using BinnedDiffusion and
            /// field responses because epsilon above or below the
            /// impact position exactly in the middle of two wires
            /// drastically different response.
            ImpactData::pointer impact_data(int bin) const;

            // test ...
            //	    void get_charge_vec(std::vector<std::vector<std::tuple<int,int, double> > >& vec_vec_charge,
            //std::vector<int>& vec_impact);

            /// Return the range of pitch containing depos out to
            /// given nsigma and without bounds checking.
            std::pair<double, double> pitch_range(double nsigma = 0.0) const;

            /// Return the half open bin range of impact bins,
            /// constrained so that either number is in [0,nimpacts].
            std::pair<int, int> impact_bin_range(double nsigma = 0.0) const;

            /// Return the range of time containing depos out to given
            /// nsigma and without bounds checking.
            std::pair<double, double> time_range(double nsigma = 0.0) const;

            /// Return the half open bin range for time bins
            /// constrained so that either number is in [0,nticks].
            std::pair<int, int> time_bin_range(double nsigma = 0.0) const;

            double get_nsigma() const { return m_nsigma; };

           private:
            const Pimpos& m_pimpos;
            const IDFT::pointer& m_dft;
            const Binning& m_tbins;

            double m_nsigma;
            IRandom::pointer m_fluctuate;
            ImpactDataCalculationStrategy m_calcstrat;

            // current window set by user.
            std::pair<int, int> m_window;
            // the content of the current window
            std::map<int, ImpactData::mutable_pointer> m_impacts;
            // std::vector<std::shared_ptr<GaussianDiffusion> > m_diffs;
            // This GausDiffTimeCompare works in BinnedDiffusion_transform but not tested
            // in BinnedDiffusion, details in issue #447
            std::set<std::shared_ptr<GaussianDiffusion>, GausDiffTimeCompare> m_diffs;
            // std::set<std::shared_ptr<GaussianDiffusion> > m_diffs;

            int m_outside_pitch;
            int m_outside_time;
        };

    }  // namespace Gen

}  // namespace WireCell

#endif
