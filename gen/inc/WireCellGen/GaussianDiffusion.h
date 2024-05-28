#ifndef WIRECELLGEN_GAUSSIANDIFFUSION
#define WIRECELLGEN_GAUSSIANDIFFUSION

#include "WireCellUtil/Array.h"
#include "WireCellUtil/Binning.h"
#include "WireCellIface/IDepo.h"
#include "WireCellIface/IRandom.h"

#include <memory>

namespace WireCell {
    namespace Gen {

        /** A GausDesc describes a Gaussian distribution.
         *
         * Two are used by GaussianDiffusion.  One describes the
         * transverse dimension along the direction of wire pitch (and
         * for a given wire plane) and one the longitudinal dimension
         * is along the drift direction as measured in time.  */
        struct GausDesc {
            /// The absolute location of the mean of the Gaussian as
            /// measured relative to some externally defined origin.
            double center;
            /// The Gaussian sigma (half) width.
            double sigma;

            GausDesc(double center, double sigma)
              : center(center)
              , sigma(sigma)
            {
            }

            /// Return the distance in number of sigma that x is from the center
            double distance(double x)
            {
                double ret = 0.0;
                if (!sigma) {
                    ret = x - center;
                }
                else {
                    ret = (x - center) / sigma;
                }
                return ret;
            }

            std::pair<double, double> sigma_range(double nsigma = 3.0)
            {
                return std::make_pair(center - sigma * nsigma, center + sigma * nsigma);
            }

            /** Sample the Gaussian at points on a uniform linear grid. */
            std::vector<double> sample(double start, double step, int nsamples) const;

            /** Integrate Gaussian across uniform bins.  Result is
             * normalized assuming integral of Gaussian over entire
             * domain is 1.0. */
            std::vector<double> binint(double start, double step, int nbins) const;

            /** Integrate Gaussian diffusion with linear weighting
             *  to redistribute the charge to the two neartest impact positions
             *  for linear interpolation of the field response */
            std::vector<double> weight(double start, double step, int nbins, std::vector<double> pvec) const;
        };

        class GaussianDiffusion {
           public:
            typedef std::shared_ptr<GaussianDiffusion> pointer;

            /// A patch is a 2D array of diffuse charge in (nimpacts X
            /// nticks) bins.  patch[0] is drifted/diffused charge
            /// waveform at impact position 0 (relative to min pitch
            /// for the patch).  patch[0][0] is charge at this impact
            /// position at time = 0 (relative to min time for the
            /// patch).  See `bin()`.
            typedef Array::array_xxf patch_t;

            /** Create a diffused deposition.
             */

            GaussianDiffusion(const IDepo::pointer& depo, const GausDesc& time_desc, const GausDesc& pitch_desc);

            /// This fills the patch once matching the given time and
            /// pitch binning. The patch is limited to the 2D sample
            /// points that cover the subdomain determined by the
            /// number of sigma.  If there should be Poisson
            /// fluctuations applied pass in an IRandom.  Total charge
            /// is charge is preserved.  Each cell of the patch
            /// represents the 2D bin-centered sampling of the
            /// Gaussian.

            void set_sampling(const Binning& tbin, const Binning& pbin, double nsigma = 3.0,
                              IRandom::pointer fluctuate = nullptr,
                              unsigned int weightstrat = 1 /*see BinnedDiffusion ImpactDataCalculationStrategy*/);
            void clear_sampling();

            /// Get the diffusion patch as an array of N_pitch rows X
            /// N_time columns.  Index as patch(i_pitch, i_time).
            /// Call set_sampling() first.
            const patch_t& patch() const;

            const std::vector<double> weights() const;

            /// Return the absolute time bin in the binning corresponding to column 0 of the patch.
            int toffset_bin() const { return m_toffset_bin; }

            /// Return the absolute impact bin in the binning corresponding to column 0 of the patch.
            int poffset_bin() const { return m_poffset_bin; }

            /// Access deposition.
            IDepo::pointer depo() const { return m_deposition; }

            double depo_time() const { return m_deposition->time(); }
            double depo_x() const { return m_deposition->pos().x(); }

            const GausDesc pitch_desc() { return m_pitch_desc; }
            const GausDesc time_desc() { return m_time_desc; }

           private:
            IDepo::pointer m_deposition;  // just for provenance

            GausDesc m_time_desc, m_pitch_desc;

            patch_t m_patch;
            std::vector<double> m_qweights;

            int m_toffset_bin;
            int m_poffset_bin;
        };
    }  // namespace Gen
}  // namespace WireCell

#endif
