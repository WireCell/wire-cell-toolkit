/** ChargeSolving solves for charge in blobs.

    This attempts to faithfully implement the algorithms from the
    original Wire-Cell Prototype version:
    
    https://github.com/BNLIF/wire-cell-2dtoy/blob/master/src/ChargeSolving.cxx

    It includes these features:

    - uses channel charge and its uncertainty.

    - configurable threshold on the amount of charge for a measurement to be considered.

    - uses b-m connectivity to separate the solving to smaller problems.

    - solving uses b-weights based on b-b connectivity.

    - a Cholsky decomposition "whitening" of the matrices to accomodate LassoModel solver issues.

    https://en.wikipedia.org/wiki/Whitening_transformation

 */
#ifndef WIRECELL_CHARGESOLVING_H
#define WIRECELL_CHARGESOLVING_H

#include "WireCellIface/IClusterFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"

namespace WireCell {

    namespace Img {

        class ChargeSolving : public Aux::Logger, public IClusterFilter, public IConfigurable {
          public:
            ChargeSolving();
            virtual ~ChargeSolving();

            virtual void configure(const WireCell::Configuration& cfg);
            virtual WireCell::Configuration default_configuration() const;

            virtual bool operator()(const input_pointer& in, output_pointer& out);

          private:

            // Used to hold measurement and blob values
            // (central+uncertainty).
            using value_t = ISlice::value_t;

            // Config: meas_{value,error}_threshold.  The value
            // threshold is a lower bound and the error threshold is
            // an upper bound to which measures are accepted for
            // solving.  The bound is inclusive.  A threshold of 0
            // will retain measures of value 0.  Note, this threshold
            // is placed on the sum of channels in a measure, not on
            // the individual channel components.  Units are numbers
            // of electrons.
            value_t m_meas_thresh{10,1e9};

            // Config: blob_{value,error}_threshold. The value
            // threshold is a lower bound on the solved blob charge
            // for that blob to be output.  The value threshold is
            // inclusive.  Blobs with signal value equal or larger
            // than the threshold are output.  Note, the uncertainty
            // is not currently considered.
            value_t m_blob_thresh{-1,0};

            // config. The tolerance passed to LassoModel solver.
            double m_lasso_tolerance{1e-3};

            // config. The minimum norm passed to LassoModel solver.
            double m_lasso_minnorm{1e-6};

            double m_lasso_lambda_scale{1.0};

            // config. The blob weighting strategies for each of a
            // number of solving rounds.  This is an arbitrary long
            // sequence of strategy names from the set: "uniform" (all
            // blob weights same), "simple" (weights based on simple
            // overlap connectivity measure), "distance" (weights
            // based on a distance between neighbors).  Each strategy
            // will correspond to one round of solving with that
            // strategy followed by a selection based on
            // m_minimum_charge.
            //std::vector<std::string> m_weighting_strategies{"uniform", "simple"};
            std::vector<std::string> m_weighting_strategies{"uboone"};
            
            // controls ress param choice and one blob graph solving
            // uboone, simple
            std::string m_solve_config{"uboone"};

            // Config: if false, any measurement uncertainty is
            // ignored.  If true, the measurement covariance matrix is
            // formed on the measurement uncertainties and is
            // "whitened" via Cholesky decomposition.
            bool m_whiten{true};

            int m_count{0};
        };

    }  // namespace Img

}  // namespace WireCell

#endif /* WIRECELL_CHARGESOLVING_H */
