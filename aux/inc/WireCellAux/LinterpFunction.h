#include "WireCellIface/IScalarFunction.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Interpolate.h"

namespace WireCell::Aux {

    /** A scalar function implemented as linear interpolation on regularly or
     * irregularly spaced points.
     *
     * Extrapolation returns end point values.
     */
    class LinterpFunction : public WireCell::IScalarFunction, WireCell::IConfigurable {
    public:
        virtual ~LinterpFunction();
        
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const {
            return Configuration{}; // user must supply all
        }

        virtual double scalar_function(double x);
    private:

        /// Config: values
        ///
        /// Array of function samples.
        ///
        /// Config: coords
        ///
        /// Array of the coordinates at which values are sampled.
        /// 
        /// Config: start
        ///
        /// The first abscissa value.
        ///
        /// Config: step
        ///
        /// Distance between regular samples
        ///
        /// Note: coords and (start,step) are mutually exclusive.  The former
        /// implies irregular sample interpolation and the latter regular.

        std::function<double(const double&)> m_terp;
    };

}
