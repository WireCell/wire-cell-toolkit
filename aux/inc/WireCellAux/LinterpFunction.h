#include "WireCellIface/IScalarFunction.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Interpolate.h"

namespace WireCell::Aux {

    /** A scalar function implemented as linear interpolation on regularly spaced points.
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
        /// Array of regularly sampled values.
        ///
        /// Config: start
        ///
        /// The first abscissa value.
        ///
        /// Config: step
        ///
        WireCell::linterp<double> m_terp;

    };

    /** A scalar function implemented as linear interpolation on irregularly spaced points.
     *
     * 
     */
    class IrrterpFunction : public WireCell::IScalarFunction, WireCell::IConfigurable {
    public:
        virtual ~IrrterpFunction();
        
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const {
            return Configuration{}; // user must supply all
        }

        virtual double scalar_function(double x);

    private:
        /// Config: values
        ///
        /// Array of regularly sampled values.
        ///
        /// Config: coords
        ///
        /// Array of the coordinates at which values are sampled.

        WireCell::irrterp<double> m_terp;

    };
}
