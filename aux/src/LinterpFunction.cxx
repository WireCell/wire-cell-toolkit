#include "WireCellAux/LinterpFunction.h"

namespace WireCell::Aux {

    /** A scalar function implemented as linear interpolation on regularly spaced points.
     *
     * Extrapolation returns end point values.
     */
    LinterpFunction::~LinterpFunction()
    {
    }
        
    void LinterpFunction::configure(const WireCell::Configuration& cfg)
    {
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

        auto values = get<std::vector<double>>(cfg, "values");
        auto start = get<double>(cfg, "start");
        auto step = get<double>(cfg, "step");
        m_terp = linterp<double>(values.begin(), values.end(), start, step);
    }

    double LinterpFunction::scalar_function(double x)
    {
        return m_terp(x);
    }

    IrrterpFunction::~IrrterpFunction()
    {
    }
        
    void IrrterpFunction:: configure(const WireCell::Configuration& cfg)
    {
        /// Config: values
        ///
        /// Array of regularly sampled values.
        ///
        /// Config: coords
        ///
        /// Array of the coordinates at which values are sampled.

        auto values = cfg["values"];
        auto coords = cfg["coords"];
        int npts = values.size();
        for (int ind=0; ind<npts; ++ind) {
            m_terp.add_point(coords[ind].asDouble(), values[ind].asDouble());
        }
    }


    double IrrterpFunction::scalar_function(double x)
    {
        return m_terp(x);
    }

}
