#include "WireCellAux/LinterpFunction.h"
#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(LinterpFunction, WireCell::Aux::LinterpFunction,
                 WireCell::IScalarFunction, WireCell::IConfigurable)

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
        auto values = get<std::vector<double>>(cfg, "values");

        if (cfg["coords"].isArray()) {
            auto coords = get<std::vector<double>>(cfg, "coords");
            irrterp<double> terp;
            int npts = values.size();
            for (int ind=0; ind<npts; ++ind) {
                terp.add_point(coords[ind], values[ind]);
            }
            m_terp = terp;
        }
        else {
            auto start = get<double>(cfg, "start");
            auto step = get<double>(cfg, "step");
            m_terp = linterp<double>(values.begin(), values.end(), start, step);
        }
    }

    double LinterpFunction::scalar_function(double x)
    {
        return m_terp(x);
    }


}
