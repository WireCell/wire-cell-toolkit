#ifndef WIRECELL_ISCALARFUNCTION
#define WIRECELL_ISCALARFUNCTION

#include "WireCellUtil/IComponent.h"

namespace WireCell {


    /** A scalar function maps R -> R.
     *
     * See also IWaveform which provides a similar concept but for a regularly
     * sampled function.
     */
    class IScalarFunction : public IComponent<IScalarFunction> {
    public:

        virtual ~IScalarFunction() {};

        /// Implementation provides the function. 
        virtual double scalar_function(double x) = 0;
        
    };
}

#endif
