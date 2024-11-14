#include "WireCellUtil/Measurement.h"
#include <cassert>
#include <iostream>

using namespace WireCell::Measurement;

float32 downacast(float64 dd)
{
    return float32(dd);
}


int main()
{
    {
        float32 f1;
        float32 f2(42);
        // float32 f3(42,6.9);
        float64 d1;
        float64 d2(42);
        float64 d3(42,6.9);

        assert(f1 == 0);
        assert(!f1);
        assert(d1 == 0);
        assert(!d1);
        assert(f2.uncertainty() == 0);
        assert(d2.uncertainty() == 0);

        f1 = 11;
        d1 = 12;
        d2 = d3;

        // No cross-type support.
        // d1 = f1;
    }
    {   // not-too-large ints should be held exact
        float  x1 [[maybe_unused]] = float32(12345678);
        double x2 [[maybe_unused]] = float64(12345678);
        assert (x1 == x2);
    }
    {   // even larger for doubles
        double y1 [[maybe_unused]] = float64(1234567890);
        double y2 [[maybe_unused]] = float64(1234567890);
        assert (y1 == y2);
    }

    {   // basic arithmetic
        float32 a(10,0.1);
        float32 b(20,0.1);
        std::cerr << "a=" << a << " b=" << b
                  << " +:" << a+b
                  << " -:" << a-b
                  << " *:" << a*b
                  << " /:" << a/b << "\n";
    }

    return 0;
}
