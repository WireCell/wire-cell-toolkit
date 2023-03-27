// These intervals are for numberic ranges, eg for error propagation.
// https://www.boost.org/doc/libs/1_80_0/libs/numeric/interval/doc/interval.htm
#include <boost/numeric/interval.hpp>
#include <iostream>

template <class T, class Policies>
std::ostream &operator<<(std::ostream &os, const boost::numeric::interval<T, Policies> &x)
{
    os << "[" << x.lower() << ", " << x.upper() << "]";
    return os;
}

using namespace std;
int main()
{
    typedef boost::numeric::interval<double> I;
    I a = 42.0, b{5.0, 7.0};
    I c = a / b;

    cout << a << " / " << b << " = " << c << endl;

    return 0;
}
