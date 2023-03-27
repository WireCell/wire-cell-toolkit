#include "WireCellUtil/Binning.h"
#include "WireCellUtil/Testing.h"

#include <vector>
#include <iostream>

using namespace WireCell;

void test_subset()
{
    Binning a(10, 0, 10);
    std::cerr << "a: " << a << "\n";
    {
        auto s = subset(a, 0.5, 10);
        std::cerr << "s1: " << s << "\n";
        Assert(s.nbins()   == a.nbins());
        Assert(s.binsize() == a.binsize());
        Assert(s.min()     == a.min());
        Assert(s.max()     == a.max());
    }
    {
        auto s = subset(a, 1, 10);
        std::cerr << "s2: " << s << "\n";
        Assert(s.nbins()   == a.nbins() - 1);
        Assert(s.binsize() == a.binsize());
        Assert(s.min()     == a.min() + 1);
        Assert(s.max()     == a.max());
    }
    {
        auto s = subset(a, -1, 10);
        std::cerr << "s3: " << s << "\n";
        Assert(s.nbins()   == a.nbins());
        Assert(s.binsize() == a.binsize());
        Assert(s.min()     == a.min());
        Assert(s.max()     == a.max());
    }
    {
        Binning b(122, 1.528e+06, 1.772e+06);
        std::cerr << "b: " << b << "\n";
        auto s = subset(b, 1791637.8683138976, 1795821.8405835067);
        std::cerr << "s: " << s << "\n";
        Assert(s.nbins() == 0);
        
    }
}

void test_gaussian(const Binning bins, double mean=0, double sigma=1, double expected=1)
{
    std::cerr << "bins: " << bins << " " << "mean="<< mean << " sigma=" << sigma << "\n";
    std::vector<double> g(bins.nbins(), 0);
    const size_t nbins = bins.nbins();
    const double total = gaussian(g.begin(), bins, mean, sigma);
    double sum = 0;
    for (size_t ind=0; ind<nbins; ++ind) {
        const double got = g[ind];
        sum += got;
        std::cerr << ind <<": [" << bins.edge(ind) << "] " << got << " [" << bins.edge(ind+1) << "]\n";
    }
    const double diff = total-sum;
    std::cerr << "total=" << total << " sum=" << sum << " diff=" << diff << "\n";
    Assert(std::abs(diff) < 1e-6);
    const double err = expected-total;
    std::cerr << "total=" << total << " expect=" << expected << " err=" << err << "\n";
    Assert(std::abs(err) < 1e-6);
}

void test_gcumulative(double L=0, double mean=0, double sigma=1, double expected=1)
{
    const double got = gcumulative(L,mean,sigma);
    const double err = expected-got;
    std::cerr << "L=" << L << " mean="<< mean << " sigma=" << sigma 
              << " got=" << got << " expect=" << expected << " err=" << err << "\n";
    Assert(std::abs(err) < 1e-6);
}

void test_gbounds(double L1=-1.0, double L2=1.0, double mean=0, double sigma=1, double expected=0.682689)
{
    const double got = gbounds(L1,L2,mean,sigma);
    const double err = expected-got;
    std::cerr << "L1=" << L1 << " L2=" << L2 << " mean="<< mean << " sigma=" << sigma 
              << " got=" << got << " expect=" << expected << " err=" << err << "\n";
    Assert(std::abs(err) < 1e-6);
}
int main()
{
    test_subset();
    test_gcumulative(0, 0, 1, 0.5);
    test_gcumulative(10, 0, 1, 1.0);
    test_gcumulative(5, 5, 30, 0.5);
    test_gcumulative(300, 5, 30, 1.0);
    test_gbounds();
    test_gbounds(1,-1);
    test_gbounds( 0, 10, 0, 1, 0.5);
    test_gbounds(-10, 0, 0, 1, 0.5);
    test_gaussian(Binning(20,-10,10), 0, 1, 1);

    return 0;
}
