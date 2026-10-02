// Tests of the FD-VD ophit10ppm port (FdvdOpHitFinder): the numpy-order sum,
// the pulse rules (MinWidth on end - start, peak threshold, pedestal from the
// first 3 samples) and the tensor-set interface (column order, ns units).

#include "WireCellUtil/doctest.h"

#include "WireCellFlash/FdvdOpHitFinder.h"
#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"

#include <cmath>
#include <cstdint>
#include <vector>

using namespace WireCell;

TEST_CASE("fdvd numpy-order sum")
{
    // a[i] = sin(i) 10^(i % 7): naive left-to-right sums differ from numpy's for n = 100, 300;
    // expected values are numpy 2.1.1 a.sum() (hex, exact).
    for (auto [n, want] : std::vector<std::pair<size_t, double>>{
             {300, 0x1.ab4159e844d66p+20}, {100, 0x1.3e62c108e40d8p+21}, {7, -0x1.75aadf1d354e2p+18}}) {
        std::vector<double> a(n);
        for (size_t i = 0; i < n; ++i) a[i] = std::sin((double) i) * std::pow(10.0, (double) (i % 7));
        CHECK(Flash::fdvd_numpy_sum(a.data(), n) == want);
    }
}

// one snippet: pedestal 500 (first 3 samples), then a flat pulse of `width + 1` samples at 500 + amp
static void add_snippet(std::vector<int>& ch, std::vector<double>& t, std::vector<int64_t>& off,
                        std::vector<double>& adc, int c, double t0, int width, double amp)
{
    if (off.empty()) off.push_back(0);
    ch.push_back(c);
    t.push_back(t0);
    for (int i = 0; i < 10; ++i) adc.push_back(500);
    for (int i = 0; i <= width; ++i) adc.push_back(500 + amp + (i == width / 2 ? 1 : 0));
    for (int i = 0; i < 10; ++i) adc.push_back(500);
    off.push_back((int64_t) adc.size());
}

TEST_CASE("fdvd ophit pulse rules")
{
    std::vector<int> ch;
    std::vector<double> t;
    std::vector<int64_t> off;
    std::vector<double> adc;
    Flash::FdvdOpHitParams par;
    add_snippet(ch, t, off, adc, 11, 100.0, 59, 20.0);   // end - start = 59 < MinWidth 60: no hit
    add_snippet(ch, t, off, adc, 12, 200.0, 60, 20.0);   // = 60: a hit
    add_snippet(ch, t, off, adc, 13, 300.0, 80, 10.0);   // peak 11 < 15: no hit
    auto hits = Flash::fdvd_find_ophits(ch, t, off, adc, par);
    REQUIRE(hits.size() == 1);
    const auto& h = hits[0];
    CHECK(h.channel == 12);
    CHECK(h.amplitude == 21.0);
    CHECK(h.start_time_us == 200.0 + par.tick_us * 10.0);
    CHECK(h.peak_time_us == 200.0 + par.tick_us * 40.0);         // first maximum, sample 10 + 30
    CHECK(h.width_us == par.tick_us * 60.0);
    CHECK(h.area == 61 * 20.0 + 1.0);
    CHECK(h.pe == h.area / 130.0 + 0.43);
}

TEST_CASE("fdvd ophit tensor interface")
{
    std::vector<int> ch;
    std::vector<double> t;
    std::vector<int64_t> off;
    std::vector<double> adc;
    add_snippet(ch, t, off, adc, 21, -5.0, 70, 30.0);
    std::vector<int32_t> ch32(ch.begin(), ch.end());
    std::vector<int16_t> adc16(adc.begin(), adc.end());
    auto ten = [](auto& v, const char* name) {
        Configuration md;
        md["name"] = name;
        return std::make_shared<Aux::SimpleTensor>(ITensor::shape_t{v.size()}, v.data(), md);
    };
    auto tv = std::make_shared<ITensor::vector>();
    tv->push_back(ten(ch32, "ch"));
    tv->push_back(ten(t, "t"));
    tv->push_back(ten(off, "off"));
    tv->push_back(ten(adc16, "adc"));
    auto in = std::make_shared<Aux::SimpleTensorSet>(7, Configuration{}, ITensor::shared_vector(tv));

    Flash::FdvdOpHitFinder hf;
    hf.configure(hf.default_configuration());
    ITensorSet::pointer out;
    REQUIRE(hf(in, out));
    REQUIRE(out);
    CHECK(out->ident() == 7);
    auto o = out->tensors()->at(0);
    CHECK(o->metadata()["name"].asString() == "ophits");
    REQUIRE(o->shape() == ITensor::shape_t{1, 9});
    const double* r = (const double*) o->data();
    CHECK(r[0] == 21.0);
    CHECK(r[1] == (-5.0 + 0.016 * 45.0) * 1000.0);   // peak [ns]
    CHECK(r[2] == 0.016 * 70.0 * 1000.0);            // width [ns]
    CHECK(r[4] == 31.0);                             // amplitude
    CHECK(r[6] == (-5.0 + 0.016 * 10.0) * 1000.0);   // start [ns]
    CHECK(r[7] == -1.0);

    ITensorSet::pointer eos;
    CHECK(hf(nullptr, eos));
    CHECK(!eos);
}
