#include "WireCellGen/Drifter.h"

#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/doctest.h"

#include <random>

using namespace WireCell;
using namespace WireCell::Gen;
using spdlog::debug;
using spdlog::info;

class DumbDepo : public WireCell::IDepo {
public:
    DumbDepo(const WireCell::Point& pos) : m_pos(pos) {};
    virtual ~DumbDepo() {};

    virtual const WireCell::Point& pos() const { return m_pos; }
    virtual double time() const { return 0; }
    virtual double charge() const { return 1000; }
    virtual double energy() const { return 1000; }
    virtual int id() const { return 0; }
    virtual int pdg() const { return 0; }
    virtual WireCell::IDepo::pointer prior() const { return nullptr; }
    virtual double extent_long() const { return 0; }
    virtual double extent_tran() const { return 0; }

private:
    WireCell::Point m_pos;
};


static void common_setup()
{
    PluginManager& pm = PluginManager::instance();
    WireCell::Plugin* pi = pm.add("WireCellGen");
    REQUIRE(pi != nullptr);
    {
        auto icfg = Factory::lookup<IConfigurable>("Random");
        auto cfg = icfg->default_configuration();
        icfg->configure(cfg);
    }

}

static void test_flat(double sign)
{
    Drifter drifter;
    Configuration cfg = drifter.default_configuration();
    const double cathode = 2*sign*units::meter;
    cfg["xregions"][0]["cathode"] = cathode;
    const double response = 10*sign*units::cm;
    cfg["xregions"][0]["response"] = response;
    const double anode = 0;
    cfg["xregions"][0]["anode"] = anode;
    drifter.configure(cfg);

    const double step = 1*units::mm;
    IDrifter::output_queue drifted;
    size_t nin=0;
    for (double x = cathode-10*step; x <= cathode+10*step; x += step) {
        IDrifter::output_queue depos;
        const Point pt(x,0,0);
        auto indepo = std::make_shared<DumbDepo>(pt);
        drifter(indepo, depos);
        drifted.insert(drifted.end(), depos.begin(), depos.begin());
        ++nin;
        // debug("nin={} x={}, drifted out={}", nin, x, depos.size());
    }
    {   // flush
        IDrifter::output_queue depos;
        drifter(nullptr, depos);
        drifted.insert(drifted.end(), depos.begin(), depos.end());
        debug("flush drifted={} depos={}", drifted.size(), depos.size());
    }
    debug("nin={} nout={}", nin, drifted.size());
    CHECK(drifted.size() > 0);

    for (const auto& idepo : drifted) {
        if (!idepo) { continue; } // skip EOS
        auto pos = idepo->prior()->pos();
        if (sign > 0)
            CHECK(pos.x() < cathode);
        else
            CHECK(pos.x() > cathode);
    }


}

TEST_CASE("drifter flat cathode")
{
    common_setup();
    test_flat(+1);
    test_flat(-1);
}

static std::vector<double> uniform_values(size_t seed, size_t number = 3, 
                                          double vmin=-1.0, double vmax=1.0)
{
    std::default_random_engine re{seed};
    std::uniform_real_distribution<double> dist(vmin, vmax);

    std::vector<double> vals(number);
    for (size_t ind=0; ind<number; ++ind) {
        vals[ind] = dist(re);
    }
    return vals;
}

static void test_bent(double sign)
{
    Drifter drifter;
    Configuration cfg = drifter.default_configuration();
    const double cathode = 2*sign*units::meter;
    const double half = units::meter;
    const double scale = 2*half;
    const double yzstep = 10*units::cm;
    const double mag = 1*units::cm;
    const double xmax = 1*units::cm;
    const double xstep = 1*units::mm;
    const double epsilon = 0.1*units::mm;

    auto catx = [&](double x, double y) -> double {
        return cathode + sign * mag * ((x*x + y*y)/(scale*scale));
    };
    const double response = 10*sign*units::cm;
    cfg["xregions"][0]["response"] = response;
    const double anode = 0;
    cfg["xregions"][0]["anode"] = anode;

    int sample_index=0;
    for (double y=-half; y<=half; y+=yzstep) {
        for (double z=-half; z<=half; z+=yzstep) {
            cfg["xregions"][0]["cathode"]["x"][sample_index] = catx(y,z);
            cfg["xregions"][0]["cathode"]["y"][sample_index] = y;
            cfg["xregions"][0]["cathode"]["z"][sample_index] = z;
            ++sample_index;
        }
    }
    drifter.configure(cfg);

    IDrifter::output_queue drifted;
    for (double y : uniform_values(12345, 100, -half, half)) {
        for (double z : uniform_values(12345, 100, -half, half)) {
            for (double x = cathode-xmax; x <= cathode+xmax; x += xstep) {
                IDrifter::output_queue depos;
                const Point pt(x,y,z);
                auto indepo = std::make_shared<DumbDepo>(pt);
                drifter(indepo, depos);
                drifted.insert(drifted.end(), depos.begin(), depos.begin());
            }
        }
    }
    {   // flush
        IDrifter::output_queue depos;
        drifter(nullptr, depos);
        drifted.insert(drifted.end(), depos.begin(), depos.end());
        debug("flush drifted={} depos={}", drifted.size(), depos.size());
    }



    for (const auto& idepo : drifted) {
        if (!idepo) { continue; } // skip EOS
        auto pos = idepo->prior()->pos();

        double xcat = catx(pos.y(), pos.z());

        if (sign > 0)
            CHECK(pos.x() < (xcat + epsilon));
        else
            CHECK(pos.x() > (xcat - epsilon));
    }

}

TEST_CASE("drifter bent cathode")
{
    common_setup();
    test_bent(+1);
    test_bent(-1);
}

TEST_CASE("drifter null face")
{
    common_setup();

    Drifter drifter;
    Configuration cfg = drifter.default_configuration();
    cfg["xregions"][0]["cathode"] = 2*units::meter;
    cfg["xregions"][0]["response"] = 10*units::cm;
    cfg["xregions"][0]["anode"] = 0;
    cfg["xregions"][1] = Json::nullValue;
    drifter.configure(cfg);

    struct Trial { double x; bool in; };

    std::vector<Trial> trials = {
        {3*units::meter, false}, // behind cathode
        {1*units::meter, true},  // inside drift region
        {5*units::cm, true},     // in response region
        {-5*units::cm, false},   // behind anode
    };

    for (auto& [x,want] : trials) {
        const Point pt(x,0,0);
        auto depo = std::make_shared<DumbDepo>(pt);
        bool got = drifter.insert(depo);
        debug("test point x={} want={} got={}", x, want, got);
        CHECK(got == want);
    }
}
