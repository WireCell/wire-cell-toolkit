#include "WireCellGen/Drifter.h"

#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/doctest.h"

#include <algorithm>
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
    std::default_random_engine re{static_cast<std::default_random_engine::result_type>(seed)};
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
    const double epsilon = 1*units::mm;

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
                drifted.insert(drifted.end(), depos.begin(), depos.end());
            }
        }
    }
    {   // flush
        IDrifter::output_queue depos;
        drifter(nullptr, depos);
        drifted.insert(drifted.end(), depos.begin(), depos.end());
        debug("flush drifted={} depos={}", drifted.size(), depos.size());
    }



    int nerrors = 0;

    for (const auto& idepo : drifted) {
        if (!idepo) { continue; } // skip EOS
        auto pos = idepo->prior()->pos();

        double xcat = catx(pos.y(), pos.z());

        if (sign > 0) {
            if (pos.x() > (xcat + epsilon)) {
                ++nerrors;
            }
        }
        else {
            if (pos.x() < (xcat - epsilon)) {
                ++nerrors;
            }
        }
    }
    CHECK(nerrors == 0);

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

class TieDepo : public WireCell::IDepo {
public:
    TieDepo(double t, const WireCell::Point& pos, int id) : m_t(t), m_pos(pos), m_id(id) {};
    virtual ~TieDepo() {};

    virtual const WireCell::Point& pos() const { return m_pos; }
    virtual double time() const { return m_t; }
    virtual double charge() const { return -1000; }
    virtual double energy() const { return 1000; }
    virtual int id() const { return m_id; }
    virtual int pdg() const { return 0; }
    virtual WireCell::IDepo::pointer prior() const { return nullptr; }
    virtual double extent_long() const { return 0; }
    virtual double extent_tran() const { return 0; }

private:
    double m_t;
    WireCell::Point m_pos;
    int m_id;
};

// Drift depos that all reach the response plane at the same time (same t and x, different y) and return the output ids in
// order.  The input order (and with it the allocation order of the drifted depos) is the permutation given.  A late depo
// ripens the tied ones through flush_ripe(), a second tied group is left for the EOS flush().
static std::vector<int> drift_tied(bool stable_order, const std::vector<int>& perm)
{
    Drifter drifter;
    Configuration cfg = drifter.default_configuration();
    cfg["xregions"][0]["cathode"] = 2*units::meter;
    cfg["xregions"][0]["response"] = 10*units::cm;
    cfg["xregions"][0]["anode"] = 0;
    cfg["fluctuate"] = false;
    if (stable_order) {
        cfg["stable_order"] = true;
    }
    drifter.configure(cfg);

    std::vector<int> ids;
    auto take = [&](const IDrifter::output_queue& q) {
        for (const auto& d : q) {
            if (d) ids.push_back(d->id());
        }
    };
    const double late = 10*units::ms;
    for (double t0 : {0.0, late}) {
        for (int k : perm) {
            IDrifter::output_queue q;
            drifter(std::make_shared<TieDepo>(t0, Point(1*units::meter, k*units::mm, 0), k + (t0 > 0 ? 100 : 0)), q);
            take(q);
        }
    }
    IDrifter::output_queue q;
    drifter(nullptr, q);
    take(q);
    return ids;
}

TEST_CASE("drifter stable order of equal-time depos")
{
    common_setup();

    {
        Drifter drifter;
        CHECK(drifter.default_configuration()["stable_order"].asBool() == false);
    }

    const int n = 20;
    std::vector<int> up(n), down(n), shuffled(n);
    for (int k = 0; k < n; ++k) {
        up[k] = k;
        down[k] = n - 1 - k;
        shuffled[k] = (7 * k) % n;  // 7 and 20 are coprime
    }

    // stable_order=true: the output order is fixed by content (here y), whatever the input / allocation order.
    std::vector<int> want;
    for (int k = 0; k < n; ++k) want.push_back(k);
    for (int k = 0; k < n; ++k) want.push_back(100 + k);
    CHECK(drift_tied(true, up) == want);
    CHECK(drift_tied(true, down) == want);
    CHECK(drift_tied(true, shuffled) == want);

    // default (legacy): the same depos come out, the order of ties is not specified.
    for (const auto& perm : {up, down, shuffled}) {
        auto got = drift_tied(false, perm);
        std::sort(got.begin(), got.end());
        CHECK(got == want);
    }
}
