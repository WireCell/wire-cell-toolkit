#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Faninout.h"

#include <boost/signals2.hpp>
#include <boost/optional.hpp>

#include <vector>

using namespace std;
using namespace WireCell;

struct Counter {
    int count;
    Counter(int c = 0)
      : count(c)
    {
    }
    int operator()() { return count++; }
};

struct Echoer {
    int operator()(int x) { return x; }
};

struct Consumer {
    typedef boost::signals2::signal<int()> signal;
    typedef signal::slot_type slot_type;
    signal sig;
    int operator()()
    {
        int val = *sig();
        spdlog::debug("Consumed: {}", val);
        return val;
    }
};

TEST_SUITE("faninout") {

TEST_CASE("plug and play") {
    Counter counter;
    Fanout<int> fanout;
    fanout.connect(counter);

    int addresses[] = {1, 2, 3, 0, 0, 0, 2, 4, 4, 4, 4, 4, 1, 1, 2, 1, 1, -1};
    int want_val[]  = {0, 1, 2, 3, 4, 5, 2, 0, 1, 2, 3, 4, 1, 2, 3, 3, 4, -1};

    fanout.address(4);

    for (int ind = 0; addresses[ind] >= 0; ++ind) {
        int addr = addresses[ind];
        int want = want_val[ind];
        int got = fanout(addr);
        spdlog::debug("addr={} got:{} want:{}", addr, got, want);
        CHECK(want == got);
    }
}

TEST_CASE("ref") {
    Echoer echo;

    Addresser<int> addr0(0), addr1(1);
    addr0.connect(echo);
    addr1.connect(echo);

    Consumer consumer0, consumer1;
    consumer0.sig.connect(boost::ref(addr0));
    consumer1.sig.connect(boost::ref(addr1));

    CHECK(consumer0() == 0);
    CHECK(consumer1() == 1);
}

TEST_CASE("fanout addressing") {
    Counter counter;
    Fanout<int> fanout;
    fanout.address(0);
    fanout.address(1);
    fanout.connect(counter);

    Addresser<int> addr0(0), addr1(1);
    addr0.connect(boost::ref(fanout));
    addr1.connect(boost::ref(fanout));

    Consumer consumer0, consumer1;
    consumer0.sig.connect(boost::ref(addr0));
    consumer1.sig.connect(boost::ref(addr1));

    spdlog::debug("{}", consumer1());
    spdlog::debug("{}", consumer0());
    spdlog::debug("{}", consumer0());
    spdlog::debug("{}", consumer0());
    spdlog::debug("{}", consumer1());
    spdlog::debug("{}", consumer0());
    // TODO: add meaningful assertions
}

TEST_CASE("fanin") {
    Counter c1(10), c2(20), c3(30);

    boost::signals2::signal<int(), Fanin<std::vector<int>>> sig;
    sig.connect(c1);
    sig.connect(c2);
    sig.connect(c3);

    for (int ind = 0; ind < 10; ++ind) {
        auto vec = sig();
        REQUIRE(vec.size() == 3);
        CHECK(vec[0] == 10 + ind);
        CHECK(vec[1] == 20 + ind);
        CHECK(vec[2] == 30 + ind);
        spdlog::debug("ind={} values: {} {} {}", ind, vec[0], vec[1], vec[2]);
    }
}

}  // TEST_SUITE("faninout")
