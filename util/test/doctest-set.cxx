#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"

#include <vector>
#include <set>
#include <deque>

#include <iterator>   // std::back_inserter
#include <algorithm>  // std::set_difference
#include <memory>
#include <random>

using namespace std;

TEST_SUITE("set") {

TEST_CASE("set difference") {
    // http://stackoverflow.com/a/10604500
    vector<int> items = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    set<int> die = {2, 4, 5, 42};
    set<int> die2;

    vector<int> result;
    set_difference(items.begin(), items.end(), die.begin(), die.end(), back_inserter(result));

    spdlog::debug("items size: {}", items.size());
    for (auto item : items) {
        spdlog::debug("  item: {}", item);
    }
    spdlog::debug("alive size: {}", result.size());
    for (auto alive : result) {
        spdlog::debug("  alive: {}", alive);
    }

    vector<int> result2;
    set_difference(items.begin(), items.end(), die2.begin(), die2.end(), back_inserter(result2));

    CHECK(items.size() == 10);
    CHECK(die.size() == 4);
    CHECK(result.size() == 7);

    CHECK(die2.size() == 0);
    CHECK(result2.size() == 10);

    // result should exclude 2, 4, 5 (42 is not in items)
    vector<int> expected = {0, 1, 3, 6, 7, 8, 9};
    CHECK(result == expected);
}

TEST_CASE("deque shared_ptr pop") {
    std::random_device rd;
    std::default_random_engine re(rd());
    std::uniform_real_distribution<> dist(0, 1000);

    typedef std::shared_ptr<int> Pint;
    deque<Pint> queue;
    const int nitems = 1000;
    for (int ind = 0; ind < nitems; ++ind) {
        queue.push_back(Pint(new int(dist(re))));
    }
    CHECK(static_cast<int>(queue.size()) == nitems);
    for (int ind = 0; ind < nitems; ++ind) {
        Pint front = queue.front();
        queue.pop_front();
        spdlog::debug("{}: popped:{} now with: {} items", ind, *front, queue.size());
    }
    CHECK(queue.empty());
}

}  // TEST_SUITE("set")
