#include "WireCellUtil/TupleHelpers.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell;

TEST_SUITE("tuple") {

TEST_CASE("type names") {
    typedef std::tuple<int, float, double, char, std::string> IFDCS;
    typedef tuple_helper<IFDCS> IFDCS_helper;
    IFDCS_helper ifdcs_helper;
    std::vector<std::string> typenames = ifdcs_helper.type_names();
    CHECK(typenames.size() == 5);
    for (auto tn : typenames) {
        spdlog::debug("{}", tn);
    }
}

TEST_CASE("as_any and from_any") {
    typedef std::tuple<int, float, double, char, std::string> IFDCS;
    typedef tuple_helper<IFDCS> IFDCS_helper;
    IFDCS_helper ifdcs_helper;

    IFDCS ifdcs{1, 2.2f, 3.0e-9, 'a', "foo"};
    std::vector<boost::any> anyvec = ifdcs_helper.as_any(ifdcs);
    CHECK(anyvec.size() == 5);
    CHECK(boost::any_cast<int>(anyvec[0]) == 1);
    CHECK(boost::any_cast<std::string>(anyvec[4]) == "foo");
    spdlog::debug("{}", boost::any_cast<int>(anyvec[0]));
    spdlog::debug("{}", boost::any_cast<std::string>(anyvec[4]));

    auto ifdcs2 = ifdcs_helper.from_any(anyvec);
    CHECK(std::get<0>(ifdcs2) == 1);
    CHECK(std::get<4>(ifdcs2) == std::string("foo"));
}

TEST_CASE("as_any_queue and from_any_queue") {
    typedef std::tuple<int, float, double, char, std::string> IFDCS;
    typedef shared_queued<IFDCS> IFDCS_shqed;
    typedef IFDCS_shqed::shared_queued_tuple_type IFDCS_queues;

    IFDCS_shqed ifdcs_shqed;

    IFDCS_queues qs;
    std::get<0>(qs).push_back(std::make_shared<int>(1));
    std::get<1>(qs).push_back(std::make_shared<float>(2.2f));
    std::get<2>(qs).push_back(std::make_shared<double>(3.0e-9));
    std::get<3>(qs).push_back(std::make_shared<char>('a'));
    std::get<4>(qs).push_back(std::make_shared<std::string>("foo"));

    auto any_q = ifdcs_shqed.as_any_queue(qs);
    CHECK(any_q.size() == 5);
    for (auto q : any_q) {
        CHECK(q.size() == 1);
    }

    auto qs2 = ifdcs_shqed.from_any_queue(any_q);
    spdlog::debug("First element from each queue:");
    spdlog::debug("{}", *std::get<0>(qs2)[0]);
    spdlog::debug("{}", *std::get<1>(qs2)[0]);
    spdlog::debug("{}", *std::get<2>(qs2)[0]);
    spdlog::debug("{}", *std::get<3>(qs2)[0]);
    spdlog::debug("{}", *std::get<4>(qs2)[0]);

    CHECK(*std::get<0>(qs2)[0] == 1);
    CHECK(*std::get<4>(qs2)[0] == std::string("foo"));
}

TEST_CASE("wrapped vector types") {
    typedef std::tuple<int, float, double, char, std::string> IFDCS;
    typedef tuple_helper<IFDCS> IFDCS_helper;
    typedef typename IFDCS_helper::Wrapped<std::vector>::type ifdcs_vectors;
    ifdcs_vectors vs;
    std::get<0>(vs).push_back(1);
    std::get<1>(vs).push_back(2.2f);
    std::get<2>(vs).push_back(3.0e-9);
    std::get<3>(vs).push_back('a');
    std::get<4>(vs).push_back(std::string("foo"));
    CHECK(std::get<0>(vs).size() == 1);
    CHECK(std::get<0>(vs)[0] == 1);
    CHECK(std::get<4>(vs)[0] == std::string("foo"));
}

TEST_CASE("type_repeater") {
    type_repeater<3, std::string>::type ahahah("one", "two", "three");
    spdlog::debug("{} {} {}", std::get<0>(ahahah), std::get<1>(ahahah), std::get<2>(ahahah));
    CHECK(std::get<0>(ahahah) == "one");
    CHECK(std::get<1>(ahahah) == "two");
    CHECK(std::get<2>(ahahah) == "three");

    type_repeater<16, std::string>::type t16("0","1","2","3","4","5","6","7","8","9","a","b","c","d","e","f");
    spdlog::debug("{} {}", std::get<0>(t16), std::get<15>(t16));
    CHECK(std::get<0>(t16) == "0");
    CHECK(std::get<15>(t16) == "f");
}

}  // TEST_SUITE("tuple")
