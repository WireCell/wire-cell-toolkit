#include "WireCellSpng/Testing.h"
#include "WireCellSpng/DFT.h"

using namespace WireCell::SPNG;

TEST_CASE("spng dft size prime factor")
{
    CHECK(2 == largest_prime_factor(2));
    CHECK(3 == largest_prime_factor(3));
    CHECK(2 == largest_prime_factor(8));
    CHECK(7 == largest_prime_factor(2*3*5*7));
}
TEST_CASE("spng dft size largest prime")
{
    CHECK(2 == faster_dft_size(2));
    CHECK(3 == faster_dft_size(3));
    CHECK(216 == faster_dft_size(2*3*5*7, 5)); // 210 -> 216

    auto fds = make_faster_dft_size_primes(5);
    REQUIRE(fds);
    CHECK(216 == (*fds)(2*3*5*7)); // 210 -> 216    

}
TEST_CASE("spng dft size largest measured")
{
    auto fds = make_faster_dft_size_measured(torch::kCPU, "this file does not exist");
    CHECK(nullptr == fds);
    
}

TEST_CASE("spng dft size cached")
{
    FasterDftSize dfs;
    CHECK(2 == dfs(2));
    CHECK(3 == dfs(3));

    dfs.reset(make_faster_dft_size_primes(5));
    CHECK(216 == dfs(2*3*5*7)); // 210 -> 216    

}
