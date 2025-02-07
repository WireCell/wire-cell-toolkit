#include "WireCellUtil/doctest.h"
#include "WireCellUtil/LMN.h"

using namespace WireCell;

TEST_CASE("lmn gcd")
{
    CHECK(4 == LMN::gcd(512.0, 500.0));
    CHECK(4 == LMN::gcd(500.0, 512.0));
    CHECK(4 == LMN::gcd(64.0, 100.0));
}

TEST_CASE("lmn rational")
{
    // upsample from 512ns to 500ns in units of ns and us
    CHECK(125 == LMN::rational(512.0, 500.0));
    CHECK(125 == LMN::rational(0.512, 0.5));

    // downsample from 500ns to 512ns in units of ns and us
    CHECK(128 == LMN::rational(500.0, 512.0));
    CHECK(128 == LMN::rational(0.5, 0.512));

    // etc between 64 and 100 ns (relevant to field response)
    CHECK(16 == LMN::rational(100.0, 64.0));
    CHECK(25 == LMN::rational(64.0, 100.0));
}

TEST_CASE("lmn nhalf")
{
    CHECK(1 == LMN::nhalf(3));
    CHECK(1 == LMN::nhalf(4));
    CHECK(2 == LMN::nhalf(5));
    CHECK(2 == LMN::nhalf(6));
}

TEST_CASE("lmn resample")
{
    using namespace std::complex_literals;
    Array::array_xxc sampled(1, 5);
    sampled << 0.0f+0.0if, 2.0f+1.0if, 1.0f+0.5if, 1.0f-0.5if, 2.0f-1.0if;

    Array::array_xxc us = LMN::resample(sampled, 7);
    CHECK(us.cols() == 7);
    CHECK(us(0,0) == 0.0f);
    CHECK(us(0,3) == 0.0f);
    Array::array_xxc ds = LMN::resample(sampled, 3);
    CHECK(ds.cols() == 3);
    CHECK(ds(0,0) == 0.0f);
    CHECK(ds(0,1) == 2.0f+1.0if);
    
}
