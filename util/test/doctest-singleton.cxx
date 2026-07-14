#include "WireCellUtil/Singleton.h"
#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"

using spdlog::info;

class Foo {
   public:
    Foo() { info("Foo() at {:p}", (void*) this); }
    virtual ~Foo() { info("~Foo() at {:p}", (void*) this); }
    virtual void chirp() { info("Foo::chirp() at {:p}", (void*) this); }
};

typedef WireCell::Singleton<Foo> OnlyFoo;

TEST_SUITE("singleton") {

TEST_CASE("singleton returns same instance") {
    Foo* foo1 = &OnlyFoo::Instance();
    Foo* foo2 = &OnlyFoo::Instance();
    REQUIRE(foo1 == foo2);

    OnlyFoo::Instance().chirp();
    OnlyFoo::Instance().chirp();

    Foo* foo3 = &WireCell::Singleton<Foo>::Instance();
    CHECK(foo3 == foo1);
}

}  // TEST_SUITE("singleton")
