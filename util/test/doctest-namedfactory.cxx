#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Exceptions.h"

using namespace WireCell;

class ISomeComponent : public virtual WireCell::IComponent<ISomeComponent> {
    Log::logptr_t l;

   public:
    ISomeComponent()
      : l(Log::logger("SomeComponent"))
    {
        l->debug("SomeComponent() at {:p}", (void*) this);
    }
    virtual ~ISomeComponent() { l->debug("~SomeComponent()"); }
    virtual void chirp() = 0;
};

class SomeConcrete : public virtual ISomeComponent {
    Log::logptr_t l;

   public:
    SomeConcrete()
      : l(Log::logger("SomeConcrete"))
    {
        l->debug("SomeConcrete() at {:p}", (void*) this);
    }
    virtual ~SomeConcrete() { l->debug("~SomeConcrete()"); }
    virtual void chirp() { l->info("SomeConcrete::chirp() at {:p}", (void*) this); }
};

WIRECELL_FACTORY(SomeConcrete, SomeConcrete, ISomeComponent)

TEST_SUITE("namedfactory") {

TEST_CASE("lookup registered component") {
    // Cheat: force factory registration since this isn't in a shared library
    make_SomeConcrete_factory();

    auto ins = WireCell::Factory::lookup<ISomeComponent>("SomeConcrete");
    REQUIRE_MESSAGE(ins, "Failed to lookup 'SomeConcrete' with interface 'ISomeComponent'");
    spdlog::info("Got SomeConcrete @ {:p}", (void*) ins.get());
    ins->chirp();
}

TEST_CASE("lookup nonexistent component throws") {
    make_SomeConcrete_factory();

    bool caught = false;
    try {
        auto should_fail = WireCell::Factory::lookup<ISomeComponent>("NothingNamedThis");
    }
    catch (WireCell::FactoryException& e) {
        spdlog::warn(errstr(e));
        spdlog::info("^^^ Successfully failed to lookup a nonexistent component");
        caught = true;
    }
    REQUIRE_MESSAGE(caught, "Failed to throw");
}

TEST_CASE("find_maybe_tn nonexistent named instance returns null") {
    make_SomeConcrete_factory();

    std::shared_ptr<ISomeComponent> ptr = WireCell::Factory::find_maybe_tn<ISomeComponent>("SomeConcrete:doesnotexist");
    REQUIRE_MESSAGE(ptr == nullptr, "Got non null for nonexistent named component");
}

TEST_CASE("find_maybe_tn existing component returns non-null") {
    make_SomeConcrete_factory();

    std::shared_ptr<ISomeComponent> ptr2 = WireCell::Factory::find_maybe_tn<ISomeComponent>("SomeConcrete");
    spdlog::debug("Got ptr2 @ {:p}", (void*) ptr2.get());
    REQUIRE_MESSAGE(ptr2 != nullptr, "Got null for existing component");
}

}  // TEST_SUITE("namedfactory")
