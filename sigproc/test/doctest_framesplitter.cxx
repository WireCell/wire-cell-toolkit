// FrameSplitter trivially broadcasts a single frame to two outputs.
// Both outputs share identity with the input.

#include "WireCellUtil/doctest.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/NamedFactory.h"

#include "WireCellIface/IFrameSplitter.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleTrace.h"

using namespace WireCell;
using namespace WireCell::Aux;

TEST_CASE("FrameSplitter broadcasts one frame to both outputs") {
    PluginManager& pm = PluginManager::instance();
    pm.add("WireCellSigProc");
    pm.add("WireCellAux");

    ITrace::vector traces { std::make_shared<SimpleTrace>(7, 0, std::vector<float>{1.0f, 2.0f}) };
    auto in = std::make_shared<SimpleFrame>(42, 0.0, traces, 0.5);

    auto fs = Factory::lookup_tn<IFrameSplitter>("FrameSplitter");
    IFrameSplitter::output_tuple_type out;
    REQUIRE(fs->operator()(in, out));

    CHECK(std::get<0>(out).get() == in.get());
    CHECK(std::get<1>(out).get() == in.get());
    CHECK(std::get<0>(out)->ident() == 42);
    CHECK(std::get<1>(out)->ident() == 42);
}

TEST_CASE("FrameSplitter passes EOS (null) on both outputs") {
    PluginManager& pm = PluginManager::instance();
    pm.add("WireCellSigProc");

    auto fs = Factory::lookup_tn<IFrameSplitter>("FrameSplitter");
    IFrameSplitter::output_tuple_type out;
    REQUIRE(fs->operator()(nullptr, out));
    CHECK(!std::get<0>(out));
    CHECK(!std::get<1>(out));
}
