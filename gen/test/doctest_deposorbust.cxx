// DeposOrBust routes an empty depo set to its frame port as an empty frame.
// Frame consumers (eg Aux::dump_frame via DumpFrames) dereference
// IFrame::traces(), so the empty frame must give an empty trace vector, not
// a null pointer.

#include "WireCellGen/DeposOrBust.h"
#include "WireCellAux/SimpleDepoSet.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell;

TEST_CASE("gen deposorbust empty frame has traces")
{
    Gen::DeposOrBust dob;

    Gen::DeposOrBust::input_queues iqs;
    Gen::DeposOrBust::output_queues oqs;
    std::get<0>(iqs).push_back(std::make_shared<Aux::SimpleDepoSet>(42, IDepo::vector{}));

    REQUIRE(dob(iqs, oqs));
    CHECK(std::get<0>(iqs).empty());
    CHECK(std::get<0>(oqs).empty());
    REQUIRE(std::get<1>(oqs).size() == 1);

    auto frame = std::get<1>(oqs).front();
    REQUIRE(frame);
    CHECK(frame->ident() == 42);
    auto traces = frame->traces();
    REQUIRE(traces);
    CHECK(traces->empty());
}
