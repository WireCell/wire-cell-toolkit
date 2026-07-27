#include "WireCellSpng/Stream.h" // for write() to ostream
#include "WireCellUtil/Stream.h" // for boost iostreams

#include "WireCellSpng/Testing.h"


using namespace WireCell;
using namespace WireCell::Stream;

TEST_CASE("spng torch stream")
{
    auto ten = torch::rand({5,2});
    auto dt = dtype_of(ten);
    REQUIRE(dt == "<f4");
    auto to = tensoroptions_from_dtype(dt);
    REQUIRE(to.dtype() == ten.dtype());

    auto shape = shape_of(ten);
    REQUIRE(shape[0] == 5);
    REQUIRE(shape[1] == 2);
    auto here = ten.to(torch::kCPU).contiguous();
    REQUIRE(here.data_ptr() != nullptr);
    REQUIRE(here.element_size() == 4);
    REQUIRE(here.numel() == 5*2);


    // Fixme: locate this in a tempdir
    std::string archive = "test-torch-stream.npz";

    boost::iostreams::filtering_ostream out;
    output_filters(out, archive);
    REQUIRE(out.size() >= 1);
    const std::string fname = "rand.npy";
    write(out, fname, ten);
    REQUIRE(out);
    out.flush();
    out.pop();

    boost::iostreams::filtering_istream in;
    input_filters(in, archive);
    REQUIRE(in.size() >= 1);

    torch::Tensor ten2;
    std::string fname2;

    read(in, fname2, ten2);
    REQUIRE(in);
    REQUIRE(fname2 == fname);
    REQUIRE(torch::all(ten == ten2).item<bool>());
}
