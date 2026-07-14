#include "WireCellUtil/Stream.h"
#include "WireCellUtil/Array.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell::Stream;
using namespace WireCell::Array;

const int nrow = 3;
const int ncol = 4;

static void do_writing(std::string fname)
{
    boost::iostreams::filtering_ostream so;
    output_filters(so, fname);

    array_xxf arr = array_xxf::Zero(3, 4);
    std::vector<int> cols;
    for (int icol = 0; icol < ncol; ++icol) {
        cols.push_back(icol);
        for (int irow = 0; irow < nrow; ++irow) {
            arr(irow, icol) = irow + icol * nrow;
        }
    }

    write(so, "twodee.npy", arr);
    REQUIRE(so);
    write(so, "onedee.npy", cols);
    REQUIRE(so);
}

static void do_reading(std::string fname)
{
    boost::iostreams::filtering_istream si;
    input_filters(si, fname);

    array_xxf arr;
    std::vector<int> cols;

    std::string aname;
    read(si, aname, arr);
    spdlog::debug("aname is {}", aname);
    REQUIRE(aname == "twodee.npy");

    read(si, aname, cols);
    REQUIRE(aname == "onedee.npy");

    for (int icol = 0; icol < ncol; ++icol) {
        CHECK(cols[icol] == icol);
        for (int irow = 0; irow < nrow; ++irow) {
            CHECK(arr(irow, icol) == irow + icol * nrow);
        }
    }
}

TEST_SUITE("util stream") {

TEST_CASE("write and read tar.gz archive") {
    std::string fname = "/tmp/doctest-util-stream.tar.gz";
    do_writing(fname);
    do_reading(fname);
}

}  // TEST_SUITE("util stream")
