#include "WireCellUtil/Stream.h"
#include "WireCellUtil/Array.h"
#include <iostream>

using namespace WireCell::Stream;
using namespace WireCell::Array;

const int nrow=3;
const int ncol=4;

void do_writing(std::string fname)
{
    boost::iostreams::filtering_ostream so;
    output_filters(so, fname);

    array_xxf arr = array_xxf::Zero(3,4);
    std::vector<int> cols;
    for (int icol=0; icol<ncol; ++icol) {
        cols.push_back(icol);
        for (int irow=0; irow<nrow; ++irow) {
            arr(irow, icol) = irow + icol*nrow;
        }
    }
        
    write(so, "twodee.npy", arr);
    assert(so);
    write(so, "onedee.npy", cols);
    assert(so);
}

void do_reading(std::string fname)
{
    boost::iostreams::filtering_istream si;
    input_filters(si, fname);

    array_xxf arr;
    std::vector<int> cols;

    std::string aname;
    read(si, aname, arr);
    std::cerr << "aname is " << aname <<"\n";
    assert(aname == "twodee.npy");

    read(si, aname, cols);
    assert(aname == "onedee.npy");
    
    for (int icol=0; icol<ncol; ++icol) {
        assert(cols[icol] == icol);
        for (int irow=0; irow<nrow; ++irow) {
            assert(arr(irow, icol) == irow + icol*nrow);
        }
    }
}

int main(int argc, char* argv[])
{
    std::string name = argv[0];
    std::string fname = name + ".tar.gz";

    do_writing(fname);
    do_reading(fname);

    return 0;
}
