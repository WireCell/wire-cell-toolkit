#include "WireCellUtil/Testing.h"
#include <boost/filesystem.hpp>
#include <iostream>
using namespace boost::filesystem;

int main(int argc, char* argv[])
{
    path me(argv[0]);
    Assert(exists(me));
    Assert(is_regular_file(me));
    auto p = me.parent_path();
    auto f = me.filename();
    auto me2 = p / f;
    Assert(me == me2);
    std::cerr << me << std::endl;
    std::cerr << absolute(me).string() << std::endl;
    std::cerr << canonical(me).string() << std::endl;
    
    path you = "path/to/you/";
    std::cerr << you << std::endl;
    std::cerr << absolute(you).string() << std::endl;
    try {
        std::cerr << canonical(you).string() << std::endl;
    }
    catch (const filesystem_error& err) {
        std::cerr << "no canonical version of " << you << std::endl;
    }

    std::cerr << you / "another" << std::endl;

    return 0;
}
