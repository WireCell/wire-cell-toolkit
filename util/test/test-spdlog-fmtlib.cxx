/**
   Test for fmtlib's change in how it handles (or doesn't handle) implicitly
   using ostream insertion operator<< starting at fmtlib 9.

   Test like:

   spack install spdlog@1.8.2+shared ^fmt@7.1.3
   spack install spdlog@1.9.2+shared ^fmt@8.1.1
   spack install spdlog@1.10.0+shared ^fmt@8.1.1
   spack install spdlog@1.11.0+shared ^fmt@9.1.0
   spack install spdlog@1.13.0+shared

   The goal is to force a version of fmtlib near to the version of the bundled
   fmtlib for a given spdlog version.

   For each,

   spack view add -i spdlog-<version>/local spdlog@<version>
   echo 'load_prefix local' > spdlog-<version>/.envrc
   cd spdlog-<version>/
   direnv allow

   g++ -Wall -std=c++17 \
        -o test-spdlog-fmtlib \
        ../test-spdlog-fmtlib.cxx \
        (pkg-config spdlog --cflags --libs | string split ' ')    

        (fish syntax)

   ./test-spdlog-fmtlib 

 */

#include <iostream>
namespace NS {

    struct S { int x; };

    inline std::ostream& operator<<(std::ostream& os, const S& s)
    {
        os << "<S.x=" << s.x << ">";
        return os;
    }
}

#include "WireCellUtil/Spdlog.h"


// The goal with the mess above is to make this "new style" requirement work with old fmtlib.
template <> struct fmt::formatter<NS::S> : fmt::ostream_formatter {};



int main()
{
    NS::S s{42};
    std::cerr << s << "\n";
    spdlog::info("{}", s);
    spdlog::info("FMT_VERSION={}", FMT_VERSION);

    return 0;
}
