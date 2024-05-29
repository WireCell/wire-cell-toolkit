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

#include <spdlog/spdlog.h>

#if defined(SPDLOG_FMT_EXTERNAL)
#include <fmt/core.h>
#else
#include <spdlog/fmt/bundled/core.h>
#endif

#if FMT_VERSION >= 90000
#include <fmt/ostream.h>
// This provides fmt::ostream_formatter.
#else
#include <spdlog/fmt/ostr.h>
// This id copied from fmt/ostream.h from fmtlib 10.2 with const_cast tweak.
namespace fmt {
    template <typename Char>
    struct basic_ostream_formatter : formatter<basic_string_view<Char>, Char> {
        void set_debug_format() = delete;

        template <typename T, typename OutputIt>
        auto format(const T& value, basic_format_context<OutputIt, Char>& ctx) const
            -> OutputIt {
            auto buffer = basic_memory_buffer<Char>();
            detail::format_value(buffer, value);
            // Must const_cast as apparently the format() method from fmtlib 7.1.3 is not const.
            return const_cast<basic_ostream_formatter<Char>*>(this)->formatter<basic_string_view<Char>, Char>::format(
                {buffer.data(), buffer.size()}, ctx);
        }
    };
    using ostream_formatter = basic_ostream_formatter<char>;
}
#endif


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
