/**
   Include this file instead of directly including spdlog or fmtlib (fmt)
   headers.

   It takes care of changing fmtlib API in the face of the fact that spdlog may
   be compiled against an external fmtlib or use its "bundled" copy.

   More info: https://github.com/WireCell/wire-cell-spack/issues/12

 */

#ifndef WIRECELLUTIL_SPDLOG
#define WIRECELLUTIL_SPDLOG

// Prefer SPDLOG_LOGGER_DEBUG() or SPDLOG_DEBUG() over bare calls to
// log->debug() or spdlog::debug().
//
// Always use SPDLOG_LOGGER_TRACE() or SPDLOG_TRACE() for trace level logs.
//
//
// To set default minmum level for these macros use, eg:
//
//   ./wcb configure --with-spdlog-active-level=trace [...]
//
// See util/docs/logging.org for more info.

// #define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG
#include "WireCellUtil/BuildConfig.h"

#ifdef __clang__
#  if defined(__has_warning)
#    define HAS_WARNING(warning) __has_warning(warning)
#  else
#    define HAS_WARNING(warning) 1
#  endif
#else
#  define HAS_WARNING(warning) 1
#endif

#if HAS_WARNING("-Wstringop-overflow")
#pragma GCC diagnostic push
#pragma GCC diagnostic warning "-Wstringop-overflow"
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif

#include <spdlog/spdlog.h>

// We need FMT version but it is found in core.h in different locations
// depending on if fmtlib is external to or bundled with spdlog.
#if defined(SPDLOG_FMT_EXTERNAL)
#include <fmt/core.h>
#else
#include <spdlog/fmt/bundled/core.h>
#endif

// Newer fmtlib no longer implicitly uses the ostream insertion operator<<.  So
// we must have some ugliness to allow us to provide new, required specialization:
//
// template <> struct fmt::formatter<TYPE> : fmt::ostream_formatter {};
//
// while backporting that base type to old fmtlib.
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


#if HAS_WARNING("-Wstringop-overflow")
#pragma GCC diagnostic pop
#endif

#endif
