/**
    Wire Cell code may throw exceptions.

    The core Wire Cell code should only throw WireCell::Exception (or
    a subclass).

    Boost exceptions are used to allow the programmer to provide
    dynamic information about an error and to collect file name and
    line number.

    Exceptions should be thrown something like:

    if (omg) {
        THROW(ValueError() << errmsg{"I didn't expect that"});
    }
 */

#ifndef WIRECELL_EXCEPTIONS
#define WIRECELL_EXCEPTIONS

#include <boost/exception/all.hpp>

// Defining this will give more info in the back trace but requires
// linking to -lbacktrace which is usually provided by the compiler.
// This library is fairly prevalent but not everywhere.  The main
// wscript auto-detects if it is available.
#include "WireCellUtil/BuildConfig.h"
#ifdef HAVE_BACKTRACE_LIB
#define BOOST_STACKTRACE_USE_BACKTRACE
#endif

#include <boost/stacktrace.hpp>
#include <exception>
#include <string>

using stack_traced_t = boost::error_info<struct tag_stacktrace, boost::stacktrace::stacktrace>;
// template <class E>
// void throw_with_trace(const E& e) {
//     BOOST_THROW_EXCEPTION(boost::enable_error_info(e) << stack_traced_t(boost::stacktrace::stacktrace()));
// }
// #define THROW(e) throw_with_trace(e)
#define THROW(e) BOOST_THROW_EXCEPTION(boost::enable_error_info(e) << stack_traced_t(boost::stacktrace::stacktrace()))
//#define THROW(e) BOOST_THROW_EXCEPTION(e)
#define errstr(e) boost::diagnostic_information(e)


namespace WireCell {

    // Get the stacktrace as an object.  You must test for non-nullptr.
    // Or, just rely on e.what().
    inline
    const boost::stacktrace::stacktrace* stacktrace(const std::exception& e) {
        return boost::get_error_info<stack_traced_t>(e);
    }


    /// The base wire cell exception.
    struct Exception : virtual public std::exception, virtual boost::exception {
        char const *what() const throw() {
            return diagnostic_information_what(*this);
        }
    };

    /// Thrown when a wrong value has been encountered.
    struct ValueError : virtual public Exception {
    };

    /// Thrown when a wrong index is used.
    struct IndexError : virtual public Exception {
    };

    /// Thrown when a wrong key or has been encountered.
    struct KeyError : virtual public Exception {
    };

    /// Thrown when an error involving accessing input or output has occurred.
    class IOError : virtual public Exception {
    };

    /// Thrown when an error occurs during the data processing
    class RuntimeError : virtual public Exception {
    };

    /// Thrown when an assertion fails
    class AssertionError : virtual public Exception {
    };

    typedef boost::error_info<struct tag_errmsg, std::string> errmsg;

}  // namespace WireCell
#endif
