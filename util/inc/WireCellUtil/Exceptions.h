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

#include "WireCellUtil/String.h"

#include <boost/stacktrace.hpp>
#include <exception>

/// We can use THROW() to add features to the exceptions.
///
/// THROW(ValueError() << errmsg{"some useful description"});
///
/// See raise<T>(...) for a less verbose way.
using stack_traced_t = boost::error_info<struct tag_stacktrace, boost::stacktrace::stacktrace>;
#define THROW(e) BOOST_THROW_EXCEPTION(boost::enable_error_info(e) << stack_traced_t(boost::stacktrace::stacktrace()))
#define errstr(e) boost::diagnostic_information(e)


namespace WireCell {

    using errmsg = boost::error_info<struct tag_errmsg, std::string>;


    /// Simpler way to throw with messages
    ///
    /// raise<ValueError>("a formatted string with %s codes", "printf")
    ///
    template<typename Err, typename... TYPES>
    void raise(const std::string& form, TYPES... objs) {
        THROW(Err() << errmsg{WireCell::String::format<TYPES...>(form, objs...)});
    }
    template<typename Err>
    void raise() {
        THROW(Err());
    }

    // Get the stacktrace as an object.  You must test for non-nullptr.
    // Or, just rely on e.what().
    inline
    const boost::stacktrace::stacktrace* stacktrace(const std::exception& e) {
        return boost::get_error_info<stack_traced_t>(e);
    }

    /// The base wire cell exception.  Do not throw directly.
    struct Exception : virtual public std::exception, virtual boost::exception {
        char const *what() const throw() {
            return diagnostic_information_what(*this);
        }
    };

    /// Thrown when a wrong value is encountered (function argument, calculated value).
    struct ValueError : virtual public Exception {
    };

    /// Thrown when a wrong type is encountered (function argument, calculated value).
    struct TypeError : virtual public Exception {
    };

    /// Thrown when a wrong value is to be used as an index.
    struct IndexError : virtual public ValueError {
    };

    /// Thrown when a wrong value is to be used as a key.
    struct KeyError : virtual public ValueError {
    };

    /// Thrown when encoutered error with external resource (filesystem, network).
    class IOError : virtual public Exception {
    };

    /// Thrown when an error is found that is due to events beyond the
    /// scope of the program and can not easily be predicted.
    class RuntimeError : virtual public Exception {
    };

    /// An error is encountered which must be due to faulty logical or
    /// is logically inconsistent with assumptions or violates logical
    /// preconditions.
    class LogicError : virtual public Exception {
    };

    /// Avoid this, but do use it instead of an assert().
    class AssertionError : virtual public Exception {
    };

}  // namespace WireCell
#endif
