/// Some wrappers around ftmlib as provided by spdlog.
///

#ifndef WIRECELL_UTIL_FMT
#define WIRECELL_UTIL_FMT

#include "WireCellUtil/Configuration.h"

#include "WireCellUtil/Spdlog.h"



namespace WireCell::Fmt {

    /**
       @brief Formats a string pattern using parameters from a Json::Value.

       This function inspects the type of the 'params' object and formats the
       'pattern' string accordingly using the fmt library.

       If 'params' is a scalar (string, number, bool), it's used as a single
       positional argument.  The 'pattern' should have a single "{}" marker or
       one or more "{0}" markers. a

       If 'params' is an array of scalars, its elements are used as multiple
       positional arguments.  This may use "{}" or "{0}, {1}, ..." markers.

       If 'params' is an object with scalar values, its key-value pairs are used as
       named arguments.  This should use "{name}" attribute name markers.

       @param pattern The fmt-style format string.

       @param params The Json::Value containing the parameters.

       @return The formatted string.

       @throws std::invalid_argument if 'params' or its elements are of an

       unsupported type (e.g., nested objects/arrays).
    */
    std::string format(const std::string& pattern, const Configuration& params);
}

#endif
