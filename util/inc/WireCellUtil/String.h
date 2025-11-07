#ifndef WIRECELLUTIL_STRING
#define WIRECELLUTIL_STRING

// fixme: watchme: Boost started to deprecate some internal header
// inclusion which is not, as best as I can tell, any of our problem.
// The message is:
//
// ../../../../../opt/boost-1-76-0/include/boost/config/pragma_message.hpp:24:34: note: ‘#pragma message: This header is deprecated. Use <iterator> instead.’
//
//  This arises from a deeply nested #include well beyond anything
//  which is obvious here.
//
//  If/when this is cleaned up in Boost, remove this comment and the
//  next line.
#define BOOST_ALLOW_DEPRECATED_HEADERS 1
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

#include <vector>
#include <string>
#include <map>
#include <sstream>

namespace WireCell {

    namespace String {

        std::vector<std::string> split(const std::string& in, const std::string& delim = ":");

        std::pair<std::string, std::string> parse_pair(const std::string& in, const std::string& delim = ":");

        // format_flatten converts from "%"-list to variadic function call.
        template <typename TYPE>
        boost::format format_flatten(boost::format f, TYPE obj)
        {
            return f % obj;
        }
        template <typename TYPE, typename... MORE>
        boost::format format_flatten(boost::format start, TYPE o, MORE... objs)
        {
            auto next = start % o;
            return format_flatten(next, objs...);
        }
        inline boost::format format_flatten(boost::format f) { return f; }

        /** The format() function provides an sprintf() like function.

            It's a wrapper on boost::format() but returns a string
            instead of a stream and printf() like semantics.

            int a=42;
            std::string foo = "bar";
            // Assuming positional order
            std::string msg = format("the answer to %s is %d", foo, a);
            // Changing positional order
            std::string msg2 = format("the answer %2% to the question %1%", foo, a);

            WARNING, the "%"-formatting markup syntax used here IS DIFFERENT
            than that what is used in SPDLOG logging functions.  It is a very
            common mistake to confuse the two.  Making this all uniform is a
            "fixme" in WCT that can be done once C++20 and std::fmt is
            available.  See #429.
        */

        template <typename... TYPES>
        std::string format(const std::string& form, TYPES... objs)
        {
            auto final = format_flatten(boost::format(form), objs...);
            std::stringstream ss;
            ss << final;
            return ss.str();
        }

        // Turn some type into a string via streaming.
        template <typename T>
        std::string stringify(T obj)
        {
            std::stringstream ss;
            ss << obj;
            return ss.str();
        }

        inline
        std::string indent(size_t num, const std::string& tab="\t")
        {
            std::stringstream ss;
            for (size_t ind=0; ind<num; ++ind) ss << tab;
            return ss.str();
        }

        bool endswith(const std::string& whole, const std::string& part);
        bool startswith(const std::string& whole, const std::string& part);

    }  // namespace String
}  // namespace WireCell

#endif
