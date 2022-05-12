/** This is a temporary compatibility between JsonCPP and
 * nlohmann::json to help convertering from the former to the latter.
 */

#ifndef NLJS2JCPP_HPP
#define NLJS2JCPP_HPP

#include "WireCellUtil/Persist.h"
#include <nlohmann/json.hpp>
#include <json/json.h>


using nljs_t = nlohmann::json;

// nljs converters need to be in in the namespace of the user data
// which here is JsonCPP so we use its ns:
namespace Json {


    // To nlohmann::json
    inline void to_json(nljs_t& nljs, const Value& jval) {
        std::string str = WireCell::Persist::dumps(jval);
        nljs = nljs_t::parse(str);
    }

    inline void from_json(const nljs_t& nljs, Value& jval) {
        std::string str = nljs.dump();
        jval = WireCell::Persist::loads(str);
    }
}

#endif
