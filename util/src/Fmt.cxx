#include "WireCellUtil/Fmt.h"

namespace WireCell::Fmt {

    std::string format(const std::string& pattern, const Configuration& params)
    {

        if (params.isObject()) {
            auto args = fmt::dynamic_format_arg_store<fmt::format_context>(); 
            const auto member_names = params.getMemberNames();
            for (const auto& key : member_names) {
                const Json::Value& value = params[key];
                if (value.isString()) {
                    args.push_back(fmt::arg(key.c_str(), value.asString()));
                } else if (value.isIntegral()) {
                    args.push_back(fmt::arg(key.c_str(), value.asInt64()));
                } else if (value.isDouble()) {
                    args.push_back(fmt::arg(key.c_str(), value.asDouble()));
                } else if (value.isBool()) {
                    args.push_back(fmt::arg(key.c_str(), value.asBool()));
                } else {
                    throw std::invalid_argument(
                        "Object values must be scalar (string, int, double, bool).");
                }
            }
            return fmt::vformat(pattern, args);
        }

        if (params.isArray()) {
            auto args = fmt::dynamic_format_arg_store<fmt::format_context>();
            for (const auto& element : params) {
                if (element.isString()) {
                    args.push_back(element.asString());
                } else if (element.isIntegral()) {
                    args.push_back(element.asInt64());
                } else if (element.isDouble()) {
                    args.push_back(element.asDouble());
                } else if (element.isBool()) {
                    args.push_back(element.asBool());
                } else {
                    throw std::invalid_argument(
                        "Array elements must be scalar (string, int, double, bool).");
                }
            }
            return fmt::vformat(pattern, args);
        }

        if (params.isString()) {
            return fmt::format(pattern, params.asString());
        }

        if (params.isIntegral()) {
            return fmt::format(pattern, params.asInt64());
        }

        if (params.isDouble()) {
            return fmt::format(pattern, params.asDouble());
        }

        if (params.isBool()) {
            return fmt::format(pattern, params.asBool());
        }

        if (params.isNull()) {
            // fmt cannot format a nullptr directly with {}, decide on a representation.
            // Here we throw, but you could return fmt::format(pattern, "null");
            throw std::invalid_argument("JSON null is not a supported scalar type for formatting.");
        }

        throw std::invalid_argument("Unsupported Json::Value type for formatting.");
    }
    
}
