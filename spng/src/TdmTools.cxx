#include "WireCellSpng/TdmTools.h"
#include "WireCellUtil/Fmt.h"
#include <regex>

namespace WireCell::SPNG::TDM {

    Configuration derive_metadata(Configuration md,
                                  Configuration from,
                                  const std::string& datapath_format,
                                  const std::string& tag)
    {
        auto df = from["datapath"];
        md = update(from, md);  // want md to take precedence
        if (! tag.empty()) {
            md["tag"] = tag;
        }
        if (! df.isNull()) {
            md["derived_from"] = df;
        }
        if (datapath_format.size()) {
            md["datapath"] = Fmt::format(datapath_format, md);
        }
        return md;
    }


    Configuration derive_metadata(Configuration md,
                                  const std::vector<Configuration>& froms,
                                  const std::string& datapath_format,
                                  const std::string& tag)
    {
        Configuration dmd;
        Configuration dfs;
        for (const auto& one : froms) {
            auto from = one["datapath"];
            if (! from.isNull()) {
                dfs.append(from);
            }
            update(dmd, one);
        }
        md = update(dmd, md);
        if (! tag.empty()) {
            md["tag"] = tag;
        }
        if (dfs.size()) {
            md["derived_from"] = dfs;
        }
        if (datapath_format.size()) {
            md["datapath"] = Fmt::format(datapath_format, md);
        }
        return md;
    }
    

    bool match_object(const Configuration& match, const Configuration& object)
    {
        Json::StreamWriterBuilder builder;

        // We must complete the loop to return true.
        for (const auto& key : match.getMemberNames()) {
            auto jtarget = object[key];

            if (jtarget.isNull()) { // must have all attributes in match
                // std::cerr << "No match attribute: " << key << "\n";
                return false;
            }

            auto jmatch = match[key];
            if (jtarget == jmatch) { // test for simple equality
                // std::cerr << "Found: " << key << " " << jtarget << "\n";
                continue;
            }

            // one last chance, regex

            if (! jmatch.isString()) { // can't be a regex
                // std::cerr << "No match: " << key << " not a string " << jmatch << "\n";
                return false;
            }

            // coerce target to string
            std::string text;
            if (jtarget.isString()) {
                text = jtarget.asString();
            }
            else if (jtarget.isInt()) {
                text = std::to_string(jtarget.asInt());
            }
            else if (jtarget.isDouble()) {
                text = std::to_string(jtarget.asDouble());
            }
            else {
                // object or array
                text = Json::writeString(builder, jtarget);
            }
            std::regex re(jmatch.asString());

            if (std::regex_match(text, re)) {
                // std::cerr << "Found: " << key << " regex " << jmatch.asString() << " against '" << text << "'\n";
                continue;
            }
            // std::cerr << "No match: " << key << " regex " << jmatch.asString() << " against '" << text << "'\n";
            return false;  // fail one, fail all
        }
        // std::cerr << "Matched: " << object << "\n";
        return true;
    }


    // Note, this is a multi-nested loop!  Probably slow especially given the
    // conversion to json that can happen when the test for simple equality
    // fails.
    ITorchTensor::vector select_tensors(const ITorchTensor::vector& tensors,
                                        const Configuration& match)
    {
        Configuration jmatches = Json::arrayValue;
        if (match.isObject()) {
            jmatches[0] = match;
        }
        else {
            jmatches = match;
        }
        ITorchTensor::vector selected;

        for (const auto& jmatch : jmatches) {
            for (const auto& ten : tensors) {
                if (match_object(jmatch, ten->metadata())) {
                    selected.push_back(ten);
                }
            }
        }
        return selected;
    }

    std::map<std::string, ITorchTensor::vector> by_datatype(const ITorchTensor::vector& tensors)
    {
        std::map<std::string, ITorchTensor::vector> ret;        
        for (const auto& ten : tensors) {
            auto dt = get<std::string>(ten->metadata(), "datatype", "");
            ret[dt].push_back(ten);
        }
        return ret;
    }

}
