/*

A tag rule is an association of a match rule and a single string or an
array of strings.  The match is in general a regular expression
(regex).  Multiple tag rules are collected into a "tag rule set" which
is expressed as an associative array aka dict aka object aka map with
the match as the key.  Tag rule sets may be further collected into
contexts.

For example, the FrameFanin/FrameFanout pair from gen allows Jsonnet
expressions like:

tags_rules : [
 {
     frame: {
         tagXY: ["tagX%d"%n, "tagY%d"%n],
         ["is-tagged%d"%n]: "was-tagged",
         "*": "fanned",
         "*": "saw-pipeline-%d"%n,
     },
     // Output traces are tagged based on a map from an
     // input to one or more output tags.
     trace: {
         intagA: "outtagA%d"%n,
         intagB: ["outtagB%d"%n, "outtagBB%d"%n],
     }
 }
} for n in std.range(0,self.multiplicity-1)],

In C++ this file provides utilities to parse and query this structure.

*/

#include "WireCellUtil/Configuration.h"

#include <unordered_set>
#include <unordered_map>
#include <string>
#include <vector>
#include <regex>

namespace WireCell {

    namespace tagrules {

        typedef std::string tag_t;
        typedef std::string match_t;
        typedef std::unordered_set<tag_t> tagset_t;
        typedef std::pair<std::regex, tagset_t> rule_t;
        typedef std::vector<rule_t> ruleset_t;

        // return rule's tagset if rule matches tag else empty set
        tagset_t match(const tag_t& tag, const rule_t& rs);

        // Compare a rule set to a tag and fill ret with first
        // matching rule value and return true if match found.  Else
        // return false.  If "all_rules" is false, then stop checking
        // rules after the first is found.
        bool match(const tag_t& tag, const ruleset_t& rs, tagset_t& ret, bool all_rules = true);

        // Given a tagset and a ruleset return a tagset which
        // transforms all input tags according to rules.
        tagset_t transform(const tagset_t& ts, const ruleset_t& rs, bool all_rules = true);

        /* A tagrule::Context represents a collection of tag
         * transforms, each at an index and with a name.
         */
        class Context {
            std::unordered_map<std::string, std::vector<ruleset_t> > m_rulesets;

          public:

            // Return the number of rulesets under a given rule name.
            size_t nrules(const std::string& name) const {
                auto it = m_rulesets.find(name);
                if (it == m_rulesets.end()) return 0;
                return it->second.size();
            }

            // This should be an array of objects each keyed by a
            // context name and with values providing a ruleset.
            void configure(const Configuration& cfg);

            // Transform tag in context and return set of produced tags.
            tagset_t transform(size_t ind, const std::string& name, const tag_t& tag);

            // Transform a collection of tags in a context.
            template <typename Tags>
            Tags transform(size_t ind, const std::string& name, const Tags& tags)
            {
                const auto& rsv = m_rulesets[name];
                if (rsv.empty() or ind >= rsv.size()) {
                    return Tags();
                }
                tagrules::tagset_t ts(tags.begin(), tags.end());
                const auto& rs = rsv[ind];
                auto out = tagrules::transform(ts, rs);
                return Tags(out.begin(), out.end());
            }
        };

    }  // namespace tagrules

    // Some Configuration cverter helpers for targules types.

    template <>
    inline tagrules::tagset_t convert<tagrules::tagset_t>(const Configuration& cfg, const tagrules::tagset_t& def)
    {
        tagrules::tagset_t ret;
        if (cfg.isString()) {
            ret.insert(cfg.asString());
            return ret;
        }
        if (cfg.isArray()) {
            for (auto one : cfg) {
                ret.insert(one.asString());
            }
            return ret;
        }
        return def;
    }

    template <>
    inline tagrules::ruleset_t convert<tagrules::ruleset_t>(const Configuration& cfg, const tagrules::ruleset_t& def)
    {
        tagrules::ruleset_t ret;
        for (auto key : cfg.getMemberNames()) {
            auto ts = convert<tagrules::tagset_t>(cfg[key]);
            if (ts.empty()) {
                continue;
            }
            ret.push_back(make_pair(std::regex(key), ts));
        }
        if (ret.empty()) {
            return def;
        }
        return ret;
    }

}  // namespace WireCell
