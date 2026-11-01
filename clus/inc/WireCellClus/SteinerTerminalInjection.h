/** Steiner terminals handed to Steiner::Grapher::create_steiner_tree from a file, in place of its Phase 1.

    pr_ml doc 36 sec 2.2 (E1): learned terminals are judged after the exact router, and the router runs only here, so the test
    needs terminals in C++ while the model stays in Python.  CreateSteinerGraph key terminals_inject_file; absent = the rule.

    The file is plain text, one record per line, '#' starts a comment:

        T <cluster key> <point index> <score>   a Phase-1 terminal, a point index of the retiled cluster
        N <cluster key> <n>                     keep the n best-scored after phases 2, 3 and 3b; absent = keep all

    A cluster with no T record keeps the rule.  n is applied after the filters and the thinning (PDVD sets terminal_min_separation) so that a matched terminal count holds for the set
    handed to the router, as doc 36 sec 2.2 registers.  Phase 4 (the extreme points) is added after that, as for the rule. */
#ifndef WIRECELLCLUS_STEINERTERMINALINJECTION
#define WIRECELLCLUS_STEINERTERMINALINJECTION

#include "WireCellUtil/Exceptions.h"

#include <algorithm>
#include <cstddef>
#include <istream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace WireCell::Clus::Steiner {

    struct InjectedTerminals {
        std::map<size_t, double> score;   // point index -> score; ordered, so iteration is deterministic
        long keep{-1};                    // < 0: keep all
    };
    using TerminalInjection = std::map<int, InjectedTerminals>;   // cluster key -> its terminals

    inline TerminalInjection parse_terminal_injection(std::istream& in)
    {
        TerminalInjection out;
        std::string line;
        size_t nline = 0;
        while (std::getline(in, line)) {
            ++nline;
            const auto hash = line.find('#');
            if (hash != std::string::npos) line.resize(hash);
            std::istringstream ss(line);
            std::string rec;
            if (!(ss >> rec)) continue;   // blank
            int key = 0;
            std::string extra;
            if (rec == "T") {
                long idx = -1;
                double score = 0;
                if (!(ss >> key >> idx >> score) || idx < 0 || (ss >> extra)) {
                    raise<ValueError>("terminal file line %d: bad T record", nline);
                }
                if (!out[key].score.emplace((size_t) idx, score).second) {
                    raise<ValueError>("terminal file line %d: point %d of cluster %d listed twice", nline, idx, key);
                }
            }
            else if (rec == "N") {
                long n = -1;
                if (!(ss >> key >> n) || n < 0 || (ss >> extra)) {
                    raise<ValueError>("terminal file line %d: bad N record", nline);
                }
                out[key].keep = n;
            }
            else {
                raise<ValueError>("terminal file line %d: unknown record '%s'", nline, rec.c_str());
            }
        }
        return out;
    }

    /// The n best-scored of terminals; equal scores by ascending index, so the result does not depend on anything but the file.
    template <typename VertexSet>
    VertexSet keep_best_terminals(const VertexSet& terminals, const InjectedTerminals& inj)
    {
        if (inj.keep < 0 || terminals.size() <= (size_t) inj.keep) return terminals;
        std::vector<std::pair<double, size_t>> ranked;
        ranked.reserve(terminals.size());
        for (const auto idx : terminals) {
            const auto it = inj.score.find((size_t) idx);
            if (it == inj.score.end()) {
                raise<ValueError>("terminal %d has no score", (size_t) idx);
            }
            ranked.emplace_back(-it->second, (size_t) idx);
        }
        std::sort(ranked.begin(), ranked.end());
        VertexSet out;
        for (size_t i = 0; i < (size_t) inj.keep; ++i) out.insert(ranked[i].second);
        return out;
    }

}  // namespace WireCell::Clus::Steiner

#endif
