#ifndef WIRECELL_CLUS_KINESEGMENTLINKS_H
#define WIRECELL_CLUS_KINESEGMENTLINKS_H

/// sbnd_xin/docs/128-129: the segment -> T_kine row map.
///
/// fill_kine_tree pushes one kine row per particle: a track row is one
/// segment, a shower row is every member segment of the shower.  This helper
/// inverts that into one link per segment, over plain ints (graph indices), so
/// the rule is unit-testable without a PR graph.
///
/// Rule (owner ruling 3, doc 128 sec 7): a segment that fed several rows
/// links to the lowest-index SHOWER row if it fed any shower row, else to its
/// lowest row, and reports how many rows it fed.  With
/// kine_mainvtx_used_guard off a shower member at the main vertex feeds its
/// own track row AND its shower's row (the pr/101 K5 double count); two
/// showers sharing a start segment can also share members.

#include <cstddef>
#include <map>
#include <set>
#include <vector>

namespace WireCell::Clus::PR {

    struct KineRowMembers {
        bool is_shower{false};
        std::vector<int> graph_indices;  // the segment(s) whose energy the row carries
    };

    struct KineSegmentLink {
        int graph_index{-1};
        int kine_index{-1};
        int n_rows{0};
    };

    /// One link per distinct graph index, sorted by graph index.
    inline std::vector<KineSegmentLink> build_kine_segment_links(const std::vector<KineRowMembers>& rows)
    {
        struct Acc {
            int first_row{-1};
            int first_shower_row{-1};
            int n_rows{0};
        };
        std::map<int, Acc> acc;
        for (std::size_t r = 0; r < rows.size(); ++r) {
            const int row = static_cast<int>(r);
            // A row listing a segment twice still feeds it once.
            std::set<int> seen;
            for (int gi : rows[r].graph_indices) {
                if (!seen.insert(gi).second) continue;
                Acc& a = acc[gi];
                ++a.n_rows;
                if (a.first_row < 0) a.first_row = row;
                if (rows[r].is_shower && a.first_shower_row < 0) a.first_shower_row = row;
            }
        }
        std::vector<KineSegmentLink> out;
        out.reserve(acc.size());
        for (const auto& [gi, a] : acc) {
            out.push_back({gi, a.first_shower_row >= 0 ? a.first_shower_row : a.first_row, a.n_rows});
        }
        return out;
    }

}  // namespace WireCell::Clus::PR

#endif
