#include "WireCellSpng/TdmFrame.h"

#include <regex>

namespace WireCell::SPNG {

    TdmFrame::TdmFrame(const TensorIndex::tree_type& frame_node)
    {
        for (const auto& node : frame_node.depth()) {
            auto ten = node.value;
            auto md = ten->metadata();
            auto dt = md["datatype"].asString();

            if (! frame) {
                // If "frame_node" is actually the root node, or future TDM
                // allows for a frame to be a parent of another tensor, this
                // will (eventually) catch the first "frame" type found in the
                // DFS descent.
                if (dt != "frame") {
                    continue;
                }
                frame = ten;
                continue;
            }
            
            if (dt == "chmasks") {
                chmasks[md["label"].asString()] = ten;
                continue;
            }
            auto tag = md["tag"].asString();
            if (dt == "traces") {
                traces[tag] = ten;
                continue;
            }
            if (dt == "summaries") {
                summaries[tag] = ten;
                continue;
            }
            if (dt == "chids") {
                chids[tag] = ten;
                continue;
            }
        }
    }

    TdmFrame frame_at_index(TensorIndex& ti, size_t index)
    {
        size_t count = 0;
        /// TDM says frame is an ultimate parent.  If that changes in the future
        /// and a "frame" node is allowed to be a child, change this from
        /// iterating on children() to depth().
        for (const auto* node : ti.tree().children()) {
            auto md = node->value->metadata();
            if (md["datatype"].asString() != "frame") {
                continue;
            }
            if (count == index) {
                return TdmFrame(*node);
            }
            ++count;
        }
        return TdmFrame();
    }

    TdmFrame frame_at_match(TensorIndex& ti, const std::string& frame_path_match)
    {
        std::regex re(frame_path_match);

        for (const auto* node : ti.tree().children()) {
            auto md = node->value->metadata();
            if (md["datatype"].asString() != "frame") {
                continue;
            }
            auto dp = md["datapath"].asString();
            if (std::regex_match(dp, re)) {
                return TdmFrame(*node);
            }
        }
        return TdmFrame();
    }

}
