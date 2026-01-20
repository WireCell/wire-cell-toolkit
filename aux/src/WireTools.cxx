#include "WireCellAux/WireTools.h"

namespace WireCell::Aux::WireTools {

    IChannel::vector wan_ordered_channels(IAnodePlane::pointer anode, const std::vector<int>& wpid_nums)
    {
        IChannel::vector all_chids;

        for (int wpid_num : wpid_nums) {

            WirePlaneId wpid(std::abs(wpid_num));
            auto face = anode->face(wpid.face());
            auto plane = face->planes()[wpid.index()];

            IChannel::vector chids;
            for (const auto& ich : plane->channels()) {
                chids.push_back(ich);
            }
            if (wpid_num < 0) {
                std::reverse(chids.begin(), chids.end());
            }
            all_chids.insert(all_chids.end(), chids.begin(), chids.end());
        }
        return all_chids;
    }
    
    std::vector<size_t> nwires_wpid(IAnodePlane::pointer anode, const std::vector<int>& wpid_nums)
    {
        std::vector<size_t> sizes;
        for (int wpid_num : wpid_nums) {
            WirePlaneId wpid(std::abs(wpid_num));
            auto face = anode->face(wpid.face());
            auto plane = face->planes()[wpid.index()];
            sizes.push_back(plane->wires().size());
        }
        return sizes;
    }

    IChannel::vector wire_channels(IWire::vector wires, const IChannel::vector& chans)
    {
        std::unordered_map<int, size_t> ich2ind;
        const size_t nchans = chans.size();
        for (size_t cind=0; cind<nchans; ++cind) {
            ich2ind[chans[cind]->ident()] = cind;
        }

        const size_t nwires = wires.size();
        IChannel::vector out(nwires, nullptr);

        for (size_t wind=0; wind<nwires; ++wind) {
            auto wire = wires[wind];
            auto it = ich2ind.find(wire->channel());
            if (it != ich2ind.end()) {
                size_t cind = it->second;
                out[wind] = chans[cind];
            }
        }
        return out;
    }
}
