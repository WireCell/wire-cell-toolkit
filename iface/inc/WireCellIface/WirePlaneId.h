#ifndef WIRECELL_WIREPLANEID
#define WIRECELL_WIREPLANEID

// fixme: should move into WirePlaneIdCfg.h or similar. (more below)
#include "WireCellUtil/Configuration.h"
#include "WireCellUtil/Spdlog.h"
#include <ostream>
#include <functional>

namespace WireCell {

    /// Enumerate layer IDs.  These are not indices!
    enum WirePlaneLayer_t { kUnknownLayer = 0, kUlayer = 1, kVlayer = 2, kWlayer = 4 };
    const WirePlaneLayer_t iplane2layer[3] = {kUlayer, kVlayer, kWlayer};

    class WirePlaneId {
       public:
        explicit WirePlaneId(WirePlaneLayer_t layer, int face = 0, int apa = 0);
        explicit WirePlaneId(int packed);

        /// Unit ID as integer
        int ident() const;

        /// Layer as enum
        WirePlaneLayer_t layer() const;

        /// Layer as integer (not index!)
        int ilayer() const;

        /// Layer as index number (0,1 or 2).  -1 if unknown
        int index() const;

        /// per-Anode face index NOT ident!
        int face() const;

        /// APA number
        int apa() const;

        /// return true if valid
        // operator bool() const;
        bool valid() const;

        bool operator==(const WirePlaneId& rhs);

        bool operator!=(const WirePlaneId& rhs);

        bool operator<(const WirePlaneId& rhs);

       private:
        int m_pack;
    };

    std::ostream& operator<<(std::ostream& os, const WireCell::WirePlaneId& wpid);
    std::ostream& operator<<(std::ostream& o, const WireCell::WirePlaneLayer_t& layer);

    // fixme: should move into WirePlaneIdCfg.h or similar.
    template <>
    inline WireCell::WirePlaneId convert<WireCell::WirePlaneId>(const Configuration& cfg,
                                                                const WireCell::WirePlaneId& def)
    {
        return WireCell::WirePlaneId(iplane2layer[convert<int>(cfg[0])], convert<int>(cfg[1], 0),
                                     convert<int>(cfg[2], 0));
    }

}  // namespace WireCell

// implement hash() so WirePlaneId an be used as a key in unordered STL containers.
namespace std {
    template <>
    struct hash<WireCell::WirePlaneId> {
        std::size_t operator()(const WireCell::WirePlaneId& wpid) const { return std::hash<int>()(wpid.ident()); }
    };
}  // namespace std

template <> struct fmt::formatter<WireCell::WirePlaneId> : fmt::ostream_formatter {};

#endif
