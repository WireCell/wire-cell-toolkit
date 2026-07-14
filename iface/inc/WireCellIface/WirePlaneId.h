#ifndef WIRECELL_WIREPLANEID
#define WIRECELL_WIREPLANEID

// fixme: should move into WirePlaneIdCfg.h or similar. (more below)
#include "WireCellUtil/Configuration.h"
#include "WireCellUtil/Spdlog.h"
#include <ostream>
#include <functional>
#include "WireCellUtil/Spdlog.h"

namespace WireCell {

    /// Enumerate layer IDs.  These are not indices but are masks!  A wpid can
    /// have a "layer" that is multiple layers.  
    enum WirePlaneLayer_t {
        kUnknownLayer = 0,
        kUlayer = 1,
        kVlayer = 2,
        kWlayer = 4,
        kAllLayers=7             // represents anode+face context
    };
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

        /// Layer as index number (0,1 or 2).  -1 is returned when the layer is not well defined.
        int index() const;

        /// per-Anode face index NOT ident!
        int face() const;

        /// APA number
        int apa() const;

        /// Return true if apa, face and layer are all valid numbers.
        bool valid() const;

        /// Return true if the wpid has only legal values for apa, face and
        /// layer.  Layer must be well defined as a single layer (u,v,w) or as
        /// "all" layers  which then represents the anode+face context.
        // operator bool() const;

        bool operator==(const WirePlaneId& rhs) const;

        bool operator!=(const WirePlaneId& rhs) const;

        bool operator<(const WirePlaneId& rhs) const;

        /// Return a new wpid defined with the given layer value but same apa/face.
        WirePlaneId to_layer(WirePlaneLayer_t layer) const;        

        /// Return a new wpid with a well defined plane but same apa/face.        
        WirePlaneId to_u() const;
        WirePlaneId to_v() const;
        WirePlaneId to_w() const;

        /// Return a new wpid brocaded to apply to all planes but same apa/face.
        WirePlaneId to_all() const;

        /// Return a standardized name of the form:
        ///
        ///     "a{apa()}f{face()}p{layer()}"
        ///
        /// Eg, "a2f1pU" is the U-plane of face index 1 on anode ident 2.  Layer
        /// letter may be "A" for all or "?" for unknown.
        std::string name() const;

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
