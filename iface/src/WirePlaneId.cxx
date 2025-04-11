#include "WireCellIface/WirePlaneId.h"

using namespace WireCell;

static const int layer_mask = 0x7;
static const int face_shift = 3;
static const int apa_shift = 4;

WireCell::WirePlaneId::WirePlaneId(WirePlaneLayer_t layer, int face, int apa)
  : m_pack((((int) layer) & layer_mask) | (face << face_shift) | (apa << apa_shift))
{
}
WireCell::WirePlaneId::WirePlaneId(int packed)
  : m_pack(packed)
{
}

int WireCell::WirePlaneId::ident() const { return m_pack; }
WireCell::WirePlaneLayer_t WireCell::WirePlaneId::layer() const { return WirePlaneLayer_t(ilayer()); }
int WireCell::WirePlaneId::ilayer() const { return m_pack & layer_mask; }

int WireCell::WirePlaneId::index() const
{
    switch (layer()) {
    case kUlayer:
        return 0;
    case kVlayer:
        return 1;
    case kWlayer:
        return 2;
    default:
        return -1;
    }
    return -1;
}
int WireCell::WirePlaneId::face() const { return (m_pack & (1 << face_shift)) >> 3; }
int WireCell::WirePlaneId::apa() const { return m_pack >> apa_shift; }

bool WireCell::WirePlaneId::valid() const
{
    if (apa() < 0) return false;
    if (face() < 0) return false;
    if (layer() == kAllLayers) return true;
    const int ind = index();
    return 0 <= ind && ind < 3;
}

WireCell::WirePlaneId::operator bool() const
{
    if (apa() < 0) return false;
    if (face() < 0) return false;
    if (layer() == kAllLayers) return true;
    const int ind = index();
    return 0 <= ind && ind < 3;
}

WirePlaneId WirePlaneId::to_layer(WirePlaneLayer_t layer) const
{
    return WirePlaneId(layer, face(), apa());
}
WirePlaneId WirePlaneId::to_u() const
{
    return to_layer(kUlayer);
}
WirePlaneId WirePlaneId::to_v() const
{
    return to_layer(kVlayer);
}
WirePlaneId WirePlaneId::to_w() const
{
    return to_layer(kWlayer);
}
WirePlaneId WirePlaneId::to_all() const
{
    return to_layer(kAllLayers);
}
std::string WirePlaneId::name() const
{
    std::stringstream ss;
    ss << "a" << apa() << "f" << face() << "p" << layer();
    return ss.str();
}


bool WireCell::WirePlaneId::operator==(const WirePlaneId& rhs) const { return m_pack == rhs.m_pack; }

bool WireCell::WirePlaneId::operator!=(const WirePlaneId& rhs) const { return !(*this == rhs); }

bool WireCell::WirePlaneId::operator<(const WirePlaneId& rhs) const
{
    return m_pack < rhs.m_pack;
    // if (!this->valid() || !rhs.valid()) {
    //     return false;
    // }

    // if (apa() == rhs.apa()) {
    //     if (face() == rhs.face()) {
    //         return index() < rhs.index();
    //     }
    //     return face() < rhs.face();
    // }
    // return apa() < rhs.apa();
}

std::ostream& WireCell::operator<<(std::ostream& o, const WireCell::WirePlaneId& wpid)
{
    o << "[WirePlaneId \"" << wpid.name() << "\" ident=" << wpid.ident()
      << " ind:" << wpid.index() << " layer:" << wpid.layer()
      << " apa:" << wpid.apa() << " face:" << wpid.face();
    if (wpid.valid()) {
        o << " valid";
    }
    else {
        o << " bogus";
    }
    if (wpid) {
        o << " true";
    }
    else {
        o << " false";
    }
    o << "]";
    return o;
}

std::ostream& WireCell::operator<<(std::ostream& o, const WireCell::WirePlaneLayer_t& layer)
{
    switch (layer) {
    case WireCell::kUlayer:
        o << "U";
        break;
    case WireCell::kVlayer:
        o << "V";
        break;
    case WireCell::kWlayer:
        o << "W";
        break;
    case WireCell::kAllLayers:
        o << "A";
        break;
    default:
        o << "?";
        break;
    }
    return o;
}
