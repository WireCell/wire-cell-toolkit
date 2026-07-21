#include "WireCellIface/IChannel.h"

using namespace WireCell;

IChannel::~IChannel() {}

WirePlaneId IChannel::planeid() const
{
    static const WirePlaneId bogus(kUnknownLayer, -1, -1);
    const IWire::vector& w = wires();
    if (w.empty()) {
        return bogus;
    }
    return w.front()->planeid();
}

size_t IChannel::global_order() const
{
    const size_t wpid = planeid().ident();
    //auto wire0 = wires()[0];

    const size_t wan = index();
    return ((size_t)wpid << 32) & (0xFFFFFFFF & wan);
    
}
