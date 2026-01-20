#pragma once

#include "WireCellIface/IAnodePlane.h"

namespace WireCell::Aux::WireTools {
    /// @brief Enumerate channel IDENT numbers (chids) in channel INDEX (WAN) order.
    ///
    /// This returns channels in wire-attachment number (WAN, aka
    /// IChannel::index()) order for the given WirePlaneId (WPID) numbers.  The
    /// WPID number may be signed and a negative WPID number is given to reverse
    /// the "natural" order as provided by IWirePlane::channels().
    ///
    /// @param anode A pointer to IAnodePlane.
    /// @param wpid_nums A vector of SIGNED WirePlaneId number.
    /// @return Vector of IChannel::pointer
    ///
    /// A WPID number (absolute value) can be produced with WirePlaneId::ident()
    /// or in Jsonnet with wirecell.WirePlaneId().  The "apa" number is ignored
    /// as the anode is given explicitly.
    ///
    /// Note, the intended use is to return channels in a common plane that are
    /// ordered "around/across" two-faced anodes.  However, caller is free to
    /// call with wpid_nums that span different planes.
    IChannel::vector wan_ordered_channels(IAnodePlane::pointer anode, const std::vector<int>& wpid_nums);
    // fixme: this is generic and should go into aux

    /// Return the ordered number of wires in the given wpids.
    ///
    /// This is useful in order to partition range that spans multiple faces or planes.
    std::vector<size_t> nwires_wpid(IAnodePlane::pointer anode, const std::vector<int>& wpid_nums);
    // fixme: this is generic and should go into aux

    /// @brief Return subset of given channels corresponding to and ordered by
    /// the wires.
    ///
    /// @param wires The ordered wires to consider.
    /// @param chans The ordered channels to consider.
    /// @return Vector of channels of same length as the vector of wires.
    ///
    /// If a wire's channel is not provided, a nullptr entry will be produced for that wire.
    ///
    /// See nwires_wpid() for the ingredients to partition your wires into per face blocks.
    IChannel::vector wire_channels(IWire::vector wires, const IChannel::vector& chans);
    // fixme: this is generic and should go into aux

}
