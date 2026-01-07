#pragma once

#include "WireCellSpng/RayGrid.h"
#include "WireCellIface/IAnodeFace.h"

namespace WireCell::SPNG::CellBasis {

    /// Return a cell basis tensor.
    ///
    /// By default this constructs in U-V-W order.
    /// FIXME: could make the order a parameter
    torch::Tensor cell_basis(IAnodeFace::pointer face);

    
    /// Get wire endpoints as (nwires, 2-endpoints, 2-dimensions) tensor
    /// suitable for RayGrid.  Besides packing to tensor, a coordinate mapping
    /// is applied to the input points to get final dim columns: (X,Y,Z)->(Z,Y).
    torch::Tensor wire_endpoints(const IWire::vector& wires_vec);
    
    /// Return a 1D tensor holding channel IDs
    torch::Tensor channel_idents(const IChannel::vector& chans);

    /// Return a 1D tensor holding the index into chans for each in wires.
    torch::Tensor wire_channel_index(IWire::vector wires, const IChannel::vector& chans);

    /// Use a cell basis index tensor to index a set of data tensors.
    ///
    /// @param basis A cell basis tensor holding indices.
    /// @param data A vector of tensors indexed by the indices.
    ///
    /// This works by indexing each vector in indices with the corresponding
    /// column in basis.
    ///
    /// Each data tensor must span the indices of the corresponding column.
    ///
    /// Examples:
    ///
    /// The basis tensor holds the wire indices for each plane in a cell and the
    /// data tensor holds channel indices for each wire.  The function returns a
    /// new basis tensor holding the channel indices for each cell.
    ///
    /// The basis tensor holds channel indices for each cell and the data
    /// tensors hold charge for each channel.  The result is a basis tensor
    /// holding a trio of charge for each cell.
    torch::Tensor index(torch::Tensor basis, std::vector<torch::Tensor>& data);

}
