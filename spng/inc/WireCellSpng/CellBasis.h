#pragma once

#include "WireCellSpng/RayGrid.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IAnodeFace.h"

namespace WireCell::SPNG::CellBasis {

    /// Return a cell basis tensor.
    ///
    /// This constructs in U-V-W order and returns shape (Ncells, 3)
    // FIXME: could make the order a parameter
    torch::Tensor cell_basis(IAnodeFace::pointer face);

    
    /// Get wire endpoints as (nwires, 2-endpoints, 2-dimensions) tensor
    /// suitable for RayGrid.  Besides packing to tensor, a coordinate mapping
    /// is applied to the input points to get final dim columns: (X,Y,Z)->(Z,Y).
    torch::Tensor wire_endpoints(const IWire::vector& wires_vec);
    
    /// Return a 1D tensor holding channel IDs
    torch::Tensor channel_idents(const IChannel::vector& chans);

    /// Return a 1D tensor holding the index into chans for each in wires.
    ///
    /// An index value in the output of -1 indicates a wire in wires had no channel in chans.
    torch::Tensor wire_channel_index(IWire::vector wires, const IChannel::vector& chans);

    /// Use a cell basis index tensor to index a set of data tensors.
    ///
    /// @param basis A cell basis tensor holding indices.  Shape (Ncells, Nviews)
    /// @param data A vector of Nviews tensors indexed by the indices.
    /// @return A tensor of shape (Ncells, Nviews) holding values from data.
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

    /// Return canonical channel ordering.  Each vector holds a channel ID
    /// number.  Indices into this canonical ordering are used for the cell
    /// views algorithm.  Input tensors MUST IMPLICITLY FOLLOW this ordering!
    std::vector<IChannel::vector> channels_per_view(IAnodePlane::pointer ianode,
                                                    const std::vector<int>& face_idents);

    /// Make a cell basis tensor holding channel index into canonical channel order.
    ///
    /// This returns a tensor of shape (Ncells, 3) holding indices into the
    /// canonical channel ordering for each view.  See channels_per_view for that ordering.
    torch::Tensor cell_channel_indices(IAnodePlane::pointer ianode,
                                       const std::vector<int>& face_idents = {0,1});


    /// The core CellViews algorithm
    ///
    /// @param uvw_roi The input ROIs as a vector of tensors.  Each tensor is shaped (nbatch, nchan, ntick)
    /// @param indices Output such as from cell_channel_indices(). A cell basis tensor shape (ncells, 3)
    /// @param out_views A vector of view indices to output cell view info.
    /// @param cell_views A vector of cell view information type.
    /// @param chunk_size The size along the tick dimension to process.  Default value of zero will use full tick domain.  Use a finite chunk size if you face memory pressure in large detectors.
    /// @return A vector spanning out_views of tensors shaped (nbatch, nfeat, nchan, ntick) where nfeat spans cell_views.
    ///
    /// The device of the cell_channel_indices is used.  The uvw_roi will be
    /// moved to that device and the result wil be on that device.
    std::vector<torch::Tensor> cell_views(const std::vector<torch::Tensor>& uvw_roi,
                                          const torch::Tensor& indices,
                                          const std::vector<int>& out_views = {0,1},
                                          const std::vector<std::string>& cell_views = {"mp2", "mp3"},
                                          int64_t chunk_size = 0);

}
