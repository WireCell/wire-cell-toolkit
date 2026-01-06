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
    
}
