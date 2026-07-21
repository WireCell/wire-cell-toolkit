#include "WireCellSpng/DeconTools.h"
#include "WireCellSpng/Torch.h"
#include "WireCellUtil/Exceptions.h"

namespace WireCell::SPNG {

    torch::Tensor fr_tensor(const Response::Schema::FieldResponse& fr_data,
                            int plane_id,
                            bool zero_centered)
    {
        for (auto & plane : fr_data.planes) {
            if (plane.planeid != plane_id) continue;
        
            const int64_t npaths = plane.paths.size();
            if (npaths == 0) {
                raise<ValueError>("FR for plane %d lacks any paths", plane_id);
            }

            const int64_t nsteps = plane.paths[0].current.size();
            if (nsteps == 0) {
                raise<ValueError>("FR for plane %d got empty path", plane_id);
            }

            std::vector<int64_t> shape = {npaths, nsteps};

            auto response = torch::zeros(shape, torch::kFloat32);
            auto accessor = response.accessor<float,2>();
        
            for (int64_t irow = 0; irow < npaths; ++irow) {
                auto& path = plane.paths[irow];
                for (int64_t icol = 0; icol < nsteps; ++icol) {
                    accessor[irow][icol] = path.current[icol];
                }
            }

            if (! zero_centered) {
                return response;
            }

            /// Move "wire of interest" row from center row to row zero in order
            /// to remove the artificial shift.
            return torch::roll(response, {(npaths+1)/2}, {0});
        }
        raise<ValueError>("no plane %d in field response", plane_id);
        return torch::zeros({0}); // avoid compiler warning.
    }

    torch::Tensor fr_average_tensor(const Response::Schema::FieldResponse& fr_data,
                                    int plane_id, bool zero_centered)
    {
        auto fr_avg = Response::wire_region_average(fr_data);
        return fr_tensor(fr_avg, plane_id, zero_centered);
    }
}
