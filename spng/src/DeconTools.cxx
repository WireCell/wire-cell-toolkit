#include "WireCellSpng/DeconTools.h"
#include "WireCellSpng/Torch.h"
#include "WireCellUtil/Exceptions.h"

namespace WireCell::SPNG {

    torch::Tensor fr_average_tensor(IFieldResponse::pointer ifr, int plane_id,
                                    bool zero_centered)
    {
        // Originally extracted from Decon.cxx

        // FIXME: this could be made more generic by moving these first to lines
        // out and making this function take a Response::Schema::FieldResponse.
        auto fr_fine = ifr->field_response();
        auto fr_avg = Response::wire_region_average(fr_fine);

        for (auto & plane : fr_avg.planes) {
            if (plane.planeid != plane_id) continue;
        
            const int64_t nchans = plane.paths.size();
            if (nchans == 0) {
                raise<ValueError>("FR for plane %d lacks any channels", plane_id);
            }

            const int64_t nticks = plane.paths[0].current.size();
            if (nticks == 0) {
                raise<ValueError>("FR for plane %d got empty path", plane_id);
            }

            std::vector<int64_t> shape = {nchans, nticks};

            auto response = torch::zeros(shape, torch::kFloat32);
            auto accessor = response.accessor<float,2>();
        
            for (int64_t irow = 0; irow < nchans; ++irow) {
                auto& path = plane.paths[irow];
                for (int64_t icol = 0; icol < nticks; ++icol) {
                    accessor[irow][icol] = path.current[icol];
                }
            }

            if (! zero_centered) {
                return response;
            }

            /// Move "wire of interest" row from center row to row zero in order
            /// to remove the artificial shift.
            return torch::roll(response, {nchans/2}, {0});
        }
        raise<ValueError>("no plane %d in field response", plane_id);
        return torch::zeros({0}); // avoid compiler warning.
    }
}
