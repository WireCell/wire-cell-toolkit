/** PhotonLibraryModel: gridded photon-visibility lookup for QLMatching.
 *
 * Loads a regular-grid visibility library sampled offline from a detector
 * optical model (e.g. the ProtoDUNE-VD PDFastSimANN computable graph, see
 * pdvd/photlib/sample_ann.py + export_wct_photlib.py) and returns per-OpDet
 * visibilities at arbitrary 3D points by trilinear interpolation.
 *
 * The library is described by a JSON meta file:
 *   { "origin_cm": [x,y,z], "step_cm": [dx,dy,dz], "n": [nx,ny,nz],
 *     "nchan": N, "vis_npy": "file.npy" }
 * with vis_npy a float32 C-order numpy array of shape [nx,ny,nz,N] holding
 * the visibility of each grid node for each OpDet channel.  vis_npy is
 * resolved relative to the meta file's directory unless absolute.
 *
 * Points outside the grid are clamped to the grid boundary.  No same-TPC
 * gating is applied: the library itself encodes the detector's optical
 * shadowing (e.g. cathode opacity and double-sided cathode PDs).
 */

#ifndef WIRECELL_MATCH_PHOTONLIBRARYMODEL
#define WIRECELL_MATCH_PHOTONLIBRARYMODEL

#include "WireCellUtil/Point.h"

#include <array>
#include <string>
#include <vector>

namespace WireCell::Match {

    class PhotonLibraryModel {
      public:
        // meta_file: JSON meta path, resolved through Persist (WIRECELL_PATH).
        explicit PhotonLibraryModel(const std::string& meta_file);

        size_t nchan() const { return m_nch; }

        // Fill vis (resized to nchan) with per-channel visibilities at the
        // point (cm).  Channels with mask[i]==0 are skipped (left 0).
        void visibilities(std::vector<double>& vis, const WireCell::Point& p_cm,
                          const std::vector<unsigned int>* mask = nullptr) const;

      private:
        std::array<double, 3> m_origin{};  // cm
        std::array<double, 3> m_step{};    // cm
        std::array<int, 3> m_n{};
        size_t m_nch{0};
        std::vector<float> m_vis;          // [nx*ny*nz*nch], channel fastest
    };

}  // namespace WireCell::Match

#endif
