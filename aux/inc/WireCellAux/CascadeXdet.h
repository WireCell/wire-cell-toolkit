/** The G + D input columns of the pooled cross-detector cascade model, read from a file (wcfm doc 53 W7 / P4).

    The doc 53 pooled model (wcp-porting-img wcfm/scripts/xdet/lean_x.py CascadeGNNLeanX) reads, beside the
    level graph of CascadeGraph, 8 blob columns, 4 wire columns and a live-wire mask that the toolkit does not compute
    yet (that is P5, behind a later `ml_xdet_columns` knob).  This test path reads them from an npz written in python by
    wcfm/scripts/d53_xdet_columns_file.py with xdet/columns.gd_columns, so that the exported model can be run in the C++
    cascade and compared with python (P4 parity).

    Every column is a function of the node alone (its shape, its own bw edges, its channel) and of the event's
    dead-channel mask, so the file is keyed, not ordered: whatever survives a level under any model is looked up.
    Per cascade level i (0, 1, ... in the order of CascadeDeghosting's `levels`), the npz holds

      L<i>_bkey   int32 [Nb, 8]   face, slice ident, u0, u1, v0, v1, w0, w1 (the half-open wire bounds of the
                                  blob's U, V, W strips, layers 2-4 of the ray grid; 0 for a missing layer)
      L<i>_bcols  float [Nb, 8]   w_mm[3], n_shared, n_dead_planes, dead_frac[3]
      L<i>_wkey   int32 [Mw, 2]   plane, first channel of the (super-)wire node
      L<i>_wcols  float [Mw, 4]   pitch_mm, k_bin, shared, dead
      L<i>_live   float [Mw]      1 = the channel measured, 0 = dead

    A key given twice with different values is refused when the file is read; a node whose key is missing is refused
    when a level is assembled.
 */
#ifndef WIRECELLAUX_CASCADEXDET
#define WIRECELLAUX_CASCADEXDET

#include "WireCellAux/CascadeGraph.h"
#include "WireCellIface/IBlob.h"

#include <array>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace WireCell::Aux::Cascade {

    constexpr size_t xdet_blob_cols = 8;
    constexpr size_t xdet_wire_cols = 4;

    using XdetBlobKey = std::array<int32_t, 8>;   // face, slice ident, u0, u1, v0, v1, w0, w1
    using XdetWireKey = std::array<int32_t, 2>;   // plane, first channel of the (super-)wire

    struct XdetLevelTable {
        std::map<XdetBlobKey, std::array<float, xdet_blob_cols>> blob;
        std::map<XdetWireKey, std::array<float, xdet_wire_cols + 1>> wire;   // the 4 columns, then live
    };

    struct XdetColumns {
        std::string path;
        std::vector<XdetLevelTable> levels;   // index = cascade level
    };

    /// Read the npz (layout above).  Levels L0, L1, ... are read until the first absent L<i>_bkey.  Throws
    /// ValueError on a malformed array or on a key given twice with different values.
    XdetColumns load_xdet_columns(const std::string& path);

    /// The blob's key: face, slice ident, the U / V / W strip bounds (the dump's `wip` row).
    XdetBlobKey xdet_blob_key(const IBlob::pointer& blob);

    /// The columns of one level graph in its node order: xb_cols [N x 8], xw_cols [M x 4], live [M].
    struct XdetLevelCols {
        std::vector<float> xb_cols, xw_cols, live;
    };

    /// Look up every node of `lev` (built from `blobs`, node order) in level `level` of the file.  Throws ValueError
    /// when the file has no such level or when a node's key is missing (the message gives the counts and the first
    /// missing key).
    XdetLevelCols xdet_level_columns(const XdetColumns& xc, size_t level, const std::vector<IBlob::pointer>& blobs,
                                     const Level& lev);

}  // namespace WireCell::Aux::Cascade

#endif
