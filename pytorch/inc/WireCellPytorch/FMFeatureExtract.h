/** FMFeatureExtract -- the frame -> foundation-model feature stage (wcfm/docs/01 sec 4.2, docs/04).

    An IFrameTensorSet: one instance per anode.  For every configured wire plane it packs the
    tagged traces of the input frame the way the FM training packs were built --

      rows     = the plane's channels in channel-ident order (row = ident - first ident; U [0,800),
                 V [800,1600), W [1600,2560) per DUNE FD-HD APA; the W image is face 1 then face 0
                 in channel order, the training layout),
      columns  = 4-tick slices k from the frame's tick 0 (tick_span = MaskSlices' production span),
      pixel    = input_scale * SUM of the tick_span ticks (+ input_offset),   active <=> pixel > active_threshold,
      value    = the plane's VIEW_NORM log map  2 (log10(x+m) - log10 m) / (log10(M+m) - log10 m) - 1,
      canvas   = [1, 2, h, w] = (value, 0/1 active mask) on the tight bounding box of the active
                 pixels, floored at min_canvas on the high side (dense_mae_adapter._rasterize_one),
                 optionally padded to a multiple of bbox_pad --

    calls the named ITensorForward (TorchService holding the scripted MBV3 student), and gathers the
    [1, C, h, w] reply at the active pixels.  Output: one ITensorSet per frame with, per plane,

      coords    (N, 2) i4   [channel ident, slice k], sorted by (channel, slice)
      feat      (N, C) f4   or, with store_half, feat_half (N, C) u2 = IEEE half bits of the same

    and the sidecar metadata (frame ident/time/tick, tick_span, input tag/scale, per-plane view_norm,
    bbox, canvas, n_active, plus a caller-supplied "provenance" object, e.g. the model file and sha).

    Large canvases (h*w > max_dense_pixels, 0 = never) are tiled along the slice axis with a halo
    of `halo` slices on each side of every tile; only the tile's core columns are gathered so a
    tile edge never enters the stored features.

    The DNN-ROI front end (DNNROIFinding.cxx, Util.cxx) is not touched: this node builds its own
    canvas and reads the forward's reply through the ITensor data pointer directly.
*/

#ifndef WIRECELLPYTORCH_FMFEATUREEXTRACT
#define WIRECELLPYTORCH_FMFEATUREEXTRACT

#include "WireCellIface/IFrameTensorSet.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/ITensorForward.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellAux/Logger.h"
#include "WireCellUtil/Configuration.h"

#include <array>
#include <string>
#include <vector>

namespace WireCell::Pytorch {

    class FMFeatureExtract : public Aux::Logger, public IFrameTensorSet, public IConfigurable {
      public:
        FMFeatureExtract();
        virtual ~FMFeatureExtract();

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        virtual bool operator()(const input_pointer& in, output_pointer& out);

      private:
        struct PlaneRows {
            int plane{0};
            int base{0};              // first channel ident (row 0)
            int nrows{0};
            std::vector<int> chlist;  // ident of row i = base + i
        };

        // The anode whose planes are packed.
        std::string m_anode_tn{"AnodePlane"};
        // The ITensorForward holding the model.
        std::string m_forward_tn{"TorchService"};
        // Tag of the input traces (e.g. "gauss10").
        std::string m_input_tag{""};
        // Wire plane indices to pack, each producing coords+feat tensors.
        std::vector<int> m_planes{0, 1, 2};
        int m_tick0{0};
        int m_nticks{6000};
        int m_tick_span{4};
        double m_input_scale{0.25};
        double m_input_offset{0.0};
        double m_active_threshold{0.0};
        // Per plane index (m, M) of the log map; the training VIEW_NORM.
        std::vector<std::array<double, 2>> m_view_norm{{2.77, 144977.0}, {2.97, 157944.0}, {3.75, 83861.2}};
        int m_min_canvas{64};
        int m_bbox_pad{1};
        int m_feature_dim{128};
        long m_max_dense_pixels{0};
        int m_halo{64};
        bool m_store_half{true};
        Configuration m_provenance;

        IAnodePlane::pointer m_anode;
        ITensorForward::pointer m_forward;
        std::vector<PlaneRows> m_rows;
        int m_count{0};
    };

}  // namespace WireCell::Pytorch

#endif
