/* Note, this file should not be #include'd into WCT components
 * directly but may be included into tests.

   Uniformly slice a frame by charge summation

   This file provides two variants.

   - MaskSlicer :: a function node which produces an IFrameSlice

   - MaskSlices :: a queuedout node which produces an ISlice

   The two have trace-offs.

   r: keeps context, makes DFP graph simpler, forces monolithic downstream/more memory/node

   s: the opposite

   Notes:

   - Any sample must be different than 0.0 to be considered for
     addition to a slice but no other thresholding is done here (do it
     in some other component).

*/

#ifndef WIRECELLIMG_MASKSLICE
#define WIRECELLIMG_MASKSLICE

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IFrameSlicer.h"
#include "WireCellIface/IFrameSlices.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellAux/Logger.h"

#include <string>
#include <unordered_map>

namespace WireCell::Img {

    namespace Data {
        // ISlice class is held temporarily as concrete.
        class Slice;
    }  // namespace Data

    class MaskSliceBase : public Aux::Logger, public IConfigurable {
    public:
        MaskSliceBase();
        virtual ~MaskSliceBase();

        // IConfigurable
        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

    protected:
        typedef std::map<size_t, Data::Slice*> slice_map_t;
        void slice(const IFrame::pointer& in, slice_map_t& sm);
        bool thresholding(const WireCell::ITrace::ChargeSequence& wiener_charge,
                          const WireCell::ITrace::ChargeSequence& gauss_charge, const size_t qind,
                          const double threshold, const int tick_span, const bool verbose = false);

        size_t m_count{0};
    private:
        IAnodePlane::pointer m_anode;
        int m_tick_span{4};
        std::string m_wiener_tag{"wiener"};
        std::string m_summary_tag{"wiener"};
        std::string m_charge_tag{"gauss"};
        std::string m_error_tag{"gauss_error"};
        std::vector<int> m_active_planes;
        std::vector<int> m_dummy_planes;
        std::vector<int> m_masked_planes;
        double m_dummy_charge{0.0}; // 0.1 in previous version
        double m_dummy_error{1e12};
        double m_masked_charge{0.0}; // 0.2 in previous version
        double m_masked_error{1e12};
        /// NOTE: if set to zero then use the default threshold!!!
        std::vector<double> m_nthreshold{3.6, 3.6, 3.6};
        // from uboone u, v, w settings
        std::vector<double> m_default_threshold{5.87819e+02 * 4.0, 8.36644e+02 * 4.0, 5.67974e+02 * 4.0};
        // If both are zero then determine the tbin range from
        // input frame.
        int m_min_tbin{0};
        int m_max_tbin{0}; // not including
        // doc sbnd_xin/113 "quiet-plane 2-view tiling".  A masked plane's channels are normally filled
        // only from the "bad" channel mask, so a track whose signal one plane does not see (a track
        // running along that plane's wires or along the drift, whose prolonged signal the SP does not
        // pass) yields NO blob at all: GridTiling needs every layer non-empty.  With
        // quiet_mask_window >= 0, a masked-plane channel is ALSO marked masked ({masked_charge,
        // masked_error}) in every slice that has no active tick on that channel within +-window
        // slices, so the 2-view branch can tile the two live planes there.  -1 (default) = off,
        // byte-identical legacy behaviour.
        int m_quiet_mask_window{-1};
        // quiet_mask_gap >= 0 restricts the quiet masking to the slices within +-gap slices of an ACTIVE slice
        // of the same channel (a gap in a channel that does carry signal nearby in time: the prolonged-signal
        // stretch of a track that plane sees at its ends), leaving channels that are silent everywhere live and
        // empty, so the 2-view tiling cannot ghost on them.  -1 (default) = no restriction (window rule only).
        int m_quiet_mask_gap{-1};
    };

    class MaskSlicer : public MaskSliceBase, public IFrameSlicer {
    public:
        virtual ~MaskSlicer();

        // IFrameSlicer
        bool operator()(const input_pointer& in, output_pointer& out);
    };

    class MaskSlices : public MaskSliceBase, public IFrameSlices {
    public:
        virtual ~MaskSlices();

        // IFrameSlices
        virtual bool operator()(const input_pointer& depo, output_queue& slices);
    };

}  // namespace WireCell::Img
#endif
