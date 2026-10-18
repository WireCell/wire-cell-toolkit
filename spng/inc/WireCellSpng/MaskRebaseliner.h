#pragma once

#include "WireCellSpng/FanBase.h"
#include "WireCellSpng/Rebaseline.h"
#include "WireCellSpng/Torch.h"

#include "WireCellUtil/HanaJsonCPP.h"
#include "WireCellIface/IConfigurable.h"

namespace WireCell::SPNG {

    /** Configuration for SPNGMaskRebaseliner.  The algorithm options are
     *  documented on MaskRebaselineOptions in Rebaseline.h.
     */
    struct MaskRebaselinerConfig {

        /// The dimension along which ROIs extend.
        int dim = -1;

        /// Samples outside the tight ROI at which the baseline is anchored.
        /// 0 reproduces SPNGRebaseliner (anchors on the ROI's own edge samples).
        int anchor_pad = 4;

        /// Samples averaged per anchor, extending away from the ROI.
        int anchor_width = 1;

        /// ROIs with (end - start) <= min_roi_size - 1 are not rebaselined.
        int min_roi_size = 1;

        /// Keep (rebaselined) samples in the pad region, i.e. output the
        /// dilated ROI instead of the tight one.
        bool keep_pad = false;

        /// Clamp the output to be non-negative.
        bool remove_negative = false;

        /// [begin,end) ranges of the channel (dim -2) index to rebaseline.
        /// Empty means all.  Channels not covered are output as waveform*mask.
        std::vector<std::vector<int>> rows = {};

        /// If false, no rebaselining is done: output = waveform*mask.
        bool enabled = true;

        /// Set the "tag" metadata attribute on the produced tensor.
        std::string tag = "signal";

        /// Format applied to the metadata to produce the output datapath.
        /// Empty keeps the waveform's datapath.
        std::string datapath_format = "";
    };
}

BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::MaskRebaselinerConfig,
                        dim, anchor_pad, anchor_width, min_roi_size,
                        keep_pad, remove_negative, rows, enabled,
                        tag, datapath_format);

namespace WireCell::SPNG {

    /** Apply an ROI mask to a waveform tensor and subtract a per-ROI linear
     *  baseline anchored outside the ROI.
     *
     *  Two inputs: port 0 the waveform (e.g. the gauss-filtered deconvolution),
     *  port 1 the ROI mask (bool, or nonzero means in-ROI), same shape.  One
     *  output tensor of the waveform's shape and dtype.  The output metadata
     *  derives from the waveform's.
     *
     *  This replaces the pair SPNGReduce(mul) + SPNGRebaseliner.  Unlike that
     *  pair it sees the unmasked waveform, so the baseline anchors can sit on
     *  samples outside the ROI, where the signal has returned to baseline,
     *  instead of on the ROI's edge samples which are still on the pulse.
     *  See rebaseline_masked() in Rebaseline.h.
     */
    struct MaskRebaseliner : FaninBase<ITorchTensor>, virtual public IConfigurable {

        MaskRebaseliner();
        virtual ~MaskRebaseliner() = default;

        /// FaninBase API
        virtual void fanin_combine(const input_vector& inv, output_pointer& out);

        /// IConfigurable
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;

    private:
        MaskRebaselinerConfig m_config;
        MaskRebaselineOptions m_opts;
        std::vector<std::pair<int64_t,int64_t>> m_rows;
    };
}
