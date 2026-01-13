#ifndef WIRECELL_SPNG_FRAMETOTORCHSET
#define WIRECELL_SPNG_FRAMETOTORCHSET

#include "WireCellSpng/IFrameToTorchSet.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IFrameFilter.h"
#include "WireCellUtil/Array.h"
#include "WireCellAux/Logger.h"

#include <unordered_set>
namespace WireCell {
    namespace SPNG {

        /**
           FrameToTorchSet converts an input IFrame into a TorchTensorSet.

           Configuration:

           - (none yet)

        */
        class FrameToTorchSet
            : public Aux::Logger,
              public IFrameToTorchSet, public IConfigurable {
           public:
            FrameToTorchSet();
            virtual ~FrameToTorchSet() {};

            // INode, override because we get multiplicity at run time.
            // virtual std::vector<std::string> output_types();

            // IFrameToTorchSet
            virtual bool operator()(const input_pointer& in, output_pointer& out);

            // IConfigurable
            virtual void configure(const WireCell::Configuration& cfg);
            virtual WireCell::Configuration default_configuration() const;
        struct DNNROIFindingCfg {

            // The anode to focus on
            std::string anode{"AnodePlane"};

            // The plane index number (in 0,1,2) to determine which
            // channels span the data.
            int plane{0};

            // If true, sort the selected channels by their channel ID
            // value.  If false, the ordering given by the channel
            // placement (ordered according to the so called "wire
            // atachement number").
            bool sort_chanids{false};

            // DNN needs consistent scaling with trained model.  This is
            // multiplied to the input charge values.
            double input_scale{1.0 / 4000};

            // Charge offset added to input charge values.
            double input_offset{0.0};

            // The new, typically tighter ROIs are followed by a local
            // rebaselining which can result in a bias of the resulting
            // charge.  This output scale is a rough correction for this bias.
            // It will be multiplied to output charge.  Caution: this is an
            // ad-hoc fix and setting it to something other than 1.0 is really
            // a sign that the model is somehow imperfect.
            double output_scale{1.0};

            // A charge offset added to output charge values.
            double output_offset{0.0};

            // The first tick and number of ticks to consider from the
            // input traces.  See also chids.
            int tick0{0};
            int nticks{6000};

            // Set the threshold on the ROI "likelihood" which is a
            // probability-like value.  It is tuned to balance efficiency and
            // noise reduction and strictly is best optimized for the given
            // model.
            double mask_thresh{0.5};

            // The IForward service to use
            std::string forward{"TorchService"};

            // Tags of sets of traces to use as input.  These are
            // usually "loose_lfN", "mp2_roiN" and "mp3_roiN" with "N"
            // replaced by the anode number.
            std::vector<std::string> intags;

            std::string summary_tag{""};

            // The tag used for the input decon charge.  This is
            // usually "decon_chargeN" with "N" replaced with the
            // anode number.
            std::string decon_charge_tag{""};

            // The model downsamples/rebins in time by a number of ticks.
            int tick_per_slice{10};

            // An output file used for special debugging.
            std::string debugfile{""};


            // The output trace tag, likely should be set to "dnnspN"
            // with "N" marking the anode number.
            std::string outtag{""};

            int nchunks{1};

            // if true, save the negative parts of the charge traces
            bool save_negative_charge{false};
        };
           private:
            DNNROIFindingCfg m_cfg;
            std::unordered_set<int> m_chset;
            std::vector<int> m_chlist;
            size_t m_nrows{0}, m_ncols{0};
            IFrame::trace_list_t m_trace_indices;
            torch::Tensor traces_to_tensor(ITrace::vector traces);
            ITrace::vector select(ITrace::vector traces);

        };
    }  // namespace SPNG
}  // namespace WireCell

#endif