/** Wrapper for Troch Script Model
 *
 */

#ifndef WIRECELLPYTORCH_TSMODEL
#define WIRECELLPYTORCH_TSMODEL

#include "WireCellIface/ITensorForward.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IFrameFilter.h"
#include "WireCellUtil/Array.h"
#include "WireCellAux/Logger.h"

#include <unordered_set>

namespace WireCell {
    namespace Pytorch {

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
            double mask_thresh{0.7};

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
        };

        class DNNROIFinding : public Aux::Logger,
                              public IFrameFilter, public IConfigurable {
          public:
            DNNROIFinding();
            virtual ~DNNROIFinding();

            /// working operation - interface from IFrameFilter
            /// executed when called by pgrapher
            virtual bool operator()(const IFrame::pointer &inframe, IFrame::pointer &outframe);

            /// interfaces from IConfigurable

            /// exeexecuted once at node creation
            virtual WireCell::Configuration default_configuration() const;

            /// executed once after node creation
            virtual void configure(const WireCell::Configuration &config);

          private:
            DNNROIFindingCfg m_cfg;

            // A set of channel IDs that we consider
            std::unordered_set<int> m_chset;
            // Ordered channel IDs defining rows of output dense array.
            std::vector<int> m_chlist;

            // size for dense trace array based on channel ID span (rows) and tick0/nticks (cols)
            size_t m_nrows{0}, m_ncols{0};
            IFrame::trace_list_t m_trace_indices;

            // The heavy lifting is done with a "forward" service
            ITensorForward::pointer m_forward{nullptr}; 

            /// Return copy of traces vector filled with those having
            /// channel IDs in our set.
            ITrace::vector select(ITrace::vector traces);
            
            // Convert traces to a dense array
            Array::array_xxf traces_to_eigen(ITrace::vector traces);

            // extract trace summary to mathch the eigen format order
            IFrame::trace_summary_t get_summary_e(const IFrame::pointer& inframe, const std::string &tag) const;

            // Convert dense array to (dense) traces
            ITrace::shared_vector eigen_to_traces(const Array::array_xxf& arr);

            int m_save_count;  // count frames saved
        };
    }  // namespace Pytorch
}  // namespace WireCell

#endif
