#ifndef WIRECELL_SPNGDNNROI
#define WIRECELL_SPNGDNNROI

#include "WireCellAux/Logger.h"
#include "WireCellSpng/ITorchToTensorSet.h"
#include "WireCellSpng/ITorchForward.h"
#include "WireCellSpng/ITorchTensorSetFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellSpng/ITorchSpectrum.h"


namespace WireCell {
    namespace SPNG {
        struct DNNROIFindingCfg {

            // The APA to focus on
            std::string apa{"AnodePlane"};

            // The plane index number (in 0,1,2) to determine which
            // channels span the data.
           // int plane{-1};
            std::string plane{"na"};
            // If true, sort the selected channels by their channel ID
            // value.  If false, the ordering given by the channel
            // placement (ordered according to the so called "wire
            // atachement number").
            bool sort_chanids{false};

            // DNN needs consistent scaling with trained model.  This is
            // multiplied to the input charge values.
            double input_scale{1.0 / 8000};

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
            //Should be some variant of TorchService (AB)
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

        //ROI data has target plane and mp2 and mp3 planes
        struct ROIData {
            std::vector<torch::Tensor> tensors; //induction planes
            std::vector<std::string> r_tags; //tags associated with each tensor
        };

        class DNNROI : public Aux::Logger, 
                      public ITorchTensorSetFilter,                      
                      public IConfigurable {
        public:
            DNNROI( );
            virtual ~DNNROI();

            virtual bool operator()(const input_pointer& in, output_pointer& out);
            virtual void configure(const Configuration& config);
            virtual WireCell::Configuration default_configuration() const {
                Configuration cfg;
                return cfg;
            };
            virtual void finalize();
        private:
            DNNROIFindingCfg m_cfg;
            
            std::unordered_set<int>m_chset; // channels to processstd
            std::vector<int>m_chlist; // channels to process in order
            size_t m_nrows{0}, m_ncols{0};
            
            //AB: Forward the pointer
            ITorchForward::pointer m_forward{nullptr};
            
            
            //std::shared_ptr<ITorchSpectrum> base_frer_spectrum, base_wire_filter;
            int m_coarse_time_offset = 0;
            int m_save_count = 0;
            bool m_is_gpu{false};

        };
    }
}
#endif