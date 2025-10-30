#ifndef WIRECELL_SPNGROIPROCESS
#define WIRECELL_SPNGROIPROCESS

#include "WireCellSpng/ITorchTensorSet.h"
#include "WireCellSpng/ITorchForward.h"
#include "WireCellSpng/ITorchTensorSetFilter.h"
#include "WireCellSpng/DNNROIPreProcess.h"
#include "WireCellSpng/DNNROIPostProcess.h"
#include "WireCellAux/Logger.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IData.h"
#include "WireCellSpng/TorchnvTools.h"
namespace WireCell {
    namespace SPNG{
        struct DNNROIFindingCfg {

            // The APA to focus on
            std::string apa{"AnodePlane"};
            std::string plane{"na"};

            std::string forward{"TorchService"};


            std::string summary_tag{""};

            std::string decon_charge_tag{""};

            // The model downsamples/rebins in time by a number of ticks.
            int tick_per_slice{10};

            // The output trace tag, likely should be set to "dnnspN"
            // with "N" marking the anode number.
            std::string outtag{""};
        };

        class DNNROIProcess: public Aux::Logger,
                            public ITorchTensorSetFilter,
                            public IConfigurable{
        public:
            DNNROIProcess();
            virtual ~DNNROIProcess();

            virtual bool operator()(const input_pointer &in, output_pointer& out);
            virtual void configure(const Configuration& cfg);
            virtual WireCell::Configuration default_configuration() const {
                Configuration cfg;
                return cfg;
            };
            virtual void finalize();

        private:
            DNNROIFindingCfg m_cfg;

            //AB Forward the Pointer
            ITorchForward::pointer m_forward;

            int m_coarse_time_offset{0};
            int m_save_count{0};
            bool m_is_gpu{true};

            //pointers for the pre and post processing
            //IDNNROIPostProcess::pointer m_postprocess;
            //IDNNROIPreProcess::pointer m_preprocess;
            std::shared_ptr<DNNROIPreProcess> m_preprocess;
            std::shared_ptr<DNNROIPostProcess> m_postprocess;

            
            //Additional (and not necessary) helper functions for debugging
            //something to save the tensors for debugging
            //Instead of putting something in config file...just add here
            bool m_is_debug{false};
            // a vector of torch::Tensor that only pushes when m_is_debug is true
            std::vector<torch::Tensor> m_debug_tensors;

            void push_debug_tensor(const torch::Tensor& tensor) {
                if (m_is_debug) {
                    m_debug_tensors.push_back(tensor.clone());
                }
            }
            void push_debug_tensor(const std::vector<torch::Tensor>& tensors) {
                if (m_is_debug) {
                    for (const auto& tensor : tensors) {
                        m_debug_tensors.push_back(tensor.clone());
                    }
                }
            }
            void save_debug_tensors(const std::string& filename) {
                if (m_is_debug && !m_debug_tensors.empty()) {
                    torch::save(m_debug_tensors, filename);
                    m_debug_tensors.clear();
                }
            }

        };
    }
}


#endif