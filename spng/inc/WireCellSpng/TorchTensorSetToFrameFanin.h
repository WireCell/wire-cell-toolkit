#ifndef WIRECELL_SPNG_TORCHTENSORSETTOFRAMEFANIN
#define WIRECELL_SPNG_TORCHTENSORSETTOFRAMEFANIN

#include "WireCellSpng/ITorchTensorSetToFrameFanin.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Logging.h"
#include "WireCellAux/Logger.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/WirePlaneId.h"

namespace WireCell {
    namespace SPNG {

        /**
           TorchSetToFrameFanin converts a non-TDM tensor set to IFrame.
        */
        class TorchTensorSetToFrameFanin
            : public Aux::Logger,
              public ITorchTensorSetToFrameFanin, public IConfigurable {
           public:
            TorchTensorSetToFrameFanin();
            virtual ~TorchTensorSetToFrameFanin() {};

            // INode, override because we get multiplicity at run time.
            virtual std::vector<std::string> input_types();

            // IFanout
            virtual bool operator()(const input_vector& inv, output_pointer& out);

            // IConfigurable
            virtual void configure(const WireCell::Configuration& cfg);
            virtual WireCell::Configuration default_configuration() const;

           private:
            int m_multiplicity;

            //Wire Planes to pack together into an output TorchTensorSet
            std::map<const WirePlaneId, int> m_input_groups;
            //How many wires in the TorchTensor in the ith output
            std::map<int, int> m_output_nchannels;

            //Channel Map
            std::unordered_map<int, int> m_channel_map;
            std::unordered_map<int, int> m_channel_to_output_group;
            std::vector<std::unordered_map<int, int>> m_per_group_channel_map;


            //Expand dimensionality on output
            bool m_unsqueeze_output{false};

            std::string m_anode_tn{"AnodePlane"};
            IAnodePlane::pointer m_anode;

            size_t m_count{0};

        };
    }  // namespace Aux
}  // namespace WireCell

#endif
