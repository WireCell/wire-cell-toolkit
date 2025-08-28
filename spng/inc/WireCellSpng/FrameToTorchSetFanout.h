/**
   FrameToTorchSetFanout converts an input IFrame into sets of tensors.

   The frame's set of traces is partitioned into multiple subsets based on a
   selection on channel IDs.  Each subset of traces is transformed into a tensor
   and sent out the corresponding output port.

   Configuration:

   - anode :: the type/name of the anode defining the scope of traces to consider.

   - expected_nticks :: number of ticks that the dense output spans (FIXME:
     should be renamed "nticks" to match global conventions).

   - output_groups :: array of group definitions, the size of which sets the
     output multiplicity.  Each group is an array of WirePlaneId values.  The
     index of the group in the array of groups determines the output port to
     which matching data will be placed.


 */

#ifndef WIRECELL_SPNG_FRAMETOTORCHSETFANOUT
#define WIRECELL_SPNG_FRAMETOTORCHSETFANOUT

#include "WireCellSpng/IFrameToTorchSetFanout.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Logging.h"
#include "WireCellAux/Logger.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/WirePlaneId.h"

namespace WireCell {
    namespace SPNG {

        // Fan out 1 frame to N set at construction or configuration time.
        class FrameToTorchSetFanout
            : public Aux::Logger,
              public IFrameToTorchSetFanout, public IConfigurable {
           public:
            FrameToTorchSetFanout();
            virtual ~FrameToTorchSetFanout() {};

            // INode, override because we get multiplicity at run time.
            virtual std::vector<std::string> output_types();

            // IFanout
            virtual bool operator()(const input_pointer& in, output_vector& outv);

            // IConfigurable
            virtual void configure(const WireCell::Configuration& cfg);
            virtual WireCell::Configuration default_configuration() const;

           private:
            int m_multiplicity;

            //Wire Planes to pack together into an output TorchTensorSet
            std::map<const WirePlaneId, int> m_output_groups;
            //How many wires in the TorchTensor in the ith output
            std::map<int, int> m_output_nchannels;

            //Channel Map
            std::unordered_map<int, int> m_channel_map;
            std::unordered_map<int, int> m_channel_to_output_group;
            std::vector<Json::Value> m_per_group_channel_map;


            //Expand dimensionality on output
            bool m_unsqueeze_output{false};

            std::string m_anode_tn{"AnodePlane"};
            IAnodePlane::pointer m_anode;

            bool m_debug_force_cpu{false};


            std::unordered_map<size_t, std::vector<std::pair<size_t,size_t>>> m_channel_ranges;
            //TODO -- possibly add configuration
            //allowing for different behavior if receive unexpected ticks

            WireCell::Configuration m_cfg;
            // Log::logptr_t log;
        };
    }  // namespace Aux
}  // namespace WireCell

#endif
