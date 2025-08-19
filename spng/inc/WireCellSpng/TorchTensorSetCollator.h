#ifndef WIRECELL_SPNG_TORCHTENSORSETCOLLATOR
#define WIRECELL_SPNG_TORCHTENSORSETCOLLATOR

#include "WireCellSpng/ITorchTensorSetFanin.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Logging.h"
#include "WireCellAux/Logger.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/WirePlaneId.h"

namespace WireCell {
    namespace SPNG {

        // Fan out 1 frame to N set at construction or configuration time.
        class TorchTensorSetCollator
            : public Aux::Logger,
              public ITorchTensorSetFanin, public IConfigurable {
           public:
            TorchTensorSetCollator();
            virtual ~TorchTensorSetCollator() {};

            // INode, override because we get multiplicity at run time.
            virtual std::vector<std::string> input_types();

            // IFanin
            virtual bool operator()(const input_vector& inv, output_pointer& out);

            // IConfigurable
            virtual void configure(const WireCell::Configuration& cfg);
            virtual WireCell::Configuration default_configuration() const;

           private:
            int m_multiplicity{0};
            Json::Value m_output_set_tag;
            WireCell::Configuration m_cfg;
        };
    }  // namespace Aux
}  // namespace WireCell

#endif