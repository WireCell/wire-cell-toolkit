#ifndef WIRECELL_SPNG_TORCHTOTENSOR
#define WIRECELL_SPNG_TORCHTOTENSOR

#include "WireCellIface/IConfigurable.h"
#include "WireCellSpng/ITorchToTensorSet.h"
#include "WireCellAux/Logger.h"

namespace WireCell {
namespace SPNG {

class TorchToTensor : public Aux::Logger,
                      public WireCell::IConfigurable,
                      public WireCell::SPNG::ITorchToTensorSet {
 public:
  TorchToTensor();
  virtual ~TorchToTensor();
  virtual WireCell::Configuration default_configuration() const;
  virtual bool operator()(const input_pointer& in, output_pointer& out);
  virtual void configure(const WireCell::Configuration& cfg);
};

}  // namespace SPNG
}  // namespace WireCell

#endif