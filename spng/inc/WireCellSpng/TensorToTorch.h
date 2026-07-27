#ifndef WIRECELL_SPNG_TENSORTOTORCH
#define WIRECELL_SPNG_TENSORTOTORCH

#include "WireCellIface/IConfigurable.h"
#include "WireCellSpng/ITensorToTorchSet.h"
#include "WireCellAux/Logger.h"

namespace WireCell {
namespace SPNG {

class TensorToTorch : public Aux::Logger,
                      public WireCell::IConfigurable,
                      public WireCell::SPNG::ITensorToTorchSet {
 public:
  TensorToTorch();
  virtual ~TensorToTorch();
  virtual WireCell::Configuration default_configuration() const;
  virtual bool operator()(const input_pointer& in, output_pointer& out);
  virtual void configure(const WireCell::Configuration& cfg);
};

}  // namespace SPNG
}  // namespace WireCell

#endif