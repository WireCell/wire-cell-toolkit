#ifndef WIRECELL_GEN_SCALER
#define WIRECELL_GEN_SCALER

#include "WireCellIface/IDrifter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IYZMap.h"
#include "WireCellAux/Logger.h"
#include "WireCellUtil/BoundingBox.h"

namespace WireCell {

  namespace Gen {

    class Scaler : public Aux::Logger,
      public IDrifter, public IConfigurable {
    public:
      Scaler();
      virtual ~Scaler();

      virtual void reset();
      virtual bool operator()(const input_pointer& depo, output_queue& outq);

      /// WireCell::IConfigurable interface.
      virtual void configure(const WireCell::Configuration& config);
      virtual WireCell::Configuration default_configuration() const;

      // Implementation methods.

    private:
      std::vector<WireCell::BoundingBox> m_boxes;

      int plane{0};
      std::string anode_name;
      IYZMap::pointer m_yzmap_svc{nullptr};


    };  // Scaler

  }  // namespace Gen

}  // namespace WireCell

#endif
