#ifndef WIRECELL_GEN_SCALER
#define WIRECELL_GEN_SCALER

#include "WireCellIface/IDrifter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IRandom.h"
#include "WireCellUtil/Units.h"
#include "WireCellAux/Logger.h"
#include "WireCellUtil/BoundingBox.h"

#include <set>

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
      std::size_t m_count{0};

      double bin_width;
      double tpc_width;
      double bin_height;
      int plane;
      std::string anode_name;
      //Json::Value map = Json::arrayValue;
      //Json::arrayValue map;
      Json::Value jmap;


    };  // Scaler

  }  // namespace Gen

}  // namespace WireCell

#endif
