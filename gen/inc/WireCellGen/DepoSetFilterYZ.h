//  Module is designed to filter depos based on provided APA dimentions.
//  It takes pointer to full stack of depos from DepoFanout and outputs
//  a poonter to only ones contained in a given volume

#ifndef WIRECELLGEN_DEPOSETFILTERYZ
#define WIRECELLGEN_DEPOSETFILTERYZ

#include "WireCellIface/IDepoSetFilter.h"
#include "WireCellIface/INamed.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IYZMap.h"
#include "WireCellAux/Logger.h"
#include "WireCellUtil/BoundingBox.h"

namespace WireCell::Gen {

  class DepoSetFilterYZ : public Aux::Logger, public IDepoSetFilter, public IConfigurable {
  public:
    DepoSetFilterYZ();
    virtual ~DepoSetFilterYZ();

    // IDepoSetFilterYZ
    virtual bool operator()(const input_pointer& in, output_pointer& out);

    /// WireCell::IConfigurable interface.
    virtual void configure(const WireCell::Configuration& config);
    virtual WireCell::Configuration default_configuration() const;

  private:
    std::vector<WireCell::BoundingBox> m_boxes;
    std::size_t m_count{0};

    int resp{0};
    int plane{0};
    std::string anode_name;
    IYZMap::pointer m_yzmap_svc{nullptr};


  };

}  // namespace WireCell::Gen

#endif
