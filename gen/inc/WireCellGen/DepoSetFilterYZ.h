//  Module is designed to filter depos based on provided APA dimentions.
//  It takes pointer to full stack of depos from DepoFanout and outputs
//  a poonter to only ones contained in a given volume

#ifndef WIRECELLGEN_DEPOSETFILTERYZ
#define WIRECELLGEN_DEPOSETFILTERYZ

#include "WireCellIface/IDepoSetFilter.h"
#include "WireCellIface/INamed.h"
#include "WireCellIface/IConfigurable.h"
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

    double bin_width;
    double tpc_width;
    double bin_height;
    double yoffset;
    double zoffset;    
    int nbinsy;
    int nbinsz;
    int resp;
    int plane;
    std::string anode_name;
    //Json::Value map = Json::arrayValue;
    //Json::arrayValue map;
    Json::Value jmap;


  };

}  // namespace WireCell::Gen

#endif
