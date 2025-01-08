#include "WireCellGen/FrameSummerYZ.h"

#include "WireCellAux/FrameTools.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellAux/SimpleFrame.h"

#include <iostream>

WIRECELL_FACTORY(FrameSummerYZ, WireCell::Gen::FrameSummerYZ, WireCell::IFrameFanin, WireCell::IConfigurable)

using namespace WireCell;
using WireCell::Aux::SimpleFrame;

Configuration Gen::FrameSummerYZ::default_configuration() const
{
    // fixme: maybe add operators, scaleing, offsets.

    Configuration cfg;

    cfg["multiplicity"] = (int) m_multiplicity;

    return cfg;
}

std::vector<std::string> Gen::FrameSummerYZ::input_types()
{
  const std::string tname = std::string(typeid(input_type).name());
  std::vector<std::string> ret(m_multiplicity,tname);    
  return ret;
}

void Gen::FrameSummerYZ::configure(const Configuration& cfg)
{
  
  int m = get<int>(cfg, "multiplicity", (int) m_multiplicity);
  if (m <= 0) {
    log->critical("illegal multiplicity: {}", m);
    THROW(ValueError() << errmsg{"FrameFanin multiplicity must be positive"});
  }
  m_multiplicity = m;

}


bool Gen::FrameSummerYZ::operator()(const input_vector& invec, output_pointer& out)
{
  //New

  out = nullptr;
  size_t neos = 0;
  
  for (const auto& fr : invec) {
    if (!fr) {
      ++neos;
    }    
  }
  
  if (neos) {
    log->debug("EOS at call={} with {}", m_count, neos);
    ++m_count;
    return true;
  }
  
  if (invec.size() != m_multiplicity) {
    log->critical("input vector size={} my multiplicity={}", invec.size(), m_multiplicity);
    THROW(ValueError() << errmsg{"input vector size mismatch"});
  }
  
  out = Aux::sum(invec, invec[0]->ident());
  return true;
}

Gen::FrameSummerYZ::FrameSummerYZ(size_t multiplicity)
  :Aux::Logger("FrameSummerYZ", "glue")
  , m_multiplicity(multiplicity)
{
}

Gen::FrameSummerYZ::~FrameSummerYZ() {}
