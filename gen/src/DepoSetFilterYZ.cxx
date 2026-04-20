#include "WireCellGen/DepoSetFilterYZ.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellAux/SimpleDepoSet.h"
#include "WireCellIface/IDepo.h"
#include "WireCellIface/IAnodePlane.h"

WIRECELL_FACTORY(DepoSetFilterYZ, WireCell::Gen::DepoSetFilterYZ, WireCell::INamed, WireCell::IDepoSetFilter,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Gen;

DepoSetFilterYZ::DepoSetFilterYZ()
  : Aux::Logger("DepoSetFilterYZ", "gen")
{
}
DepoSetFilterYZ::~DepoSetFilterYZ() {}

WireCell::Configuration DepoSetFilterYZ::default_configuration() const
{
  Configuration cfg;
  cfg["yzmap"]  = "";    // required: type:name of an IYZMap service
  cfg["resp"]   = resp;
  cfg["anode"]  = "";    // required: type:name of an IAnodePlane
  cfg["plane"]  = plane;
  return cfg;
}

void DepoSetFilterYZ::configure(const WireCell::Configuration& cfg)
{
  const std::string yzmap_tn = cfg["yzmap"].asString();
  if (yzmap_tn.empty()) {
    THROW(ValueError() << errmsg{"DepoSetFilterYZ: yzmap service is required"});
  }
  m_yzmap_svc = Factory::find_tn<IYZMap>(yzmap_tn);

  const std::string anode_tn = cfg["anode"].asString();
  if (anode_tn.empty()) {
    THROW(ValueError() << errmsg{"DepoSetFilterYZ: anode is required"});
  }
  WireCell::IAnodePlane::pointer anode = Factory::find_tn<IAnodePlane>(anode_tn);
  if (anode == nullptr) {
    THROW(ValueError() << errmsg{"DepoSetFilterYZ: anode is a nullptr"});
  }
  for (auto face : anode->faces()) {
    m_boxes.push_back(face->sensitive());
  }

  resp       = get<int>        (cfg, "resp",  resp);
  plane      = get<int>        (cfg, "plane", plane);
  anode_name = get<std::string>(cfg, "anode", anode_name);
}

bool DepoSetFilterYZ::operator()(const input_pointer& in, output_pointer& out)
{
  out = nullptr;
  if (!in) {
    log->debug("DepoSetFilterYZ fail with no input on call = {}", m_count);
    return true;
  }
  IDepo::vector output_depos;

  for (auto idepo : *(in->depos())) {
    bool pass_resp = false;
    bool pass_anod = false;

    for (auto box : m_boxes) {
      if (box.inside(idepo->pos())) {
	pass_anod = true;
	break;
      }
    }

    if(pass_anod == false){
      continue;}

    double depo_y = idepo->pos().y();
    double depo_z = idepo->pos().z();

    if ((int)m_yzmap_svc->value(anode_name, plane, depo_y, depo_z) == resp+1) {
      pass_resp = true;
    }

    if (pass_resp && pass_anod) {
      //log->debug(" Passed! Resp {} at Y : {} and Z : {} on Plane {} and Anode {} ", resp+1, depo_y, depo_z, plane, anode_name);
      output_depos.push_back(idepo);
    }
  }

  log->debug("call={} Number of Depos for a give APA={}", m_count, output_depos.size());
  out = std::make_shared<WireCell::Aux::SimpleDepoSet>(m_count, output_depos);
  ++m_count;
  return true;
}
