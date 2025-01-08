#include "WireCellGen/DepoSetFilterYZ.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellAux/SimpleDepoSet.h"
#include "WireCellIface/IDepo.h"
#include "WireCellUtil/Persist.h"
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
  //../../util/src/Response.cxx

  cfg["yzmap_filename"] = "YZMap_filename";
  cfg["bin_width"]      = "BinWidth";
  cfg["tpc_width"]      = "TPCWidth";
  cfg["bin_height"]     = "BinHeight";
  cfg["yoffset"]        = "YOffset";
  cfg["zoffset"]        = "ZOffset";
  cfg["nbinsy"]         = "NbinsY";
  cfg["nbinsz"]         = "NbinsZ";
  cfg["resp"]           = "Response";
  cfg["anode"]          = "AnodePlane";
  cfg["plane"]          = "WirePlane";
  return cfg;
}

void DepoSetFilterYZ::configure(const WireCell::Configuration& cfg)
{

  //  log->debug("Conifg File Name");
  const std::string filename = cfg["yzmap_filename"].asString();
  if (filename.empty()) {
    THROW(ValueError() << errmsg{"DepoSetFilterYZ requires an YZ region"});
  }

  //    log->debug("Conifg Anode TN");    
  const std::string anode_tn = cfg["anode"].asString();
  if (anode_tn.empty()) {
    THROW(ValueError() << errmsg{"DepoSetFilterYZ requires an anode plane"});
  }
  //    log->debug("Config Andoe");
  WireCell::IAnodePlane::pointer anode = Factory::find_tn<IAnodePlane>(anode_tn);
  if (anode == nullptr) {
    THROW(ValueError() << errmsg{"Input anode is a nullptr"});
  }
  //    log->debug("Conifg Andode Face");
  IAnodeFace::vector abode_faces = anode->faces();
  for (auto face : abode_faces) {
    m_boxes.push_back(face->sensitive());
  }
  //    log->debug("Rest...");
  bin_width =  get<double>(cfg, "bin_width");
  tpc_width =  get<double>(cfg, "tpc_width");
  bin_height = get<double>(cfg, "bin_height");
  yoffset =    get<double>(cfg, "yoffset");
  zoffset =    get<double>(cfg, "zoffset");
  nbinsy  =    get<int>   (cfg, "nbinsy");
  nbinsz  =    get<int>   (cfg, "nbinsz");
  resp =       get<int>   (cfg, "resp");
  plane =      get<int>   (cfg, "plane");
  anode_name = get<std::string>(cfg, "anode");
  jmap = WireCell::Persist::load(filename);

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
      WireCell::Ray r = box.bounds();

      if (box.inside(idepo->pos())) {
	pass_anod = true;
	break;
      }
    }

    if(pass_anod == false){
      continue;}

    double depo_y = idepo->pos().y()*units::mm;
    double depo_z = idepo->pos().z()*units::mm;
    //double yoffset = 180*units::cm;
    //double zoffset = 900*units::cm;


    int depo_bin_y = std::floor((depo_y+yoffset)/bin_height);
    int depo_bin_z = std::floor((depo_z+zoffset)/bin_width);

    if(depo_bin_y < 0)
      depo_bin_y = 0;

    if(depo_bin_z < 0)
      depo_bin_z = 0;

    if(depo_bin_y > nbinsy)
      depo_bin_y = nbinsy;

    if(depo_bin_z > nbinsz)
      depo_bin_z = nbinsz;

    if (jmap[anode_name][std::to_string(plane)][depo_bin_z][depo_bin_y].asInt() == resp+1) {pass_resp = true;}

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
