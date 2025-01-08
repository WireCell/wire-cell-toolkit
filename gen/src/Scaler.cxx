#include "WireCellGen/Scaler.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/String.h"

#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IAnodeFace.h"
#include "WireCellIface/IWirePlane.h"
#include "WireCellAux/SimpleDepo.h"

#include "WireCellUtil/Persist.h"

#include <boost/range.hpp>

#include <sstream>

WIRECELL_FACTORY(Scaler, WireCell::Gen::Scaler,
                 WireCell::INamed,
		 WireCell::IDrifter,
                 WireCell::IConfigurable)

using namespace std;
using namespace WireCell;

Gen::Scaler::Scaler()
  : Aux::Logger("Scaler", "gen")
{
}

Gen::Scaler::~Scaler() {}

WireCell::Configuration Gen::Scaler::default_configuration() const
{
  Configuration cfg;

  cfg["yzmap_scale_filename"] = "YZMap_Scale_filename";
  cfg["bin_width"]      = "BinWidth";
  cfg["tpc_width"]      = "TPCWidth";
  cfg["bin_height"]     = "BinHeight";
  cfg["anode"]          = "AnodePlane";
  cfg["plane"]          = "WirePlane";
  return cfg;
}

void Gen::Scaler::configure(const WireCell::Configuration& cfg)
{
  reset();

  const std::string filename = cfg["yzmap_scale_filename"].asString();
  if (filename.empty()) {
    THROW(ValueError() << errmsg{"Scaler requires an YZ region"});
  }

  //    log->debug("Conifg Anode TN");    
  const std::string anode_tn = cfg["anode"].asString();
  if (anode_tn.empty()) {
    THROW(ValueError() << errmsg{"Scaler requires an anode plane"});
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
  plane =      get<int>   (cfg, "plane");

  anode_name = get<std::string>(cfg, "anode");

  jmap = WireCell::Persist::load(filename);

}

void Gen::Scaler::reset() { }

bool scaler_by_time(const IDepo::pointer& lhs, const IDepo::pointer& rhs) { return lhs->time() < rhs->time(); }

// always returns true because by hook or crook we consume the input.
bool Gen::Scaler::operator()(const input_pointer& depo, output_queue& outq)
{


  const double Qi = depo->charge();

  if (Qi == 0.0) {
    // Yes, some silly depo sources ask us to drift nothing....
    return false;
  }

  for (auto box : m_boxes) {
    WireCell::Ray r = box.bounds();

    if (box.inside(depo->pos()) == false) {
      return false;
    }
  }

  double depo_y = depo->pos().y()*units::mm;
  double depo_z = depo->pos().z()*units::mm;
  double yoffset = 180*units::cm;
  double zoffset = 900*units::cm;

  int depo_bin_y = std::floor((depo_y+yoffset)/bin_height);
  int depo_bin_z = std::floor((depo_z+zoffset)/bin_width);

  if(depo_bin_y < 0){
    depo_bin_y = 0;
  }

  if(depo_bin_z < 0){
    depo_bin_z = 0;
  }

  if(depo_bin_y > 31){
    depo_bin_y = 31;
  }

  if(depo_bin_z > 180){
    depo_bin_z = 180;
  }

  //  double scale = (jmap[anode_name][std::to_string(plane)][depo_bin_y][depo_bin_z].asDouble());
  double scale = (jmap[anode_name][std::to_string(plane)][depo_bin_z][depo_bin_y].asDouble());

  auto newdepo = make_shared<Aux::SimpleDepo>(depo->time(), depo->pos(), Qi*scale, depo, depo->extent_long(), depo->extent_tran());

  outq.push_back(newdepo);

  std::sort(outq.begin(), outq.end(), scaler_by_time);

  return true;
}
