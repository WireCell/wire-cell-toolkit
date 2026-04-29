#include "WireCellGen/Scaler.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"

#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IAnodeFace.h"
#include "WireCellAux/SimpleDepo.h"

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
  cfg["yzmap"] = "";    // required: type:name of an IYZMap service
  cfg["anode"] = "";    // required: type:name of an IAnodePlane
  cfg["plane"] = plane;
  return cfg;
}

void Gen::Scaler::configure(const WireCell::Configuration& cfg)
{
  reset();

  const std::string yzmap_tn = cfg["yzmap"].asString();
  if (yzmap_tn.empty()) {
    THROW(ValueError() << errmsg{"Scaler: yzmap service is required"});
  }
  m_yzmap_svc = Factory::find_tn<IYZMap>(yzmap_tn);

  const std::string anode_tn = cfg["anode"].asString();
  if (anode_tn.empty()) {
    THROW(ValueError() << errmsg{"Scaler: anode is required"});
  }
  WireCell::IAnodePlane::pointer anode = Factory::find_tn<IAnodePlane>(anode_tn);
  if (anode == nullptr) {
    THROW(ValueError() << errmsg{"Scaler: anode is a nullptr"});
  }
  for (auto face : anode->faces()) {
    m_boxes.push_back(face->sensitive());
  }

  plane      = get<int>        (cfg, "plane", plane);
  anode_name = get<std::string>(cfg, "anode", anode_name);
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
    if (box.inside(depo->pos()) == false) {
      return false;
    }
  }

  double depo_y = depo->pos().y();
  double depo_z = depo->pos().z();

  double scale = m_yzmap_svc->value(anode_name, plane, depo_y, depo_z);

  auto newdepo = make_shared<Aux::SimpleDepo>(depo->time(), depo->pos(), Qi*scale, depo->prior(), depo->extent_long(), depo->extent_tran(), depo->prior()->id(), depo->prior()->pdg(), depo->prior()->energy());

  outq.push_back(newdepo);

  std::sort(outq.begin(), outq.end(), scaler_by_time);

  return true;
}
