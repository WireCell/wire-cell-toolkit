#include "WireCellGen/DepoSetScaler.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellAux/SimpleDepoSet.h"

WIRECELL_FACTORY(DepoSetScaler, WireCell::Gen::DepoSetScaler,
                 WireCell::INamed,
                 WireCell::IDepoSetFilter, WireCell::IConfigurable)


using namespace WireCell;
using namespace WireCell::Gen;

DepoSetScaler::DepoSetScaler()
  : Aux::Logger("DepoSetScaler", "gen")
{
}
DepoSetScaler::~DepoSetScaler()
{
}

WireCell::Configuration DepoSetScaler::default_configuration() const
{
  Configuration cfg;
  // The typename of the scaler to do the real work.
  cfg["scaler"] = "Scaler";
  return cfg;
}

void DepoSetScaler::configure(const WireCell::Configuration& cfg)
{
  auto name = get<std::string>(cfg, "scaler", "Scaler");
  m_scaler = Factory::find_tn<IDrifter>(name);
}

bool DepoSetScaler::operator()(const input_pointer& in, output_pointer& out)
{
  out = nullptr;
  if (!in) {                  // EOS
    log->debug("EOS at call={}", m_count);
    return true;
  }

  log->debug("DepoSetScaler about to scale our depos! ");

  // make a copy so we can append an EOS to flush the per depo
  // drifter.
  IDepo::vector in_depos(in->depos()->begin(), in->depos()->end());
  //    in_depos.push_back(nullptr); // input EOS

  log->debug("DepoSetScaler gonna scale {} depos! ", in_depos.size());

  double charge_in = 0, charge_out=0;
  IDepo::vector all_depos;
  for (auto idepo : in_depos) {
    IDrifter::output_queue more;        
    (*m_scaler)(idepo, more);
    all_depos.insert(all_depos.end(), more.begin(), more.end());



    //if (idepo) {
    //    charge_in += idepo->charge();
    // }
    //for (const auto& d : more) {
    //    if (d) {
    //        charge_out += d->charge();
    //    }
    //  }
  }
  // The EOS comes through
  //all_depos.pop_back();

  log->debug("call={} drifted ndepos={} in={} out={}", m_count, all_depos.size(), charge_in, charge_out);
  out = std::make_shared<Aux::SimpleDepoSet>(m_count, all_depos);
  ++m_count;

  return true;
}
