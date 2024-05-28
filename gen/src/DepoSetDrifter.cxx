#include "WireCellGen/DepoSetDrifter.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellAux/SimpleDepoSet.h"

WIRECELL_FACTORY(DepoSetDrifter, WireCell::Gen::DepoSetDrifter,
                 WireCell::INamed,
                 WireCell::IDepoSetFilter, WireCell::IConfigurable)


using namespace WireCell;
using namespace WireCell::Gen;

DepoSetDrifter::DepoSetDrifter()
    : Aux::Logger("DepoSetDrifter", "gen")
{
}
DepoSetDrifter::~DepoSetDrifter()
{
}

WireCell::Configuration DepoSetDrifter::default_configuration() const
{
    Configuration cfg;
    // The typename of the drifter to do the real work.
    cfg["drifter"] = "Drifter";
    return cfg;
}

void DepoSetDrifter::configure(const WireCell::Configuration& cfg)
{
    auto name = get<std::string>(cfg, "drifter", "Drifter");
    m_drifter = Factory::find_tn<IDrifter>(name);
}

bool DepoSetDrifter::operator()(const input_pointer& in, output_pointer& out)
{
    out = nullptr;
    if (!in) {                  // EOS
        log->debug("EOS at call={}", m_count);
        return true;
    }

    // make a copy so we can append an EOS to flush the per depo
    // drifter.
    IDepo::vector in_depos(in->depos()->begin(), in->depos()->end());
    in_depos.push_back(nullptr); // input EOS

    double charge_in = 0, charge_out=0;
    IDepo::vector all_depos;
    for (auto idepo : in_depos) {
        IDrifter::output_queue more;        
        (*m_drifter)(idepo, more);
        all_depos.insert(all_depos.end(), more.begin(), more.end());

        if (idepo) {
            charge_in += idepo->charge();
        }
        for (const auto& d : more) {
            if (d) {
                charge_out += d->charge();
            }
        }
    }
    // The EOS comes through
    all_depos.pop_back();
        
    log->debug("call={} drifted ndepos={} Qout={} ({}%)", m_count, all_depos.size(), charge_out, 100.0*charge_out/charge_in);
    out = std::make_shared<Aux::SimpleDepoSet>(m_count, all_depos);
    ++m_count;
    return true;
}

