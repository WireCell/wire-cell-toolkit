#include "WireCellGen/DepoSetDrifter.h"
#include "WireCellUtil/BoundingBox.h"
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

    BoundingBox bb_in, bb_out, tt_in, tt_out;

    double charge_in = 0, charge_out=0;
    IDepo::vector all_depos;
    for (auto idepo : in_depos) {

        if (idepo) {            // could be eos
            charge_in += idepo->charge();
            bb_in(idepo->pos());
            tt_in(Point(idepo->time(), idepo->extent_long(), idepo->extent_tran()));
        }

        IDrifter::output_queue more;        
        (*m_drifter)(idepo, more);

        all_depos.insert(all_depos.end(), more.begin(), more.end());

        for (const auto& d : more) {
            if (!d) {            // EOS gets forwarded out
                continue;
            }
            charge_out += d->charge();
            bb_out(d->pos());
            tt_out(Point(d->time(), d->extent_long(), d->extent_tran()));
        }
    }
    // The EOS comes through
    all_depos.pop_back();
        
    log->debug("call={} drifted ndepos: {}->{}, Q: {}->{} ({}%)", m_count,
               in_depos.size(), all_depos.size(),
               charge_in, charge_out, 100.0*charge_out/charge_in);

    log->debug("depos in,  pos: {} -> {}", bb_in.bounds().first, bb_in.bounds().second);
    log->debug("depos in,  t/D: {} -> {}", tt_in.bounds().first, tt_in.bounds().second);
    log->debug("depos out, pos: {} -> {}", bb_out.bounds().first, bb_out.bounds().second);
    log->debug("depos out, t/D: {} -> {}", tt_out.bounds().first, tt_out.bounds().second);

    out = std::make_shared<Aux::SimpleDepoSet>(m_count, all_depos);
    ++m_count;
    return true;
}

