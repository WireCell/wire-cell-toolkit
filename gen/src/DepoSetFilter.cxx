#include "WireCellGen/DepoSetFilter.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellAux/SimpleDepoSet.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IDepo.h"

WIRECELL_FACTORY(DepoSetFilter, WireCell::Gen::DepoSetFilter, WireCell::INamed, WireCell::IDepoSetFilter,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Gen;

DepoSetFilter::DepoSetFilter()
  : Aux::Logger("DepoSetFilter", "gen")
{
}
DepoSetFilter::~DepoSetFilter() {}

WireCell::Configuration DepoSetFilter::default_configuration() const
{
    Configuration cfg;
    cfg["anode"] = "AnodePlane";
    return cfg;
}

void DepoSetFilter::configure(const WireCell::Configuration& cfg)
{
    const std::string anode_tn = cfg["anode"].asString();
    if (anode_tn.empty()) {
        THROW(ValueError() << errmsg{"DepoSetFilter requires an anode plane"});
    }
    WireCell::IAnodePlane::pointer anode = Factory::find_tn<IAnodePlane>(anode_tn);
    if (anode == nullptr) {
        THROW(ValueError() << errmsg{"Input anode is a nullptr"});
    }
    IAnodeFace::vector anode_faces = anode->faces();
    for (auto face : anode_faces) {
        auto bb = face->sensitive();
        m_boxes.push_back(bb);
        log->debug("face={} bounds={}", face->ident(), bb.bounds());
    }
}

bool DepoSetFilter::operator()(const input_pointer& in, output_pointer& out)
{
    out = nullptr;
    if (!in) {
        log->debug("EOS at call={}", m_count);
        return true;
    }
    IDepo::vector output_depos;
    auto indepos = in->depos();
    for (auto idepo : *indepos) {
        for (auto box : m_boxes) {
            WireCell::Ray r = box.bounds();
            if (! box.inside(idepo->pos())) {
                continue;
            }
            output_depos.push_back(idepo);
            break;
        }
    }
    const int nin = indepos->size();
    const int nout = output_depos.size();
    log->debug("call={} depos: input={} output={} dropped={}",
               m_count, nin, nout, nout-nin);
    out = std::make_shared<WireCell::Aux::SimpleDepoSet>(m_count, output_depos);
    ++m_count;
    return true;
}
