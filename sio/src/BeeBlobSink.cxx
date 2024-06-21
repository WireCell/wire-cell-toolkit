#include "WireCellSio/BeeBlobSink.h"

#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"

#include "WireCellAux/Bee.h"

WIRECELL_FACTORY(BeeBlobSink,
                 WireCell::Sio::BeeBlobSink,
                 WireCell::INamed,
                 WireCell::IBlobSetSink,
                 WireCell::ITerminal,
                 WireCell::IConfigurable)


using namespace WireCell;

Sio::BeeBlobSink::BeeBlobSink()
    : Aux::Logger("BeeBlobSink", "sio")
{
}

Sio::BeeBlobSink::~BeeBlobSink()
{
}


void Sio::BeeBlobSink::finalize()
{
    flush();
    m_sink.close();
}


void Sio::BeeBlobSink::configure(const WireCell::Configuration& cfg)
{
    m_outname = get(cfg, "outname", m_outname);
    m_sink.reset(m_outname);

    // required
    if (cfg["geom"].isNull()) {
        raise<ValueError>("BeeBlobSink requires a 'geom' to give canonical Bee detector name");
    }
    m_obj.detector(cfg["geom"].asString());
    // optional
    m_obj.algorithm(get<std::string>(cfg, "type", "unknown"));

    m_ident_offset = get(cfg,"evt",0);
    m_obj.rse(get(cfg,"run",0), get(cfg,"sub",0), m_ident_offset);

    // String or array of string, can not be empty
    m_samplers.clear();
    auto samplers = cfg["samplers"];
    if (samplers.isNull() || samplers.empty()) {
        raise<ValueError>("BeeBlobSink requires at least one entry in the \"samplers\" configuration parameter");
    }
    if (samplers.isString()) {
        m_samplers.push_back(Factory::find_tn<IBlobSampler>(samplers.asString()));
    }
    else {
        for (auto jtn : samplers) {
            m_samplers.push_back(Factory::find_tn<IBlobSampler>(jtn.asString())); 
        }
    }

}

WireCell::Configuration Sio::BeeBlobSink::default_configuration() const
{
    Configuration cfg;
    // required: "geom"

    cfg["samplers"] = Json::arrayValue;
    cfg["outname"] = m_outname;
    cfg["type"] = "unknown";
    cfg["run"] = cfg["sub"] = cfg["evt"] = 0;
    return cfg;
}

void Sio::BeeBlobSink::flush(int ident)
{
    if (m_obj.empty()) {
        return;
    }
    log->debug("writing to {} at call={}", m_outname, m_calls);
    m_sink.write(m_obj);
    m_obj.reset(ident);
    m_last_ident = ident;
}

bool Sio::BeeBlobSink::operator()(const IBlobSet::pointer& bs)
{
    if (!bs) {
        flush();
        log->debug("EOS at call={}", m_calls++);
        return true;       // EOS
    }

    auto blobs = bs->blobs();
    if (blobs.empty()) {
        log->trace("have no blobs at call={}", m_calls++);
        return true;
    }
    log->trace("input {} blobs at call={}", blobs.size(), m_calls);

    // Get a new "event number".  Ideally this is frame number but if there is
    // no frame settle for slice or blob set.  The latter two will make each Bee
    // "event identifier" span one time slice.
    int ident = m_ident_offset;
    if (bs->slice()) {
        if (bs->slice()->frame()) {
            ident += bs->slice()->frame()->ident();
        }
        else {
            ident += bs->slice()->ident();
        }
    }
    else {
        ident += bs->ident();
    }

    if (m_last_ident < 0) {     // first time.
        m_obj.reset(ident);
        m_last_ident = ident;
    }
    else if (m_last_ident != ident) {
        flush(ident);
    }

    for (auto sampler : m_samplers) {
        auto obj = Aux::Bee::dump(bs, sampler);
        log->trace("sampled {} blobs getting {} points at call={}",
                   blobs.size(), obj.size(), m_calls);
        m_obj.append(obj);
    }

    log->trace("have {} points at call={}", m_obj.size(), m_calls);
    ++m_calls;
    return true;
}
