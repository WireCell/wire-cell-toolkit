#include "WireCellClus/BeeSink.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/String.h"

WIRECELL_FACTORY(BeeSink, WireCell::Clus::BeeSink,
                 WireCell::INamed, WireCell::Clus::IBeeSink, WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Clus;

BeeSink::BeeSink()
  : Aux::Logger("BeeSink", "clus")
{
}

BeeSink::~BeeSink()
{
}

Configuration BeeSink::default_configuration() const
{
    Configuration cfg;
    cfg["outname"] = m_outname;
    cfg["initial_index"] = m_initial_index;
    return cfg;
}

void BeeSink::configure(const Configuration& cfg)
{
    m_outname = get<std::string>(cfg, "outname", m_outname);
    m_initial_index = get<int>(cfg, "initial_index", m_initial_index);
    m_templated = m_outname.find('%') != std::string::npos;
    if (m_templated) {
        // Opened lazily in write(), one zip per event (see the header).
        log->debug("shared Bee zip, one per event, name template {}", m_outname);
    }
    else {
        m_sink.reset(m_outname, m_initial_index);
        log->debug("shared Bee zip {} (initial index {})", m_outname, m_initial_index);
    }
}

void BeeSink::acquire()
{
    std::lock_guard<std::mutex> lk(m_mtx);
    ++m_refs;
}

void BeeSink::release()
{
    std::lock_guard<std::mutex> lk(m_mtx);
    if (m_refs > 0) { --m_refs; }
    if (m_refs == 0) {
        if (!m_templated || m_open_evt >= 0) {
            m_sink.close();
        }
        m_open_evt = -1;
        log->debug("closed shared Bee zip {}", m_outname);
    }
}

size_t BeeSink::write(const Bee::Object& obj, size_t index, int run, int sub, int evt)
{
    std::lock_guard<std::mutex> lk(m_mtx);

    if (m_templated) {
        if (evt != m_open_evt) {
            if (m_open_evt >= 0) {
                m_sink.close();
            }
            const std::string name = String::format(m_outname, evt);
            m_sink.reset(name, m_initial_index);
            m_open_evt = evt;
            log->debug("shared Bee zip {} (event {})", name, evt);
        }
        // Each per-event zip starts its own numbering, so the caller's
        // running event index (which counts events in the process) must not
        // be used as the slot inside this file.
        index = m_initial_index;
    }

    // Explicit index: set_index() clears the "seen" set so Bee::Sink's
    // name-collision auto-increment never fires; every object lands at the
    // caller-supplied event index regardless of inter-node write ordering.
    m_sink.set_rse(run, sub, evt);
    m_sink.set_index(index);
    return m_sink.write(obj);
}
