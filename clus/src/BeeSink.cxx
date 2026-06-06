#include "WireCellClus/BeeSink.h"

#include "WireCellUtil/NamedFactory.h"

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
    m_sink.reset(m_outname, m_initial_index);
    log->debug("shared Bee zip {} (initial index {})", m_outname, m_initial_index);
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
        m_sink.close();
        log->debug("closed shared Bee zip {}", m_outname);
    }
}

size_t BeeSink::write(const Bee::Object& obj, size_t index, int run, int sub, int evt)
{
    std::lock_guard<std::mutex> lk(m_mtx);
    // Explicit index: set_index() clears the "seen" set so Bee::Sink's
    // name-collision auto-increment never fires; every object lands at the
    // caller-supplied event index regardless of inter-node write ordering.
    m_sink.set_rse(run, sub, evt);
    m_sink.set_index(index);
    return m_sink.write(obj);
}
