#include "WireCellSio/TensorFileSink.h"
#include "WireCellUtil/Stream.h"
#include "WireCellUtil/String.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(TensorFileSink, WireCell::Sio::TensorFileSink,
                 WireCell::INamed,
                 WireCell::ITensorSetSink,
                 WireCell::ITerminal,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Sio;

TensorFileSink::TensorFileSink()
    : Aux::Logger("TensorFileSink", "sio")
{
}

TensorFileSink::~TensorFileSink()
{
}

WireCell::Configuration TensorFileSink::default_configuration() const
{
    Configuration cfg;
    cfg["outname"] = m_outname;
    cfg["prefix"] = m_prefix;
    return cfg;
}

void TensorFileSink::configure(const WireCell::Configuration& cfg)
{
    m_outname = get(cfg, "outname", m_outname);
    m_prefix = get<std::string>(cfg, "prefix", m_prefix);

    // A '%' in the name means one container per ITensorSet ident, opened
    // lazily in operator() (see the header).  Without one, everything below is
    // exactly as it was: a single container, opened here.
    m_templated = m_outname.find('%') != std::string::npos;
    if (m_templated) {
        log->debug("sink one container per ident, name template {} with prefix \"{}\"",
                   m_outname, m_prefix);
    }
    else {
        m_out.clear();
        custard::output_filters(m_out, m_outname);
        if (m_out.empty()) {
            const std::string msg = "ClusterFileSink: unsupported outname: " + m_outname;
            log->critical(msg);
            THROW(ValueError() << errmsg{msg});
        }
        m_open = true;
        log->debug("sink through {} filters to {} with prefix \"{}\"",
                   m_out.size(), m_outname, m_prefix);
    }
    m_dump_mode = get<bool>(cfg, "dump_mode", m_dump_mode);
    log->debug("dump_mode={}", m_dump_mode);
}

void TensorFileSink::closeout()
{
    if (!m_open) return;
    // reset() pops and closes every filter and device in the chain.  The
    // single-container path below keeps its historical pop() so its bytes are
    // untouched; a reopen must close the whole chain or the container is
    // truncated.
    m_out.reset();
    m_out.clear();
    m_open = false;
    m_open_ident = -1;
}

void TensorFileSink::reopen(int ident)
{
    closeout();
    const std::string name = String::format(m_outname, ident);
    custard::output_filters(m_out, name);
    if (m_out.empty()) {
        const std::string msg = "TensorFileSink: unsupported outname: " + name;
        log->critical(msg);
        THROW(ValueError() << errmsg{msg});
    }
    m_open = true;
    m_open_ident = ident;
    log->debug("sink through {} filters to {} with prefix \"{}\"",
               m_out.size(), name, m_prefix);
}

void TensorFileSink::finalize()
{
    log->debug("closing {} after {} calls", m_outname, m_count);
    if (m_templated) {
        closeout();
        return;
    }
    m_out.pop();
}

void TensorFileSink::numpyify(ITensor::pointer ten, const std::string& fname)
{
    const auto shape = ten->shape();
    if (shape.empty()) {
        return;
    }
    Stream::write(m_out, fname, ten->data(), shape, ten->dtype());
    m_out.flush();
}

void TensorFileSink::jsonify(const Configuration& md, const std::string& fname)
{
    // Stringify json
    std::stringstream ss;
    ss << md;
    auto mdstr = ss.str();

    // Custard stream protocol.
    m_out << "name " << fname << "\n"
          << "body " << mdstr.size() << "\n"
          << mdstr.data();

    m_out.flush();
}

bool TensorFileSink::operator()(const ITensorSet::pointer &in)
{
    if (!in) {             // EOS
        log->debug("see EOS at call={}", m_count++);
        return true;
    }

    if(m_dump_mode) {
        log->debug("dumping tensor set ident={} at call {}",
                   in->ident(), m_count);
        return true;
    }

    if (m_templated && in->ident() != m_open_ident) {
        reopen(in->ident());
    }

    const std::string pre = m_prefix + "tensor";
    const std::string sident = std::to_string(in->ident());
    auto tens = in->tensors();
    const size_t ntens = tens->size();

    jsonify(in->metadata(), pre + "set_" + sident + "_metadata.json");
    for (size_t ind=0; ind<ntens; ++ind) {
        auto ten = tens->at(ind);
        const std::string ppre = pre + "_" + sident + "_" + std::to_string(ind);
        jsonify(ten->metadata(), ppre + "_metadata.json");
        numpyify(ten, ppre + "_array.npy");
    }

    log->debug("write tensor set ident={} ntensors={} at call {}",
               sident, ntens, m_count);
    ++m_count;
    return true;
}

