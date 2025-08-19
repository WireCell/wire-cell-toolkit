#include "WireCellSpng/TorchTensorFileSink.h"
#include "WireCellUtil/Stream.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(SPNGTorchTensorFileSink, WireCell::SPNG::TorchTensorFileSink,
                 WireCell::INamed,
                 WireCell::ITorchTensorSetSink,
                 WireCell::ITerminal,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::SPNG;

TorchTensorFileSink::TorchTensorFileSink()
    : Aux::Logger("SPNGTorchTensorFileSink", "spng")
{
}

TorchTensorFileSink::~TorchTensorFileSink()
{
}

WireCell::Configuration TorchTensorFileSink::default_configuration() const
{
    Configuration cfg;
    cfg["outname"] = m_outname;
    cfg["prefix"] = m_prefix;
    return cfg;
}

void TorchTensorFileSink::configure(const WireCell::Configuration& cfg)
{
    m_outname = get(cfg, "outname", m_outname);
    m_out.clear();
    custard::output_filters(m_out, m_outname);
    if (m_out.empty()) {
        const std::string msg = "ClusterFileSink: unsupported outname: " + m_outname;
        log->critical(msg);
        THROW(ValueError() << errmsg{msg});
    }
    m_prefix = get<std::string>(cfg, "prefix", m_prefix);
    log->debug("sink through {} filters to {} with prefix \"{}\"",
               m_out.size(), m_outname, m_prefix);
    m_dump_mode = get<bool>(cfg, "dump_mode", m_dump_mode);
    log->debug("dump_mode={}", m_dump_mode);
}

void TorchTensorFileSink::finalize()
{
    log->debug("closing {} after {} calls", m_outname, m_count);
    m_out.pop();
}

void TorchTensorFileSink::numpyify(ITorchTensor::pointer ten, const std::string& fname)
{
    //Clone and make contiguous in memory
    auto tensor_clone = ten->tensor().clone().contiguous();
    const auto shape = ten->shape();
    if (shape.empty()) {
        return;
    }
    
    std::vector<float> as_vec(
            tensor_clone.data_ptr<float>(),
            tensor_clone.data_ptr<float>() + tensor_clone.numel());
    Stream::write<float>(
        m_out, fname,
        as_vec);
        //shape, ten->dtype());
    m_out.flush();
}

void TorchTensorFileSink::jsonify(const Configuration& md, const std::string& fname)
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

bool TorchTensorFileSink::operator()(const ITorchTensorSet::pointer &in)
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

