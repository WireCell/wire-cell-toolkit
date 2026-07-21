#include "WireCellHio/TensorSink.h"
#include "WireCellHio/HIO.h"
#include "WireCellHio/Tensors.h"
#include "WireCellIface/INamed.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Fmt.h"

WIRECELL_FACTORY(HioTensorSink,
                 WireCell::Hio::TensorSink,
                 WireCell::ITensorSetSink,
                 WireCell::IConfigurable,
                 WireCell::INamed);

using WireCell::HanaJsonCPP::to_json;
using WireCell::HanaJsonCPP::from_json;

namespace WireCell::Hio {

    TensorSink::TensorSink()
        : Aux::Logger("HioTensorSink", "hio")
    {
    }

    TensorSink::~TensorSink()
    {
        if (m_file_id >= 0) {
            close(m_file_id);
            m_file_id = -1;
        }
    }

    void TensorSink::configure(const WireCell::Configuration& jconfig)
    {
        from_json(m_config, jconfig);

        if (m_config.filename.empty()) {
            raise<ValueError>("TensorSink: filename must not be empty");
        }

        // Open file in truncate mode
        m_file_id = open(m_config.filename, FileMode::trunc);
    }

    WireCell::Configuration TensorSink::default_configuration() const
    {
        return to_json(m_config);
    }

    bool TensorSink::operator()(const std::shared_ptr<const ITensorSet>& in)
    {
        // Handle EOS (end of stream)
        if (!in) {
            // Close file on EOS
            if (m_file_id >= 0) {
                close(m_file_id);
                m_file_id = -1;
            }
            return true;
        }

        // Get metadata and add ident
        Configuration params = in->metadata();
        params["ident"] = in->ident();

        // Compute datapath using pattern
        std::string datapath = Fmt::format(m_config.datapath_pattern, params);

        // Write the tensor set to HDF5
        write_itensorset(m_file_id, in, datapath, m_config.gzip, m_config.chunks);
        log->debug("wrote tensor set {} to {}", in->ident(), datapath);

        // Create links if configured
        if (m_config.links) {
            link_itensorset(m_file_id, in, datapath);
        }

        return true;
    }

}  // namespace WireCell::Hio
