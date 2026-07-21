#ifndef WIRECELLHIO_TENSORSINK
#define WIRECELLHIO_TENSORSINK

#include "WireCellIface/ITensorSetSink.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/INamed.h"
#include "WireCellAux/Logger.h"
#include "WireCellUtil/HanaJsonCPP.h"

#include <hdf5.h>
#include <string>
#include <vector>

namespace WireCell::Hio {

    /// Configuration for HioTensorSink
    struct TensorSinkConfig {

        /// The file name to write.  Required.
        std::string filename = "";

        /// Gzip compression level.  0 is no compression (default) with higher
        /// compression and longer running time up to an including level 9.  The
        /// trade-off of compression factor and speed varies greatly by the
        /// content of the data.
        int gzip = 0;

        /// Chunk sizes for each dimension for compression.  Generally, total
        /// size should be around or less 1MB.  This can be left empty for a
        /// heuristic choice that is probably good enough.
        std::vector<int> chunks = {};

        /// The pattern to determine the datapath of tensor sets.  The pattern
        /// may include format variables including the "ident" of the set and
        /// any that are be included in the set's metadta.
        std::string datapath_pattern = "tensorsets/{ident}";

        /// Whether to create HDF5 links from metadata "datapath" entries.
        /// Default is true.
        bool links = true;
    };

}  // namespace WireCell::Hio

BOOST_HANA_ADAPT_STRUCT(WireCell::Hio::TensorSinkConfig,
                         filename, gzip, chunks, datapath_pattern, links);

namespace WireCell::Hio {

    class TensorSink : public Aux::Logger,
                       public virtual IConfigurable,
                       public virtual ITensorSetSink
    {
    public:
        TensorSink();
        virtual ~TensorSink();

        // IConfigurable API
        virtual void configure(const WireCell::Configuration& jconfig);
        virtual WireCell::Configuration default_configuration() const;

        // ITensorSetSink API
        virtual bool operator()(const std::shared_ptr<const ITensorSet>& in);

        // INamed API
        virtual void set_name(const std::string& name) { m_name = name; }
        virtual std::string get_name() const { return m_name; }

    private:
        TensorSinkConfig m_config;
        hid_t m_file_id{-1};
        std::string m_name;
    };

}  // namespace WireCell::Hio

#endif
