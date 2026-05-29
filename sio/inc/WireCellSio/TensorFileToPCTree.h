#ifndef WIRECELL_SIO_TENSORFILETOPCTREE
#define WIRECELL_SIO_TENSORFILETOPCTREE

#include "WireCellAux/Logger.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/ITensorSetFilter.h"
#include "WireCellIface/ITensorSetSource.h"

#include <string>

namespace WireCell::Sio {

    /// Read a tensor archive from disk and attach its tensor as a named point
    /// cloud onto the live root node of an incoming point-cloud tree.
    ///
    /// This is a generic "side-data into the pctree root node" primitive: it
    /// does not interpret the tensor, it just stores it verbatim so a
    /// downstream component can read it back from the tree (mirroring the
    /// MicroBooNE root/UbooneClusterSource pattern of placing auxiliary data on
    /// the root node). One use is SBND charge-light matching, where the optical
    /// flash matrix is attached as a "flash" PC for QLMatching to consume.
    ///
    /// Input  (ITensorSetFilter): a pctree tensor set (live + dead).
    /// Output: the same tensor set with a "<pcname>" point cloud added to the
    ///         live root node, holding the file's first tensor verbatim as a
    ///         2-D "value" array.
    ///
    /// The file is read by a composed Sio::TensorFileSource (one tensor set per
    /// incoming pctree; the file and the pctree stream must be event-aligned).
    class TensorFileToPCTree : public Aux::Logger,
                               public ITensorSetFilter,
                               public IConfigurable {
    public:
        TensorFileToPCTree();
        virtual ~TensorFileToPCTree();

        bool operator()(const input_pointer& in, output_pointer& out) override;

        void configure(const WireCell::Configuration& cfg) override;
        WireCell::Configuration default_configuration() const override;

    private:
        // Tensor archive to read (e.g. "opflash_apa0.tar.gz") and the in-archive
        // tensor prefix, passed through to the composed source.
        std::string m_input{""};
        std::string m_prefix{""};
        // pctree datapath template and the name of the point cloud attached to
        // the live root node. pcname is required (no default).
        std::string m_inpath{"pointtrees/%d"};
        std::string m_pcname{""};

        std::size_t m_count{0};

        // Composed reader for the raw tensor archive.
        ITensorSetSource::pointer m_filesrc;
    };

} // namespace WireCell::Sio

#endif
