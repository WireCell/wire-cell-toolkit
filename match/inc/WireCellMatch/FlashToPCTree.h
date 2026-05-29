#ifndef WIRECELL_MATCH_FLASHTOPCTREE
#define WIRECELL_MATCH_FLASHTOPCTREE

#include "WireCellAux/Logger.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/ITensorSetFilter.h"
#include "WireCellIface/ITensorSetSource.h"

#include <string>

namespace WireCell::Match {

    /// Read SBND optical flashes from a tensor-archive file and attach them as
    /// a named point cloud onto the live root node of the incoming cluster
    /// point-cloud tree.
    ///
    /// This mirrors the MicroBooNE design (root/UbooneClusterSource), where the
    /// optical/light data is placed on the live root node and the downstream
    /// charge-light matcher reads it from there rather than from a separate
    /// input port. SBND's dedicated light I/O is this component.
    ///
    /// Input  (ITensorSetFilter): cluster pctree tensor set (live + dead).
    /// Output: the same tensor set with a "<pcname>" point cloud added to the
    ///         live root node holding the per-event flash matrix
    ///         [nflash, 1+nchan] (col 0 = time, cols 1..nchan = per-channel PE).
    ///
    /// The raw flash bytes are read via a composed Sio::TensorFileSource so the
    /// parsed values are byte-identical to the previous direct-tensor input.
    class FlashToPCTree : public Aux::Logger,
                          public ITensorSetFilter,
                          public IConfigurable {
    public:
        FlashToPCTree();
        virtual ~FlashToPCTree();

        bool operator()(const input_pointer& in, output_pointer& out) override;

        void configure(const WireCell::Configuration& cfg) override;
        WireCell::Configuration default_configuration() const override;

    private:
        // Flash tensor-archive file (e.g. "opflash_apa0.tar.gz") and the
        // in-archive tensor prefix passed through to the composed source.
        std::string m_input{""};
        std::string m_prefix{"opflash_"};
        // pctree datapath template (matches QLMatching's inpath) and the name
        // of the point cloud attached to the live root node.
        std::string m_inpath{"pointtrees/%d"};
        std::string m_pcname{"flash"};

        std::size_t m_count{0};

        // Composed reader for the raw flash tensor archive.
        ITensorSetSource::pointer m_flashsrc;
    };

} // namespace WireCell::Match

#endif
