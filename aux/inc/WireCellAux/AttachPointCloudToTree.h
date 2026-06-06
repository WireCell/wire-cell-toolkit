/** Attach a tensor set as a named point cloud onto a point-cloud tree's root.

    A 2-to-1 fan-in: port 0 carries a point-cloud-tree tensor set (live + dead);
    port 1 carries a tensor set whose first 2-D tensor is stored verbatim, as a
    named "value" array, in a point cloud on the LIVE root node of the tree.  The
    augmented tree (live with the new PC + dead unchanged) is emitted.

    This is the generic "park external side-data on the pctree root node"
    primitive (mirroring the MicroBooNE root/UbooneClusterSource pattern of
    placing optical data on the root node).  Pair it with any ITensorSetSource
    (e.g. Sio::TensorFileSource) feeding port 1.  It does not interpret the
    tensor; a downstream consumer reads it back from the named point cloud.
*/
#ifndef WIRECELLAUX_ATTACHPOINTCLOUDTOTREE
#define WIRECELLAUX_ATTACHPOINTCLOUDTOTREE

#include "WireCellAux/Logger.h"

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/ITensorSetFanin.h"

#include <string>
#include <vector>

namespace WireCell::Aux {

    class AttachPointCloudToTree : public Aux::Logger,
                                   public ITensorSetFanin,
                                   public IConfigurable {
    public:
        AttachPointCloudToTree();
        virtual ~AttachPointCloudToTree();

        // IConfigurable
        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        // IFanin (port 0 = pctree, port 1 = tensor set to attach)
        virtual std::vector<std::string> input_types();
        virtual bool operator()(const input_vector& invec, output_pointer& out);

    private:
        // pctree datapath template (formatted with the port-0 ident).
        std::string m_inpath{"pointtrees/%d"};
        // Name of the point cloud attached to the live root node (required).
        std::string m_pcname{""};

        const int m_multiplicity{2};
        std::size_t m_count{0};
    };

} // namespace WireCell::Aux

#endif
