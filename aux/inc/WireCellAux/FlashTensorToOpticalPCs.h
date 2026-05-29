/** Convert a flash tensor (the SBND opflash matrix) into the toolkit's canonical
    optical point clouds and attach them to a cluster point-cloud tree's live root node.

    A 2-to-1 fan-in: port 0 carries a point-cloud-tree tensor set (live + dead);
    port 1 carries a tensor set whose first 2-D tensor is the SBND opflash matrix
    [nflash, 1+nchan] (column 0 = flash time, columns 1..nchan = per-channel PE).

    The matrix is expanded into the SAME three point clouds the MicroBooNE
    root/UbooneClusterSource writes onto the root "grouping" node, so SBND flashes
    interoperate with all clus tooling (Clus::Facade::Cluster::get_flash(),
    ClusterFlashDump, retile, BEE):
      - "flash"      : time, tmin, tmax, value(=total PE), ident(=row), type
      - "light"      : ident(=channel), time, value(=PE), error
      - "flashlight" : flash(=flash row), light(=light row)  (the join table)

    Units: the stored flash/light "time" is in WCT time units (ns) just like
    UbooneTTrees::load_optical, but it is taken from the matrix verbatim (NOT
    multiplied by units::us) because the SBND opflash dump is already in ns
    (WCT ns == 1), whereas Uboone's tree is in us and is scaled. Same canonical
    units, different native input unit. "error" is 0 (the SBND dump carries no
    per-channel error; any PE_err convention is applied downstream by the matcher).

    This converter is detector-agnostic (no SBND physics): it only reshapes a
    [nflash, 1+nchan] matrix into the canonical PCs, so it lives in aux alongside
    the generic Aux::AttachPointCloudToTree.
*/
#ifndef WIRECELLAUX_FLASHTENSORTOOPTICALPCS
#define WIRECELLAUX_FLASHTENSORTOOPTICALPCS

#include "WireCellAux/Logger.h"

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/ITensorSetFanin.h"

#include <string>
#include <vector>

namespace WireCell::Aux {

    class FlashTensorToOpticalPCs : public Aux::Logger,
                                    public ITensorSetFanin,
                                    public IConfigurable {
    public:
        FlashTensorToOpticalPCs();
        virtual ~FlashTensorToOpticalPCs();

        // IConfigurable
        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        // IFanin (port 0 = pctree, port 1 = opflash matrix tensor set)
        virtual std::vector<std::string> input_types();
        virtual bool operator()(const input_vector& invec, output_pointer& out);

    private:
        // pctree datapath template (formatted with the port-0 ident).
        std::string m_inpath{"pointtrees/%d"};
        // Number of optical-detector channels; the matrix must have ncol == nchan+1.
        int m_nchan{312};

        const int m_multiplicity{2};
        std::size_t m_count{0};
    };

} // namespace WireCell::Aux

#endif
