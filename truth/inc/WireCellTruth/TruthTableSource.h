/** TruthTableSource: serve the truth_nu / truth_pf tensors of one or
    more truth files (truth/python/dump_truth_table.py) as ITensorSets.

    One tensor set per configured ident, then EOS.  The set carries the
    two tensors of TruthTable::as_tensors() at datapaths "<prefix>/nu"
    and "<prefix>/pf" with prefix = "truth/%d" (the SBND attacher's
    layout), so that a downstream MultiAlgBlobClustering with
    aux_datatypes [truth_nu, truth_pf] and the SbndPrMagnifyTrackingVisitor
    consume them unchanged.

    wcp-porting-img pr_ml/docs/19 sec 3.4 (T4).
 */

#ifndef WIRECELLTRUTH_TRUTHTABLESOURCE
#define WIRECELLTRUTH_TRUTHTABLESOURCE

#include "WireCellAux/Logger.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/ITensorSetSource.h"

#include <string>
#include <vector>

namespace WireCell::Truth {

    class TruthTableSource : public Aux::Logger,
                             public ITensorSetSource,
                             public IConfigurable
    {
      public:
        TruthTableSource();
        virtual ~TruthTableSource();

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        virtual bool operator()(ITensorSet::pointer& out);

      private:
        /// File name pattern; "%d" is replaced by the ident.
        std::string m_inname{"truth_%d.npz"};
        /// Datapath prefix pattern; "%d" is replaced by the ident.
        std::string m_prefix{"truth/%d"};
        /// The idents to serve, in order.
        std::vector<int> m_idents{0};
        /// When true a missing file ends the stream instead of throwing.
        bool m_stop_on_missing{false};
        size_t m_count{0};
        bool m_eos{false};
    };
}

#endif
