/** Row-concatenate the "ophits" tensors of N input tensor-sets into one.
 *
 * For the PDHD all-PD light reconstruction (single processing, no post-hoc
 * file merge): the self-trigger snippet branch (opch 0-119) and the
 * full-stream branch (opch 120-159) each end in an OpHitFinder "ophits"
 * tensor f8 [nhit,9] (see OpHitFinder.h), both on the SAME trigger-relative
 * clock.  This fans them in -- vstacking the hit rows -- so a single
 * OpFlashFinder builds flashes spanning all 160 PDs (in particular -x flashes
 * over the full 80-159 wall, not just one half).  OpFlashFinder's own 1 us
 * accumulator does the cross-stream grouping, so there is no arbitrary dt
 * window to choose.
 *
 * Input:  N ITensorSets, each with an "ophits" tensor (by metadata name, else
 *         tensor 0); all must share the column count (9).
 * Output: one ITensorSet, ident + metadata copied from `meta_port`, carrying a
 *         single tensor named "ophits" = the row-concatenation in port order.
 */
#ifndef WIRECELLFLASH_OPHITMERGE
#define WIRECELLFLASH_OPHITMERGE

#include "WireCellIface/ITensorSetFanin.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"

#include <vector>
#include <string>

namespace WireCell {
    namespace Flash {
        class OpHitMerge : public Aux::Logger, public ITensorSetFanin, public IConfigurable {
          public:
            OpHitMerge();
            virtual ~OpHitMerge();

            // IFanin
            virtual std::vector<std::string> input_types();
            virtual bool operator()(const input_vector& invec, output_pointer& out);

            // IConfigurable
            virtual WireCell::Configuration default_configuration() const;
            virtual void configure(const WireCell::Configuration& config);

          private:
            int m_multiplicity{2};  // number of input ports (required)
            int m_meta_port{0};     // output ident/metadata taken from this port
        };
    }  // namespace Flash
}  // namespace WireCell

#endif
