/** Read LArSoft-reconstructed PDHD optical flashes from the temporary
 * light-data ROOT files (opflashana + flashopdet + trigoff trees) and
 * emit one opflash ITensorSet per configured event.
 *
 * Output tensor-set schema (see flash/docs/design.md §3.4):
 *  - tensor 0 "opflash": f8 [nflash, 1+nchan], col 0 = flash time
 *    (ns, trigger-relative), cols 1..nchan = PE per OpDet (OpChannel
 *    == OpDet for PDHD, nchan = 160).
 *  - tensor 1 "flash_summary" (optional): f8 [nflash, 8].
 *  - tensor 2 "ophits" (optional): f8 [nhit, 9].
 *
 * The ROOT layout is a temporary exchange format; this component is
 * the only place (with PDHDOpWaveformSource) that knows about it.
 */
#ifndef WIRECELLROOT_PDHDOPFLASHSOURCE
#define WIRECELLROOT_PDHDOPFLASHSOURCE

#include "WireCellIface/ITensorSetSource.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"

namespace WireCell {
    namespace Root {
        class PDHDOpFlashSource : public Aux::Logger, public ITensorSetSource, public IConfigurable {
          public:
            PDHDOpFlashSource();
            virtual ~PDHDOpFlashSource();

            virtual bool operator()(ITensorSet::pointer& out);

            virtual WireCell::Configuration default_configuration() const;
            virtual void configure(const WireCell::Configuration& config);

          private:
            std::string m_filename{""};
            int m_run{-1};
            int m_subrun{-1};  // <0: any
            int m_event{-1};
            int m_nchan{160};
            double m_nominal_offset_us{250.0};
            bool m_include_summary{true};
            bool m_include_ophits{true};

            int m_calls{0};
        };
    }  // namespace Root
}  // namespace WireCell

#endif
