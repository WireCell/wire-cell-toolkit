/** Read SBND reco1 art/LArSoft ROOT files directly (bare ROOT, no
 * LArSoft) and emit one "opflash" ITensorSet per configured event, per
 * TPC, drop-in compatible with the tensors that larwirecell's
 * wclsOpFlashSource writes for the SBND standalone chain:
 *
 *  - tensor 0 "opflash": f8 [nflash, 1+npmts], col 0 = flash time
 *    (recob::OpFlash::Time() in us x units::microsecond, i.e. ns,
 *    trigger-relative), cols 1..npmts = PE per OpDet (padded with 0).
 *  - tensor-set metadata key "frame_apply_at_caf" (ns): per-event flash
 *    time re-reference to the CAF/trigger frame, derived by the ported
 *    sbndcode FrameShift computation (see SBNDReco1FrameShift.h).
 *    Controlled by "caf_offset_mode":
 *      "none"     -- key omitted (FlashTensorToOpticalPCs no-op; default)
 *      "auto"     -- frame_default - frame_etrig from the FrameShift port
 *      "override" -- use the "frame_apply_at_caf_override" value
 *
 * The art ROOT layout is an exchange format we do not own; this
 * component (with SBNDReco1FrameSource) is the only place that knows
 * about it.  See root/docs/sbnd-reco1-source.md.
 */
#ifndef WIRECELLROOT_SBNDRECO1OPFLASHSOURCE
#define WIRECELLROOT_SBNDRECO1OPFLASHSOURCE

#include "WireCellIface/ITensorSetSource.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"

namespace WireCell {
    namespace Root {
        class SBNDReco1OpFlashSource : public Aux::Logger,
                                       public ITensorSetSource,
                                       public IConfigurable {
          public:
            SBNDReco1OpFlashSource();
            virtual ~SBNDReco1OpFlashSource();

            virtual bool operator()(ITensorSet::pointer& out);

            virtual WireCell::Configuration default_configuration() const;
            virtual void configure(const WireCell::Configuration& config);

          private:
            std::string m_filename{""};
            // Event selection: non-negative "entry" wins, else run/event
            // (subrun < 0: any) are looked up in the Events tree.  With
            // neither set, ALL entries are streamed (one set per call).
            int m_entry{-1};
            int m_run{-1};
            int m_subrun{-1};
            int m_event{-1};

            // Art product branch names (data reco1 defaults; use the
            // opflashtpc1 branch for APA/TPC 1).
            std::string m_flash_product{"recob::OpFlashs_opflashtpc0__Reco1."};
            std::string m_ptb_product{"raw::ptb::sbndptbs_ptbdecoder__DECODE."};
            std::string m_tdc_product{"sbnd::timing::DAQTimestamps_tdcdecoder__DECODE."};
            std::string m_rawheader_prefix{"artdaq::detail::RawEventHeader_daq_RawEventHeader_"};

            int m_npmts{312};

            std::string m_caf_offset_mode{"none"};  // none | auto | override
            double m_caf_offset_override{0.0};      // ns, used by "override"

            int m_calls{0};
        };
    }  // namespace Root
}  // namespace WireCell

#endif
