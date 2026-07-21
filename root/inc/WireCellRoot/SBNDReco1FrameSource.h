/** Read SBND reco1 art/LArSoft ROOT files directly (bare ROOT, no
 * LArSoft) and emit the post-SP recob::Wire waveforms as one tagged
 * IFrame per configured event.
 *
 * Functional counterpart of larwirecell's wclsCookedFrameSource as used
 * by the SBND standalone-sample dumps (wcls-frame-dump.jsonnet): dense
 * float traces per channel (charge = frame_scale x wire value), frame
 * tag "dnnsp", per-channel summary from the wienersummary product
 * (x summary_scale), channel mask map "bad" from the badmasks product
 * ([channel, first, last] triplets).
 *
 * The art ROOT layout is an exchange format we do not own; this
 * component (with SBNDReco1OpFlashSource) is the only place that knows
 * about it.  See root/docs/sbnd-reco1-source.md.
 */
#ifndef WIRECELLROOT_SBNDRECO1FRAMESOURCE
#define WIRECELLROOT_SBNDRECO1FRAMESOURCE

#include "WireCellIface/IFrameSource.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"

namespace WireCell {
    namespace Root {
        class SBNDReco1FrameSource : public Aux::Logger, public IFrameSource, public IConfigurable {
          public:
            SBNDReco1FrameSource();
            virtual ~SBNDReco1FrameSource();

            virtual bool operator()(IFrame::pointer& out);

            virtual WireCell::Configuration default_configuration() const;
            virtual void configure(const WireCell::Configuration& config);

          private:
            std::string m_filename{""};
            // Event selection: non-negative "entry" wins, else run/event
            // (subrun < 0: any) are looked up in the Events tree.  With
            // neither set, ALL entries are streamed (one frame per call).
            int m_entry{-1};
            int m_run{-1};
            int m_subrun{-1};
            int m_event{-1};

            // Art product branch names (data reco1 defaults).
            std::string m_wire_product{"recob::Wires_sptpc2d_dnnsp_Reco1."};
            std::string m_badmask_product{"ints_sptpc2d_badmasks_Reco1."};
            std::string m_summary_product{"doubles_sptpc2d_wienersummary_Reco1."};

            std::string m_frame_tag{"dnnsp"};
            double m_frame_scale{50.0};
            double m_summary_scale{50.0};
            int m_nticks{3427};
            double m_tick{500.0};  // ns; interpreted with units at use site

            int m_calls{0};
        };
    }  // namespace Root
}  // namespace WireCell

#endif
