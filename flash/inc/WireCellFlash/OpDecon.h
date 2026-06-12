/** Optical waveform deconvolution for PDHD light reconstruction.
 *
 * A dependency-free port of the duneopdet `Deconvolution` module
 * (Wiener filter built from a single-p.e. template, optional Gauss
 * post-filter and post baseline correction), operating on tagged
 * optical-snippet traces of an IFrame.  See
 * flash/docs/stage2-reconstruction.md for the fcl correspondence.
 *
 * Input: traces tagged `intag` (raw ADC snippets, channel = OpChannel).
 * Output: the input frame plus SPE-normalized traces tagged `outtag`.
 */
#ifndef WIRECELLFLASH_OPDECON
#define WIRECELLFLASH_OPDECON

#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IDFT.h"
#include "WireCellAux/Logger.h"

#include <map>
#include <vector>

namespace WireCell {
    namespace Flash {
        class OpDecon : public Aux::Logger, public IFrameFilter, public IConfigurable {
          public:
            OpDecon();
            virtual ~OpDecon();

            virtual bool operator()(const IFrame::pointer& in, IFrame::pointer& out);

            virtual WireCell::Configuration default_configuration() const;
            virtual void configure(const WireCell::Configuration& config);

          private:
            std::string m_intag{"raw"};
            std::string m_outtag{"decon"};
            // JSON file with SPE templates and the channel -> template map.
            std::string m_spe_file{"pgrapher/experiment/pdhd/pdhd-spe-templates.json"};
            // protodunehd_deconvolution fcl values.
            int m_samples{1024};
            int m_pre_trigger{50};
            int m_pedestal_buffer{30};
            double m_line_noise_rms{4.5};
            int m_input_polarity{-1};
            bool m_auto_scale{true};       // Wiener: scale from filtered SPE response
            double m_scale{1.0};           // used when auto_scale is false
            bool m_apply_postfilter{true};
            double m_postfilter_cutoff{1.5};  // MHz, Gauss post-filter
            bool m_apply_post_blcorr{true};

            IDFT::pointer m_dft;

            // Per-template precomputed state.
            struct SPETemplate {
                std::vector<float> wave;            // time domain, m_samples
                std::vector<std::complex<float>> fft;  // full-size spectrum
                double amplitude;                   // max(wave), >= 1
            };
            std::vector<SPETemplate> m_templates;
            std::map<int, size_t> m_chan2tmpl;
            std::vector<std::complex<float>> m_postfilter;  // full-size spectrum

            std::vector<float> deconvolve(const std::vector<float>& adc, const SPETemplate& spe) const;
        };
    }  // namespace Flash
}  // namespace WireCell

#endif
