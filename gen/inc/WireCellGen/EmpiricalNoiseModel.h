/**
    A noise model based on empirically measured noise spectra.

    It requires configuration file holding a list of dictionary which
    provide association between wire length and the noise spectrum.

    TBD: document JSON file format for providing spectra and any other parameters.

*/

#ifndef WIRECELLGEN_EMPERICALNOISEMODEL
#define WIRECELLGEN_EMPERICALNOISEMODEL

#include "WireCellIface/IChannelSpectrum.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IDFT.h"
#include "WireCellIface/IChannelStatus.h"

#include "WireCellUtil/Units.h"
#include "WireCellUtil/Waveform.h"
#include "WireCellAux/Logger.h"

#include <string>
#include <vector>
#include <unordered_map>

namespace WireCell {
    namespace Gen {

        class EmpiricalNoiseModel : public Aux::Logger,
                                    public IChannelSpectrum, public IConfigurable {
          public:
            EmpiricalNoiseModel(const std::string& spectra_file = "",
                                const int nsamples = 10000,  // assuming 10k samples
                                const double period = 0.5 * units::us, const double wire_length_scale = 1.0 * units::cm,
                                /* const double time_scale = 1.0*units::ns, */
                                /* const double gain_scale = 1.0*units::volt/units::eplus, */
                                /* const double freq_scale = 1.0*units::megahertz, */
                                const std::string anode_tn = "AnodePlane",
                                const std::string chanstat_tn = "StaticChannelStatus");

            virtual ~EmpiricalNoiseModel();

          public: // interface

            /// IChannelSpectrum
            virtual const amplitude_t& channel_spectrum(int chid) const;

            // get constant term
            virtual const std::vector<float>& freq() const;
            // get json file gain
            virtual const double gain(int chid) const;
            // get json file shaping time
            virtual const double shaping_time(int chid) const;

            /// IConfigurable
            virtual void configure(const WireCell::Configuration& config);
            virtual WireCell::Configuration default_configuration() const;

          public: // Local methods, exposed to enable testing

            // fixme: this should be factored out into wire-cell-util.
            struct NoiseSpectrum {
                int plane;                 // plane identifier number
                int nsamples;              // number of samples used in preparing spectrum
                double period;             // sample period [time] used in preparing spectrum
                double gain;               // amplifier gain [voltage/charge]
                double shaping;            // amplifier shaping time [time]
                double wirelen;            // total length of wire conductor [length]
                double constant;           // amplifier independent constant noise component [voltage/frequency]
                std::vector<float> freqs;  // the frequencies at which the spectrum is sampled
                std::vector<float> amps;   // the amplitude [voltage/frequency] of the spectrum.
            };

            // Resample a NoiseSpectrum to match what the model was
            // configured to provide.  This method modifies in place.
            void resample(NoiseSpectrum& spectrum) const;

            // generate the default electronics response function at 1 mV/fC gain
            void gen_elec_resp_default();

            // Return new amplitude spectrum and constant term that
            // has been interpolated by wire_length between those
            // given in the spectra file.
            using spectrum_data_t = std::pair<double, amplitude_t>;
            spectrum_data_t interpolate_wire_length(int plane, double wire_length) const;
            const spectrum_data_t& get_spectrum_data(int iplane, int ilen) const;

            int get_nsamples() { return m_nsamples; };

           private:
            IAnodePlane::pointer m_anode;
            IChannelStatus::pointer m_chanstat;
            IDFT::pointer m_dft;

            std::string m_spectra_file;
            int m_nsamples;
            int m_fft_length;
            double m_period;
            double m_wlres{1.0 * units::cm}; // wire length scale -> ~10 - 800 cm
            double m_gres{0.1 * units::mV / units::fC}; // gain resolution, 0.1 mV/fC, 4.7, 7.8, 14, 25 -> 47 - 250
            double m_sres{0.1 * units::microsecond}; // shaping resolution, 0.5, 1, 2, 3. x1.1 when cold -> 10 - 30

            // double m_tres, m_gres, m_fres;
            std::string m_anode_tn, m_chanstat_tn;

            std::map<int, std::vector<NoiseSpectrum*> > m_spectral_data;

            // cache amplitudes to the nearest integral distance.
            mutable std::unordered_map<int, int> m_chid_to_intlen;

            typedef std::unordered_map<int, spectrum_data_t> len_amp_cache_t;
            mutable std::vector<len_amp_cache_t> m_len_amp_cache;

            // two layered caching to avoid per-channel cache
            // chid -> pack-key
            mutable std::unordered_map<int, unsigned int> m_channel_packkey_cache;
            using channel_amp_cache_t = std::unordered_map<unsigned int, amplitude_t>;
            mutable channel_amp_cache_t m_channel_amp_cache;
            
            // need to convert the electronics response in here ...
            Waveform::realseq_t m_elec_resp_freq;
            mutable std::unordered_map<int, Waveform::realseq_t> m_elec_resp_cache;

        };

    }  // namespace Gen
}  // namespace WireCell

#endif
