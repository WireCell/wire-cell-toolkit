#ifndef WIRECELLSIGPROC_OMNICHANNELNOISEDB
#define WIRECELLSIGPROC_OMNICHANNELNOISEDB

#include "WireCellIface/IChannelNoiseDatabase.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IDFT.h"
#include "WireCellIface/IFieldResponse.h"
#include "WireCellIface/WirePlaneId.h"

#include "WireCellUtil/Waveform.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Logging.h"

#include <vector>
#include <tuple>
#include <unordered_map>
#include <memory>

namespace WireCell {
    namespace SigProc {

        class OmniChannelNoiseDB : public WireCell::IChannelNoiseDatabase, public WireCell::IConfigurable {
           public:
            /// Create a configurable channel noise DB for digitized
            /// waveforms with the given size and number of samples.
            /// Default is for microphone.
            OmniChannelNoiseDB();
            virtual ~OmniChannelNoiseDB();

            // IConfigurable
            virtual void configure(const WireCell::Configuration& config);
            virtual WireCell::Configuration default_configuration() const;

            // Push the actual frame size from the first frame into the DB so that
            // frequency-domain spectra (freqmasks, rcrc, …) are always sized to
            // match the real FFT length.  Only OmniChannelNoiseDB supports this;
            // OmnibusNoiseFilter calls it via dynamic_cast.
            int nsamples() const { return m_nsamples; }
            void set_nsamples(int n);

            // IChannelNoiseDatabase
            virtual double sample_time() const;

            virtual double nominal_baseline(int channel) const;
            virtual double gain_correction(int channel) const;
            virtual double response_offset(int channel) const;

            virtual double min_rms_cut(int channel) const;
            virtual double max_rms_cut(int channel) const;

            virtual int pad_window_front(int channel) const;
            virtual int pad_window_back(int channel) const;

            virtual float coherent_nf_decon_limit(int channel) const;
            virtual float coherent_nf_decon_lf_cutoff(int channel) const;
            virtual float coherent_nf_decon_limit1(int channel) const;
            virtual float coherent_nf_adc_limit(int channel) const;
            virtual float coherent_nf_protection_factor(int channel) const;
            virtual float coherent_nf_min_adc_limit(int channel) const;
            virtual float coherent_nf_roi_min_max_ratio(int channel) const;

            virtual const filter_t& rcrc(int channel) const;
            virtual const filter_t& config(int channel) const;
            virtual const filter_t& noise(int channel) const;
            virtual const filter_t& response(int channel) const;

            // todo:

            virtual std::vector<channel_group_t> coherent_channels() const { return m_channel_groups; }
            virtual channel_group_t bad_channels() const { return m_bad_channels; }
            virtual channel_group_t miscfg_channels() const { return m_miscfg_channels; }

           protected:
            // Allow subclasses some access so that they may leverage
            // all the configuration code this class provides while
            // supplying some subset from an external source.  These
            // setters should be called before each use of
            // IChannelNoiseDatabase interface and may be called
            // multiple times

            // Override the bad channels.
            virtual void set_bad_channels(const channel_group_t& bad) { m_bad_channels = bad; }

            // Override the reconfigured spectrum for a set of
            // channels.  If `reset` is true then all channels
            // reconfig spectrum is set to the default before applying
            // the reconfig spectrum associated with the given from/to
            // gain/shaping.  If false, then other channels are not
            // touched.
            virtual void set_misconfigured(const std::vector<int>& channels, double from_gain, double from_shaping,
                                           double to_gain, double to_shaping, bool reset = false);

           private:
            double m_tick;
            int m_nsamples;
            Json::Value m_cfg_channel_info;  // stashed for set_nsamples() rebuild
            IAnodePlane::pointer m_anode;
            IFieldResponse::pointer m_fr;
            int m_rc_layers;

            // per-channel summed wire length in cm; built lazily on first use
            std::unordered_map<int, double> m_wire_length_cm;
            bool m_wire_length_cached{false};
            void cache_wire_lengths();

            typedef std::shared_ptr<filter_t> shared_filter_t;
            typedef std::vector<shared_filter_t> filter_vector_t;

            // Embody the "database" entry for one channel.
            struct ChannelInfo {
                int chid;

                // direct scalar values
                double nominal_baseline, gain_correction, response_offset, min_rms_cut, max_rms_cut;
                int pad_window_front, pad_window_back;

                float decon_limit;
                float decon_lf_cutoff;
                float adc_limit;
                float decon_limit1;
                float protection_factor;
                float min_adc_limit;
                float roi_min_max_ratio;

                // parameters

                // frequency space filters
                shared_filter_t rcrc, config, noise, response;

                ChannelInfo();
            };

            // std::vector<ChannelInfo> m_db;
            // std::unordered_map<int, ChannelInfo*> m_db;
            std::unordered_map<int, ChannelInfo> m_db;

            const ChannelInfo& dbget(int ch) const
            {
                auto it = m_db.find(ch);
                if (it == m_db.end()) {
                    // it = m_db.find(defch);
                    // return *(it->second);
                    // m_db.insert(std::make_pair(ch, new ChannelInfo));
                    // return *(m_db[ch]);
                    THROW(KeyError() << errmsg{String::format("no db info for channel %d", ch)});
                }
                return it->second;
                // return m_db.at(ch);
            }

            std::vector<channel_group_t> m_channel_groups;
            channel_group_t m_bad_channels;
            channel_group_t m_miscfg_channels;

            // JSON parsing.  Exhausting.
            std::vector<int> parse_channels(const Json::Value& jchannels);
            shared_filter_t make_filter(std::complex<float> defval = std::complex<float>(1, 0));
            shared_filter_t default_filter();
            shared_filter_t parse_freqmasks(Json::Value jfm);
            shared_filter_t parse_rcrc(Json::Value jrcrc, int nrc);
            double parse_gain(Json::Value jreconfig);
            shared_filter_t parse_reconfig(Json::Value jreconfig);
            shared_filter_t get_reconfig(double from_gain, double from_shaping, double to_gain, double to_shaping);
            shared_filter_t parse_response(Json::Value jreconfig);
            // ChannelInfo* make_ci(int chid, Json::Value jci);
            void update_channels(Json::Value cfg);
            ChannelInfo& get_ci(int chid);

            // Reuse the same filter spectra for matching input parameters.

            // lookup by truncated rcrc value
            std::unordered_map<int, shared_filter_t> m_rcrc_cache;
            // lookup by OR of the four truncated values
            std::unordered_map<int, shared_filter_t> m_reconfig_cache;
            // lookup by explicit waveform id
            std::unordered_map<int, shared_filter_t> m_waveform_cache;
            // lookup by WirePlaneId::ident()
            std::unordered_map<int, shared_filter_t> m_response_cache;

            Log::logptr_t log;
            IDFT::pointer m_dft;
        };
    }  // namespace SigProc

}  // namespace WireCell

#endif

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
