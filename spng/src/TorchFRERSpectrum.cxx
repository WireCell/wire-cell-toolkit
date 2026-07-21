#include "WireCellSpng/TorchFRERSpectrum.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellIface/IFieldResponse.h"
#include "WireCellIface/IWaveform.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/FFTBestLength.h"
#include "WireCellUtil/NumpyHelper.h"
// #include <torch/torch.h>

WIRECELL_FACTORY(TorchFRERSpectrum,
                 WireCell::SPNG::TorchFRERSpectrum,
                 WireCell::INamed,
                 WireCell::ITorchSpectrum,
                 WireCell::IConfigurable)

using namespace torch::indexing;

namespace WireCell::SPNG {

    TorchFRERSpectrum::TorchFRERSpectrum()
        : Logger("TorchFRERSpectrum", "spng")
    {
    }

    SPNG::TorchFRERSpectrum::~TorchFRERSpectrum() {}

    WireCell::Configuration SPNG::TorchFRERSpectrum::default_configuration() const
    {
        Configuration cfg = this->ContextBase::default_configuration();
        cfg["elec_response"] = m_elec_response_name;
        cfg["default_nticks"] = m_default_nticks;
        cfg["readout_period"] = m_readout_period;
        cfg["default_nchans"] = m_default_nchans;

        cfg["gain"] = m_gain;
        cfg["ADC_mV"] = m_ADC_mV;

        cfg["field_response"] = m_field_response_name;
        cfg["fr_plane_id"] = m_plane_id;

        return cfg;
    }

    void SPNG::TorchFRERSpectrum::configure(const WireCell::Configuration& cfg)
    {
        this->ContextBase::configure(cfg);

        m_field_response_name = get(cfg, "field_response", m_field_response_name);
        m_elec_response_name = get(cfg, "elec_response", m_elec_response_name);
        m_plane_id = get(cfg, "fr_plane_id", m_plane_id);

        m_readout_period = get(cfg, "readout_period", m_readout_period);

        if (! cfg["debug_force_cpu"].isNull()) {
            log->warn("the debug_force_cpu option is ignored, use 'device'");
        }

        m_gain = get(cfg, "gain", m_gain);
        m_ADC_mV = get(cfg, "ADC_mV", m_ADC_mV);

        auto ifr = Factory::find_tn<IFieldResponse>(m_field_response_name);
        m_field_response = ifr->field_response();
        m_field_response_avg = Response::wire_region_average(
            m_field_response
            );

        auto tfloat32 = tensor_options();

        bool found_plane = false;
        for (auto & plane : m_field_response_avg.planes) {
            if (plane.planeid != m_plane_id) continue;
        
            found_plane = true;
            m_fravg_nticks = fft_best_length(
                plane.paths[0].current.size()
                );
            m_fravg_period = m_field_response_avg.period;
            log->debug(String::format("Got %f as avg period", m_fravg_period));
            m_fravg_nchans = plane.paths.size();
            if (m_fravg_nchans == 0) {
                THROW(ValueError() <<
                      errmsg{String::format("TorchFRERSpectrum::%s: ", m_field_response_name)} <<
                      errmsg{"Got 0 nrows (electron paths)"});
            }

            if (m_fravg_nticks == 0) {
                THROW(ValueError() <<
                      errmsg{String::format("TorchFRERSpectrum::%s: ", m_field_response_name)} <<
                      errmsg{"Got 0 ncols"});
            }

            m_shape = {m_fravg_nchans, m_fravg_nticks};

            m_total_response = torch::zeros(m_shape, tfloat32);
            log->debug("Got {} {}", m_fravg_nchans, m_fravg_nticks);
            auto accessor = m_total_response.accessor<float,2>();
        
            for (int irow = 0; irow < m_fravg_nchans; ++irow) {
                auto& path = plane.paths[irow];
                for (size_t icol = 0; icol < plane.paths[0].current.size(); ++icol) {
                    accessor[irow][icol] = path.current[icol];
                }
            }

            log->debug("Doing FFT field resp");
            m_total_response = torch::fft::rfft(m_total_response, std::nullopt, 1);

        }
        if (!found_plane) {
            THROW(ValueError() <<
                  errmsg{String::format("TorchFRERSpectrum::%s: ", m_field_response_name)} <<
                  errmsg{String::format("Could not find plane %d", m_plane_id)});
        }


        m_elec_response = Factory::find_tn<IWaveform>(m_elec_response_name);
        torch::Tensor elec_response_tensor = torch::zeros(m_fravg_nticks, tfloat32);
        WireCell::Binning tbins(m_fravg_nticks, 0, m_fravg_nticks*m_fravg_period);
        auto ewave = m_elec_response->waveform_samples(tbins);
        auto accessor = elec_response_tensor.accessor<float,1>();
        for (int i = 0; i < m_fravg_nticks; ++i) {
            accessor[i] = ewave[i];
        }
        elec_response_tensor = (elec_response_tensor * m_gain * m_ADC_mV * (-1));

        log->debug("Doing FFT elec resp");
        elec_response_tensor = torch::fft::rfft(elec_response_tensor);

        log->debug("Multiplying FFT'd elec and field resp");
        m_total_response = (m_total_response*elec_response_tensor)*m_fravg_period;
    
        log->debug("IFFTing FRER");
        m_total_response = torch::fft::irfft(m_total_response, std::nullopt, 1);
        log->debug("Done");


        // TODO -- Improve redigitization. Brett suggests -- if larger sizes is multiple of other --
        //         integer downsampling in Freq. domain -- chop 1-N/M bits of freq.
        //          Or do LMN resampling technique
        int redigitized_nticks = std::ceil(m_total_response.size(1) * m_fravg_period / m_readout_period);
        m_redigitized_response = torch::zeros( {m_total_response.size(0), redigitized_nticks}, tfloat32);


        for (int i = 0; i < redigitized_nticks; i++) {
            double ctime = m_readout_period*i;
            int fcount = std::ceil(ctime/m_fravg_period);

            if (fcount < m_fravg_nticks) {
            
                auto vals_fcount = m_total_response.index({Slice(), fcount});
                auto vals_fcount_m1 = m_total_response.index({Slice(), fcount-1});

                m_redigitized_response.index_put_(
                    {Slice(None, m_fravg_nchans), i},
                    vals_fcount_m1 * ((ctime - m_fravg_period*(fcount - 1)) / m_fravg_period) +
                    vals_fcount * ((m_fravg_period*fcount - ctime) / m_fravg_period)
                    );
            }
            else {
                // std::cerr << "ERROR PASSED AVG NTICKS" << std::endl;
                THROW(RuntimeError() << errmsg {"Passed average nticks: " + std::to_string(m_fravg_nticks)});
            }
        }

    }

    torch::Tensor SPNG::TorchFRERSpectrum::spectrum(const std::vector<int64_t> & shape) const
    {

        auto the_tensor = torch::zeros(shape, tensor_options());
        the_tensor.index_put_(
            {Slice(None, m_redigitized_response.size(0)), Slice(None, m_redigitized_response.size(1))},
            m_redigitized_response
            );

        return the_tensor;
    }

    torch::Tensor SPNG::TorchFRERSpectrum::spectrum() const {
        return m_redigitized_response.clone();
    }

    // Get any shifts of the response
    std::vector<int64_t> SPNG::TorchFRERSpectrum::shifts() const {
        return {
            (m_fravg_nchans - 1)/2,
            int(m_field_response.origin/m_field_response.speed)//Maybe std::round here?
        };
    };

}
