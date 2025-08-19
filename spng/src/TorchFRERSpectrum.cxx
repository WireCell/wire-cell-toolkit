#include "WireCellSpng/TorchFRERSpectrum.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellIface/IFieldResponse.h"
#include "WireCellIface/IWaveform.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/FFTBestLength.h"
#include "WireCellUtil/NumpyHelper.h"
// #include <torch/torch.h>

WIRECELL_FACTORY(TorchFRERSpectrum, WireCell::SPNG::TorchFRERSpectrum, WireCell::ITorchSpectrum, WireCell::IConfigurable)

using namespace WireCell;

SPNG::TorchFRERSpectrum::TorchFRERSpectrum() : Aux::Logger("TorchFRERSpectrum", "spng"), m_cache(5) {}

SPNG::TorchFRERSpectrum::~TorchFRERSpectrum() {}

WireCell::Configuration SPNG::TorchFRERSpectrum::default_configuration() const
{
    Configuration cfg;
    cfg["anode_num"] = m_anode_num;
    cfg["elec_response"] = m_elec_response_name;
    cfg["extra_scale"] = m_extra_scale;
    // cfg["do_fft"] = m_do_fft;
    cfg["default_nticks"] = m_default_nticks;
    cfg["default_period"] = m_default_period;
    cfg["debug_force_cpu"] = m_debug_force_cpu;
    cfg["default_nchans"] = m_default_nchans;
    // cfg["tick_period"] = m_tick_period;
    // cfg["shaping"] = m_shaping;
    // cfg["gain"] = m_gain;
    cfg["inter_gain"] = m_inter_gain;
    cfg["ADC_mV"] = m_ADC_mV;

    cfg["field_response"] = m_field_response_name;
    cfg["fr_plane_id"] = m_plane_id;
    // cfg["do_fft"] = m_do_fft;
    // cfg["do_average"] = m_do_average;

    return cfg;
}

void SPNG::TorchFRERSpectrum::configure(const WireCell::Configuration& cfg)
{
    m_field_response_name = get(cfg, "field_response", m_field_response_name);
    m_elec_response_name = get(cfg, "elec_response", m_elec_response_name);
    m_plane_id = get(cfg, "fr_plane_id", m_plane_id);
    m_extra_scale = get(cfg, "extra_scale", m_extra_scale);

    m_default_nticks = get(cfg, "default_nticks", m_default_nticks);
    m_default_nchans = get(cfg, "default_nchans", m_default_nchans);
    m_default_period = get(cfg, "default_period", m_default_period);

    m_debug_force_cpu = get(cfg, "debug_force_cpu", m_debug_force_cpu);

    m_inter_gain = get(cfg, "inter_gain", m_inter_gain);
    m_ADC_mV = get(cfg, "ADC_mV", m_ADC_mV);
    m_shape = {m_default_nchans, m_default_nticks};

    m_anode_num = get(cfg, "anode_num", m_anode_num);

    auto ifr = Factory::find_tn<IFieldResponse>(m_field_response_name);
    m_field_response = ifr->field_response();
    m_field_response_avg = Response::wire_region_average(
        m_field_response
    );

    {
        int nticks = 0;
        int nchans = 0;
        for (auto & plane : m_field_response.planes) {
            if (plane.planeid != m_plane_id) continue;
            
            //Metadata of avg FR
            nticks = plane.paths[0].current.size();
            nchans = plane.paths.size();
            break;
        }
        const std::string fname = "test_output_fr_full.npz";
        const std::string mode = "a";

        Array::array_xxf arr = Array::array_xxf::Zero(nchans, nticks);
        for (int irow = 0; irow < m_fravg_nchans; ++irow) {
            for (int icol = 0; icol < m_fravg_nticks; ++icol) {
                arr(irow, icol) = m_field_response.planes[m_plane_id].paths[irow].current[icol];
            }
        }
        const std::string aname = String::format(
            "apa_%i_plane_%i", m_anode_num, m_plane_id
        );
        WireCell::Numpy::save2d(arr, aname, fname, mode);
    }

    bool found_plane = false;
    for (auto & plane : m_field_response_avg.planes) {
        if (plane.planeid != m_plane_id) continue;
        
        found_plane = true;
        //Metadata of avg FR
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

        m_total_response = torch::zeros(m_shape, torch::TensorOptions().dtype(torch::kFloat64));
        log->debug("Got {} {}", m_fravg_nchans, m_fravg_nticks);
        auto accessor = m_total_response.accessor<double,2>();
        
        Array::array_xxf arr_alt = Array::array_xxf::Zero(m_fravg_nchans, m_fravg_nticks);

        for (int irow = 0; irow < m_fravg_nchans; ++irow) {
            auto& path = plane.paths[irow];
            for (int icol = 0; icol < plane.paths[0].current.size(); ++icol) {


                // if (std::abs(path.current[icol]) > 1.) {
                //     log->warn(String::format("Setting large value %f to 0.", path.current[icol]));
                //     accessor[irow][icol] = 0.;
                // }
                // else {
                    accessor[irow][icol] = path.current[icol];
                // }

                arr_alt(irow, icol) = path.current[icol];
            }
        }

        //Saving FR -- with & without large val protection
        {
            const std::string fname = "test_output_fr_alt.npz";
            const std::string mode = "a";
    
            const std::string aname = String::format(
                "apa_%i_plane_%i", m_anode_num, m_plane_id
            );
            WireCell::Numpy::save2d(arr_alt, aname, fname, mode);
        }

        {
            const std::string fname = "test_output_fr.npz";
            const std::string mode = "a";
    
            Array::array_xxf arr = Array::array_xxf::Zero(m_fravg_nchans, m_fravg_nticks);
            for (int irow = 0; irow < m_fravg_nchans; ++irow) {
                for (int icol = 0; icol < m_fravg_nticks; ++icol) {
                    arr(irow, icol) = accessor[irow][icol];
                }
            }
            const std::string aname = String::format(
                "apa_%i_plane_%i", m_anode_num, m_plane_id
            );
            WireCell::Numpy::save2d(arr, aname, fname, mode);
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
    torch::Tensor elec_response_tensor = torch::zeros(m_fravg_nticks, torch::TensorOptions().dtype(torch::kFloat64));
    WireCell::Binning tbins(m_fravg_nticks, 0, m_fravg_nticks*m_fravg_period);
    auto ewave = m_elec_response->waveform_samples(tbins);
    auto accessor = elec_response_tensor.accessor<double,1>();
    for (int i = 0; i < m_fravg_nticks; ++i) {
        accessor[i] = ewave[i]/* * m_inter_gain * m_ADC_mV * (-1)*/;
    }
    // elec_response_tensor *= m_inter_gain * m_ADC_mV * (-1);
    elec_response_tensor = (elec_response_tensor * m_inter_gain * m_ADC_mV * (-1));
    // std::cout << elec_response_tensor << std::endl;

    {
        const std::string fname = "test_output_er.npz";
        const std::string mode = "a";
        const std::string aname = String::format(
            "apa_%i_plane_%i", m_anode_num, m_plane_id
        );
        cnpy::npz_save(fname, aname, ewave.data(), {m_fravg_nticks}, mode);
    }
    
    log->debug("Doing FFT elec resp");
    elec_response_tensor = torch::fft::rfft(elec_response_tensor);

    log->debug("Multiplying FFT'd elec and field resp");
    m_total_response = (m_total_response*elec_response_tensor)*m_fravg_period;
    
    std::cout << "sizes" << std::endl;
    for (auto & s : m_total_response.sizes()) std::cout << s << std::endl;

    log->debug("IFFTing FRER");
    m_total_response = torch::fft::irfft(m_total_response, std::nullopt, 1);
    for (auto & s : m_total_response.sizes()) std::cout << s << std::endl;
    log->debug("Done");

    //Redigitize according to default nchans, nticks.
    std::vector<int64_t> default_shape = {m_default_nchans, m_default_nticks};
    redigitize(default_shape);
    
    const std::string fname = "test_output_frer.npz";
    Array::array_xxf arr = Array::array_xxf::Zero(default_shape[0], default_shape[1]);
    auto output = m_cache.get(default_shape).value().to(torch::kCPU);
    auto output_accessor = output.accessor<double,2>();
    for (int irow = 0; irow < default_shape[0]; ++irow) {
        for (int icol = 0; icol < default_shape[1]; ++icol) {
            arr(irow, icol) = output_accessor[irow][icol];
        }
    }
    const std::string mode = "a";
    const std::string aname = String::format(
        "apa_%i_plane_%i_%i_%i", m_anode_num, m_plane_id, default_shape[0], default_shape[1]
    );
    WireCell::Numpy::save2d(arr, aname, fname, mode);

}


// TODO -- Improve redigitization. Brett suggests -- if larger sizes is multiple of other --
//         integer downsampling in Freq. domain -- chop 1-N/M bits of freq.
//          Or do LMN resampling technique
void SPNG::TorchFRERSpectrum::redigitize(
    const std::vector<int64_t> & input_shape
) {

    //TODO -- check that input shape is 2-dim
    // and anything else i.e. > fravg shape?
    // auto nchans = input_shape[0];
    auto nticks = input_shape[1];

    //This means we've gotten here by mistake -- TODO consider throwing
    if (m_cache.contains(input_shape)) {
        return;
    }

    auto the_tensor = torch::zeros(input_shape, torch::TensorOptions().dtype(torch::kFloat64));
    //If not, we need to create a new version
    // auto result_accessor = m_cache.get(input_shape).value().accessor<float,2>();
    auto result_accessor = the_tensor.accessor<double,2>();
    auto total_response_accessor = m_total_response.accessor<double,2>();
    for (int irow = 0; irow < m_fravg_nchans; ++irow) {
        // std::cout << "Redigitizing " << irow << std::endl;
        int fcount = 1;
        for (int i = 0; i < nticks; i++) {
            double ctime = m_default_period*i;
            // if (irow == 0) 
            //         std::cout << i << " Ctime: " << ctime << std::endl;
            if (fcount < m_fravg_nticks) {
                while (ctime > fcount*m_fravg_period) {
                    // if (irow == 0) std::cout << "\tftime: " << fcount*m_fravg_period << std::endl;
                    fcount++;
                    if (fcount >= m_fravg_nticks) break;
                }
            }

            if (fcount < m_fravg_nticks) {
                result_accessor[irow][i] = (
                    (ctime - m_fravg_period*(fcount - 1)) / m_fravg_period * total_response_accessor[irow][fcount - 1] +
                    (m_fravg_period*fcount - ctime) / m_fravg_period * total_response_accessor[irow][fcount]
                );
                // if (irow == 0) {
                //     std::cout << "\t" << ctime << " " << fcount << " " << std::setprecision(10) << m_fravg_period*(fcount - 1) << " " << std::setprecision(10) << m_fravg_period*fcount << std::endl;
                //     std::cout << "\tdiff: " << (ctime - m_fravg_period*(fcount - 1)) << " " << (m_fravg_period*fcount - ctime) << std::endl;
                //     std::cout << "\t" << m_fravg_period << " " << total_response_accessor[irow][fcount - 1] << " " << total_response_accessor[irow][fcount] << std::endl;
                //     std::cout << "\t" << result_accessor[irow][i] << std::endl;
                // }
            }

        }
    }
    bool has_cuda = torch::cuda::is_available();
    torch::Device device((
        (has_cuda && !m_debug_force_cpu) ? 
        torch::kCUDA : 
        torch::kCPU));
    m_cache.insert(input_shape, the_tensor.to(device));
}

torch::Tensor SPNG::TorchFRERSpectrum::spectrum(const std::vector<int64_t> & shape) {

    //TODO -- check that input shape is 2-dim

    //If this shape has been used recenty, return the cached value
    if (m_cache.contains(shape)) {
        return m_cache.get(shape).value();
    }

    //If not, we need to create a new version -- this is done in redigitize
    redigitize(shape);

    const std::string fname = "test_output_frer.npz";
    Array::array_xxf arr = Array::array_xxf::Zero(shape[0], shape[1]);
    auto output = m_cache.get(shape).value().to(torch::kCPU);
    auto output_accessor = output.accessor<double,2>();
    for (int irow = 0; irow < shape[0]; ++irow) {
        for (int icol = 0; icol < shape[1]; ++icol) {
            arr(irow, icol) = output_accessor[irow][icol];
        }
    }
    const std::string mode = "a";
    const std::string aname = String::format(
        "apa_%i_plane_%i_%i_%i", m_anode_num, m_plane_id, shape[0], shape[1]
    );
    WireCell::Numpy::save2d(arr, aname, fname, mode);
    //Finally, return the result
    return m_cache.get(shape).value();
}

torch::Tensor SPNG::TorchFRERSpectrum::spectrum() const {
    return torch::zeros({m_default_nchans, m_default_nticks});
}

/// Get any shifts of the response
std::vector<int64_t> SPNG::TorchFRERSpectrum::shifts() const {
    std::cout << "origin " << m_field_response.origin << " speed " <<
                 m_field_response.speed << " ratio " <<
                 m_field_response.origin/m_field_response.speed << std::endl;
    return {
        (m_fravg_nchans - 1)/2,
        int(m_field_response.origin/m_field_response.speed)//Maybe std::round here?
    };
};