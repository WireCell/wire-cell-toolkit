#include "WireCellSpng/Torch1DSpectrum.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"
// #include <torch/torch.h>

WIRECELL_FACTORY(Torch1DSpectrum, WireCell::SPNG::Torch1DSpectrum, WireCell::ITorchSpectrum, WireCell::IConfigurable)

using namespace WireCell;

SPNG::Torch1DSpectrum::Torch1DSpectrum() : Aux::Logger("Torch1DSpectrum", "spng"), m_cache(5) {}

SPNG::Torch1DSpectrum::~Torch1DSpectrum() {}

WireCell::Configuration SPNG::Torch1DSpectrum::default_configuration() const
{
    Configuration cfg;
    cfg["default_length"] = m_default_length;
    cfg["debug_force_cpu"] = m_debug_force_cpu;
    return cfg;
}

void SPNG::Torch1DSpectrum::configure(const WireCell::Configuration& cfg)
{

    m_default_length = get(cfg, "default_length", m_default_length);
    m_debug_force_cpu = get(cfg, "debug_force_cpu", m_debug_force_cpu);
    m_total_spectrum = torch::ones(m_default_length);
    /*auto accessor = m_total_spectrum.accessor<float,1>();*/

    if (cfg.isMember("spectra")) {
        m_spectra.clear();
        m_spectra_tns.clear();
        for (auto tn: cfg["spectra"]) {
            m_spectra_tns.push_back(tn.asString());
            m_spectra.push_back(Factory::find_tn<IFilterWaveform>(tn.asString()));
            // m_spectrum_tensors.push_back(torch.zeros({m_default_length}));
            // auto accessor = m_spectrum_tensors.back().accessor<float,1>();
            // log->debug("Adding {} to list of spectra", m_spectra_tns.back());

            // auto vals = m_spectra.back()->filter_waveform(m_default_length);
            // for (size_t i = 0; i != vals.size(); i++) {
            //     accessor[i] = vals.at(i);
            // }

            // //FFT this spectrum
            // m_spectrum_tensors.back() = torch::fft::rfft(m_spectrum_tensors.back());
        }
    }

}

torch::Tensor SPNG::Torch1DSpectrum::spectrum() const {
    return torch::zeros({m_default_length});
}

torch::Tensor SPNG::Torch1DSpectrum::spectrum(const std::vector<int64_t> & shape) {
    //TODO -- add in throw if not 1D shape

    //If this shape has been used recenty, return the cached value
    if (m_cache.contains(shape)) {
        return m_cache.get(shape).value();
    }

    //If not, we need to create a new version
    // m_cache.insert(shape, torch::ones(shape));

    //Loop over the spectra
    bool has_cuda = torch::cuda::is_available();

    torch::Device device((
        (has_cuda && !m_debug_force_cpu) ? torch::kCUDA : torch::kCPU
    ));
    for (size_t i = 0; i < m_spectra.size(); ++i) {
        const auto & spectrum = m_spectra[i];

        //Get the waveform for this size
        auto vals = spectrum->filter_waveform(shape[0]);
        
        //Fill a tensor from it
        auto this_tensor = torch::ones(shape);
        auto accessor = this_tensor.accessor<float,1>();
        for (size_t i = 0; i < vals.size(); i++) {
            accessor[i] = vals.at(i);
        }

        //Multiply the cached tensor by the FFT of this tensor
        if (i == 0) {
            m_cache.insert(shape, /*torch::fft::fft*/(this_tensor.to(device)));
        }
        else {
            m_cache.get(shape).value() *= /*torch::fft::fft*/(this_tensor.to(device));
        }
    }
    return m_cache.get(shape).value();
}