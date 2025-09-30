#include "WireCellSpng/Decon.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"
// #include "WireCellSpng/ITorchFieldResponse.h"
#include "WireCellSpng/ITorchSpectrum.h"
#include "WireCellUtil/FFTBestLength.h"
#include "WireCellSpng/Util.h"
#include "WireCellSpng/TensorTools.h"

// #include "WireCellSpng/ITorchColdElecResponse.h"

WIRECELL_FACTORY(SPNGDecon,
                 WireCell::SPNG::Decon,
                 WireCell::SPNG::ITorchTensorSetFilter,
                 WireCell::INamed,
                 WireCell::IConfigurable)

using namespace std::placeholders; // For _1, _2, _3

namespace WireCell::SPNG {

    Decon::Decon()
        : Logger("SPNGDecon", "spng")
    {
    }


    Decon::~Decon() {};

    WireCell::Configuration Decon::default_configuration() const
    {
        // We do not currently set any default values, so just forward our base classes'
        auto cfg = this->ContextBase::default_configuration();
        auto cfg2 = this->Logger::default_configuration();
        return update(cfg, cfg2);
    };


    void Decon::configure(const WireCell::Configuration& config)
    {
        this->ContextBase::configure(config);
        this->Logger::configure(config);

        m_passthrough.append("channel_map");

        m_frer_spectrum = get(config, "frer_spectrum", m_frer_spectrum);
        m_base_frer_spectrum = Factory::find_tn<ITorchSpectrum>(m_frer_spectrum);

        m_wire_filter = get(config, "wire_filter", m_wire_filter);
        m_base_wire_filter = Factory::find_tn<ITorchSpectrum>(m_wire_filter);

        m_coarse_time_offset = get(config, "coarse_time_offset", m_coarse_time_offset);
        m_pad_wire_domain = get(config, "pad_wire_domain", m_pad_wire_domain);
        m_use_fft_best_length = get(config, "use_fft_best_length", m_use_fft_best_length);

        m_debug_no_frer = get(config, "debug_no_frer", m_debug_no_frer);
        m_debug_no_wire_filter = get(config, "debug_no_wire_filter", m_debug_no_wire_filter);
        m_debug_no_roll = get(config, "debug_no_roll", m_debug_no_roll);

        m_unsqueeze_input = get(config, "unsqueeze_input", m_unsqueeze_input);

        m_output_set_tag = get(config, "output_set_tag", m_output_set_tag);
        m_output_tensor_tag = get(config, "output_tensor_tag", m_output_tensor_tag);
        log->debug("Will tag with Set:{} Tensor:{}", m_output_set_tag.asString(),
                   m_output_tensor_tag.asString());

        m_tensor_index = get<int>(config, "tensor_index", m_tensor_index);

        // auto wt = config["waveforms_tensor"];
        // if (wt.isString()) {
        //     const std::string datapath = wt.asString();
        //     get_input_itensor = [this, datapath](const input_pointer& in) {
        //         return this->get_input_tdm(in, datapath);
        //     };
        //     make_output_tensor = std::bind(Decon::make_output_old, this, _1, _2, _3)
        // }
        // else {
        //     size_t tensor_index = 0;
        //     if (wt.isInt()) {
        //         tensor_index = wt.asInt();
        //     }
        //     get_input_itensor = [this, tensor_index](const input_pointer& in) {
        //         return this->get_input_old(in, tensor_index);
        //     };
        //     make_output_tensor = std::bind(Decon::make_output_tdm, this, _1, _2, _3)
        // }
    }

    ITorchTensor::pointer Decon::get_input(const ITorchTensorSet::pointer& in, size_t index) const
    {
        auto iten = in->tensors()->at(index);
        if (iten == nullptr) {
            log->critical("decon finds no tensor at set index {}, call={}.  Check your config?", index, m_count);
            raise<ValueError>("no tensor input to decon, fix configuration?");
        }
        return iten;
    }


    ITorchTensorSet::pointer Decon::make_output(const ITorchTensorSet::pointer& intenset,
                                                const ITorchTensor::pointer& inten,
                                                torch::Tensor outten) const
    {
        Configuration set_md, tensor_md;

        tensor_md["tag"] = m_output_tensor_tag;
        metadata_passthrough(inten->metadata(), tensor_md, m_passthrough);

        set_md["tag"] = m_output_set_tag;
        auto output_tensor = std::make_shared<SimpleTorchTensor>(outten, tensor_md);
        std::vector<ITorchTensor::pointer> itv{
            output_tensor
        };
        return std::make_shared<SimpleTorchTensorSet>(
            intenset->ident(), set_md,
            std::make_shared<std::vector<ITorchTensor::pointer>>(itv)
            );
    }
    
    bool Decon::operator()(const input_pointer& in, output_pointer& out)
    {
        out = nullptr;
        if (!in) {
            logit("EOS");
            ++m_count;
            return true;
        }
        logit(in, "decon");

        // //Get the cloned tensor from the input
        // size_t ntensors = in->tensors()->size();
        // log->debug("Got {} tensors", ntensors);
        // std::vector<int64_t> prev_shape = {-999, -999};
        // for (size_t i = 0; i < ntensors; ++i) {
        //     auto this_sizes = in->tensors()->at(i)->tensor().sizes();
        //     if (i == 0) {
        //         prev_shape = this_sizes;
        //     }
        //     // if (prev_shape[0] != this_sizes[0] || prev_shape[1] != this_sizes[1]) {
        //     //     log->debug("Size mismatch in tensors");
        //     //     //TODO -- throw error
        //     // }
        //     log->debug("Tensor {} has shape {} {}", i, this_sizes[0], this_sizes[1]);
        // }

        // Get input tensor.
        auto input_itensor = get_input(in, m_tensor_index);
        const double period = get<double>(in->metadata(), "period", 0.0);

        auto orig_tensor = input_itensor->tensor();
        if (! orig_tensor.numel()) {
            log->critical("empty torch tensor on input at call={}, passing through the tensor set. Fix your config?",
                          m_count);
            raise<ValueError>("decon gets itensor with empty torch tensor, fix configuration?");
        }

        // Why do we clone?  Is there an early operation that is in-place?  Just to be safe?
        auto tensor_clone = to(orig_tensor.clone());
    
        if (m_unsqueeze_input) {
            tensor_clone = torch::unsqueeze(tensor_clone, 0);
        }

        tensor_clone = decon(tensor_clone, period);

        if (m_unsqueeze_input) {
            torch::squeeze(tensor_clone, 0);
        }

        out = make_output(in, input_itensor, tensor_clone);

        logit(out, "deconed");

        // log->debug("Tagged output with Set:{} Tensor:{}",
        //            out->metadata()["tag"],
        //            out->tensors()->at(0)->metadata()["tag"]);

        ++m_count;
        return true;
    }

    // The actual deconvolution.  The input tensor provides ADC waveforms and is
    // assumed to be batched.
    torch::Tensor Decon::decon(torch::Tensor waveforms, double period) const
    {
        auto sizes = waveforms.sizes();
        std::vector<int64_t> shape;
        for (const auto & s : sizes) {
            shape.push_back(s);
        }
        std::vector<int64_t> response_shape(shape.begin()+1, shape.end());

        int64_t original_nchans = shape[1];

    
        // TODO -- Padding in time domain then trim
        //      -- Later down the line overlap/add for extended readout
        // Always Pad in time because of non-periodicity in this domain
    
        if (m_pad_wire_domain) {
            // auto input_length = shape[0];
            auto response_length = m_base_frer_spectrum->shape()[0];

            auto padding_size = response_length - 1;

            if (m_use_fft_best_length) {
                auto suggested_length = fft_best_length(original_nchans + response_length - 1);
                log->debug("FFT best length suggested {}", suggested_length);
                padding_size = suggested_length - original_nchans;
            }

            //Pad to M+N-1 size 
            waveforms = torch::nn::functional::pad(
                waveforms,
                //The padding function option argument is the padding applied to the various sides of the tensor.
                //It's arranged last_dim_"left", last_dim_right, second_to_last_left, etc...
                //Later on, the FRER spectrum we get will be padded automatically on the wires dimension on the right
                //hence the order given here
                // Also we only give the additional size it must be padded with
                torch::nn::functional::PadFuncOptions({0, 0, 0, padding_size}).mode(torch::kConstant)

                );

            log->debug("Padded to {} ", waveforms.sizes()[1]);
            response_shape[0] = waveforms.sizes()[1];
        }

        //FFT on time dim
        waveforms = torch::fft::rfft(waveforms, std::nullopt, 2);

        //FFT on chan dim
        waveforms = torch::fft::fft(waveforms, std::nullopt, 1);

        //Get the Field x Elec. Response and do FFT in both dimensions.  Assure it is on our device.
        auto frer_spectrum_tensor = to(m_base_frer_spectrum->spectrum(response_shape).clone());

        frer_spectrum_tensor = torch::fft::rfft2(frer_spectrum_tensor);


        //Apply to input data
        if (!m_debug_no_frer)
            waveforms = waveforms / frer_spectrum_tensor;

        //Get the Wire filter -- already FFT'd
        //TODO -- fix the log here because of HfFilter weirdness
        auto wire_filter_tensor = m_base_wire_filter->spectrum({response_shape[0]});

        //Multiply along the wire dimension
        if (!m_debug_no_wire_filter)
            waveforms = waveforms * wire_filter_tensor.view({-1,1});

        //Inverse FFT in both dimensions
        waveforms = torch::fft::irfft2(waveforms);
    
        //Get the wire shift
        int wire_shift = m_base_frer_spectrum->shifts()[0];
        log->debug("Preparing to shift by {} wires in", wire_shift);

        //Shift along time dimension
        int time_shift = (int) (
            (m_coarse_time_offset + m_base_frer_spectrum->shifts()[1]) / period
            );

        if (!m_debug_no_roll) {
            //Shift along wire dimension
            waveforms = waveforms.roll(wire_shift, 1);
            waveforms = waveforms.roll(time_shift, 2);
        }

        if (m_pad_wire_domain) {
            //This selects the post-roll block
            using namespace torch::indexing;
            // auto unpadded_width = shape[0] - (m_base_frer_spectrum->shape()[0]-1);
            // auto unpadded_width = shape[0] - (m_base_frer_spectrum->shape()[0]-1);
            waveforms = waveforms.index({
                    Slice(), //Over all of the planes
                    Slice(None, original_nchans) //0 to N orig channels
                });
        }
        return waveforms;
    }

}
