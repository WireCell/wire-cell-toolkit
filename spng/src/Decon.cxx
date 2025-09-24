#include "WireCellSpng/Decon.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"
// #include "WireCellSpng/ITorchFieldResponse.h"
#include "WireCellSpng/ITorchSpectrum.h"
#include "WireCellUtil/FFTBestLength.h"
#include "WireCellSpng/Util.h"

// #include "WireCellSpng/ITorchColdElecResponse.h"

WIRECELL_FACTORY(SPNGDecon,
                 WireCell::SPNG::Decon,
                 WireCell::SPNG::ITorchTensorSetFilter,
                 WireCell::INamed,
                 WireCell::IConfigurable)

WireCell::SPNG::Decon::Decon()
  : Aux::Logger("SPNGDecon", "spng") {

}

WireCell::SPNG::Decon::~Decon() {};

WireCell::Configuration WireCell::SPNG::Decon::default_configuration() const
{
    Configuration cfg;
    cfg["tensor_index"] = m_tensor_index;
    return cfg;
};


void WireCell::SPNG::Decon::configure(const WireCell::Configuration& config) {

    m_passthrough.append("channel_map");

    m_frer_spectrum = get(config, "frer_spectrum", m_frer_spectrum);
    base_frer_spectrum = Factory::find_tn<ITorchSpectrum>(m_frer_spectrum);

    m_wire_filter = get(config, "wire_filter", m_wire_filter);
    base_wire_filter = Factory::find_tn<ITorchSpectrum>(m_wire_filter);

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

    m_tensor_index = get(config, "tensor_index", m_tensor_index);
}

bool WireCell::SPNG::Decon::operator()(const input_pointer& in, output_pointer& out) {
    out = nullptr;
    if (!in) {
        log->debug("EOS ");
        return true;
    }
    log->debug("Running Decon");

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

    auto tensor_clone = in->tensors()->at(m_tensor_index)->tensor().clone();
    
    if (m_unsqueeze_input)
        tensor_clone = torch::unsqueeze(tensor_clone, 0);

    auto sizes = tensor_clone.sizes();
    std::vector<int64_t> shape;
    for (const auto & s : sizes) {
        shape.push_back(s);
    }
    std::vector<int64_t> response_shape(shape.begin()+1, shape.end());

    int64_t original_nchans = shape[1];

    
    // TODO -- Padding in time domain then trim
    //      -- Later down the line overlap/add for extended readout
    //Always Pad in time because of non-periodicity in this domain

    if (m_pad_wire_domain) {
        // auto input_length = shape[0];
        auto response_length = base_frer_spectrum->shape()[0];

        auto padding_size = response_length - 1;

        if (m_use_fft_best_length) {
            auto suggested_length = fft_best_length(original_nchans + response_length - 1);
            log->debug("FFT best length suggested {}", suggested_length);
            padding_size = suggested_length - original_nchans;
        }

        //Pad to M+N-1 size 
        tensor_clone = torch::nn::functional::pad(
            tensor_clone,
            //The padding function option argument is the padding applied to the various sides of the tensor.
            //It's arranged last_dim_"left", last_dim_right, second_to_last_left, etc...
            //Later on, the FRER spectrum we get will be padded automatically on the wires dimension on the right
            //hence the order given here
            // Also we only give the additional size it must be padded with
            torch::nn::functional::PadFuncOptions({0, 0, 0, padding_size}).mode(torch::kConstant)

        );

        log->debug("Padded to {} ", tensor_clone.sizes()[1]);
        response_shape[0] = tensor_clone.sizes()[1];
    }



    //FFT on time dim
    tensor_clone = torch::fft::rfft(tensor_clone, std::nullopt, 2);

    //FFT on chan dim
    tensor_clone = torch::fft::fft(tensor_clone, std::nullopt, 1);

    //Get the Field x Elec. Response and do FFT in both dimensons
    auto frer_spectrum_tensor = base_frer_spectrum->spectrum(response_shape).clone();

    

    frer_spectrum_tensor = torch::fft::rfft2(frer_spectrum_tensor);

    //Get the wire shift
    int wire_shift = base_frer_spectrum->shifts()[0];
    log->debug("Preparing to shift by {} wires in", wire_shift);

    //Apply to input data
    if (!m_debug_no_frer)
        tensor_clone = tensor_clone / frer_spectrum_tensor;

    //Get the Wire filter -- already FFT'd
    //TODO -- fix the log here because of HfFilter weirdness
    auto wire_filter_tensor = base_wire_filter->spectrum({response_shape[0]});

    //Multiply along the wire dimension
    if (!m_debug_no_wire_filter)
        tensor_clone = tensor_clone * wire_filter_tensor.view({-1,1});

    //Inverse FFT in both dimensions
    tensor_clone = torch::fft::irfft2(tensor_clone);
    
    //Shift along time dimension
    int time_shift = (int) (
        (m_coarse_time_offset + base_frer_spectrum->shifts()[1]) /
        in->metadata()["period"].asDouble()
    );

    if (!m_debug_no_roll) {
        //Shift along wire dimension
        tensor_clone = tensor_clone.roll(wire_shift, 1);
        tensor_clone = tensor_clone.roll(time_shift, 2);
    }

    if (m_pad_wire_domain) {
        //This selects the post-roll block
        using namespace torch::indexing;
        // auto unpadded_width = shape[0] - (base_frer_spectrum->shape()[0]-1);
        // auto unpadded_width = shape[0] - (base_frer_spectrum->shape()[0]-1);
        tensor_clone = tensor_clone.index({
            Slice(), //Over all of the planes
            Slice(None, original_nchans) //0 to N orig channels
        });
        log->debug("Unpadded to {}", tensor_clone.sizes()[1]);
    }

    if (m_unsqueeze_input)
        torch::squeeze(tensor_clone, 0);
    
    Configuration set_md, tensor_md;

    tensor_md["tag"] = m_output_tensor_tag;
    metadata_passthrough(in->tensors()->at(0)->metadata(), tensor_md, m_passthrough);
    // log->debug("Passed channel_map\n{}", tensor_md["channel_map"]);
    set_md["tag"] = m_output_set_tag;
    auto output_tensor = std::make_shared<SimpleTorchTensor>(tensor_clone, tensor_md);
    std::vector<ITorchTensor::pointer> itv{
        output_tensor
    };
    out = std::make_shared<SimpleTorchTensorSet>(
        in->ident(), set_md,
        std::make_shared<std::vector<ITorchTensor::pointer>>(itv)
    );
    log->debug("Tagged output with Set:{} Tensor:{}",
               out->metadata()["tag"],
               out->tensors()->at(0)->metadata()["tag"]);

    return true;
}
