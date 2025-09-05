#include "WireCellSpng/NoTileMPCoincidence.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"
#include "WireCellSpng/ITorchSpectrum.h"
#include "WireCellSpng/RayTest.h"
#include "WireCellSpng/RayTiling.h"
#include <fstream>

using tensor_map = torch::Dict<std::string, torch::Tensor>;
using namespace torch::indexing;
WIRECELL_FACTORY(SPNGNoTileMPCoincidence, WireCell::SPNG::NoTileMPCoincidence,
                 WireCell::INamed,
                 WireCell::ITorchTensorSetFilter, WireCell::IConfigurable)

WireCell::SPNG::NoTileMPCoincidence::NoTileMPCoincidence()
  : Aux::Logger("SPNGNoTileMPCoincidence", "spng") {

}

WireCell::SPNG::NoTileMPCoincidence::~NoTileMPCoincidence() {};


void WireCell::SPNG::NoTileMPCoincidence::configure(const WireCell::Configuration& config) {
    m_rebin_val = get(config, "rebin_val", m_rebin_val);

    m_debug_force_cpu = get(config, "debug_force_cpu", m_debug_force_cpu);

    //Get the indices of the planes we're working with.
    //We apply MP2/MP3 finding to some target plane n
    //And we need planes l & m to determine those.
    m_face_index = get(config, "face_index", m_face_index);
    m_target_plane_index = get(config, "target_plane_index", m_target_plane_index);
    m_aux_plane_l_index = get(config, "aux_plane_l_index", m_aux_plane_l_index);
    m_aux_plane_m_index = get(config, "aux_plane_m_index", m_aux_plane_m_index);
    m_output_torch_name = get(config, "output_torch_name", m_output_torch_name);
    m_debug_output = get(config, "debug_output", m_debug_output);
    m_test_style = get(config, "test_style", m_test_style);
    //Check that we aren't requesting any of the same 2 planes
    if ((m_target_plane_index == m_aux_plane_l_index) ||
        (m_target_plane_index == m_aux_plane_m_index) ||
        (m_aux_plane_m_index == m_aux_plane_l_index)) {
        THROW(ValueError() <<
            errmsg{"Must request unqiue indices for the target and auxiliary planes. Provided:\n"} <<
            errmsg{String::format("\tTarget (n): %d\n", m_target_plane_index)} <<
            errmsg{String::format("\tAux (l): %d\n", m_aux_plane_l_index)} <<
            errmsg{String::format("\tAux (m): %d\n", m_aux_plane_m_index)}
        );
    }

    m_readout_plane_width = get(config, "readout_plane_width", m_readout_plane_width); //Unused
    m_readout_plane_height = get(config, "readout_plane_height", m_readout_plane_height); //Unused
    m_pitch = get(config, "pitch", m_pitch); //Unused
    m_angle_in_radians = get(config, "angle_in_radians", m_angle_in_radians); //Unused 

    //Get trivial blobs
    m_trivial_blobs = WireCell::SPNG::RayGrid::trivial_blobs(); //unused
    //Create the views & coordinates used in RayGrid


    m_anode_tn = get(config, "anode", m_anode_tn);
    m_anode = Factory::find_tn<IAnodePlane>(m_anode_tn);

    // std::cout << "Nfaces: " << m_anode->faces().size() << std::endl;

    // std::vector<torch::Tensor> m_pitch_tensors;
    // for (const auto& face : m_anode->faces()) {
    const auto & face = m_anode->face(m_face_index);
    m_raygrid_views = torch::zeros({5, 2, 2}/*, options*/);

    // std::cout << "Face: " << face->ident() << std::endl;
    const auto & coords = face->raygrid();
    const auto & centers = coords.centers();
    const auto & pitch_dirs = coords.pitch_dirs();
    const auto & pitch_mags = coords.pitch_mags();
    auto next_rays = centers;
    for (int ilayer = 0; ilayer < coords.nlayers(); ++ilayer)
        next_rays[ilayer] += pitch_dirs[ilayer]*pitch_mags[ilayer];

    //Construct views by hand. We need to do this based off of our target plane.
    //We create the 2 trivial blobs first, then 
    std::vector<int> layers = {0, 1, m_aux_plane_l_index+2, m_aux_plane_m_index+2, m_target_plane_index+2};
    int layer_count = 0;
    for (const auto & ilayer : layers) {
        // std::cout << "\tCenter: " << centers[ilayer] <<
        // " Pitch Dir (Mag): " <<
        // pitch_dirs[ilayer] << 
        // " (" << pitch_mags[ilayer] << ")" << std::endl;

        // next_rays[ilayer] += pitch_dirs[ilayer]*pitch_mags[ilayer];

        // std::cout << "\t\tNext ray: " << next_rays[ilayer] << std::endl;

        //Set the values in the tensor
        m_raygrid_views.index_put_({layer_count, 0, 0}, centers[ilayer][2]);
        m_raygrid_views.index_put_({layer_count, 0, 1}, centers[ilayer][1]);
        m_raygrid_views.index_put_({layer_count, 1, 0}, next_rays[ilayer][2]);
        m_raygrid_views.index_put_({layer_count, 1, 1}, next_rays[ilayer][1]);
        ++layer_count;
    }
    // }

    //Build up map to/from wires to channels for fast lookup
    for (const auto & plane : face->planes()) {
        if (face->ident() == m_face_index) {//Only do this once
            // std::cout << "Plane wires: " << plane->wires().size() << std::endl;
            m_plane_nwires[plane->ident()] = (int)plane->wires().size();
            m_plane_wires_to_channels[plane->ident()] = torch::zeros({(int)plane->wires().size()}, torch::TensorOptions().dtype(torch::kInt32));
        }
    }

    std::vector<int> plane_to_nchans(3);
    std::vector<size_t> max_nwires(3, 0);
    for (const auto & iface : m_anode->faces()) {
        //Hardcoding this until I figure out a better solution
        //Reset the collection plane
        plane_to_nchans[2] = 0;
        for (const auto & plane : iface->planes()) {

            // if (iface->ident() == m_face_index) {//Only do this once
            //     std::cout << "Plane wires: " << plane->wires().size() << std::endl;
            //     m_plane_nwires[plane->ident()] = plane->wires().size();
            //     m_plane_wires_to_channels[plane->ident()] = torch::zeros({plane->wires().size()}, torch::TensorOptions().dtype(torch::kInt32));
            // }
            auto wires_to_chans_accessor = m_plane_wires_to_channels[plane->ident()].accessor<int, 1>();

            auto & map = m_chan_index_to_wires[plane->ident()];
            // int ichan = 0;
            for (const auto & plane_chan : plane->channels()) {
                int & ichan = plane_to_nchans[plane->ident()];
                // std::cout << "Plane & Chan: " << plane->ident() << " " << ichan << " " << plane_chan->ident() << " wires" << std::endl;
                for (const auto & w : plane_chan->wires()) {
                    // std::cout << "\t" << w->index() << " " <<  w->planeid().face() << std::endl;
                    if (w->planeid().face() == face->ident()) { //Have to check against the target face
                        auto & temp_chans_to_wires = map[ichan];
                        // long these_nwires = static_cast<long>(temp_chans_to_wires.size());
                        
                        temp_chans_to_wires.push_back(w->index());
                        if (temp_chans_to_wires.size() > max_nwires[plane->ident()])
                            max_nwires[plane->ident()] = temp_chans_to_wires.size();
                        // if ((these_nwires+1) > max_nwires[plane->ident()]) {
                        //     max_nwires[plane->ident()] = these_nwires;
                        //     std::cout << "Increasing max_nwires " << plane->ident() << " " << these_nwires << std::endl;
                        // }
                        
                        wires_to_chans_accessor[w->index()] = ichan;
                    }
                }
                ++ichan;
            }
        }
    }

    for (size_t iplane = 0; iplane < plane_to_nchans.size(); ++iplane) {
        auto nchans = plane_to_nchans[iplane];
        // std::cout << iplane << " Making nchans/wires " << nchans << " " << max_nwires[iplane] << std::endl;
        auto & chan_to_wires_tensor = m_plane_channels_to_wires[iplane];
        chan_to_wires_tensor = torch::full(
            {nchans, static_cast<long>(max_nwires[iplane])},
            m_plane_nwires[iplane], //for default channel
            torch::TensorOptions().dtype(torch::kInt64));
            
        
        for (auto & [chan, wires] : m_chan_index_to_wires[iplane]) {
            // std::cout << "Chan: " << chan << " wires: " << wires.size() << std::endl;
            for (size_t iw = 0; iw < wires.size(); iw++) {
                // chan_to_wires_tensor[chan][iw] = wires[iw];
                chan_to_wires_tensor.index_put_({chan, static_cast<long>(iw)}, wires[iw]);
            }
        }

        // std::cout << "Mapped channels to wires\n" << chan_to_wires_tensor << std::endl;
    }

    
    // for (const auto & plane : face->planes()) {
    //     std::cout << "Wires to chans\n" << m_plane_wires_to_channels[plane->ident()] << std::endl;
    // }

}

void WireCell::SPNG::NoTileMPCoincidence::convert_wires_to_channels(
        torch::Tensor & input, torch::Tensor & indices) {

    //Add a zero which acts as the default value when accessed by the indices later
    input = torch::cat(
        torch::TensorList({
            input,
             torch::zeros({input.size(0), 1}, torch::TensorOptions(m_device))}),
        1//Wire dimension
    ).unsqueeze(1)//Add a dimension so we can broadcast within gather
    .expand({-1, indices.size(0), -1});//And repeat it according to the number of channels
    
    input = torch::gather(
        input,
        2,//New wire dimension
        indices.unsqueeze(0).expand({input.size(0), -1, -1})
    );

    input = torch::any(//Check if any wire segments according to this channel are active
        input,
        2
    )//.permute({1, 0})//Put channels first again 
    // .unsqueeze(0)//and add another 'batch' dimension
    .to(torch::kFloat64);//And make it double
}

torch::Tensor WireCell::SPNG::NoTileMPCoincidence::make_indices(torch::Tensor & frame) {
    auto reshaped = frame.reshape({-1, frame.size(-1)});
    auto counts = reshaped.sum(-1);
    auto max_size = torch::max(counts).item<long>();
    auto results = torch::full({reshaped.size(0), max_size}, -1);

    auto col_indices = reshaped.nonzero().index({"...", -1});
    int64_t current_offset = 0;
    for (int64_t i = 0; i < reshaped.size(0); ++i) {
        auto count = counts.index({i}).item<int64_t>();
        // std::cout << i << " " << count << std::endl;
        if (count > 0) {
            // Get the indices for the current row
            auto current_col_indices = col_indices.slice(0, current_offset, current_offset + count);

            // Fill the corresponding slice in the result tensor
            results.index_put_({i, torch::indexing::Slice(0, count)}, current_col_indices.to(torch::kInt));
        }
        current_offset += count;
    }

    return results;
}

bool WireCell::SPNG::NoTileMPCoincidence::operator()(const input_pointer& in, output_pointer& out) {
    out = nullptr;
    if (!in) {
        log->debug("EOS ");
        return true;
    }
    log->debug("Running NoTileMPCoincidence");

    m_device = ((
        (torch::cuda::is_available() && !m_debug_force_cpu) ? torch::kCUDA : torch::kCPU
    ));

    m_raygrid_views = m_raygrid_views.to(m_device);

    // std::cout << m_raygrid_views[0] << std::endl;

    tensor_map to_save;

    WireCell::SPNG::RayGrid::Coordinates m_raygrid_coords =
            WireCell::SPNG::RayGrid::Coordinates(m_raygrid_views);
    m_raygrid_coords.to(m_device);

    m_plane_channels_to_wires[m_target_plane_index] = m_plane_channels_to_wires[m_target_plane_index].to(m_device);

    //Clone the inputs  
    auto target_tensor_n = (*in->tensors())[m_target_plane_index]->tensor().clone().to(m_device);
    // target_tensor_n.index_put_({0, Slice(), Slice()}, 0.);
    // target_tensor_n.index_put_({0, 317, 4887}, 1.);
    // {
    //     auto name = "target_" + m_output_torch_name;
    //     std::cerr << "writing " << name << "\n";
    //     std::ofstream output_file(name, std::ios::out | std::ios::binary);
    //     auto data = torch::pickle_save(target_tensor_n.to(torch::kCPU));
    //     output_file.write(data.data(), data.size());
    //     output_file.close();
    // }


    auto aux_tensor_l = (*in->tensors())[m_aux_plane_l_index]->tensor().clone().to(m_device);
    // aux_tensor_l.index_put_({0, Slice(), Slice()}, 0.);
    // aux_tensor_l.index_put_({0, 104, 4887}, 1.);
    // {
    //     auto name = "aux_tensor_l_" + m_output_torch_name;
    //     std::cerr << "writing " << name << "\n";
    //     std::ofstream output_file(name, std::ios::out | std::ios::binary);
    //     auto data = torch::pickle_save(aux_tensor_l.to(torch::kCPU));
    //     output_file.write(data.data(), data.size());
    //     output_file.close();
    // }

    auto aux_tensor_m = (*in->tensors())[m_aux_plane_m_index]->tensor().clone().to(m_device);
    // aux_tensor_m.index_put_({0, Slice(), Slice()}, 0.);
    // aux_tensor_m.index_put_({0, 544, 4887}, 1.);
    // {
    //     auto name = "aux_tensor_m_" + m_output_torch_name;
    //     std::cerr << "writing " << name << "\n";
    //     std::ofstream output_file(name, std::ios::out | std::ios::binary);
    //     auto data = torch::pickle_save(aux_tensor_m.to(torch::kCPU));
    //     output_file.write(data.data(), data.size());
    //     output_file.close();
    // }

    //Transform into bool tensors (activities)
    // auto tester = torch::zeros({1}).to(m_device);
    torch::nn::MaxPool1d pool(torch::nn::MaxPool1dOptions(4));
    aux_tensor_l = (pool(aux_tensor_l) > 0);
    aux_tensor_m = (pool(aux_tensor_m) > 0);

    target_tensor_n = pool(target_tensor_n);

    if (m_debug_output) {
        to_save.insert("aux_tensor_l", aux_tensor_l);
        to_save.insert("aux_tensor_m", aux_tensor_m);
        to_save.insert("target_tensor_n", target_tensor_n);
    }

    auto nbatch = target_tensor_n.size(0);
    auto nsamples = target_tensor_n.size(-1);

    log->debug("Running on nbatches:{} nsamples:{}", nbatch, nsamples);

    torch::Tensor output_tensor_active = torch::zeros({
        nbatch,
        m_plane_nwires[m_target_plane_index],
        nsamples},
        torch::TensorOptions(m_device).dtype(torch::kFloat64));
    torch::Tensor output_tensor_inactive = torch::zeros_like(output_tensor_active, torch::TensorOptions(m_device).dtype(torch::kFloat64));

    
    auto l_rows = aux_tensor_l.index({"...", m_plane_wires_to_channels[m_aux_plane_l_index], torch::indexing::Slice()});
    auto m_rows = aux_tensor_m.index({"...", m_plane_wires_to_channels[m_aux_plane_m_index], torch::indexing::Slice()});
    auto target_rows = target_tensor_n.index({"...", m_plane_wires_to_channels[m_target_plane_index], torch::indexing::Slice()});
    
    l_rows = l_rows.transpose(-1, -2);
    l_rows = l_rows.reshape({-1, l_rows.size(-1)});
    
    m_rows = m_rows.transpose(-1, -2);
    m_rows = m_rows.reshape({-1, m_rows.size(-1)});
    
    target_rows = target_rows.transpose(-1, -2);
    target_rows = target_rows.reshape({-1, target_rows.size(-1)});
    
    output_tensor_active = output_tensor_active.transpose(-1, -2).reshape({-1, m_plane_nwires[m_target_plane_index]});
    output_tensor_inactive = output_tensor_inactive.transpose(-1, -2).reshape({-1, m_plane_nwires[m_target_plane_index]});

    for (long int irow = 0; irow < l_rows.size(-2); ++irow) {

        auto l_row = l_rows.index({irow});
        auto m_row = m_rows.index({irow});
        auto l_hi = torch::where(l_row > 0)[0];
        auto m_hi = torch::where(m_row > 0)[0];

        //Make Nl x Nm pairs of indices of the active wires (rays) in each plane
        auto cross = torch::zeros({l_hi.size(0), m_hi.size(0), 2}, torch::TensorOptions(m_device).dtype(l_hi.dtype()));
        cross.index_put_({"...", 1}, m_hi);
        cross = cross.permute({1, 0, 2});
        cross.index_put_({"...", 0}, l_hi);
        cross = cross.reshape({-1, 2});

        //Get the ray indices for each plane/view
        auto r1 = cross.index({Slice(), 0});
        auto r2 = cross.index({Slice(), 1});

        //Make view tensors in shapes of the ray pairs
        auto view1 = torch::full_like(r1, 2, torch::TensorOptions(m_device).dtype(l_hi.dtype()));
        auto view2 = torch::full_like(r1, 3, torch::TensorOptions(m_device).dtype(l_hi.dtype()));
        auto view3 = torch::full_like(r1, 4, torch::TensorOptions(m_device).dtype(l_hi.dtype()));
        //Also just a scalar
        auto view3_short = torch::tensor(4, torch::TensorOptions(m_device).dtype(l_hi.dtype()));

        //Get the locations of crossing wires within view3
        auto view3_locs = m_raygrid_coords.pitch_location(view1, r1, view2, r2, view3);

        //Convert the locations to pitch indices
        auto results = m_raygrid_coords.pitch_index(
            view3_locs,
            view3_short
        );
        
        if (m_debug_output) {
            to_save.insert("index" + std::to_string(irow), results.to(torch::kCPU));
            to_save.insert("locs" + std::to_string(irow), view3_locs.to(torch::kCPU));
            auto view3_zy = torch::zeros({view3_locs.size(0), 2}, torch::TensorOptions(m_device).dtype(view3_locs.dtype()));
            view3_zy.index_put_({Slice(), 0}, view3_locs*m_raygrid_coords.pitch_dir.index({4, 0}));
            view3_zy.index_put_({Slice(), 1}, view3_locs*m_raygrid_coords.pitch_dir.index({4, 1}));
            to_save.insert("locs_zy_" + std::to_string(irow), view3_zy.to(torch::kCPU));
        }

        //Make sure the returned indices are within the active volume
        results = std::get<0>(at::_unique(results.index(
            {torch::where((results >= 0) & (results < m_plane_nwires[m_target_plane_index]))[0]}
        )));

        //Get the corresponding statuses in the target plane
        // auto target_row = target_rows.index({results, irow});
        auto target_row = target_rows.index({irow, results});

        //MP3 means a high wire in the target plane
        //overlaps with a crossing pair in the other two planes
        auto mp3_indices = torch::nonzero(target_row);
        // output_tensor_active.index_put_({0, results.index({mp3_indices}), irow}, 1.);
        output_tensor_active.index_put_({irow, results.index({mp3_indices})}, 1.);
        
        //MP2 means a low wire in the target plane
        //overlaps with a crossing pair in the other two planes, is low
        auto mp2_indices = torch::nonzero(target_row == 0);
        output_tensor_inactive.index_put_({irow, results.index({mp2_indices})}, 1.);

    }

    convert_wires_to_channels(output_tensor_active, m_plane_channels_to_wires[m_target_plane_index]);
    output_tensor_active = output_tensor_active.reshape({nbatch, nsamples, -1}).transpose(-1, -2);
    convert_wires_to_channels(output_tensor_inactive, m_plane_channels_to_wires[m_target_plane_index]);
    output_tensor_inactive = output_tensor_inactive.reshape({nbatch, nsamples, -1}).transpose(-1, -2);


    // TODO: set md?
    Configuration set_md, mp2_md, mp3_md;
    set_md["tag"] = "";//m_output_set_tag;
    mp2_md["tag"] = "mp2";
    mp3_md["tag"] = "mp3";

    std::vector<ITorchTensor::pointer> itv{
        std::make_shared<SimpleTorchTensor>(output_tensor_active.clone(), mp3_md),
        std::make_shared<SimpleTorchTensor>(output_tensor_inactive.clone(), mp2_md),
    };
    out = std::make_shared<SimpleTorchTensorSet>(
        in->ident(), set_md,
        std::make_shared<std::vector<ITorchTensor::pointer>>(itv)
    );

    // {//Writing l_rows
    //     auto name = "rows_l_" + m_output_torch_name;
    //     std::cerr << "writing " << name << "\n";
    //     std::ofstream output_file(name, std::ios::out | std::ios::binary);
    //     auto data = torch::pickle_save(l_rows.to(torch::kCPU));
    //     output_file.write(data.data(), data.size());
    //     output_file.close();
    // }
    // {//Writing m_rows
    //     auto name = "rows_m_" + m_output_torch_name;
    //     std::cerr << "writing " << name << "\n";
    //     std::ofstream output_file(name, std::ios::out | std::ios::binary);
    //     auto data = torch::pickle_save(m_rows.to(torch::kCPU));
    //     output_file.write(data.data(), data.size());
    //     output_file.close();
    // }
    // {//Writing target_rows
    //     auto name = "rows_target_" + m_output_torch_name;
    //     std::cerr << "writing " << name << "\n";
    //     std::ofstream output_file(name, std::ios::out | std::ios::binary);
    //     auto data = torch::pickle_save(target_rows.to(torch::kCPU));
    //     output_file.write(data.data(), data.size());
    //     output_file.close();
    // }

    // {//Writing views
    //     auto name = "views_" + m_output_torch_name;
    //     std::cerr << "writing " << name << "\n";
    //     std::ofstream output_file(name, std::ios::out | std::ios::binary);
    //     auto data = torch::pickle_save(m_raygrid_views.to(torch::kCPU));
    //     output_file.write(data.data(), data.size());
    //     output_file.close();
    // }

    if (m_debug_output) {
        to_save.insert("coords", m_raygrid_views.to(torch::kCPU));
        to_save.insert("target_rows", target_rows.to(torch::kCPU));
        to_save.insert("m_rows", m_rows.to(torch::kCPU));
        to_save.insert("l_rows", l_rows.to(torch::kCPU));


        std::ofstream output_file(m_output_torch_name, std::ios::out | std::ios::binary);
        auto data = torch::pickle_save(to_save);
        output_file.write(data.data(), data.size());
        output_file.close();
    }
    return true;
}
