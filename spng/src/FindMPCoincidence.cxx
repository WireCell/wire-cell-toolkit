#include "WireCellSpng/FindMPCoincidence.h"
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
WIRECELL_FACTORY(SPNGFindMPCoincidence, WireCell::SPNG::FindMPCoincidence,
                 WireCell::INamed,
                 WireCell::ITorchTensorSetFilter, WireCell::IConfigurable)

WireCell::SPNG::FindMPCoincidence::FindMPCoincidence()
  : Aux::Logger("SPNGFindMPCoincidence", "spng") {

}

WireCell::SPNG::FindMPCoincidence::~FindMPCoincidence() {};


void WireCell::SPNG::FindMPCoincidence::configure(const WireCell::Configuration& config) {
    m_rebin_val = get(config, "rebin_val", m_rebin_val);

    m_debug_force_cpu = get(config, "debug_force_cpu", m_debug_force_cpu);

    //Get the indices of the planes we're working with.
    //We apply MP2/MP3 finding to some target plane n
    //And we need planes l & m to determine those.
    m_target_plane_index = get(config, "target_plane_index", m_target_plane_index);
    m_aux_plane_l_index = get(config, "aux_plane_l_index", m_aux_plane_l_index);
    m_aux_plane_m_index = get(config, "aux_plane_m_index", m_aux_plane_m_index);
    m_output_torch_name = get(config, "output_torch_name", m_output_torch_name);
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

    m_readout_plane_width = get(config, "readout_plane_width", m_readout_plane_width);
    m_readout_plane_height = get(config, "readout_plane_height", m_readout_plane_height);
    m_pitch = get(config, "pitch", m_pitch);
    m_angle_in_radians = get(config, "angle_in_radians", m_angle_in_radians);

    //Get trivial blobs
    m_trivial_blobs = WireCell::SPNG::RayGrid::trivial_blobs();
    //Create the views & coordinates used in RayGrid


    m_anode_tn = get(config, "anode", m_anode_tn);
    m_anode = Factory::find_tn<IAnodePlane>(m_anode_tn);

    std::cout << "Nfaces: " << m_anode->faces().size() << std::endl;

    // std::vector<torch::Tensor> m_pitch_tensors;
    // for (const auto& face : m_anode->faces()) {
    const auto & face = m_anode->face(m_face_index);
    m_raygrid_views = torch::zeros({5, 2, 2}/*, options*/);

    std::cout << "Face: " << face->ident() << std::endl;
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
        std::cout << "\tCenter: " << centers[ilayer] <<
        " Pitch Dir (Mag): " <<
        pitch_dirs[ilayer] << 
        " (" << pitch_mags[ilayer] << ")" << std::endl;

        // next_rays[ilayer] += pitch_dirs[ilayer]*pitch_mags[ilayer];

        std::cout << "\t\tNext ray: " << next_rays[ilayer] << std::endl;

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
            std::cout << "Plane wires: " << plane->wires().size() << std::endl;
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
                        wires_to_chans_accessor[w->index()] = ichan;
                    }
                }
                ++ichan;
            }
        }
    }

    for (size_t iplane = 0; iplane < plane_to_nchans.size(); ++iplane) {
        auto nchans = plane_to_nchans[iplane];
        auto & chan_to_wires_tensor = m_plane_channels_to_wires[iplane];
        chan_to_wires_tensor = torch::full(
            {nchans, static_cast<long>(max_nwires[iplane])},
            m_plane_nwires[iplane], //for default channel
            torch::TensorOptions().dtype(torch::kInt64));
            
        
        for (auto & [chan, wires] : m_chan_index_to_wires[iplane]) {
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

void WireCell::SPNG::FindMPCoincidence::convert_wires_to_channels(
        torch::Tensor & input, torch::Tensor & indices) {
    input = torch::cat(
        torch::TensorList({
            input.squeeze(0),
             torch::zeros({1, input.size(-1)}, torch::TensorOptions(m_device))}),
        0//Wire dimension
    ).permute({1, 0})//Put wire dimension last
    .unsqueeze(1)//Add a dimension so we can broadcast within gather
    .expand({-1, indices.size(0), -1});//And repeat it according to the number of channels
    
    input = torch::gather(
        input,
        2,//New wire dimension
        indices.unsqueeze(0).expand({input.size(0), -1, -1})
    );

    input = torch::any(//Check if any wire segments according to this channel are active
        input,
        2
    ).permute({1, 0})//Put channels first again 
    .unsqueeze(0)//and add another 'batch' dimension
    .to(torch::kFloat64);//And make it double
}

bool WireCell::SPNG::FindMPCoincidence::operator()(const input_pointer& in, output_pointer& out) {
    out = nullptr;
    if (!in) {
        log->debug("EOS ");
        return true;
    }
    log->debug("Running FindMPCoincidence");

    m_device = ((
        (torch::cuda::is_available() && !m_debug_force_cpu) ? torch::kCUDA : torch::kCPU
    ));
    m_trivial_blobs = m_trivial_blobs.to(m_device);

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
    auto target_tensor_map = (*in->tensors())[m_target_plane_index]->metadata()["channel_map"];
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
    auto aux_tensor_l_map = (*in->tensors())[m_aux_plane_l_index]->metadata()["channel_map"];
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
    auto aux_tensor_m_map = (*in->tensors())[m_aux_plane_m_index]->metadata()["channel_map"];
    // {
    //     auto name = "aux_tensor_m_" + m_output_torch_name;
    //     std::cerr << "writing " << name << "\n";
    //     std::ofstream output_file(name, std::ios::out | std::ios::binary);
    //     auto data = torch::pickle_save(aux_tensor_m.to(torch::kCPU));
    //     output_file.write(data.data(), data.size());
    //     output_file.close();
    // }

    //Transform into bool tensors (activities)
    auto tester = torch::zeros({1}).to(m_device);
    torch::nn::MaxPool1d pool(torch::nn::MaxPool1dOptions(4));
    aux_tensor_l = (pool(aux_tensor_l) > tester);
    aux_tensor_m = (pool(aux_tensor_m) > tester);

    target_tensor_n = pool(target_tensor_n);

    torch::Tensor output_tensor_active = torch::zeros({
        target_tensor_n.size(0),
        m_plane_nwires[m_target_plane_index],
        target_tensor_n.size(-1)},
        torch::TensorOptions(m_device).dtype(torch::kFloat64));
    torch::Tensor output_tensor_inactive = torch::zeros_like(output_tensor_active, torch::TensorOptions(m_device).dtype(torch::kFloat64));


    auto l_rows = aux_tensor_l.index({0, m_plane_wires_to_channels[m_aux_plane_l_index], torch::indexing::Slice()});
    auto m_rows = aux_tensor_m.index({0, m_plane_wires_to_channels[m_aux_plane_m_index], torch::indexing::Slice()});
    auto target_rows = target_tensor_n.index({0, m_plane_wires_to_channels[m_target_plane_index], torch::indexing::Slice()});

    // torch::Tensor output_l_wires = torch::zeros({
    //     aux_tensor_l.size(0),
    //     m_plane_nwires[m_aux_plane_l_index],
    //     aux_tensor_l.size(-1)},
    //     torch::TensorOptions(m_device).dtype(torch::kFloat64));
    // torch::Tensor output_m_wires = torch::zeros({
    //     aux_tensor_m.size(0),
    //     m_plane_nwires[m_aux_plane_m_index],
    //     aux_tensor_m.size(-1)},
    //     torch::TensorOptions(m_device).dtype(torch::kFloat64));

    torch::Tensor output_l_blobs = torch::zeros({
        aux_tensor_l.size(0),
        m_plane_nwires[m_aux_plane_l_index],
        aux_tensor_l.size(-1)},
        torch::TensorOptions(m_device).dtype(torch::kFloat64));
    torch::Tensor output_m_blobs = torch::zeros({
        aux_tensor_m.size(0),
        m_plane_nwires[m_aux_plane_m_index],
        aux_tensor_m.size(-1)},
        torch::TensorOptions(m_device).dtype(torch::kFloat64));

    //Apply the first two 'real' layers -- the order doesn't matter?
    // auto coords = m_raygrid_coords[0]; // For testing -- just one side of the APA

    // std::cout << "Trivial Blobs" << m_trivial_blobs << std::endl;
    for (long int irow = 0; irow < aux_tensor_l.sizes().back(); ++irow) {

        auto l_row = l_rows.index({Slice(), irow});
        auto m_row = m_rows.index({Slice(), irow});
        auto target_row = target_rows.index({Slice(), irow});

        // auto raygrid_row_l = element_tensor_l.to(m_device);
        auto blobs = WireCell::SPNG::RayGrid::apply_activity(m_raygrid_coords, m_trivial_blobs, l_row);

        // std::cout << "First layer done" << std::endl;
        // std::cout << blobs.sizes() << std::endl;
        if (blobs.size(0) == 0) {
            // std::cout << "Found no blobs. Moving on" << std::endl;
            continue;
        }
        // for (int iblob = 0; iblob < blobs.size(0); ++iblob) {
        //     // std::cout << mp3_accessor[iblob][4][0] << " " << mp3_accessor[iblob][4][1] << std::endl;
        //     output_l_blobs.index_put_(
        //         {
        //             0,
        //             torch::indexing::Slice(blobs.index({iblob, -1, 0}).item<long>(), blobs.index({iblob, -1, 1}).item<long>()),
        //             torch::indexing::Slice(irow, (irow+1))
        //         },
        //         1.
        //     );
        // }
        // auto raygrid_row_m = element_tensor_m.to(m_device);
        blobs = WireCell::SPNG::RayGrid::apply_activity(m_raygrid_coords, blobs, m_row/*raygrid_row_m*/);
        // {//Writing blobs with l & m
        //     std::cerr << "writing " << m_output_torch_name << "\n";
        //     std::ofstream output_file("blobs_m_" + m_output_torch_name, std::ios::out | std::ios::binary);
        //     auto data = torch::pickle_save(blobs.to(torch::kCPU));
        //     output_file.write(data.data(), data.size());
        //     output_file.close();
        // }
        // // std::cout << "Second layer done" << std::endl;
        if (blobs.size(0) == 0) {
            // std::cout << "Found no blobs. Moving on" << std::endl;
            continue;
        }

        // for (int iblob = 0; iblob < blobs.size(0); ++iblob) {
        //     // std::cout << mp3_accessor[iblob][4][0] << " " << mp3_accessor[iblob][4][1] << std::endl;
        //     output_m_blobs.index_put_(
        //         {
        //             0,
        //             torch::indexing::Slice(blobs.index({iblob, -1, 0}).item<long>(), blobs.index({iblob, -1, 1}).item<long>()),
        //             torch::indexing::Slice(irow, (irow+1))
        //         },
        //         1.
        //     );
        // }

        //For the last layer, get the bounds of the would-be created blobs
        //MP3 means our target plane has activity overlapping with blobs from the first 2 layers
        tester = tester.to(m_device);
        auto target_active = (target_row/*raygrid_row_n*/ > tester);
        if (target_active.any().item<bool>()) {
            auto mp3_blobs = WireCell::SPNG::RayGrid::apply_activity(
                m_raygrid_coords, blobs, target_active
            );

            // {//Writing blobs with l & m
            //     std::cerr << "writing " << m_output_torch_name << "\n";
            //     std::ofstream output_file("blobs_mp3_" + m_output_torch_name, std::ios::out | std::ios::binary);
            //     auto data = torch::pickle_save(mp3_blobs.to(torch::kCPU));
            //     output_file.write(data.data(), data.size());
            //     output_file.close();
            // }

            if (mp3_blobs.size(0) > 0) {
                std::vector<torch::Tensor> indices;
                for (int iblob = 0; iblob < mp3_blobs.size(0); ++iblob) {
                    indices.push_back(
                        torch::arange(mp3_blobs.index({iblob, -1, 0}).item<long>(), mp3_blobs.index({iblob, -1, 1}).item<long>())
                    );
                }
                auto indices_tensor = std::get<0>(at::_unique(torch::cat(indices)));
                output_tensor_active.index_put_({0, indices_tensor, irow}, 1.);
            }
        }
        
        //MP2 means our target plane does not have activity overlapping with blobs from the first 2 layers
        auto target_inactive = (target_row == tester);
        if (target_inactive.any().item<bool>()) {
            auto mp2_blobs = WireCell::SPNG::RayGrid::apply_activity(
                m_raygrid_coords, blobs, target_inactive
            );
            
            if (mp2_blobs.size(0) > 0) {
                std::vector<torch::Tensor> indices;
                for (int iblob = 0; iblob < mp2_blobs.size(0); ++iblob) {
                    indices.push_back(
                        torch::arange(mp2_blobs.index({iblob, -1, 0}).item<long>(), mp2_blobs.index({iblob, -1, 1}).item<long>())
                    );
                }
                auto indices_tensor = std::get<0>(at::_unique(torch::cat(indices)));
                output_tensor_inactive.index_put_({0, indices_tensor, irow}, 1.);
            }
        }
    }

    //Have to transform back into channel basis
    convert_wires_to_channels(output_tensor_active, m_plane_channels_to_wires[m_target_plane_index]);
    convert_wires_to_channels(output_tensor_inactive, m_plane_channels_to_wires[m_target_plane_index]);


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
    // {//Writing blobs with l
    //     auto name = "blobs_l_" + m_output_torch_name;
    //     std::cerr << "writing " << name << "\n";
    //     std::ofstream output_file(name, std::ios::out | std::ios::binary);
    //     auto data = torch::pickle_save(output_l_blobs.to(torch::kCPU));
    //     output_file.write(data.data(), data.size());
    //     output_file.close();
    // }
    // {//Writing blobs with m
    //     auto name = "blobs_m_" + m_output_torch_name;
    //     std::cerr << "writing " << name << "\n";
    //     std::ofstream output_file(name, std::ios::out | std::ios::binary);
    //     auto data = torch::pickle_save(output_m_blobs.to(torch::kCPU));
    //     output_file.write(data.data(), data.size());
    //     output_file.close();
    // }

    // std::ofstream output_file(m_output_torch_name, std::ios::out | std::ios::binary);
    // to_save.insert("coords", m_raygrid_views.to(torch::kCPU));
    // auto data = torch::pickle_save(to_save);
    // output_file.write(data.data(), data.size());
    // output_file.close();
    return true;
}
