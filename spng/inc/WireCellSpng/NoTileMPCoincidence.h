#ifndef WIRECELL_SPNGNOTILEMPCOINCIDENCE
#define WIRECELL_SPNGNOTILEMPCOINCIDENCE

#include "WireCellAux/Logger.h"

#include "WireCellSpng/ITorchTensorSetFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellSpng/ITorchSpectrum.h"
#include "WireCellSpng/RayGrid.h"
#include "WireCellIface/IAnodePlane.h"


namespace WireCell {
namespace SPNG {
    class NoTileMPCoincidence : public Aux::Logger,
                    public WireCell::ITorchTensorSetFilter, public WireCell::IConfigurable {
    public:
        NoTileMPCoincidence( );
        virtual ~NoTileMPCoincidence();

        virtual bool operator()(const input_pointer& in, output_pointer& out);
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const {
            Configuration cfg;
            return cfg;
        };
    private:

        void convert_wires_to_channels(torch::Tensor & input, torch::Tensor & indices);

        int m_rebin_val{-1};
        int m_target_plane_index{0};
        int m_aux_plane_l_index{1};
        int m_aux_plane_m_index{2};
        bool m_debug_force_cpu{false};
        double m_readout_plane_width{100.},
               m_readout_plane_height{100.},
               m_pitch{5.},
               m_angle_in_radians{0.6230825}; //35.7deg
        torch::Device m_device{torch::kCPU};
        torch::Tensor m_trivial_blobs;
        torch::Tensor m_raygrid_views;
        std::string m_anode_tn{"AnodePlane"};
        IAnodePlane::pointer m_anode;
        int m_face_index{0};
        std::vector<std::unordered_map<int, std::vector<int>>> m_chan_index_to_wires{3};
        std::map<int, int> m_plane_nwires;
        std::map<int, torch::Tensor> m_plane_wires_to_channels;
        std::map<int, torch::Tensor> m_plane_channels_to_wires;
        std::string m_output_torch_name;
        bool m_debug_output{false};

        
    };
}
}

#endif