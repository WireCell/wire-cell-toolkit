#ifndef WIRECELL_CLUS_NEUTRINOTAGGERINFO_H
#define WIRECELL_CLUS_NEUTRINOTAGGERINFO_H

/// Data structures produced by the NeutrinoID tagger chain.
///
/// KineInfo  — kinematic output (particle energies, pi0 info, reco Enu).
///             Mirrors WCPPID::NeutrinoID::KineInfo.
///             Filled by PatternAlgorithms::fill_kine_tree().
///
/// TaggerInfo — BDT input features for all sub-taggers (cosmic, numu, nue,
///              singlephoton, ssm).  Mirrors WCPPID::NeutrinoID::TaggerInfo.
///              Reset to defaults by PatternAlgorithms::init_tagger_info().
///
/// Both structs use C++ default-member-initializers so that value-
/// initialisation (e.g. KineInfo{}) yields the same defaults as the
/// prototype's init_tagger_info() call.

#include <vector>

namespace WireCell::Clus::PR {

    // ------------------------------------------------------------------ //
    //  KineInfo                                                            //
    // ------------------------------------------------------------------ //
    struct KineInfo {
        // neutrino vertex (SCE-corrected when geom_helper is available)
        float kine_nu_x_corr{0};
        float kine_nu_y_corr{0};
        float kine_nu_z_corr{0};

        // reconstructed neutrino energy
        float kine_reco_Enu{0};       // sum of particle kinetic energies + add_energy
        float kine_reco_add_energy{0}; // rest masses + binding energies

        // per-particle arrays (one entry per track/shower from vertex)
        std::vector<float> kine_energy_particle;   // kinetic energy [MeV]
        std::vector<int>   kine_energy_info;       // 0=dQdx, 1=range, 2=charge
        std::vector<int>   kine_particle_type;     // PDG code
        std::vector<int>   kine_energy_included;   // 1 = included in Enu sum

        // pi0 kinematics (filled by id_pi0_with/without_vertex)
        float kine_pio_mass{0};
        int   kine_pio_flag{0};    // 0=not found, 1=with vertex, 2=without vertex
        float kine_pio_vtx_dis{0};

        float kine_pio_energy_1{0};
        float kine_pio_theta_1{0};
        float kine_pio_phi_1{0};
        float kine_pio_dis_1{0};

        float kine_pio_energy_2{0};
        float kine_pio_theta_2{0};
        float kine_pio_phi_2{0};
        float kine_pio_dis_2{0};

        float kine_pio_angle{0};
    };


    // ------------------------------------------------------------------ //
    //  TaggerInfo                                                          //
    // ------------------------------------------------------------------ //
    /// BDT input features for all NeutrinoID sub-taggers.
    ///
    /// Fields with a non-zero default carry an explicit in-class initializer
    /// matching the prototype's init_tagger_info().  All other fields
    /// default-initialize to zero (value-init: TaggerInfo{}).
    struct TaggerInfo {

        // ---- cosmic tagger (top-level flag) ----------------------------- //
        float cosmic_flag{1};
        float cosmic_filled{0};
        float cosmic_n_solid_tracks{0};
        float cosmic_energy_main_showers{0};
        float cosmic_energy_indirect_showers{0};
        float cosmic_energy_direct_showers{0};
        float cosmic_n_direct_showers{0};
        float cosmic_n_indirect_showers{0};
        float cosmic_n_main_showers{0};

        // ---- shower gap identification ---------------------------------- //
        float gap_flag{1};
        float gap_filled{0};
        float gap_n_bad{0};
        float gap_n_points{0};
        float gap_energy{0};
        float gap_flag_single_shower{0};
        float gap_flag_parallel{0};
        float gap_flag_prolong_u{0};
        float gap_flag_prolong_v{0};
        float gap_flag_prolong_w{0};
        float gap_num_valid_tracks{0};

        // ---- MIP quality ----------------------------------------------- //
        float mip_quality_flag{1};
        float mip_quality_filled{0};
        float mip_quality_energy{0};
        float mip_quality_overlap{0};
        float mip_quality_n_showers{0};
        float mip_quality_n_tracks{0};
        float mip_quality_flag_inside_pi0{0};
        float mip_quality_n_pi0_showers{0};
        float mip_quality_shortest_length{0};
        float mip_quality_acc_length{0};
        float mip_quality_shortest_angle{0};
        float mip_quality_flag_proton{0};

        // ---- MIP identification ---------------------------------------- //
        float mip_flag{1};
        float mip_filled{0};
        float mip_energy{0};
        float mip_n_end_reduction{0};
        float mip_n_first_mip{0};
        float mip_n_first_non_mip{19};
        float mip_n_first_non_mip_1{19};
        float mip_n_first_non_mip_2{19};
        float mip_vec_dQ_dx_0{1};
        float mip_vec_dQ_dx_1{1};
        float mip_max_dQ_dx_sample{1};
        float mip_n_below_threshold{19};
        float mip_n_below_zero{0};
        float mip_n_lowest{1};
        float mip_n_highest{1};
        float mip_lowest_dQ_dx{1};
        float mip_highest_dQ_dx{1};
        float mip_medium_dQ_dx{1};
        float mip_stem_length{1};
        float mip_length_main{1};
        float mip_length_total{1};
        float mip_angle_beam{0};
        float mip_iso_angle{0};
        float mip_n_vertex{1};
        float mip_n_good_tracks{0};
        float mip_E_indirect_max_energy{0};
        float mip_flag_all_above{0};
        float mip_min_dQ_dx_5{1};
        float mip_n_other_vertex{2};
        float mip_n_stem_size{20};
        float mip_flag_stem_trajectory{0};
        float mip_min_dis{0};

        // extra dQ/dx samples
        float mip_vec_dQ_dx_2{0};
        float mip_vec_dQ_dx_3{0};
        float mip_vec_dQ_dx_4{0};
        float mip_vec_dQ_dx_5{0};
        float mip_vec_dQ_dx_6{0};
        float mip_vec_dQ_dx_7{0};
        float mip_vec_dQ_dx_8{0};
        float mip_vec_dQ_dx_9{0};
        float mip_vec_dQ_dx_10{0};
        float mip_vec_dQ_dx_11{0};
        float mip_vec_dQ_dx_12{0};
        float mip_vec_dQ_dx_13{0};
        float mip_vec_dQ_dx_14{0};
        float mip_vec_dQ_dx_15{0};
        float mip_vec_dQ_dx_16{0};
        float mip_vec_dQ_dx_17{0};
        float mip_vec_dQ_dx_18{0};
        float mip_vec_dQ_dx_19{0};

        // ---- SSM (short straight muon / KDAR) tagger ------------------- //
        float ssm_flag_st_kdar{0};
        float ssm_Nsm{-999};
        float ssm_Nsm_wivtx{-999};

        float ssm_dq_dx_fwd_1{-999};
        float ssm_dq_dx_fwd_2{-999};
        float ssm_dq_dx_fwd_3{-999};
        float ssm_dq_dx_fwd_4{-999};
        float ssm_dq_dx_fwd_5{-999};
        float ssm_dq_dx_bck_1{-999};
        float ssm_dq_dx_bck_2{-999};
        float ssm_dq_dx_bck_3{-999};
        float ssm_dq_dx_bck_4{-999};
        float ssm_dq_dx_bck_5{-999};
        float ssm_d_dq_dx_fwd_12{-999};
        float ssm_d_dq_dx_fwd_23{-999};
        float ssm_d_dq_dx_fwd_34{-999};
        float ssm_d_dq_dx_fwd_45{-999};
        float ssm_d_dq_dx_bck_12{-999};
        float ssm_d_dq_dx_bck_23{-999};
        float ssm_d_dq_dx_bck_34{-999};
        float ssm_d_dq_dx_bck_45{-999};
        float ssm_max_dq_dx_fwd_3{-999};
        float ssm_max_dq_dx_fwd_5{-999};
        float ssm_max_dq_dx_bck_3{-999};
        float ssm_max_dq_dx_bck_5{-999};
        float ssm_max_d_dq_dx_fwd_3{-999};
        float ssm_max_d_dq_dx_fwd_5{-999};
        float ssm_max_d_dq_dx_bck_3{-999};
        float ssm_max_d_dq_dx_bck_5{-999};
        float ssm_medium_dq_dx{-999};
        float ssm_medium_dq_dx_bp{-999};
        float ssm_angle_to_z{-999};
        float ssm_angle_to_target{-999};
        float ssm_angle_to_absorber{-999};
        float ssm_angle_to_vertical{-999};
        float ssm_x_dir{-999};
        float ssm_y_dir{-999};
        float ssm_z_dir{-999};
        float ssm_kine_energy{-999};
        float ssm_kine_energy_reduced{-999};
        float ssm_vtx_activity{-999};
        float ssm_pdg{-999};
        float ssm_dQ_dx_cut{-999};
        float ssm_score_mu_fwd{-999};
        float ssm_score_p_fwd{-999};
        float ssm_score_e_fwd{-999};
        float ssm_score_mu_bck{-999};
        float ssm_score_p_bck{-999};
        float ssm_score_e_bck{-999};
        float ssm_score_mu_fwd_bp{-999};
        float ssm_score_p_fwd_bp{-999};
        float ssm_score_e_fwd_bp{-999};
        float ssm_length{-999};
        float ssm_direct_length{-999};
        float ssm_length_ratio{-999};
        float ssm_max_dev{-999};
        float ssm_n_prim_tracks_1{-999};
        float ssm_n_prim_tracks_3{-999};
        float ssm_n_prim_tracks_5{-999};
        float ssm_n_prim_tracks_8{-999};
        float ssm_n_prim_tracks_11{-999};
        float ssm_n_all_tracks_1{-999};
        float ssm_n_all_tracks_3{-999};
        float ssm_n_all_tracks_5{-999};
        float ssm_n_all_tracks_8{-999};
        float ssm_n_all_tracks_11{-999};
        float ssm_n_daughter_tracks_1{-999};
        float ssm_n_daughter_tracks_3{-999};
        float ssm_n_daughter_tracks_5{-999};
        float ssm_n_daughter_tracks_8{-999};
        float ssm_n_daughter_tracks_11{-999};
        float ssm_n_daughter_all_1{-999};
        float ssm_n_daughter_all_3{-999};
        float ssm_n_daughter_all_5{-999};
        float ssm_n_daughter_all_8{-999};
        float ssm_n_daughter_all_11{-999};

        // leading other primary track
        float ssm_prim_track1_pdg{-999};
        float ssm_prim_track1_score_mu_fwd{-999};
        float ssm_prim_track1_score_p_fwd{-999};
        float ssm_prim_track1_score_e_fwd{-999};
        float ssm_prim_track1_score_mu_bck{-999};
        float ssm_prim_track1_score_p_bck{-999};
        float ssm_prim_track1_score_e_bck{-999};
        float ssm_prim_track1_length{-999};
        float ssm_prim_track1_direct_length{-999};
        float ssm_prim_track1_length_ratio{-999};
        float ssm_prim_track1_max_dev{-999};
        float ssm_prim_track1_kine_energy_range{-999};
        float ssm_prim_track1_kine_energy_range_mu{-999};
        float ssm_prim_track1_kine_energy_range_p{-999};
        float ssm_prim_track1_kine_energy_range_e{-999};
        float ssm_prim_track1_kine_energy_cal{-999};
        float ssm_prim_track1_medium_dq_dx{-999};
        float ssm_prim_track1_x_dir{-999};
        float ssm_prim_track1_y_dir{-999};
        float ssm_prim_track1_z_dir{-999};
        float ssm_prim_track1_add_daught_track_counts_1{-999};
        float ssm_prim_track1_add_daught_all_counts_1{-999};
        float ssm_prim_track1_add_daught_track_counts_5{-999};
        float ssm_prim_track1_add_daught_all_counts_5{-999};
        float ssm_prim_track1_add_daught_track_counts_11{-999};
        float ssm_prim_track1_add_daught_all_counts_11{-999};

        // sub-leading other primary track
        float ssm_prim_track2_pdg{-999};
        float ssm_prim_track2_score_mu_fwd{-999};
        float ssm_prim_track2_score_p_fwd{-999};
        float ssm_prim_track2_score_e_fwd{-999};
        float ssm_prim_track2_score_mu_bck{-999};
        float ssm_prim_track2_score_p_bck{-999};
        float ssm_prim_track2_score_e_bck{-999};
        float ssm_prim_track2_length{-999};
        float ssm_prim_track2_direct_length{-999};
        float ssm_prim_track2_length_ratio{-999};
        float ssm_prim_track2_max_dev{-999};
        float ssm_prim_track2_kine_energy_range{-999};
        float ssm_prim_track2_kine_energy_range_mu{-999};
        float ssm_prim_track2_kine_energy_range_p{-999};
        float ssm_prim_track2_kine_energy_range_e{-999};
        float ssm_prim_track2_kine_energy_cal{-999};
        float ssm_prim_track2_medium_dq_dx{-999};
        float ssm_prim_track2_x_dir{-999};
        float ssm_prim_track2_y_dir{-999};
        float ssm_prim_track2_z_dir{-999};
        float ssm_prim_track2_add_daught_track_counts_1{-999};
        float ssm_prim_track2_add_daught_all_counts_1{-999};
        float ssm_prim_track2_add_daught_track_counts_5{-999};
        float ssm_prim_track2_add_daught_all_counts_5{-999};
        float ssm_prim_track2_add_daught_track_counts_11{-999};
        float ssm_prim_track2_add_daught_all_counts_11{-999};

        // leading daughter track
        float ssm_daught_track1_pdg{-999};
        float ssm_daught_track1_score_mu_fwd{-999};
        float ssm_daught_track1_score_p_fwd{-999};
        float ssm_daught_track1_score_e_fwd{-999};
        float ssm_daught_track1_score_mu_bck{-999};
        float ssm_daught_track1_score_p_bck{-999};
        float ssm_daught_track1_score_e_bck{-999};
        float ssm_daught_track1_length{-999};
        float ssm_daught_track1_direct_length{-999};
        float ssm_daught_track1_length_ratio{-999};
        float ssm_daught_track1_max_dev{-999};
        float ssm_daught_track1_kine_energy_range{-999};
        float ssm_daught_track1_kine_energy_range_mu{-999};
        float ssm_daught_track1_kine_energy_range_p{-999};
        float ssm_daught_track1_kine_energy_range_e{-999};
        float ssm_daught_track1_kine_energy_cal{-999};
        float ssm_daught_track1_medium_dq_dx{-999};
        float ssm_daught_track1_x_dir{-999};
        float ssm_daught_track1_y_dir{-999};
        float ssm_daught_track1_z_dir{-999};
        float ssm_daught_track1_add_daught_track_counts_1{-999};
        float ssm_daught_track1_add_daught_all_counts_1{-999};
        float ssm_daught_track1_add_daught_track_counts_5{-999};
        float ssm_daught_track1_add_daught_all_counts_5{-999};
        float ssm_daught_track1_add_daught_track_counts_11{-999};
        float ssm_daught_track1_add_daught_all_counts_11{-999};

        // sub-leading daughter track
        float ssm_daught_track2_pdg{-999};
        float ssm_daught_track2_score_mu_fwd{-999};
        float ssm_daught_track2_score_p_fwd{-999};
        float ssm_daught_track2_score_e_fwd{-999};
        float ssm_daught_track2_score_mu_bck{-999};
        float ssm_daught_track2_score_p_bck{-999};
        float ssm_daught_track2_score_e_bck{-999};
        float ssm_daught_track2_length{-999};
        float ssm_daught_track2_direct_length{-999};
        float ssm_daught_track2_length_ratio{-999};
        float ssm_daught_track2_max_dev{-999};
        float ssm_daught_track2_kine_energy_range{-999};
        float ssm_daught_track2_kine_energy_range_mu{-999};
        float ssm_daught_track2_kine_energy_range_p{-999};
        float ssm_daught_track2_kine_energy_range_e{-999};
        float ssm_daught_track2_kine_energy_cal{-999};
        float ssm_daught_track2_medium_dq_dx{-999};
        float ssm_daught_track2_x_dir{-999};
        float ssm_daught_track2_y_dir{-999};
        float ssm_daught_track2_z_dir{-999};
        float ssm_daught_track2_add_daught_track_counts_1{-999};
        float ssm_daught_track2_add_daught_all_counts_1{-999};
        float ssm_daught_track2_add_daught_track_counts_5{-999};
        float ssm_daught_track2_add_daught_all_counts_5{-999};
        float ssm_daught_track2_add_daught_track_counts_11{-999};
        float ssm_daught_track2_add_daught_all_counts_11{-999};

        // leading other primary shower
        float ssm_prim_shw1_pdg{-999};
        float ssm_prim_shw1_score_mu_fwd{-999};
        float ssm_prim_shw1_score_p_fwd{-999};
        float ssm_prim_shw1_score_e_fwd{-999};
        float ssm_prim_shw1_score_mu_bck{-999};
        float ssm_prim_shw1_score_p_bck{-999};
        float ssm_prim_shw1_score_e_bck{-999};
        float ssm_prim_shw1_length{-999};
        float ssm_prim_shw1_direct_length{-999};
        float ssm_prim_shw1_length_ratio{-999};
        float ssm_prim_shw1_max_dev{-999};
        float ssm_prim_shw1_kine_energy_range{-999};
        float ssm_prim_shw1_kine_energy_range_mu{-999};
        float ssm_prim_shw1_kine_energy_range_p{-999};
        float ssm_prim_shw1_kine_energy_range_e{-999};
        float ssm_prim_shw1_kine_energy_cal{-999};
        float ssm_prim_shw1_kine_energy_best{-999};
        float ssm_prim_shw1_medium_dq_dx{-999};
        float ssm_prim_shw1_x_dir{-999};
        float ssm_prim_shw1_y_dir{-999};
        float ssm_prim_shw1_z_dir{-999};
        float ssm_prim_shw1_add_daught_track_counts_1{-999};
        float ssm_prim_shw1_add_daught_all_counts_1{-999};
        float ssm_prim_shw1_add_daught_track_counts_5{-999};
        float ssm_prim_shw1_add_daught_all_counts_5{-999};
        float ssm_prim_shw1_add_daught_track_counts_11{-999};
        float ssm_prim_shw1_add_daught_all_counts_11{-999};

        // sub-leading other primary shower
        float ssm_prim_shw2_pdg{-999};
        float ssm_prim_shw2_score_mu_fwd{-999};
        float ssm_prim_shw2_score_p_fwd{-999};
        float ssm_prim_shw2_score_e_fwd{-999};
        float ssm_prim_shw2_score_mu_bck{-999};
        float ssm_prim_shw2_score_p_bck{-999};
        float ssm_prim_shw2_score_e_bck{-999};
        float ssm_prim_shw2_length{-999};
        float ssm_prim_shw2_direct_length{-999};
        float ssm_prim_shw2_length_ratio{-999};
        float ssm_prim_shw2_max_dev{-999};
        float ssm_prim_shw2_kine_energy_range{-999};
        float ssm_prim_shw2_kine_energy_range_mu{-999};
        float ssm_prim_shw2_kine_energy_range_p{-999};
        float ssm_prim_shw2_kine_energy_range_e{-999};
        float ssm_prim_shw2_kine_energy_cal{-999};
        float ssm_prim_shw2_kine_energy_best{-999};
        float ssm_prim_shw2_medium_dq_dx{-999};
        float ssm_prim_shw2_x_dir{-999};
        float ssm_prim_shw2_y_dir{-999};
        float ssm_prim_shw2_z_dir{-999};
        float ssm_prim_shw2_add_daught_track_counts_1{-999};
        float ssm_prim_shw2_add_daught_all_counts_1{-999};
        float ssm_prim_shw2_add_daught_track_counts_5{-999};
        float ssm_prim_shw2_add_daught_all_counts_5{-999};
        float ssm_prim_shw2_add_daught_track_counts_11{-999};
        float ssm_prim_shw2_add_daught_all_counts_11{-999};

        // leading daughter shower
        float ssm_daught_shw1_pdg{-999};
        float ssm_daught_shw1_score_mu_fwd{-999};
        float ssm_daught_shw1_score_p_fwd{-999};
        float ssm_daught_shw1_score_e_fwd{-999};
        float ssm_daught_shw1_score_mu_bck{-999};
        float ssm_daught_shw1_score_p_bck{-999};
        float ssm_daught_shw1_score_e_bck{-999};
        float ssm_daught_shw1_length{-999};
        float ssm_daught_shw1_direct_length{-999};
        float ssm_daught_shw1_length_ratio{-999};
        float ssm_daught_shw1_max_dev{-999};
        float ssm_daught_shw1_kine_energy_range{-999};
        float ssm_daught_shw1_kine_energy_range_mu{-999};
        float ssm_daught_shw1_kine_energy_range_p{-999};
        float ssm_daught_shw1_kine_energy_range_e{-999};
        float ssm_daught_shw1_kine_energy_cal{-999};
        float ssm_daught_shw1_kine_energy_best{-999};
        float ssm_daught_shw1_medium_dq_dx{-999};
        float ssm_daught_shw1_x_dir{-999};
        float ssm_daught_shw1_y_dir{-999};
        float ssm_daught_shw1_z_dir{-999};
        float ssm_daught_shw1_add_daught_track_counts_1{-999};
        float ssm_daught_shw1_add_daught_all_counts_1{-999};
        float ssm_daught_shw1_add_daught_track_counts_5{-999};
        float ssm_daught_shw1_add_daught_all_counts_5{-999};
        float ssm_daught_shw1_add_daught_track_counts_11{-999};
        float ssm_daught_shw1_add_daught_all_counts_11{-999};

        // sub-leading daughter shower
        float ssm_daught_shw2_pdg{-999};
        float ssm_daught_shw2_score_mu_fwd{-999};
        float ssm_daught_shw2_score_p_fwd{-999};
        float ssm_daught_shw2_score_e_fwd{-999};
        float ssm_daught_shw2_score_mu_bck{-999};
        float ssm_daught_shw2_score_p_bck{-999};
        float ssm_daught_shw2_score_e_bck{-999};
        float ssm_daught_shw2_length{-999};
        float ssm_daught_shw2_direct_length{-999};
        float ssm_daught_shw2_length_ratio{-999};
        float ssm_daught_shw2_max_dev{-999};
        float ssm_daught_shw2_kine_energy_range{-999};
        float ssm_daught_shw2_kine_energy_range_mu{-999};
        float ssm_daught_shw2_kine_energy_range_p{-999};
        float ssm_daught_shw2_kine_energy_range_e{-999};
        float ssm_daught_shw2_kine_energy_cal{-999};
        float ssm_daught_shw2_kine_energy_best{-999};
        float ssm_daught_shw2_medium_dq_dx{-999};
        float ssm_daught_shw2_x_dir{-999};
        float ssm_daught_shw2_y_dir{-999};
        float ssm_daught_shw2_z_dir{-999};
        float ssm_daught_shw2_add_daught_track_counts_1{-999};
        float ssm_daught_shw2_add_daught_all_counts_1{-999};
        float ssm_daught_shw2_add_daught_track_counts_5{-999};
        float ssm_daught_shw2_add_daught_all_counts_5{-999};
        float ssm_daught_shw2_add_daught_track_counts_11{-999};
        float ssm_daught_shw2_add_daught_all_counts_11{-999};

        // event level angles and vertex position
        float ssm_nu_angle_z{-999};
        float ssm_nu_angle_target{-999};
        float ssm_nu_angle_absorber{-999};
        float ssm_nu_angle_vertical{-999};
        float ssm_con_nu_angle_z{-999};
        float ssm_con_nu_angle_target{-999};
        float ssm_con_nu_angle_absorber{-999};
        float ssm_con_nu_angle_vertical{-999};
        float ssm_prim_nu_angle_z{-999};
        float ssm_prim_nu_angle_target{-999};
        float ssm_prim_nu_angle_absorber{-999};
        float ssm_prim_nu_angle_vertical{-999};
        float ssm_track_angle_z{-999};
        float ssm_track_angle_target{-999};
        float ssm_track_angle_absorber{-999};
        float ssm_track_angle_vertical{-999};
        float ssm_vtxX{-999};
        float ssm_vtxY{-999};
        float ssm_vtxZ{-999};

        // off-vertex activity
        float ssm_offvtx_length{-999};
        float ssm_offvtx_energy{-999};
        float ssm_n_offvtx_tracks_1{-999};
        float ssm_n_offvtx_tracks_3{-999};
        float ssm_n_offvtx_tracks_5{-999};
        float ssm_n_offvtx_tracks_8{-999};
        float ssm_n_offvtx_tracks_11{-999};
        float ssm_n_offvtx_showers_1{-999};
        float ssm_n_offvtx_showers_3{-999};
        float ssm_n_offvtx_showers_5{-999};
        float ssm_n_offvtx_showers_8{-999};
        float ssm_n_offvtx_showers_11{-999};

        // leading off-vertex track
        float ssm_offvtx_track1_pdg{-999};
        float ssm_offvtx_track1_score_mu_fwd{-999};
        float ssm_offvtx_track1_score_p_fwd{-999};
        float ssm_offvtx_track1_score_e_fwd{-999};
        float ssm_offvtx_track1_score_mu_bck{-999};
        float ssm_offvtx_track1_score_p_bck{-999};
        float ssm_offvtx_track1_score_e_bck{-999};
        float ssm_offvtx_track1_length{-999};
        float ssm_offvtx_track1_direct_length{-999};
        float ssm_offvtx_track1_max_dev{-999};
        float ssm_offvtx_track1_kine_energy_range{-999};
        float ssm_offvtx_track1_kine_energy_range_mu{-999};
        float ssm_offvtx_track1_kine_energy_range_p{-999};
        float ssm_offvtx_track1_kine_energy_range_e{-999};
        float ssm_offvtx_track1_kine_energy_cal{-999};
        float ssm_offvtx_track1_medium_dq_dx{-999};
        float ssm_offvtx_track1_x_dir{-999};
        float ssm_offvtx_track1_y_dir{-999};
        float ssm_offvtx_track1_z_dir{-999};
        float ssm_offvtx_track1_dist_mainvtx{-999};

        // leading off-vertex shower
        float ssm_offvtx_shw1_pdg_offvtx{-999};
        float ssm_offvtx_shw1_score_mu_fwd{-999};
        float ssm_offvtx_shw1_score_p_fwd{-999};
        float ssm_offvtx_shw1_score_e_fwd{-999};
        float ssm_offvtx_shw1_score_mu_bck{-999};
        float ssm_offvtx_shw1_score_p_bck{-999};
        float ssm_offvtx_shw1_score_e_bck{-999};
        float ssm_offvtx_shw1_length{-999};
        float ssm_offvtx_shw1_direct_length{-999};
        float ssm_offvtx_shw1_max_dev{-999};
        float ssm_offvtx_shw1_kine_energy_best{-999};
        float ssm_offvtx_shw1_kine_energy_range{-999};
        float ssm_offvtx_shw1_kine_energy_range_mu{-999};
        float ssm_offvtx_shw1_kine_energy_range_p{-999};
        float ssm_offvtx_shw1_kine_energy_range_e{-999};
        float ssm_offvtx_shw1_kine_energy_cal{-999};
        float ssm_offvtx_shw1_medium_dq_dx{-999};
        float ssm_offvtx_shw1_x_dir{-999};
        float ssm_offvtx_shw1_y_dir{-999};
        float ssm_offvtx_shw1_z_dir{-999};
        float ssm_offvtx_shw1_dist_mainvtx{-999};

        // SSM spacepoints
        int   ssmsp_Ntrack{0};
        int   ssmsp_Nsp_tot{0};
        std::vector<int>   ssmsp_Nsp;
        std::vector<int>   ssmsp_pdg;
        std::vector<int>   ssmsp_id;
        std::vector<int>   ssmsp_mother;
        std::vector<float> ssmsp_x;
        std::vector<float> ssmsp_y;
        std::vector<float> ssmsp_z;
        std::vector<float> ssmsp_dx;
        std::vector<float> ssmsp_dQ;
        std::vector<float> ssmsp_KE;
        std::vector<float> ssmsp_containing_shower_id;
        std::vector<float> ssmsp_containing_shower_ke;
        std::vector<float> ssmsp_containing_shower_flag;

        // SSM kinematic variables
        float ssm_kine_reco_Enu{-999};
        float ssm_kine_reco_add_energy{-999};
        std::vector<float> ssm_kine_energy_particle;
        std::vector<int>   ssm_kine_energy_info;
        std::vector<int>   ssm_kine_particle_type;
        std::vector<int>   ssm_kine_energy_included;
        float ssm_kine_pio_mass{-999};
        int   ssm_kine_pio_flag{-999};
        float ssm_kine_pio_vtx_dis{-999};
        float ssm_kine_pio_energy_1{-999};
        float ssm_kine_pio_theta_1{-999};
        float ssm_kine_pio_phi_1{-999};
        float ssm_kine_pio_dis_1{-999};
        float ssm_kine_pio_energy_2{-999};
        float ssm_kine_pio_theta_2{-999};
        float ssm_kine_pio_phi_2{-999};
        float ssm_kine_pio_dis_2{-999};
        float ssm_kine_pio_angle{-999};

        // ---- single-photon shower identification ----------------------- //
        float shw_sp_flag{1};
        float shw_sp_filled{0};
        float shw_sp_num_mip_tracks{0};
        float shw_sp_num_muons{0};
        float shw_sp_num_pions{0};
        float shw_sp_num_protons{0};
        float shw_sp_proton_length_1{-1};
        float shw_sp_proton_dqdx_1{-1};
        float shw_sp_proton_energy_1{-1};
        float shw_sp_proton_length_2{-1};
        float shw_sp_proton_dqdx_2{-1};
        float shw_sp_proton_energy_2{-1};
        float shw_sp_n_good_showers{0};
        float shw_sp_n_20mev_showers{0};
        float shw_sp_n_br1_showers{0};
        float shw_sp_n_br2_showers{0};
        float shw_sp_n_br3_showers{0};
        float shw_sp_n_br4_showers{0};
        float shw_sp_n_20br1_showers{0};
        std::vector<int> shw_sp_20mev_showers;
        std::vector<int> shw_sp_br1_showers;
        std::vector<int> shw_sp_br2_showers;
        std::vector<int> shw_sp_br3_showers;
        std::vector<int> shw_sp_br4_showers;
        float shw_sp_shw_vtx_dis{-1};
        float shw_sp_max_shw_dis{-1};
        float shw_sp_energy{0};
        float shw_sp_vec_dQ_dx_0{1};
        float shw_sp_vec_dQ_dx_1{1};
        float shw_sp_max_dQ_dx_sample{1};
        float shw_sp_n_below_threshold{19};
        float shw_sp_n_below_zero{0};
        float shw_sp_n_lowest{1};
        float shw_sp_n_highest{1};
        float shw_sp_lowest_dQ_dx{1};
        float shw_sp_highest_dQ_dx{1};
        float shw_sp_medium_dQ_dx{1};
        float shw_sp_stem_length{1};
        float shw_sp_length_main{1};
        float shw_sp_length_total{1};
        float shw_sp_angle_beam{0};
        float shw_sp_iso_angle{0};
        float shw_sp_n_vertex{1};
        float shw_sp_n_good_tracks{0};
        float shw_sp_E_indirect_max_energy{0};
        float shw_sp_flag_all_above{0};
        float shw_sp_min_dQ_dx_5{1};
        float shw_sp_n_other_vertex{2};
        float shw_sp_n_stem_size{20};
        float shw_sp_flag_stem_trajectory{0};
        float shw_sp_min_dis{0};

        // extra dQ/dx samples
        float shw_sp_vec_dQ_dx_2{0};
        float shw_sp_vec_dQ_dx_3{0};
        float shw_sp_vec_dQ_dx_4{0};
        float shw_sp_vec_dQ_dx_5{0};
        float shw_sp_vec_dQ_dx_6{0};
        float shw_sp_vec_dQ_dx_7{0};
        float shw_sp_vec_dQ_dx_8{0};
        float shw_sp_vec_dQ_dx_9{0};
        float shw_sp_vec_dQ_dx_10{0};
        float shw_sp_vec_dQ_dx_11{0};
        float shw_sp_vec_dQ_dx_12{0};
        float shw_sp_vec_dQ_dx_13{0};
        float shw_sp_vec_dQ_dx_14{0};
        float shw_sp_vec_dQ_dx_15{0};
        float shw_sp_vec_dQ_dx_16{0};
        float shw_sp_vec_dQ_dx_17{0};
        float shw_sp_vec_dQ_dx_18{0};
        float shw_sp_vec_dQ_dx_19{0};
        float shw_sp_vec_median_dedx{0};
        float shw_sp_vec_mean_dedx{0};

        // photon shower pi0
        float shw_sp_pio_flag{1};
        float shw_sp_pio_mip_id{0};
        float shw_sp_pio_filled{0};
        float shw_sp_pio_flag_pio{0};
        float shw_sp_pio_1_flag{1};
        float shw_sp_pio_1_mass{0};
        float shw_sp_pio_1_pio_type{0};
        float shw_sp_pio_1_energy_1{0};
        float shw_sp_pio_1_energy_2{0};
        float shw_sp_pio_1_dis_1{0};
        float shw_sp_pio_1_dis_2{0};
        std::vector<float> shw_sp_pio_2_v_dis2;
        std::vector<float> shw_sp_pio_2_v_angle2;
        std::vector<float> shw_sp_pio_2_v_acc_length;
        std::vector<float> shw_sp_pio_2_v_flag;

        // low-energy michel (single photon)
        float shw_sp_lem_shower_total_length{0};
        float shw_sp_lem_shower_main_length{0};
        float shw_sp_lem_n_3seg{0};
        float shw_sp_lem_e_charge{0};
        float shw_sp_lem_e_dQdx{0};
        float shw_sp_lem_shower_num_segs{0};
        float shw_sp_lem_shower_num_main_segs{0};
        float shw_sp_lem_flag{1};

        // bad reconstruction
        float shw_sp_br_filled{0};
        float shw_sp_br1_flag{1};
        float shw_sp_br1_1_flag{1};
        float shw_sp_br1_1_shower_type{0};
        float shw_sp_br1_1_vtx_n_segs{0};
        float shw_sp_br1_1_energy{0};
        float shw_sp_br1_1_n_segs{0};
        float shw_sp_br1_1_flag_sg_topology{0};
        float shw_sp_br1_1_flag_sg_trajectory{0};
        float shw_sp_br1_1_sg_length{0};
        float shw_sp_br1_2_flag{1};
        float shw_sp_br1_2_energy{0};
        float shw_sp_br1_2_n_connected{0};
        float shw_sp_br1_2_max_length{0};
        float shw_sp_br1_2_n_connected_1{0};
        float shw_sp_br1_2_vtx_n_segs{0};
        float shw_sp_br1_2_n_shower_segs{0};
        float shw_sp_br1_2_max_length_ratio{1};
        float shw_sp_br1_2_shower_length{0};
        float shw_sp_br1_3_flag{1};
        float shw_sp_br1_3_energy{0};
        float shw_sp_br1_3_n_connected_p{0};
        float shw_sp_br1_3_max_length_p{0};
        float shw_sp_br1_3_n_shower_segs{0};
        float shw_sp_br1_3_flag_sg_topology{0};
        float shw_sp_br1_3_flag_sg_trajectory{0};
        float shw_sp_br1_3_n_shower_main_segs{0};
        float shw_sp_br1_3_sg_length{0};

        float shw_sp_br2_flag{1};
        float shw_sp_br2_flag_single_shower{0};
        float shw_sp_br2_num_valid_tracks{0};
        float shw_sp_br2_energy{0};
        float shw_sp_br2_angle1{0};
        float shw_sp_br2_angle2{0};
        float shw_sp_br2_angle{0};
        float shw_sp_br2_angle3{0};
        float shw_sp_br2_n_shower_main_segs{0};
        float shw_sp_br2_max_angle{0};
        float shw_sp_br2_sg_length{0};
        float shw_sp_br2_flag_sg_trajectory{0};

        float shw_sp_lol_flag{1};
        float shw_sp_lol_3_flag{1};
        float shw_sp_lol_3_angle_beam{0};
        float shw_sp_lol_3_min_angle{0};
        float shw_sp_lol_3_n_valid_tracks{0};
        float shw_sp_lol_3_vtx_n_segs{0};
        float shw_sp_lol_3_energy{0};
        float shw_sp_lol_3_shower_main_length{0};
        float shw_sp_lol_3_n_sum{0};
        float shw_sp_lol_3_n_out{0};
        std::vector<float> shw_sp_lol_1_v_energy;
        std::vector<float> shw_sp_lol_1_v_vtx_n_segs;
        std::vector<float> shw_sp_lol_1_v_nseg;
        std::vector<float> shw_sp_lol_1_v_angle;
        std::vector<float> shw_sp_lol_1_v_flag;
        std::vector<float> shw_sp_lol_2_v_length;
        std::vector<float> shw_sp_lol_2_v_angle;
        std::vector<float> shw_sp_lol_2_v_type;
        std::vector<float> shw_sp_lol_2_v_vtx_n_segs;
        std::vector<float> shw_sp_lol_2_v_energy;
        std::vector<float> shw_sp_lol_2_v_shower_main_length;
        std::vector<float> shw_sp_lol_2_v_flag_dir_weak;
        std::vector<float> shw_sp_lol_2_v_flag;

        float shw_sp_br3_1_energy{0};
        float shw_sp_br3_1_n_shower_segments{0};
        float shw_sp_br3_1_sg_flag_trajectory{0};
        float shw_sp_br3_1_sg_direct_length{0};
        float shw_sp_br3_1_sg_length{0};
        float shw_sp_br3_1_total_main_length{0};
        float shw_sp_br3_1_total_length{0};
        float shw_sp_br3_1_iso_angle{0};
        float shw_sp_br3_1_sg_flag_topology{0};
        float shw_sp_br3_1_flag{1};
        float shw_sp_br3_2_n_ele{0};
        float shw_sp_br3_2_n_other{0};
        float shw_sp_br3_2_energy{0};
        float shw_sp_br3_2_total_main_length{0};
        float shw_sp_br3_2_total_length{0};
        float shw_sp_br3_2_other_fid{0};
        float shw_sp_br3_2_flag{1};
        std::vector<float> shw_sp_br3_3_v_energy;
        std::vector<float> shw_sp_br3_3_v_angle;
        std::vector<float> shw_sp_br3_3_v_dir_length;
        std::vector<float> shw_sp_br3_3_v_length;
        std::vector<float> shw_sp_br3_3_v_flag;
        float shw_sp_br3_4_acc_length{0};
        float shw_sp_br3_4_total_length{0};
        float shw_sp_br3_4_energy{0};
        float shw_sp_br3_4_flag{1};
        std::vector<float> shw_sp_br3_5_v_dir_length;
        std::vector<float> shw_sp_br3_5_v_total_length;
        std::vector<float> shw_sp_br3_5_v_flag_avoid_muon_check;
        std::vector<float> shw_sp_br3_5_v_n_seg;
        std::vector<float> shw_sp_br3_5_v_angle;
        std::vector<float> shw_sp_br3_5_v_sg_length;
        std::vector<float> shw_sp_br3_5_v_energy;
        std::vector<float> shw_sp_br3_5_v_n_main_segs;
        std::vector<float> shw_sp_br3_5_v_n_segs;
        std::vector<float> shw_sp_br3_5_v_shower_main_length;
        std::vector<float> shw_sp_br3_5_v_shower_total_length;
        std::vector<float> shw_sp_br3_5_v_flag;
        std::vector<float> shw_sp_br3_6_v_angle;
        std::vector<float> shw_sp_br3_6_v_angle1;
        std::vector<float> shw_sp_br3_6_v_flag_shower_trajectory;
        std::vector<float> shw_sp_br3_6_v_direct_length;
        std::vector<float> shw_sp_br3_6_v_length;
        std::vector<float> shw_sp_br3_6_v_n_other_vtx_segs;
        std::vector<float> shw_sp_br3_6_v_energy;
        std::vector<float> shw_sp_br3_6_v_flag;
        float shw_sp_br3_7_energy{0};
        float shw_sp_br3_7_min_angle{0};
        float shw_sp_br3_7_sg_length{0};
        float shw_sp_br3_7_shower_main_length{0};
        float shw_sp_br3_7_flag{1};
        float shw_sp_br3_8_max_dQ_dx{0};
        float shw_sp_br3_8_energy{0};
        float shw_sp_br3_8_n_main_segs{0};
        float shw_sp_br3_8_shower_main_length{0};
        float shw_sp_br3_8_shower_length{0};
        float shw_sp_br3_8_flag{1};
        float shw_sp_br3_flag{1};

        float shw_sp_br4_1_shower_main_length{0};
        float shw_sp_br4_1_shower_total_length{0};
        float shw_sp_br4_1_min_dis{0};
        float shw_sp_br4_1_energy{0};
        float shw_sp_br4_1_flag_avoid_muon_check{0};
        float shw_sp_br4_1_n_vtx_segs{0};
        float shw_sp_br4_1_n_main_segs{0};
        float shw_sp_br4_1_flag{1};
        float shw_sp_br4_2_ratio_45{1};
        float shw_sp_br4_2_ratio_35{1};
        float shw_sp_br4_2_ratio_25{1};
        float shw_sp_br4_2_ratio_15{1};
        float shw_sp_br4_2_energy{0};
        float shw_sp_br4_2_ratio1_45{1};
        float shw_sp_br4_2_ratio1_35{1};
        float shw_sp_br4_2_ratio1_25{1};
        float shw_sp_br4_2_ratio1_15{1};
        float shw_sp_br4_2_iso_angle{0};
        float shw_sp_br4_2_iso_angle1{0};
        float shw_sp_br4_2_angle{0};
        float shw_sp_br4_2_flag{1};
        float shw_sp_br4_flag{1};

        float shw_sp_hol_1_n_valid_tracks{0};
        float shw_sp_hol_1_min_angle{0};
        float shw_sp_hol_1_energy{0};
        float shw_sp_hol_1_flag_all_shower{0};
        float shw_sp_hol_1_min_length{0};
        float shw_sp_hol_1_flag{1};
        float shw_sp_hol_2_min_angle{0};
        float shw_sp_hol_2_medium_dQ_dx{1};
        float shw_sp_hol_2_ncount{0};
        float shw_sp_hol_2_energy{0};
        float shw_sp_hol_2_flag{1};
        float shw_sp_hol_flag{1};

        // ---- shower pi0 identification --------------------------------- //
        float pio_flag{1};
        float pio_mip_id{0};
        float pio_filled{0};
        float pio_flag_pio{0};
        float pio_1_flag{1};
        float pio_1_mass{0};
        float pio_1_pio_type{0};
        float pio_1_energy_1{0};
        float pio_1_energy_2{0};
        float pio_1_dis_1{0};
        float pio_1_dis_2{0};
        std::vector<float> pio_2_v_dis2;
        std::vector<float> pio_2_v_angle2;
        std::vector<float> pio_2_v_acc_length;
        std::vector<float> pio_2_v_flag;

        // ---- single shower pi0 case ------------------------------------ //
        std::vector<float> sig_1_v_angle;
        std::vector<float> sig_1_v_flag_single_shower;
        std::vector<float> sig_1_v_energy;
        std::vector<float> sig_1_v_energy_1;
        std::vector<float> sig_1_v_flag;
        std::vector<float> sig_2_v_energy;
        std::vector<float> sig_2_v_shower_angle;
        std::vector<float> sig_2_v_flag_single_shower;
        std::vector<float> sig_2_v_medium_dQ_dx;
        std::vector<float> sig_2_v_start_dQ_dx;
        std::vector<float> sig_2_v_flag;
        float sig_flag{1};

        // ---- multiple gamma -------------------------------------------- //
        float mgo_energy{0};
        float mgo_max_energy{0};
        float mgo_total_energy{0};
        float mgo_n_showers{0};
        float mgo_max_energy_1{0};
        float mgo_max_energy_2{0};
        float mgo_total_other_energy{0};
        float mgo_n_total_showers{0};
        float mgo_total_other_energy_1{0};
        float mgo_flag{1};
        float mgt_flag_single_shower{0};
        float mgt_max_energy{0};
        float mgt_energy{0};
        float mgt_total_other_energy{0};
        float mgt_max_energy_1{0};
        float mgt_e_indirect_max_energy{0};
        float mgt_e_direct_max_energy{0};
        float mgt_n_direct_showers{0};
        float mgt_e_direct_total_energy{0};
        float mgt_flag_indirect_max_pio{0};
        float mgt_e_indirect_total_energy{0};
        float mgt_flag{1};

        // ---- shower to wall -------------------------------------------- //
        float stw_1_energy{0};
        float stw_1_dis{0};
        float stw_1_dQ_dx{1};
        float stw_1_flag_single_shower{0};
        float stw_1_n_pi0{0};
        float stw_1_num_valid_tracks{0};
        float stw_1_flag{1};
        std::vector<float> stw_2_v_medium_dQ_dx;
        std::vector<float> stw_2_v_energy;
        std::vector<float> stw_2_v_angle;
        std::vector<float> stw_2_v_dir_length;
        std::vector<float> stw_2_v_max_dQ_dx;
        std::vector<float> stw_2_v_flag;
        std::vector<float> stw_3_v_angle;
        std::vector<float> stw_3_v_dir_length;
        std::vector<float> stw_3_v_energy;
        std::vector<float> stw_3_v_medium_dQ_dx;
        std::vector<float> stw_3_v_flag;
        std::vector<float> stw_4_v_angle;
        std::vector<float> stw_4_v_dis;
        std::vector<float> stw_4_v_energy;
        std::vector<float> stw_4_v_flag;
        float stw_flag{1};

        // ---- single photon case ---------------------------------------- //
        float spt_flag_single_shower{0};
        float spt_energy{0};
        float spt_shower_main_length{0};
        float spt_shower_total_length{0};
        float spt_angle_beam{0};
        float spt_angle_vertical{0};
        float spt_max_dQ_dx{1};
        float spt_angle_beam_1{0};
        float spt_angle_drift{0};
        float spt_angle_drift_1{0};
        float spt_num_valid_tracks{0};
        float spt_n_vtx_segs{0};
        float spt_max_length{0};
        float spt_flag{1};

        // ---- stem length ----------------------------------------------- //
        float stem_len_energy{0};
        float stem_len_length{0};
        float stem_len_flag_avoid_muon_check{0};
        float stem_len_num_daughters{0};
        float stem_len_daughter_length{0};
        float stem_len_flag{1};

        // ---- low-energy michel ----------------------------------------- //
        float lem_shower_total_length{0};
        float lem_shower_main_length{0};
        float lem_n_3seg{0};
        float lem_e_charge{0};
        float lem_e_dQdx{0};
        float lem_shower_num_segs{0};
        float lem_shower_num_main_segs{0};
        float lem_flag{1};

        // ---- broken muon ----------------------------------------------- //
        float brm_n_mu_segs{0};
        float brm_Ep{0};
        float brm_energy{0};
        float brm_acc_length{0};
        float brm_shower_total_length{0};
        float brm_connected_length{0};
        float brm_n_size{0};
        float brm_acc_direct_length{0};
        float brm_n_shower_main_segs{0};
        float brm_n_mu_main{0};
        float brm_flag{1};

        // ---- compare muon energy --------------------------------------- //
        float cme_mu_energy{0};
        float cme_energy{0};
        float cme_mu_length{0};
        float cme_length{0};
        float cme_angle_beam{0};
        float cme_flag{1};

        // ---- angular cut ----------------------------------------------- //
        float anc_energy{0};
        float anc_angle{0};
        float anc_max_angle{0};
        float anc_max_length{0};
        float anc_acc_forward_length{0};
        float anc_acc_backward_length{0};
        float anc_acc_forward_length1{0};
        float anc_shower_main_length{0};
        float anc_shower_total_length{0};
        float anc_flag_main_outside{0};
        float anc_flag{1};

        // ---- stem direction -------------------------------------------- //
        float stem_dir_flag{1};
        float stem_dir_flag_single_shower{0};
        float stem_dir_filled{0};
        float stem_dir_angle{0};
        float stem_dir_energy{0};
        float stem_dir_angle1{0};
        float stem_dir_angle2{0};
        float stem_dir_angle3{0};
        float stem_dir_ratio{1};

        // ---- vertex inside shower -------------------------------------- //
        float vis_1_filled{0};
        float vis_1_n_vtx_segs{0};
        float vis_1_energy{0};
        float vis_1_num_good_tracks{0};
        float vis_1_max_angle{0};
        float vis_1_max_shower_angle{0};
        float vis_1_tmp_length1{0};
        float vis_1_tmp_length2{0};
        float vis_1_particle_type{0};
        float vis_1_flag{1};
        float vis_2_filled{0};
        float vis_2_n_vtx_segs{0};
        float vis_2_min_angle{0};
        float vis_2_min_weak_track{0};
        float vis_2_angle_beam{0};
        float vis_2_min_angle1{0};
        float vis_2_iso_angle1{0};
        float vis_2_min_medium_dQ_dx{0};
        float vis_2_min_length{0};
        float vis_2_sg_length{0};
        float vis_2_max_angle{0};
        float vis_2_max_weak_track{0};
        float vis_2_flag{1};
        float vis_flag{1};

        // ---- bad reconstruction ---------------------------------------- //
        float br_filled{0};
        float br1_flag{1};
        float br1_1_flag{1};
        float br1_1_shower_type{0};
        float br1_1_vtx_n_segs{0};
        float br1_1_energy{0};
        float br1_1_n_segs{0};
        float br1_1_flag_sg_topology{0};
        float br1_1_flag_sg_trajectory{0};
        float br1_1_sg_length{0};
        float br1_2_flag{1};
        float br1_2_energy{0};
        float br1_2_n_connected{0};
        float br1_2_max_length{0};
        float br1_2_n_connected_1{0};
        float br1_2_vtx_n_segs{0};
        float br1_2_n_shower_segs{0};
        float br1_2_max_length_ratio{1};
        float br1_2_shower_length{0};
        float br1_3_flag{1};
        float br1_3_energy{0};
        float br1_3_n_connected_p{0};
        float br1_3_max_length_p{0};
        float br1_3_n_shower_segs{0};
        float br1_3_flag_sg_topology{0};
        float br1_3_flag_sg_trajectory{0};
        float br1_3_n_shower_main_segs{0};
        float br1_3_sg_length{0};

        float br2_flag{1};
        float br2_flag_single_shower{0};
        float br2_num_valid_tracks{0};
        float br2_energy{0};
        float br2_angle1{0};
        float br2_angle2{0};
        float br2_angle{0};
        float br2_angle3{0};
        float br2_n_shower_main_segs{0};
        float br2_max_angle{0};
        float br2_sg_length{0};
        float br2_flag_sg_trajectory{0};

        float br3_1_energy{0};
        float br3_1_n_shower_segments{0};
        float br3_1_sg_flag_trajectory{0};
        float br3_1_sg_direct_length{0};
        float br3_1_sg_length{0};
        float br3_1_total_main_length{0};
        float br3_1_total_length{0};
        float br3_1_iso_angle{0};
        float br3_1_sg_flag_topology{0};
        float br3_1_flag{1};
        float br3_2_n_ele{0};
        float br3_2_n_other{0};
        float br3_2_energy{0};
        float br3_2_total_main_length{0};
        float br3_2_total_length{0};
        float br3_2_other_fid{0};
        float br3_2_flag{1};
        std::vector<float> br3_3_v_energy;
        std::vector<float> br3_3_v_angle;
        std::vector<float> br3_3_v_dir_length;
        std::vector<float> br3_3_v_length;
        std::vector<float> br3_3_v_flag;
        float br3_4_acc_length{0};
        float br3_4_total_length{0};
        float br3_4_energy{0};
        float br3_4_flag{1};
        std::vector<float> br3_5_v_dir_length;
        std::vector<float> br3_5_v_total_length;
        std::vector<float> br3_5_v_flag_avoid_muon_check;
        std::vector<float> br3_5_v_n_seg;
        std::vector<float> br3_5_v_angle;
        std::vector<float> br3_5_v_sg_length;
        std::vector<float> br3_5_v_energy;
        std::vector<float> br3_5_v_n_main_segs;
        std::vector<float> br3_5_v_n_segs;
        std::vector<float> br3_5_v_shower_main_length;
        std::vector<float> br3_5_v_shower_total_length;
        std::vector<float> br3_5_v_flag;
        std::vector<float> br3_6_v_angle;
        std::vector<float> br3_6_v_angle1;
        std::vector<float> br3_6_v_flag_shower_trajectory;
        std::vector<float> br3_6_v_direct_length;
        std::vector<float> br3_6_v_length;
        std::vector<float> br3_6_v_n_other_vtx_segs;
        std::vector<float> br3_6_v_energy;
        std::vector<float> br3_6_v_flag;
        float br3_7_energy{0};
        float br3_7_min_angle{0};
        float br3_7_sg_length{0};
        float br3_7_shower_main_length{0};
        float br3_7_flag{1};
        float br3_8_max_dQ_dx{0};
        float br3_8_energy{0};
        float br3_8_n_main_segs{0};
        float br3_8_shower_main_length{0};
        float br3_8_shower_length{0};
        float br3_8_flag{1};
        float br3_flag{1};

        float br4_1_shower_main_length{0};
        float br4_1_shower_total_length{0};
        float br4_1_min_dis{0};
        float br4_1_energy{0};
        float br4_1_flag_avoid_muon_check{0};
        float br4_1_n_vtx_segs{0};
        float br4_1_n_main_segs{0};
        float br4_1_flag{1};
        float br4_2_ratio_45{1};
        float br4_2_ratio_35{1};
        float br4_2_ratio_25{1};
        float br4_2_ratio_15{1};
        float br4_2_energy{0};
        float br4_2_ratio1_45{1};
        float br4_2_ratio1_35{1};
        float br4_2_ratio1_25{1};
        float br4_2_ratio1_15{1};
        float br4_2_iso_angle{0};
        float br4_2_iso_angle1{0};
        float br4_2_angle{0};
        float br4_2_flag{1};
        float br4_flag{1};

        // ---- track overclustering -------------------------------------- //
        std::vector<float> tro_1_v_particle_type;
        std::vector<float> tro_1_v_flag_dir_weak;
        std::vector<float> tro_1_v_min_dis;
        std::vector<float> tro_1_v_sg1_length;
        std::vector<float> tro_1_v_shower_main_length;
        std::vector<float> tro_1_v_max_n_vtx_segs;
        std::vector<float> tro_1_v_tmp_length;
        std::vector<float> tro_1_v_medium_dQ_dx;
        std::vector<float> tro_1_v_dQ_dx_cut;
        std::vector<float> tro_1_v_flag_shower_topology;
        std::vector<float> tro_1_v_flag;
        std::vector<float> tro_2_v_energy;
        std::vector<float> tro_2_v_stem_length;
        std::vector<float> tro_2_v_iso_angle;
        std::vector<float> tro_2_v_max_length;
        std::vector<float> tro_2_v_angle;
        std::vector<float> tro_2_v_flag;
        float tro_3_stem_length{0};
        float tro_3_n_muon_segs{0};
        float tro_3_energy{0};
        float tro_3_flag{1};
        std::vector<float> tro_4_v_dir2_mag;
        std::vector<float> tro_4_v_angle;
        std::vector<float> tro_4_v_angle1;
        std::vector<float> tro_4_v_angle2;
        std::vector<float> tro_4_v_length;
        std::vector<float> tro_4_v_length1;
        std::vector<float> tro_4_v_medium_dQ_dx;
        std::vector<float> tro_4_v_end_dQ_dx;
        std::vector<float> tro_4_v_energy;
        std::vector<float> tro_4_v_shower_main_length;
        std::vector<float> tro_4_v_flag_shower_trajectory;
        std::vector<float> tro_4_v_flag;
        std::vector<float> tro_5_v_max_angle;
        std::vector<float> tro_5_v_min_angle;
        std::vector<float> tro_5_v_max_length;
        std::vector<float> tro_5_v_iso_angle;
        std::vector<float> tro_5_v_n_vtx_segs;
        std::vector<float> tro_5_v_min_count;
        std::vector<float> tro_5_v_max_count;
        std::vector<float> tro_5_v_energy;
        std::vector<float> tro_5_v_flag;
        float tro_flag{1};

        // ---- high/low energy overlap ----------------------------------- //
        float hol_1_n_valid_tracks{0};
        float hol_1_min_angle{0};
        float hol_1_energy{0};
        float hol_1_flag_all_shower{0};
        float hol_1_min_length{0};
        float hol_1_flag{1};
        float hol_2_min_angle{0};
        float hol_2_medium_dQ_dx{1};
        float hol_2_ncount{0};
        float hol_2_energy{0};
        float hol_2_flag{1};
        float hol_flag{1};

        float lol_flag{1};
        std::vector<float> lol_1_v_energy;
        std::vector<float> lol_1_v_vtx_n_segs;
        std::vector<float> lol_1_v_nseg;
        std::vector<float> lol_1_v_angle;
        std::vector<float> lol_1_v_flag;
        std::vector<float> lol_2_v_length;
        std::vector<float> lol_2_v_angle;
        std::vector<float> lol_2_v_type;
        std::vector<float> lol_2_v_vtx_n_segs;
        std::vector<float> lol_2_v_energy;
        std::vector<float> lol_2_v_shower_main_length;
        std::vector<float> lol_2_v_flag_dir_weak;
        std::vector<float> lol_2_v_flag;
        float lol_3_angle_beam{0};
        float lol_3_n_valid_tracks{0};
        float lol_3_min_angle{0};
        float lol_3_vtx_n_segs{0};
        float lol_3_energy{0};
        float lol_3_shower_main_length{0};
        float lol_3_n_out{0};
        float lol_3_n_sum{0};
        float lol_3_flag{1};

        // ---- cosmic tagger flags --------------------------------------- //
        float cosmict_flag_1{0};
        float cosmict_flag_2{0};
        float cosmict_flag_3{0};
        float cosmict_flag_4{0};
        float cosmict_flag_5{0};
        float cosmict_flag_6{0};
        float cosmict_flag_7{0};
        float cosmict_flag_8{0};
        float cosmict_flag_9{0};
        std::vector<float> cosmict_flag_10;
        float cosmict_flag{0};

        // single muon (cosmict_2)
        float cosmict_2_filled{0};
        float cosmict_2_particle_type{0};
        float cosmict_2_n_muon_tracks{0};
        float cosmict_2_total_shower_length{0};
        float cosmict_2_flag_inside{1};
        float cosmict_2_angle_beam{0};
        float cosmict_2_flag_dir_weak{0};
        float cosmict_2_dQ_dx_end{0};
        float cosmict_2_dQ_dx_front{0};
        float cosmict_2_theta{0};
        float cosmict_2_phi{0};
        float cosmict_2_valid_tracks{0};

        // single muon long (cosmict_3)
        float cosmict_3_filled{0};
        float cosmict_3_flag_inside{0};
        float cosmict_3_angle_beam{0};
        float cosmict_3_flag_dir_weak{0};
        float cosmict_3_dQ_dx_end{0};
        float cosmict_3_dQ_dx_front{0};
        float cosmict_3_theta{0};
        float cosmict_3_phi{0};
        float cosmict_3_valid_tracks{0};

        // kinematics muon (cosmict_4)
        float cosmict_4_filled{0};
        float cosmict_4_flag_inside{0};
        float cosmict_4_angle_beam{0};
        float cosmict_4_connected_showers{0};

        // kinematics muon long (cosmict_5)
        float cosmict_5_filled{0};
        float cosmict_5_flag_inside{0};
        float cosmict_5_angle_beam{0};
        float cosmict_5_connected_showers{0};

        // special (cosmict_6)
        float cosmict_6_filled{0};
        float cosmict_6_flag_dir_weak{0};
        float cosmict_6_flag_inside{1};
        float cosmict_6_angle{0};

        // muon + michel (cosmict_7)
        float cosmict_7_filled{0};
        float cosmict_7_flag_sec{0};
        float cosmict_7_n_muon_tracks{0};
        float cosmict_7_total_shower_length{0};
        float cosmict_7_flag_inside{1};
        float cosmict_7_angle_beam{0};
        float cosmict_7_flag_dir_weak{0};
        float cosmict_7_dQ_dx_end{0};
        float cosmict_7_dQ_dx_front{0};
        float cosmict_7_theta{0};
        float cosmict_7_phi{0};

        // muon + michel + special (cosmict_8)
        float cosmict_8_filled{0};
        float cosmict_8_flag_out{0};
        float cosmict_8_muon_length{0};
        float cosmict_8_acc_length{0};

        // front upstream dirt (cosmict_10 — vector)
        std::vector<float> cosmict_10_flag_inside;
        std::vector<float> cosmict_10_vtx_z;
        std::vector<float> cosmict_10_flag_shower;
        std::vector<float> cosmict_10_flag_dir_weak;
        std::vector<float> cosmict_10_angle_beam;
        std::vector<float> cosmict_10_length;

        // ---- numu vs NC tagger ----------------------------------------- //
        float numu_cc_flag{0};
        std::vector<float> numu_cc_flag_1;
        std::vector<float> numu_cc_1_particle_type;
        std::vector<float> numu_cc_1_length;
        std::vector<float> numu_cc_1_medium_dQ_dx;
        std::vector<float> numu_cc_1_dQ_dx_cut;
        std::vector<float> numu_cc_1_direct_length;
        std::vector<float> numu_cc_1_n_daughter_tracks;
        std::vector<float> numu_cc_1_n_daughter_all;
        std::vector<float> numu_cc_flag_2;
        std::vector<float> numu_cc_2_length;
        std::vector<float> numu_cc_2_total_length;
        std::vector<float> numu_cc_2_n_daughter_tracks;
        std::vector<float> numu_cc_2_n_daughter_all;
        float numu_cc_flag_3{0};
        float numu_cc_3_particle_type{0};
        float numu_cc_3_max_length{0};
        float numu_cc_3_acc_track_length{0};
        float numu_cc_3_max_length_all{0};
        float numu_cc_3_max_muon_length{0};
        float numu_cc_3_n_daughter_tracks{0};
        float numu_cc_3_n_daughter_all{0};

        // ---- fiducial-volume flag (placeholder: filled by TaggerCheckNeutrino) //
        float match_isFC{0};

        // ---- numu BDT scores ------------------------------------------- //
        float cosmict_2_4_score{0};
        float cosmict_3_5_score{0};
        float cosmict_6_score{0};
        float cosmict_7_score{0};
        float cosmict_8_score{0};
        float cosmict_10_score{0};
        float numu_1_score{0};
        float numu_2_score{0};
        float numu_3_score{0};
        float cosmict_score{0};
        float numu_score{0};

        // ---- nue BDT scores -------------------------------------------- //
        float mipid_score{0};
        float gap_score{0};
        float hol_lol_score{0};
        float cme_anc_score{0};
        float mgo_mgt_score{0};
        float br1_score{0};
        float br3_score{0};
        float br3_3_score{0};
        float br3_5_score{0};
        float br3_6_score{0};
        float stemdir_br2_score{0};
        float trimuon_score{0};
        float br4_tro_score{0};
        float mipquality_score{0};
        float pio_1_score{0};
        float pio_2_score{0};
        float stw_spt_score{0};
        float vis_1_score{0};
        float vis_2_score{0};
        float stw_2_score{0};
        float stw_3_score{0};
        float stw_4_score{0};
        float sig_1_score{0};
        float sig_2_score{0};
        float lol_1_score{0};
        float lol_2_score{0};
        float tro_1_score{0};
        float tro_2_score{0};
        float tro_4_score{0};
        float tro_5_score{0};
        float nue_score{0};

        float photon_flag{0};
    };

} // namespace WireCell::Clus::PR

#endif // WIRECELL_CLUS_NEUTRINOTAGGERINFO_H
