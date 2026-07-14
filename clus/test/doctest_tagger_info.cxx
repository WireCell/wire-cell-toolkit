#include "WireCellClus/NeutrinoTaggerInfo.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell::Clus::PR;

// Usage:
//   wcdoctest-clus -tc="TaggerInfo default initialization"
//   wcdoctest-clus -tc="KineInfo default initialization"

TEST_CASE("TaggerInfo default initialization")
{
    TaggerInfo ti{};

    // Top-level flags default to 1 (= "not yet classified as neutrino")
    CHECK(ti.cosmic_flag      == doctest::Approx(1.0f));
    CHECK(ti.gap_flag         == doctest::Approx(1.0f));
    CHECK(ti.mip_quality_flag == doctest::Approx(1.0f));
    CHECK(ti.mip_flag         == doctest::Approx(1.0f));

    // Filled flags start at 0
    CHECK(ti.cosmic_filled      == doctest::Approx(0.0f));
    CHECK(ti.gap_filled         == doctest::Approx(0.0f));
    CHECK(ti.mip_quality_filled == doctest::Approx(0.0f));
    CHECK(ti.mip_filled         == doctest::Approx(0.0f));

    // MIP counters with non-zero defaults
    CHECK(ti.mip_n_first_non_mip   == doctest::Approx(19.0f));
    CHECK(ti.mip_n_first_non_mip_1 == doctest::Approx(19.0f));
    CHECK(ti.mip_n_first_non_mip_2 == doctest::Approx(19.0f));
    CHECK(ti.mip_n_below_threshold == doctest::Approx(19.0f));
    CHECK(ti.mip_n_stem_size       == doctest::Approx(20.0f));
    CHECK(ti.mip_n_other_vertex    == doctest::Approx(2.0f));
    CHECK(ti.mip_n_vertex          == doctest::Approx(1.0f));
    CHECK(ti.mip_n_lowest          == doctest::Approx(1.0f));
    CHECK(ti.mip_n_highest         == doctest::Approx(1.0f));

    // MIP dQ/dx values with non-zero defaults
    CHECK(ti.mip_vec_dQ_dx_0       == doctest::Approx(1.0f));
    CHECK(ti.mip_vec_dQ_dx_1       == doctest::Approx(1.0f));
    CHECK(ti.mip_max_dQ_dx_sample  == doctest::Approx(1.0f));
    CHECK(ti.mip_lowest_dQ_dx      == doctest::Approx(1.0f));
    CHECK(ti.mip_highest_dQ_dx     == doctest::Approx(1.0f));
    CHECK(ti.mip_medium_dQ_dx      == doctest::Approx(1.0f));
    CHECK(ti.mip_stem_length       == doctest::Approx(1.0f));
    CHECK(ti.mip_length_main       == doctest::Approx(1.0f));
    CHECK(ti.mip_length_total      == doctest::Approx(1.0f));
    CHECK(ti.mip_min_dQ_dx_5       == doctest::Approx(1.0f));

    // MIP angles/positions default to 0
    CHECK(ti.mip_angle_beam == doctest::Approx(0.0f));
    CHECK(ti.mip_iso_angle  == doctest::Approx(0.0f));
    CHECK(ti.mip_min_dis    == doctest::Approx(0.0f));

    // MIP vec samples 2-19 default to 0
    CHECK(ti.mip_vec_dQ_dx_2  == doctest::Approx(0.0f));
    CHECK(ti.mip_vec_dQ_dx_10 == doctest::Approx(0.0f));
    CHECK(ti.mip_vec_dQ_dx_19 == doctest::Approx(0.0f));

    // SSM top-level counters are -999 (no SSM candidate found)
    CHECK(ti.ssm_Nsm         == doctest::Approx(-999.0f));
    CHECK(ti.ssm_Nsm_wivtx   == doctest::Approx(-999.0f));
    CHECK(ti.ssm_kine_energy == doctest::Approx(-999.0f));

    // SSM dQ/dx fields are -999
    CHECK(ti.ssm_dq_dx_fwd_1 == doctest::Approx(-999.0f));
    CHECK(ti.ssm_dq_dx_bck_1 == doctest::Approx(-999.0f));
    CHECK(ti.ssm_medium_dq_dx == doctest::Approx(-999.0f));

    // SSM flag_st_kdar defaults to 0
    CHECK(ti.ssm_flag_st_kdar == doctest::Approx(0.0f));

    // numu fields default to 0
    CHECK(ti.numu_cc_flag == doctest::Approx(0.0f));

    // match_isFC defaults to 0
    CHECK(ti.match_isFC == doctest::Approx(0.0f));
}


TEST_CASE("KineInfo default initialization")
{
    KineInfo ki{};

    // All scalar fields default to 0
    CHECK(ki.kine_nu_x_corr      == doctest::Approx(0.0f));
    CHECK(ki.kine_nu_y_corr      == doctest::Approx(0.0f));
    CHECK(ki.kine_nu_z_corr      == doctest::Approx(0.0f));
    CHECK(ki.kine_reco_Enu       == doctest::Approx(0.0f));
    CHECK(ki.kine_reco_add_energy == doctest::Approx(0.0f));
    CHECK(ki.kine_pio_mass       == doctest::Approx(0.0f));
    CHECK(ki.kine_pio_flag       == 0);
    CHECK(ki.kine_pio_vtx_dis    == doctest::Approx(0.0f));
    CHECK(ki.kine_pio_energy_1   == doctest::Approx(0.0f));
    CHECK(ki.kine_pio_energy_2   == doctest::Approx(0.0f));
    CHECK(ki.kine_pio_angle      == doctest::Approx(0.0f));

    // Per-particle vectors are empty
    CHECK(ki.kine_energy_particle.empty());
    CHECK(ki.kine_energy_info.empty());
    CHECK(ki.kine_particle_type.empty());
    CHECK(ki.kine_energy_included.empty());
}
