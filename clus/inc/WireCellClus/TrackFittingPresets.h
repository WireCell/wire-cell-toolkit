// TrackFittingPresets.h - Create this as a new header file

#ifndef WIRECELLCLUS_TRACKFITTING_PRESETS_H
#define WIRECELLCLUS_TRACKFITTING_PRESETS_H

#include "WireCellClus/TrackFitting.h"

namespace WireCell::Clus {

    /**
     * Factory class for creating TrackFitting instances with preset configurations
     */
    class TrackFittingPresets {
    public:
        
        /**
         * Create TrackFitting with your current hard-coded values
         * This gives you exactly the same behavior as before, but centralized
         */
        static TrackFitting create_with_current_values() {
            TrackFitting fitter;
            
            TrackFitting::Parameters params;
            
            // Set to exactly your current hard-coded values
            params.DL = 6.4* pow(units::cm,2)/units::second;                    // m²/s, longitudinal diffusion
            params.DT = 9.8* pow(units::cm,2)/units::second;                    // m²/s, transverse diffusion
            params.col_sigma_w_T = 0.188060 * 3*units::mm * 0.2; // Collection plane
            params.ind_sigma_u_T = 0.402993 * 3*units::mm * 0.3; // U induction plane
            params.ind_sigma_v_T = 0.402993 * 3*units::mm * 0.5; // V induction plane
            params.rel_uncer_ind = 0.075;          // Relative uncertainty for induction
            params.rel_uncer_col = 0.05;           // Relative uncertainty for collection
            params.add_uncer_ind = 0.0;            // Additional uncertainty for induction
            params.add_uncer_col = 300.0;          // Additional uncertainty for collection
            params.add_sigma_L = 1.428249 *0.5505*units::mm/ 0.5;   // Longitudinal filter effects
            
            // Additional useful parameters for charge err estimation ...
            params.rel_charge_uncer = 0.1; // 10%
            params.add_charge_uncer = 600; // electrons

            params.default_charge_th = 100;
            params.default_charge_err = 1000;

            params.scaling_quality_th = 0.5;
            params.scaling_ratio = 0.05;

            params.area_ratio1 = 1.8*units::mm;
            params.area_ratio2 = 1.7;

            params.skip_default_ratio_1 = 0.25;
            params.skip_ratio_cut = 0.97;
            params.skip_ratio_1_cut = 0.75;

            params.skip_angle_cut_1 = 160;
            params.skip_angle_cut_2 = 90;
            params.skip_angle_cut_3 = 45;
            params.skip_dis_cut = 0.5*units::cm;

            params.default_dQ_dx = 5000; // electrons

            params.end_point_factor = 0.6;
            params.mid_point_factor = 0.9;
            params.nlevel = 3;
            params.charge_cut = 2000;

            // Distance parameters (add these if you use them in your methods)
            params.low_dis_limit = 1.2*units::cm;            // cm
            params.end_point_limit = 0.6*units::cm;          // cm
            params.time_tick_cut = 20;            // time tick
            
            fitter.set_parameters(params);
            return fitter;
        }
        
        /**
         * Create TrackFitting with MicroBooNE-like parameters
         */
        static TrackFitting create_microboone() {
            return create_with_current_values(); // Same as your current values for now
        }
        
        /**
         * Create TrackFitting with alternative parameter set (for experimentation)
         */
        static TrackFitting create_alternative() {
            TrackFitting fitter = create_with_current_values();
            
            // Modify a few parameters for testing
            fitter.set_parameter("DL", 7.0e-6);           // Slightly higher diffusion
            fitter.set_parameter("min_track_length", 3.0); // Shorter minimum track
            fitter.set_parameter("rel_uncer_ind", 0.08);   // Higher uncertainty
            
            return fitter;
        }
        
        /**
         * Create TrackFitting with specific fitting type and current parameters
         */
        static TrackFitting create_with_fitting_type(TrackFitting::FittingType type) {
            TrackFitting fitter = create_with_current_values();
            fitter.set_fitting_type(type);
            return fitter;
        }
    };

} // namespace WireCell::Clus

#endif // WIRECELLCLUS_TRACKFITTING_PRESETS_H