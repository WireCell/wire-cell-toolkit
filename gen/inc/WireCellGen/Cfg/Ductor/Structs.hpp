/*
 * This file is 100% generated.  Any manual edits will likely be lost.
 *
 * This contains struct and other type definitions for shema in 
 * namespace WireCellGen::Cfg::Ductor.
 */
#ifndef WIRECELLGEN_CFG_DUCTOR_STRUCTS_HPP
#define WIRECELLGEN_CFG_DUCTOR_STRUCTS_HPP

#include <cstdint>
#include "WireCellUtil/Cfg/Base/Structs.hpp"


namespace WireCellGen::Cfg::Ductor {

    // @brief Ductor configuration
    struct Config {

        // @brief Depo truncation bound in number of Gaussian standard deviations
        WireCellUtil::Cfg::Base::Count nsigma = 3;

        // @brief Whether to fluctuate the final Gaussian deposition
        WireCellUtil::Cfg::Base::Flag fluctuate = true;

        // @brief The initial time for this ductor
        WireCellUtil::Cfg::Base::Time start_time = 0.0;

        // @brief The time span for each readout.
        WireCellUtil::Cfg::Base::Time readout_time = 5000000.0;

        // @brief The sampling period
        WireCellUtil::Cfg::Base::Tick tick = 500.0;

        // @brief Operate in streaming mode
        WireCellUtil::Cfg::Base::Flag continuous = true;

        // @brief Operate in fixed time window mode
        WireCellUtil::Cfg::Base::Flag fixed = false;

        // @brief The nominal speed of drifting electrons
        WireCellUtil::Cfg::Base::Speed drift_speed = 0.0;

        // @brief Initial value for frame count sequence
        WireCellUtil::Cfg::Base::Count first_frame_number = 0;

        // @brief Type name of IAnodePlane component
        WireCellUtil::Cfg::Base::TypeName anode = "AnodePlane";

        // @brief Type name of IRandom component
        WireCellUtil::Cfg::Base::TypeName rng = "Random";

        // @brief Sequence of type names for IPlaneImpactResponse components
        WireCellUtil::Cfg::Base::TypeNames pirs = {};
    };

} // namespace WireCellGen::Cfg::Ductor

#endif // WIRECELLGEN_CFG_DUCTOR_STRUCTS_HPP