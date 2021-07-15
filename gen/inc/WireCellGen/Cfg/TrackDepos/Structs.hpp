/*
 * This file is 100% generated.  Any manual edits will likely be lost.
 *
 * This contains struct and other type definitions for shema in 
 * namespace WireCellGen::Cfg::TrackDepos.
 */
#ifndef WIRECELLGEN_CFG_TRACKDEPOS_STRUCTS_HPP
#define WIRECELLGEN_CFG_TRACKDEPOS_STRUCTS_HPP

#include <cstdint>
#include "WireCellUtil/Cfg/Base/Structs.hpp"

#include <vector>

namespace WireCellGen::Cfg::TrackDepos {

    // @brief 
    struct Track {

        // @brief Absolute time at the start of the track
        WireCellUtil::Cfg::Base::Time time = 0.0;

        // @brief If negative, number of electrons per depo, else electrons per track
        WireCellUtil::Cfg::Base::Charge charge = -1.0;

        // @brief The ray defining the track endpoints
        WireCellUtil::Cfg::Base::Ray ray = {};
    };

    // @brief A sequence of tracks
    using Tracks = std::vector<WireCellGen::Cfg::TrackDepos::Track>;

    // @brief Configuration for TrackDepos component
    struct Config {

        // @brief Distance along track between two neighboring depos.
        WireCellUtil::Cfg::Base::Distance step_size = 1.0;

        // @brief Fraction of speed of light at which track progresses
        WireCellUtil::Cfg::Base::Normalized clight = 1.0;

        // @brief If positive, chunk the depos into groups spaning this amount of time with an EOS delimiting each group.  O.w. all depos are sent out as a stream.
        WireCellUtil::Cfg::Base::Time group_time = -1.0;

        // @brief Description of tracks on which to generate depos.
        Tracks tracks = {};
    };

} // namespace WireCellGen::Cfg::TrackDepos

#endif // WIRECELLGEN_CFG_TRACKDEPOS_STRUCTS_HPP