/*
 * This file is 100% generated.  Any manual edits will likely be lost.
 *
 * This contains configuration related types and methods for schema: 
 *
 *     WireCell.Gen.Cfg.TrackDepos
 *
 */
#ifndef WIRECELL_GEN_CFG_TRACKDEPOS
#define WIRECELL_GEN_CFG_TRACKDEPOS

#include "WireCellUtil/Cfg/Base.hpp"

#include <cstdint>
// for sequence
#include <vector>

#include <nlohmann/json.hpp>

namespace WireCell::Gen::Cfg::TrackDepos {

    //
    // Type definitions.
    //

    // @brief 
    struct Track 
    {

        // @brief Absolute time at the start of the track
        Util::Cfg::Base::Time time = 0.0;

        // @brief If negative, number of electrons per depo, else electrons per track
        Util::Cfg::Base::Charge charge = -1.0;

        // @brief The ray defining the track endpoints
        Util::Cfg::Base::Ray ray = {};
    };

    // @brief A sequence of tracks
    using Tracks = std::vector<WireCell::Gen::Cfg::TrackDepos::Track>;

    // @brief Configuration for TrackDepos component
    struct Config 
    {

        // @brief Distance along track between two neighboring depos.
        Util::Cfg::Base::Distance step_size = 1.0;

        // @brief Fraction of speed of light at which track progresses
        Util::Cfg::Base::Normalized clight = 1.0;

        // @brief If positive, chunk the depos into groups spaning this amount of time with an EOS delimiting each group.  O.w. all depos are sent out as a stream.
        Util::Cfg::Base::Time group_time = -1.0;

        // @brief Description of tracks on which to generate depos.
        Tracks tracks = {};
    };



    //
    // JSON serialization methods
    //

    using data_t = nlohmann::json;
    
    inline void to_json(data_t& j, const Track& obj) {
        j["time"] = obj.time;
        j["charge"] = obj.charge;
        j["ray"] = obj.ray;
    }
    
    inline void from_json(const data_t& j, Track& obj) {
        if (j.contains("time"))
            j.at("time").get_to(obj.time);    
        if (j.contains("charge"))
            j.at("charge").get_to(obj.charge);    
        if (j.contains("ray"))
            j.at("ray").get_to(obj.ray);    
    }
    
    inline void to_json(data_t& j, const Config& obj) {
        j["step_size"] = obj.step_size;
        j["clight"] = obj.clight;
        j["group_time"] = obj.group_time;
        j["tracks"] = obj.tracks;
    }
    
    inline void from_json(const data_t& j, Config& obj) {
        if (j.contains("step_size"))
            j.at("step_size").get_to(obj.step_size);    
        if (j.contains("clight"))
            j.at("clight").get_to(obj.clight);    
        if (j.contains("group_time"))
            j.at("group_time").get_to(obj.group_time);    
        if (j.contains("tracks"))
            j.at("tracks").get_to(obj.tracks);    
    }
} // WireCell::Gen::Cfg::TrackDepos

#endif // WIRECELL_GEN_CFG_TRACKDEPOS