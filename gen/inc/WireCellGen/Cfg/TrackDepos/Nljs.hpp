/*
 * This file is 100% generated.  Any manual edits will likely be lost.
 *
 * This contains functions struct and other type definitions for shema in 
 * namespace WireCellGen::Cfg::TrackDepos to be serialized via nlohmann::json.
 */
#ifndef WIRECELLGEN_CFG_TRACKDEPOS_NLJS_HPP
#define WIRECELLGEN_CFG_TRACKDEPOS_NLJS_HPP

// My structs
#include "WireCellGen/Cfg/TrackDepos/Structs.hpp"

// Nljs for externally referenced schema
#include "WireCellUtil/Cfg/Base/Nljs.hpp"

#include <nlohmann/json.hpp>

namespace WireCellGen::Cfg::TrackDepos {

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
    
} // namespace WireCellGen::Cfg::TrackDepos

#endif // WIRECELLGEN_CFG_TRACKDEPOS_NLJS_HPP