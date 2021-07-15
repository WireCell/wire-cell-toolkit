/*
 * This file is 100% generated.  Any manual edits will likely be lost.
 *
 * This contains functions struct and other type definitions for shema in 
 * namespace WireCellGen::Cfg::Ductor to be serialized via nlohmann::json.
 */
#ifndef WIRECELLGEN_CFG_DUCTOR_NLJS_HPP
#define WIRECELLGEN_CFG_DUCTOR_NLJS_HPP

// My structs
#include "WireCellGen/Cfg/Ductor/Structs.hpp"

// Nljs for externally referenced schema
#include "WireCellUtil/Cfg/Base/Nljs.hpp"

#include <nlohmann/json.hpp>

namespace WireCellGen::Cfg::Ductor {

    using data_t = nlohmann::json;
    
    inline void to_json(data_t& j, const Config& obj) {
        j["nsigma"] = obj.nsigma;
        j["fluctuate"] = obj.fluctuate;
        j["start_time"] = obj.start_time;
        j["readout_time"] = obj.readout_time;
        j["tick"] = obj.tick;
        j["continuous"] = obj.continuous;
        j["fixed"] = obj.fixed;
        j["drift_speed"] = obj.drift_speed;
        j["first_frame_number"] = obj.first_frame_number;
        j["anode"] = obj.anode;
        j["rng"] = obj.rng;
        j["pirs"] = obj.pirs;
    }
    
    inline void from_json(const data_t& j, Config& obj) {
        if (j.contains("nsigma"))
            j.at("nsigma").get_to(obj.nsigma);    
        if (j.contains("fluctuate"))
            j.at("fluctuate").get_to(obj.fluctuate);    
        if (j.contains("start_time"))
            j.at("start_time").get_to(obj.start_time);    
        if (j.contains("readout_time"))
            j.at("readout_time").get_to(obj.readout_time);    
        if (j.contains("tick"))
            j.at("tick").get_to(obj.tick);    
        if (j.contains("continuous"))
            j.at("continuous").get_to(obj.continuous);    
        if (j.contains("fixed"))
            j.at("fixed").get_to(obj.fixed);    
        if (j.contains("drift_speed"))
            j.at("drift_speed").get_to(obj.drift_speed);    
        if (j.contains("first_frame_number"))
            j.at("first_frame_number").get_to(obj.first_frame_number);    
        if (j.contains("anode"))
            j.at("anode").get_to(obj.anode);    
        if (j.contains("rng"))
            j.at("rng").get_to(obj.rng);    
        if (j.contains("pirs"))
            j.at("pirs").get_to(obj.pirs);    
    }
    
} // namespace WireCellGen::Cfg::Ductor

#endif // WIRECELLGEN_CFG_DUCTOR_NLJS_HPP