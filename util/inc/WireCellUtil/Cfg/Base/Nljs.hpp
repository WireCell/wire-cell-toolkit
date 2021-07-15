/*
 * This file is 100% generated.  Any manual edits will likely be lost.
 *
 * This contains functions struct and other type definitions for shema in 
 * namespace WireCellUtil::Cfg::Base to be serialized via nlohmann::json.
 */
#ifndef WIRECELLUTIL_CFG_BASE_NLJS_HPP
#define WIRECELLUTIL_CFG_BASE_NLJS_HPP

// My structs
#include "WireCellUtil/Cfg/Base/Structs.hpp"


#include <nlohmann/json.hpp>

namespace WireCellUtil::Cfg::Base {

    using data_t = nlohmann::json;
    
    inline void to_json(data_t& j, const Point& obj) {
        j["x"] = obj.x;
        j["y"] = obj.y;
        j["z"] = obj.z;
    }
    
    inline void from_json(const data_t& j, Point& obj) {
        if (j.contains("x"))
            j.at("x").get_to(obj.x);    
        if (j.contains("y"))
            j.at("y").get_to(obj.y);    
        if (j.contains("z"))
            j.at("z").get_to(obj.z);    
    }
    
    inline void to_json(data_t& j, const Ray& obj) {
        j["tail"] = obj.tail;
        j["head"] = obj.head;
    }
    
    inline void from_json(const data_t& j, Ray& obj) {
        if (j.contains("tail"))
            j.at("tail").get_to(obj.tail);    
        if (j.contains("head"))
            j.at("head").get_to(obj.head);    
    }
    
} // namespace WireCellUtil::Cfg::Base

#endif // WIRECELLUTIL_CFG_BASE_NLJS_HPP