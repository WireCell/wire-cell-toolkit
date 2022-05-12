/*
 * This file is 100% generated.  Any manual edits will likely be lost.
 *
 * This contains configuration related types and methods for schema: 
 *
 *     WireCell.Util.Cfg.Base
 *
 */
#ifndef WIRECELL_UTIL_CFG_BASE
#define WIRECELL_UTIL_CFG_BASE


#include <cstdint>
// for sequence
#include <vector>
// for string
#include <string>

#include <nlohmann/json.hpp>

namespace WireCell::Util::Cfg::Base {

    //
    // Type definitions.
    //

    // @brief Amount of charge per some unit
    using Charge = double;


    // @brief A simple counting number
    using Count = int32_t;


    // @brief A spacial distance
    using Distance = double;


    // @brief Something resembling a file system tree path name
    using Filename = std::string;

    // @brief A filename with NPZ extension
    using FilenameNPZ = std::string;

    // @brief A general boolean flag
    using Flag = bool;

    // @brief A code-friendly identifier
    using Ident = std::string;

    // @brief A real number in [0,1]
    using Normalized = double;


    // @brief A Cartesian point in 3-space.
    struct Point 
    {

        // @brief X coordinate
        Distance x = 0.0;

        // @brief Y coordinate
        Distance y = 0.0;

        // @brief Z coordinate
        Distance z = 0.0;
    };

    // @brief A directed line segment in 3-space
    struct Ray 
    {

        // @brief Start point
        Point tail = {};

        // @brief End point
        Point head = {};
    };

    // @brief A multiplicative scaling factor
    using Scaling = double;


    // @brief A speed in units distance per time
    using Speed = double;


    // @brief A simple identifying tag associated with some data
    using Tag = std::string;

    // @brief A sequence of tags
    using Tags = std::vector<WireCell::Util::Cfg::Base::Tag>;

    // @brief A sampling period
    using Tick = double;


    // @brief A temporal duration
    using Time = double;


    // @brief A component type name identifier
    using TypeName = std::string;

    // @brief A sequence of component type name identifiers
    using TypeNames = std::vector<WireCell::Util::Cfg::Base::TypeName>;



    //
    // JSON serialization methods
    //

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
} // WireCell::Util::Cfg::Base

#endif // WIRECELL_UTIL_CFG_BASE