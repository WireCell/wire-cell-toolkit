/*
 * This file is 100% generated.  Any manual edits will likely be lost.
 *
 * This contains struct and other type definitions for shema in 
 * namespace WireCellUtil::Cfg::Base.
 */
#ifndef WIRECELLUTIL_CFG_BASE_STRUCTS_HPP
#define WIRECELLUTIL_CFG_BASE_STRUCTS_HPP

#include <cstdint>

#include <vector>
#include <string>

namespace WireCellUtil::Cfg::Base {

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
    struct Point {

        // @brief X coordinate
        Distance x = 0.0;

        // @brief Y coordinate
        Distance y = 0.0;

        // @brief Z coordinate
        Distance z = 0.0;
    };

    // @brief A directed line segment in 3-space
    struct Ray {

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
    using Tags = std::vector<WireCellUtil::Cfg::Base::Tag>;

    // @brief A sampling period
    using Tick = double;


    // @brief A temporal duration
    using Time = double;


    // @brief A component type name identifier
    using TypeName = std::string;

    // @brief A sequence of component type name identifiers
    using TypeNames = std::vector<WireCellUtil::Cfg::Base::TypeName>;

} // namespace WireCellUtil::Cfg::Base

#endif // WIRECELLUTIL_CFG_BASE_STRUCTS_HPP