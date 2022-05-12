/*
 * This file is 100% generated.  Any manual edits will likely be lost.
 *
 * This contains configuration related types and methods for schema: 
 *
 *     WireCell.Sio.Cfg.NumpyFrameSaver
 *
 */
#ifndef WIRECELL_SIO_CFG_NUMPYFRAMESAVER
#define WIRECELL_SIO_CFG_NUMPYFRAMESAVER

#include "WireCellUtil/Cfg/Base.hpp"

#include <cstdint>

#include <nlohmann/json.hpp>

namespace WireCell::Sio::Cfg::NumpyFrameSaver {

    //
    // Type definitions.
    //

    // @brief Configuration for NumpyFrameSaver component
    struct Config 
    {

        // @brief If true then truncate to 16 bit short int, o.w. save as 32 bit floats
        Util::Cfg::Base::Flag digitize = false;

        // @brief Initial value to which frame samples add
        Util::Cfg::Base::Charge baseline = 0.0;

        // @brief Multiplicative scale applied to each frame sample
        Util::Cfg::Base::Scaling scale = 1.0;

        // @brief Additive offset applied to scaled frame samples
        Util::Cfg::Base::Charge offset = 0.0;

        // @brief Frame tags to consider saving, assume all frames if empty
        Util::Cfg::Base::Tags frame_tags = {};

        // @brief Name of NPZ file to which frames are saved
        Util::Cfg::Base::FilenameNPZ filename = "wct-frame.npz";
    };



    //
    // JSON serialization methods
    //

    using data_t = nlohmann::json;
    
    inline void to_json(data_t& j, const Config& obj) {
        j["digitize"] = obj.digitize;
        j["baseline"] = obj.baseline;
        j["scale"] = obj.scale;
        j["offset"] = obj.offset;
        j["frame_tags"] = obj.frame_tags;
        j["filename"] = obj.filename;
    }
    
    inline void from_json(const data_t& j, Config& obj) {
        if (j.contains("digitize"))
            j.at("digitize").get_to(obj.digitize);    
        if (j.contains("baseline"))
            j.at("baseline").get_to(obj.baseline);    
        if (j.contains("scale"))
            j.at("scale").get_to(obj.scale);    
        if (j.contains("offset"))
            j.at("offset").get_to(obj.offset);    
        if (j.contains("frame_tags"))
            j.at("frame_tags").get_to(obj.frame_tags);    
        if (j.contains("filename"))
            j.at("filename").get_to(obj.filename);    
    }
} // WireCell::Sio::Cfg::NumpyFrameSaver

#endif // WIRECELL_SIO_CFG_NUMPYFRAMESAVER