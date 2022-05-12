/*
 * This file is 100% generated.  Any manual edits will likely be lost.
 *
 * This contains configuration related types and methods for schema: 
 *
 *     WireCell.Gen.Cfg.AnodePlane
 *
 */
#ifndef WIRECELL_GEN_CFG_ANODEPLANE
#define WIRECELL_GEN_CFG_ANODEPLANE

#include "WireCellUtil/Cfg/Base.hpp"

#include <cstdint>
// for sequence
#include <vector>

#include <nlohmann/json.hpp>

namespace WireCell::Gen::Cfg::AnodePlane {

    //
    // Type definitions.
    //

    // @brief Interesting drift positions on one anode face
    struct Face 
    {

        // @brief Active volume boundary along drift direction near wires
        Util::Cfg::Base::Distance anode = 0.0;

        // @brief Position where drift simulation ends and field response calcualtion begins
        Util::Cfg::Base::Distance response = 0.0;

        // @brief Active volume boundary along drift direction near cathode
        Util::Cfg::Base::Distance cathode = 0.0;
    };

    // @brief One or two faces
    using Faces = std::vector<WireCell::Gen::Cfg::AnodePlane::Face>;

    // @brief Configuration for AnodePlane
    struct Config 
    {

        // @brief Identifier of wire plane as used in wire file
        Util::Cfg::Base::Count ident = 0;

        // @brief Type name of IWireSchema component
        Util::Cfg::Base::TypeName wire_schema = "";

        // @brief Number of impact positions per wire
        Util::Cfg::Base::Count nimpacts = 10;

        // @brief Description of each face
        Faces faces = {};
    };



    //
    // JSON serialization methods
    //

    using data_t = nlohmann::json;
    
    inline void to_json(data_t& j, const Face& obj) {
        j["anode"] = obj.anode;
        j["response"] = obj.response;
        j["cathode"] = obj.cathode;
    }
    
    inline void from_json(const data_t& j, Face& obj) {
        if (j.contains("anode"))
            j.at("anode").get_to(obj.anode);    
        if (j.contains("response"))
            j.at("response").get_to(obj.response);    
        if (j.contains("cathode"))
            j.at("cathode").get_to(obj.cathode);    
    }
    
    inline void to_json(data_t& j, const Config& obj) {
        j["ident"] = obj.ident;
        j["wire_schema"] = obj.wire_schema;
        j["nimpacts"] = obj.nimpacts;
        j["faces"] = obj.faces;
    }
    
    inline void from_json(const data_t& j, Config& obj) {
        if (j.contains("ident"))
            j.at("ident").get_to(obj.ident);    
        if (j.contains("wire_schema"))
            j.at("wire_schema").get_to(obj.wire_schema);    
        if (j.contains("nimpacts"))
            j.at("nimpacts").get_to(obj.nimpacts);    
        if (j.contains("faces"))
            j.at("faces").get_to(obj.faces);    
    }
} // WireCell::Gen::Cfg::AnodePlane

#endif // WIRECELL_GEN_CFG_ANODEPLANE