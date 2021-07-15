/*
 * This file is 100% generated.  Any manual edits will likely be lost.
 *
 * This contains functions struct and other type definitions for shema in 
 * namespace WireCellGen::Cfg::AnodePlane to be serialized via nlohmann::json.
 */
#ifndef WIRECELLGEN_CFG_ANODEPLANE_NLJS_HPP
#define WIRECELLGEN_CFG_ANODEPLANE_NLJS_HPP

// My structs
#include "WireCellGen/Cfg/AnodePlane/Structs.hpp"

// Nljs for externally referenced schema
#include "WireCellUtil/Cfg/Base/Nljs.hpp"

#include <nlohmann/json.hpp>

namespace WireCellGen::Cfg::AnodePlane {

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
    
} // namespace WireCellGen::Cfg::AnodePlane

#endif // WIRECELLGEN_CFG_ANODEPLANE_NLJS_HPP