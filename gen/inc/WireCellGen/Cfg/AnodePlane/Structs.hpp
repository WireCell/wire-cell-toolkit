/*
 * This file is 100% generated.  Any manual edits will likely be lost.
 *
 * This contains struct and other type definitions for shema in 
 * namespace WireCellGen::Cfg::AnodePlane.
 */
#ifndef WIRECELLGEN_CFG_ANODEPLANE_STRUCTS_HPP
#define WIRECELLGEN_CFG_ANODEPLANE_STRUCTS_HPP

#include <cstdint>
#include "WireCellUtil/Cfg/Base/Structs.hpp"

#include <vector>

namespace WireCellGen::Cfg::AnodePlane {

    // @brief Interesting drift positions on one anode face
    struct Face {

        // @brief Active volume boundary along drift direction near wires
        WireCellUtil::Cfg::Base::Distance anode = 0.0;

        // @brief Position where drift simulation ends and field response calcualtion begins
        WireCellUtil::Cfg::Base::Distance response = 0.0;

        // @brief Active volume boundary along drift direction near cathode
        WireCellUtil::Cfg::Base::Distance cathode = 0.0;
    };

    // @brief One or two faces
    using Faces = std::vector<WireCellGen::Cfg::AnodePlane::Face>;

    // @brief Configuration for AnodePlane
    struct Config {

        // @brief Identifier of wire plane as used in wire file
        WireCellUtil::Cfg::Base::Count ident = 0;

        // @brief Type name of IWireSchema component
        WireCellUtil::Cfg::Base::TypeName wire_schema = "";

        // @brief Number of impact positions per wire
        WireCellUtil::Cfg::Base::Count nimpacts = 10;

        // @brief Description of each face
        Faces faces = {};
    };

} // namespace WireCellGen::Cfg::AnodePlane

#endif // WIRECELLGEN_CFG_ANODEPLANE_STRUCTS_HPP