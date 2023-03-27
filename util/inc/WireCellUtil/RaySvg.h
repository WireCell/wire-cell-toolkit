/*
  RaySvg provides functions and objects to convert C++ objects related
  to RayGrid into SVG.  RayGrid objects include includes 3D objects of
  type Point, Ray (wire center or region or Activity) and Blob
  (corners, bounding rays).

 */

#ifndef WIRECELLUTIL_RAYSVG
#define WIRECELLUTIL_RAYSVG

#include "WireCellUtil/RayTiling.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/WireSchema.h"
#include "WireCellUtil/Units.h"

#include "WireCellUtil/svggpp.hpp"

#include <string>
#include <vector>

namespace WireCell::RaySvg {
        
    // Define the projection from 3D space to 2D space.  Avoid
    // performing the projection inline and instead call this.
    using svggpp::point_t;
    point_t project(const Point& pt);

    // Return a view from a BoundingBox.bounds() Ray.
    svggpp::view_t project(const Ray& bbray);

    using wire_vector_t = std::vector<WireSchema::Wire>;

    struct Geom {
        RayGrid::Coordinates coords;
        std::vector<wire_vector_t> wires;
    };

    /// The following functions render their arguments to an <svg>
    /// represented as an object.  

    /// A full rendering returning an <svg>. 
    svggpp::xml_t svg_dual(const Geom& geom, const RayGrid::activities_t& acts, const RayGrid::blobs_t& blobs);
    svggpp::xml_t svg_full(const Geom& geom, const RayGrid::activities_t& acts, const RayGrid::blobs_t& blobs);

    // Return <svg> of active area.  The view port is empty.  The
    // viewbox reflects the user coordinates.
    svggpp::xml_t g_active_area(const Geom& geom);

    // Return <g> of wire plane activities as strips.
    svggpp::xml_t g_active_strips(const Geom& geom, const RayGrid::activities_t& acts);

    // Return <g> of activities as individual wires.
    svggpp::xml_t g_active_wires(const Geom& geom, const RayGrid::activities_t& acts);

    // Return <g> of blobs represented by their corners
    svggpp::xml_t g_blob_corners(const Geom& geom, const RayGrid::blobs_t& blobs);

    // Return <g> of blobs represented by areas
    svggpp::xml_t g_blob_areas(const Geom& geom, const RayGrid::blobs_t& blobs);

    // Return <g> of activities associated to blobs
    svggpp::xml_t g_blob_activites(const Geom& geom, const RayGrid::blobs_t& blobs);

};

#endif
