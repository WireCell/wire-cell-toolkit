#include "WireCellClus/Facade_Summary.h"
#include "WireCellUtil/PointSummary.h"

using namespace WireCell;
using WireCell::PointCloud::json_summary;

Configuration
PointCloud::Facade::json_summary(const PointCloud::Facade::Grouping& grp)
{
    Configuration ret;
    ret["type"] = "Grouping";
    ret["hash"] = grp.hash();
    ret["nproj_centers"] = (int)grp.proj_centers().size();
    ret["npitch_mags"] = (int)grp.pitch_mags().size();
    ret["ndead_winds"] = (int)grp.all_dead_winds().size();
    ret["value"] = WireCell::PointCloud::json_summary(grp.value(), false);
    for (const auto* cf : grp.children()) {
        ret["clusters"].append(json_summary(*cf));
    }
    return ret;
}

Configuration PointCloud::Facade::json_summary(const PointCloud::Facade::Cluster& cls) 
{
    Configuration ret;
    // this is too huge to be exhaustive
    ret["type"] = "Cluster";
    ret["hash"] = cls.hash();
    ret["length"] = cls.get_length();
    // ret["num_slices"] = cls.get_num_time_slices();
    ret["value"] = WireCell::PointCloud::json_summary(cls.value(), false);
    for (const auto* cf : cls.children()) {
        ret["clusters"].append(json_summary(*cf));
    }
    return ret;
}

Configuration PointCloud::Facade::json_summary(const PointCloud::Facade::Blob& blb) 
{
    Configuration ret;
    ret["type"] = "Blob";
    ret["hash"] = blb.hash();
    ret["face"] = blb.wpid().face();
    ret["npoints"] = blb.npoints();
    ret["charge"] = blb.charge();
    ret["center_x"] = blb.center_x();
    ret["center_y"] = blb.center_y();
    ret["center_z"] = blb.center_z();
    ret["slice_index_min"] = blb.slice_index_min();
    ret["slice_index_max"] = blb.slice_index_max();
    ret["u_wire_index_min"] = blb.u_wire_index_min();
    ret["u_wire_index_max"] = blb.u_wire_index_max();
    ret["v_wire_index_min"] = blb.v_wire_index_min();
    ret["v_wire_index_max"] = blb.v_wire_index_max();
    ret["w_wire_index_min"] = blb.w_wire_index_min();
    ret["w_wire_index_max"] = blb.w_wire_index_max();
    ret["min_wire_interval"] = blb.get_min_wire_interval();
    ret["max_wire_interval"] = blb.get_max_wire_interval();
    ret["min_wire_type"] = blb.get_min_wire_type();
    ret["max_wire_type"] = blb.get_max_wire_type();
    ret["value"] = WireCell::PointCloud::json_summary(blb.value());
    return ret;
}



