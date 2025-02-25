#include "WireCellUtil/PointSummary.h"

using namespace WireCell;

Configuration PointCloud::json_summary(const PointCloud::Tree::Points& value, bool recur)
{
    Configuration ret;
    ret["type"] = "Points";
    for (const auto& [name, pc] : value.local_pcs()) {
        ret["pcs"][name] = json_summary(pc);
    }
    if (!recur) return ret;

    for (const auto* child : value.node()->children()) {
        const auto& cval = child->value;
        ret["children"].append(json_summary(cval));
    }
    return ret;    
}
Configuration PointCloud::json_summary(const PointCloud::Dataset& ds)
{
    Configuration ret;
    ret["type"] = "Dataset";
    ret["size_major"] = (int)ds.size_major();
    auto keys = ds.keys();
    ret["metadata"] = ds.metadata();

    for (auto key : ds.keys()) {
        auto arr = ds.get(key);
        ret["arrays"][key] = json_summary(*(arr.get()));
    }
    return ret;
}
Configuration PointCloud::json_summary(const PointCloud::Array& arr)
{
    Configuration ret;
    ret["type"] = "Array";
    ret["dtype"] = arr.dtype();
    ret["size_major"] = (int)arr.size_major();
    ret["num_elements"] = (int) arr.num_elements();
    for (const auto& s : arr.shape()) {
        ret["shape"].append((int) s);
    }
    ret["metadata"] = arr.metadata();
    return ret;
}
