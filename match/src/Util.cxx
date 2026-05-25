#include "WireCellMatch/Util.h"

#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Logging.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <boost/multi_array.hpp>
#pragma GCC diagnostic pop

#include <fstream>
#include <memory>

using namespace WireCell;
using namespace WireCell::Match;
using namespace WireCell::Clus::Facade;
using namespace WireCell::PointCloud::Tree;

void WireCell::Match::dump_bee_3d(const Points::node_t& root, const std::string& fn)
{
    using spdlog::debug;
    using WireCell::Clus::Facade::float_t;
    using WireCell::Clus::Facade::int_t;

    Json::Value bee;
    bee["runNo"] = 0;
    bee["subRunNo"] = 0;
    bee["eventNo"] = 0;
    bee["geom"] = "uboone";
    bee["type"] = "cluster";

    std::vector<float_t> x, y, z, q;
    std::vector<int_t> cluster_id;
    int_t cid = 0;

    auto grouping = root.value.facade<Grouping>();
    std::vector<Cluster*> clusters = grouping->children();
    std::sort(clusters.begin(), clusters.end(),
              [](const Cluster* a, const Cluster* b) {
                  return a->get_length() > b->get_length();
              });

    for (const auto cluster : clusters) {
        const auto& sv = cluster->sv3d();
        const auto& spcs = sv.pcs();
        for (const auto& spc : spcs) {
            if (spc.get().get("x") == nullptr) {
                debug("No x in point cloud, skip");
                continue;
            }
            const auto& x_ = spc.get().get("x")->elements<float_t>();
            const auto& y_ = spc.get().get("y")->elements<float_t>();
            const auto& z_ = spc.get().get("z")->elements<float_t>();
            const std::size_t n = x_.size();
            x.insert(x.end(), x_.begin(), x_.end());
            y.insert(y.end(), y_.begin(), y_.end());
            z.insert(z.end(), z_.begin(), z_.end());
            q.insert(q.end(), n, 1.0);
            cluster_id.insert(cluster_id.end(), n, cid);
        }
        ++cid;
    }

    Json::Value json_x(Json::arrayValue), json_y(Json::arrayValue),
                json_z(Json::arrayValue), json_q(Json::arrayValue),
                json_cluster_id(Json::arrayValue);
    for (auto v : x) json_x.append(v / units::cm);
    for (auto v : y) json_y.append(v / units::cm);
    for (auto v : z) json_z.append(v / units::cm);
    for (auto v : q) json_q.append(v);
    for (auto v : cluster_id) json_cluster_id.append(v);
    bee["x"] = json_x;
    bee["y"] = json_y;
    bee["z"] = json_z;
    bee["q"] = json_q;
    bee["cluster_id"] = json_cluster_id;

    std::ofstream file(fn);
    if (!file.is_open()) raise<ValueError>("Failed to open file: " + fn);
    Json::StreamWriterBuilder writer;
    writer["indentation"] = "    ";
    writer["precision"] = 6;
    std::unique_ptr<Json::StreamWriter> jw(writer.newStreamWriter());
    jw->write(bee, &file);
}

void WireCell::Match::dump_bee_flash(const ITensorSet::pointer& ts, const std::string& fn)
{
    using spdlog::debug;
    const auto& tens = ts->tensors();
    if (tens->size() != 1) raise<ValueError>("Expected 1 tensor, got %d", tens->size());
    const auto& ten = tens->at(0);
    using MultiArray = boost::multi_array<double, 2>;
    boost::array<MultiArray::index, 2> shape = {(int)ten->shape()[0], (int)ten->shape()[1]};
    boost::multi_array_ref<double, 2> mar((double*)ten->data(), shape);

    Json::Value data;
    data["runNo"] = 0;
    data["subRunNo"] = 0;
    data["eventNo"] = 0;
    data["geom"] = "sbnd";
    data["op_t"] = Json::Value(Json::arrayValue);
    data["op_pes"] = Json::Value(Json::arrayValue);
    data["op_pes_pred"] = Json::Value(Json::arrayValue);
    data["op_peTotal"] = Json::Value(Json::arrayValue);
    data["cluster_id"] = Json::Value(Json::arrayValue);
    data["op_nomatching_cluster_ids"] = Json::Value(Json::arrayValue);

    debug("shape: {} {}", shape[0], shape[1]);
    for (std::size_t i = 0; i < (std::size_t)shape[0]; ++i) {
        data["op_t"].append(mar[i][0]);
        data["cluster_id"].append((int)i);
        double op_peTotal = 0;
        auto op_pes = Json::Value(Json::arrayValue);
        for (std::size_t j = 1; j < (std::size_t)shape[1]; ++j) {
            op_peTotal += mar[i][j];
            op_pes.append(mar[i][j]);
        }
        data["op_peTotal"].append(op_peTotal);
        data["op_pes"].append(op_pes);
        data["op_pes_pred"].append(op_pes);
    }

    std::ofstream file(fn);
    if (!file.is_open()) raise<ValueError>("Failed to open file: " + fn);
    Json::StreamWriterBuilder writer;
    writer["indentation"] = "    ";
    writer["precision"] = 6;
    std::unique_ptr<Json::StreamWriter> jw(writer.newStreamWriter());
    jw->write(data, &file);
}

void WireCell::Match::dump_bee_bundle(const FlashBundlesMap& f2bundle,
                                      const std::map<Cluster*, int>& cluster_idx_map,
                                      const std::string& fn)
{
    Json::Value data;
    data["runNo"] = 0;
    data["subRunNo"] = 0;
    data["eventNo"] = 0;
    data["geom"] = "sbnd";
    data["op_t"] = Json::Value(Json::arrayValue);
    data["op_pes"] = Json::Value(Json::arrayValue);
    data["op_pes_pred"] = Json::Value(Json::arrayValue);
    data["op_peTotal"] = Json::Value(Json::arrayValue);
    data["cluster_id"] = Json::Value(Json::arrayValue);
    data["op_nomatching_cluster_ids"] = Json::Value(Json::arrayValue);

    for (auto it = f2bundle.begin(); it != f2bundle.end(); ++it) {
        auto flash = it->first;
        auto bundles = it->second;
        auto op_pes = Json::Value(Json::arrayValue);
        double op_peTotal = 0;
        for (auto pe : flash->get_PEs()) {
            op_pes.append(pe);
            op_peTotal += pe;
        }
        for (std::size_t i = 0; i < bundles.size(); ++i) {
            auto bundle = bundles.at(i);
            const auto pred_pe_tot = bundle->get_total_pred_light();
            if (pred_pe_tot < 100) continue;
            auto cluster = bundle->get_main_cluster();
            auto cluster_id = cluster_idx_map.at(cluster);
            auto op_cluster_id = Json::Value(Json::arrayValue);
            op_cluster_id.append(cluster_id);
            auto op_pes_pred = Json::Value(Json::arrayValue);
            for (auto pe : bundle->get_pred_flash()) op_pes_pred.append(pe);
            data["op_t"].append(flash->get_time() * 1e-3); // ns -> us
            data["op_pes"].append(op_pes);
            data["op_peTotal"].append(op_peTotal);
            data["cluster_id"].append(op_cluster_id);
            data["op_pes_pred"].append(op_pes_pred);
        }
    }

    std::ofstream file(fn);
    if (!file.is_open()) raise<ValueError>("Failed to open file: " + fn);
    Json::StreamWriterBuilder writer;
    writer["indentation"] = "    ";
    writer["precision"] = 6;
    std::unique_ptr<Json::StreamWriter> jw(writer.newStreamWriter());
    jw->write(data, &file);
}
