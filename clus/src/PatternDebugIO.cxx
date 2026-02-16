#include "WireCellClus/PatternDebugIO.h"
#include "WireCellClus/ClusteringFuncs.h"  // Flags::main_cluster
#include "WireCellClus/Graphs.h"
#include "WireCellIface/WirePlaneId.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/PointCloudDataset.h"
#include "WireCellUtil/Logging.h"

#include <boost/graph/adjacency_list.hpp>
#include <json/json.h>
#include <fstream>

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::PointCloud;

static spdlog::logger& logger() {
    static auto log = Log::logger("debugio");
    return *log;
}

// Helper: serialize a coordinate array from a Dataset to JSON
static Json::Value pc_coords_to_json(const Dataset& ds, const std::vector<std::string>& coord_names)
{
    Json::Value jpc(Json::objectValue);
    for (const auto& name : coord_names) {
        auto arr = ds.get(name);
        if (!arr) continue;
        auto elems = arr->elements<double>();
        Json::Value jarr(Json::arrayValue);
        for (double v : elems) {
            jarr.append(v);
        }
        jpc[name] = jarr;
    }
    return jpc;
}

// Helper: serialize a steiner graph to JSON
static Json::Value graph_to_json(const Graphs::Weighted::Graph& graph)
{
    Json::Value jg(Json::objectValue);
    jg["num_vertices"] = static_cast<Json::UInt>(boost::num_vertices(graph));

    auto weight_map = boost::get(boost::edge_weight, graph);
    Json::Value jedges(Json::arrayValue);
    auto [ebegin, eend] = boost::edges(graph);
    for (auto eit = ebegin; eit != eend; ++eit) {
        Json::Value edge(Json::arrayValue);
        edge.append(static_cast<Json::UInt>(boost::source(*eit, graph)));
        edge.append(static_cast<Json::UInt>(boost::target(*eit, graph)));
        edge.append(boost::get(weight_map, *eit));
        jedges.append(edge);
    }
    jg["edges"] = jedges;
    return jg;
}

// Helper: serialize TrackFitting::Parameters to JSON
static Json::Value params_to_json(const TrackFitting::Parameters& p)
{
    Json::Value j(Json::objectValue);
    j["DL"] = p.DL;
    j["DT"] = p.DT;
    j["col_sigma_w_T"] = p.col_sigma_w_T;
    j["ind_sigma_u_T"] = p.ind_sigma_u_T;
    j["ind_sigma_v_T"] = p.ind_sigma_v_T;
    j["rel_uncer_ind"] = p.rel_uncer_ind;
    j["rel_uncer_col"] = p.rel_uncer_col;
    j["add_uncer_ind"] = p.add_uncer_ind;
    j["add_uncer_col"] = p.add_uncer_col;
    j["add_sigma_L"] = p.add_sigma_L;
    j["rel_charge_uncer"] = p.rel_charge_uncer;
    j["add_charge_uncer"] = p.add_charge_uncer;
    j["default_charge_th"] = p.default_charge_th;
    j["default_charge_err"] = p.default_charge_err;
    j["scaling_quality_th"] = p.scaling_quality_th;
    j["scaling_ratio"] = p.scaling_ratio;
    j["area_ratio1"] = p.area_ratio1;
    j["area_ratio2"] = p.area_ratio2;
    j["skip_default_ratio_1"] = p.skip_default_ratio_1;
    j["skip_ratio_cut"] = p.skip_ratio_cut;
    j["skip_ratio_1_cut"] = p.skip_ratio_1_cut;
    j["skip_angle_cut_1"] = p.skip_angle_cut_1;
    j["skip_angle_cut_2"] = p.skip_angle_cut_2;
    j["skip_angle_cut_3"] = p.skip_angle_cut_3;
    j["skip_dis_cut"] = p.skip_dis_cut;
    j["default_dQ_dx"] = p.default_dQ_dx;
    j["end_point_factor"] = p.end_point_factor;
    j["mid_point_factor"] = p.mid_point_factor;
    j["nlevel"] = p.nlevel;
    j["charge_cut"] = p.charge_cut;
    j["low_dis_limit"] = p.low_dis_limit;
    j["end_point_limit"] = p.end_point_limit;
    j["time_tick_cut"] = p.time_tick_cut;
    j["share_charge_err"] = p.share_charge_err;
    j["min_drift_time"] = p.min_drift_time;
    j["search_range"] = p.search_range;
    j["dead_ind_weight"] = p.dead_ind_weight;
    j["dead_col_weight"] = p.dead_col_weight;
    j["close_ind_weight"] = p.close_ind_weight;
    j["close_col_weight"] = p.close_col_weight;
    j["overlap_th"] = p.overlap_th;
    j["dx_norm_length"] = p.dx_norm_length;
    j["lambda"] = p.lambda;
    j["div_sigma"] = p.div_sigma;
    return j;
}

// Helper: deserialize TrackFitting::Parameters from JSON
static TrackFitting::Parameters params_from_json(const Json::Value& j)
{
    TrackFitting::Parameters p;
    if (j.isMember("DL")) p.DL = j["DL"].asDouble();
    if (j.isMember("DT")) p.DT = j["DT"].asDouble();
    if (j.isMember("col_sigma_w_T")) p.col_sigma_w_T = j["col_sigma_w_T"].asDouble();
    if (j.isMember("ind_sigma_u_T")) p.ind_sigma_u_T = j["ind_sigma_u_T"].asDouble();
    if (j.isMember("ind_sigma_v_T")) p.ind_sigma_v_T = j["ind_sigma_v_T"].asDouble();
    if (j.isMember("rel_uncer_ind")) p.rel_uncer_ind = j["rel_uncer_ind"].asDouble();
    if (j.isMember("rel_uncer_col")) p.rel_uncer_col = j["rel_uncer_col"].asDouble();
    if (j.isMember("add_uncer_ind")) p.add_uncer_ind = j["add_uncer_ind"].asDouble();
    if (j.isMember("add_uncer_col")) p.add_uncer_col = j["add_uncer_col"].asDouble();
    if (j.isMember("add_sigma_L")) p.add_sigma_L = j["add_sigma_L"].asDouble();
    if (j.isMember("rel_charge_uncer")) p.rel_charge_uncer = j["rel_charge_uncer"].asDouble();
    if (j.isMember("add_charge_uncer")) p.add_charge_uncer = j["add_charge_uncer"].asDouble();
    if (j.isMember("default_charge_th")) p.default_charge_th = j["default_charge_th"].asDouble();
    if (j.isMember("default_charge_err")) p.default_charge_err = j["default_charge_err"].asDouble();
    if (j.isMember("scaling_quality_th")) p.scaling_quality_th = j["scaling_quality_th"].asDouble();
    if (j.isMember("scaling_ratio")) p.scaling_ratio = j["scaling_ratio"].asDouble();
    if (j.isMember("area_ratio1")) p.area_ratio1 = j["area_ratio1"].asDouble();
    if (j.isMember("area_ratio2")) p.area_ratio2 = j["area_ratio2"].asDouble();
    if (j.isMember("skip_default_ratio_1")) p.skip_default_ratio_1 = j["skip_default_ratio_1"].asDouble();
    if (j.isMember("skip_ratio_cut")) p.skip_ratio_cut = j["skip_ratio_cut"].asDouble();
    if (j.isMember("skip_ratio_1_cut")) p.skip_ratio_1_cut = j["skip_ratio_1_cut"].asDouble();
    if (j.isMember("skip_angle_cut_1")) p.skip_angle_cut_1 = j["skip_angle_cut_1"].asDouble();
    if (j.isMember("skip_angle_cut_2")) p.skip_angle_cut_2 = j["skip_angle_cut_2"].asDouble();
    if (j.isMember("skip_angle_cut_3")) p.skip_angle_cut_3 = j["skip_angle_cut_3"].asDouble();
    if (j.isMember("skip_dis_cut")) p.skip_dis_cut = j["skip_dis_cut"].asDouble();
    if (j.isMember("default_dQ_dx")) p.default_dQ_dx = j["default_dQ_dx"].asDouble();
    if (j.isMember("end_point_factor")) p.end_point_factor = j["end_point_factor"].asDouble();
    if (j.isMember("mid_point_factor")) p.mid_point_factor = j["mid_point_factor"].asDouble();
    if (j.isMember("nlevel")) p.nlevel = j["nlevel"].asInt();
    if (j.isMember("charge_cut")) p.charge_cut = j["charge_cut"].asDouble();
    if (j.isMember("low_dis_limit")) p.low_dis_limit = j["low_dis_limit"].asDouble();
    if (j.isMember("end_point_limit")) p.end_point_limit = j["end_point_limit"].asDouble();
    if (j.isMember("time_tick_cut")) p.time_tick_cut = j["time_tick_cut"].asDouble();
    if (j.isMember("share_charge_err")) p.share_charge_err = j["share_charge_err"].asDouble();
    if (j.isMember("min_drift_time")) p.min_drift_time = j["min_drift_time"].asDouble();
    if (j.isMember("search_range")) p.search_range = j["search_range"].asDouble();
    if (j.isMember("dead_ind_weight")) p.dead_ind_weight = j["dead_ind_weight"].asDouble();
    if (j.isMember("dead_col_weight")) p.dead_col_weight = j["dead_col_weight"].asDouble();
    if (j.isMember("close_ind_weight")) p.close_ind_weight = j["close_ind_weight"].asDouble();
    if (j.isMember("close_col_weight")) p.close_col_weight = j["close_col_weight"].asDouble();
    if (j.isMember("overlap_th")) p.overlap_th = j["overlap_th"].asDouble();
    if (j.isMember("dx_norm_length")) p.dx_norm_length = j["dx_norm_length"].asDouble();
    if (j.isMember("lambda")) p.lambda = j["lambda"].asDouble();
    if (j.isMember("div_sigma")) p.div_sigma = j["div_sigma"].asDouble();
    return p;
}

// Helper: dump one cluster's steiner data to JSON
static Json::Value dump_cluster_data(const Facade::Cluster& cluster)
{
    Json::Value jcluster(Json::objectValue);

    // Scope coordinate names
    const auto& scope = cluster.get_default_scope();
    Json::Value jcoords(Json::arrayValue);
    for (const auto& c : scope.coords) {
        jcoords.append(c);
    }
    jcluster["scope_coords"] = jcoords;

    // Steiner point cloud
    if (cluster.has_pc("steiner_pc")) {
        const auto& spc = cluster.get_pc("steiner_pc");
        jcluster["steiner_pc"] = pc_coords_to_json(spc, scope.coords);
        jcluster["steiner_pc"]["npoints"] = static_cast<Json::UInt>(spc.size_major());

        // Also dump flag_steiner_terminal if present
        auto terminal_arr = spc.get("flag_steiner_terminal");
        if (terminal_arr) {
            auto elems = terminal_arr->elements<int>();
            Json::Value jterm(Json::arrayValue);
            for (int v : elems) {
                jterm.append(v);
            }
            jcluster["steiner_pc"]["flag_steiner_terminal"] = jterm;
        }

        // Dump wpid array (WirePlaneId stored as int)
        auto wpid_arr = spc.get("wpid");
        if (wpid_arr) {
            auto elems = wpid_arr->elements<WirePlaneId>();
            Json::Value jwpid(Json::arrayValue);
            for (const auto& wpid : elems) {
                jwpid.append(wpid.ident());
            }
            jcluster["steiner_pc"]["wpid"] = jwpid;
        }
    }

    // Steiner graph
    if (cluster.has_graph("steiner_graph")) {
        const auto& graph = cluster.get_graph("steiner_graph");
        jcluster["steiner_graph"] = graph_to_json(graph);
    }

    // Flags
    Json::Value jflags(Json::objectValue);
    jflags["main_cluster"] = cluster.get_flag(Facade::Flags::main_cluster);
    jcluster["flags"] = jflags;

    // Scope transform name (e.g. "Unity", "T0Correction")
    jcluster["scope_transform"] = cluster.get_scope_transform();

    // Blob 3d point clouds: dump ALL arrays from each blob's "3d" PC
    // (not just coordinate arrays - we need corrected arrays like x_t0cor too)
    Json::Value jblobs(Json::arrayValue);
    try {
        for (const auto* blob : cluster.children()) {
            const auto& lpcs = blob->node()->value.local_pcs();
            auto it3d = lpcs.find("3d");
            if (it3d == lpcs.end()) continue;

            const auto& ds = it3d->second;
            Json::Value jblob_3d(Json::objectValue);
            Json::Value jarray_types(Json::objectValue);
            for (const auto& key : ds.keys()) {
                auto arr = ds.get(key);
                if (!arr) continue;
                jarray_types[key] = arr->dtype();
                Json::Value jarr(Json::arrayValue);
                if (arr->is_type<double>()) {
                    for (double v : arr->elements<double>()) jarr.append(v);
                } else if (arr->is_type<float>()) {
                    for (float v : arr->elements<float>()) jarr.append(v);
                } else if (arr->is_type<int>()) {
                    for (int v : arr->elements<int>()) jarr.append(v);
                } else {
                    // Fallback: try double
                    for (double v : arr->elements<double>()) jarr.append(v);
                }
                jblob_3d[key] = jarr;
            }
            jblob_3d["_array_types"] = jarray_types;

            // Also dump scalar PC for this blob
            auto it_sc = lpcs.find("scalar");
            if (it_sc != lpcs.end()) {
                const auto& sc_ds = it_sc->second;
                Json::Value jscalar(Json::objectValue);
                Json::Value jsc_types(Json::objectValue);
                for (const auto& key : sc_ds.keys()) {
                    auto arr = sc_ds.get(key);
                    if (!arr) continue;
                    jsc_types[key] = arr->dtype();
                    Json::Value jarr(Json::arrayValue);
                    if (arr->is_type<double>()) {
                        for (double v : arr->elements<double>()) jarr.append(v);
                    } else if (arr->is_type<float>()) {
                        for (float v : arr->elements<float>()) jarr.append(v);
                    } else if (arr->is_type<int>()) {
                        for (int v : arr->elements<int>()) jarr.append(v);
                    } else {
                        for (double v : arr->elements<double>()) jarr.append(v);
                    }
                    jscalar[key] = jarr;
                }
                jscalar["_array_types"] = jsc_types;
                jblob_3d["_scalar"] = jscalar;
            }

            jblobs.append(jblob_3d);
        }
    } catch (const std::exception& e) {
        logger().warn("Could not dump blob 3d PCs: {}", e.what());
    }
    jcluster["blobs"] = jblobs;

    return jcluster;
}


void PR::DebugIO::dump_init_first_segment_inputs(
    const std::string& output_path,
    const Facade::Cluster& cluster,
    const Facade::Cluster* main_cluster,
    bool flag_back_search,
    const TrackFitting& track_fitter)
{
    Json::Value root(Json::objectValue);

    root["cluster"] = dump_cluster_data(cluster);
    root["flag_back_search"] = flag_back_search;
    root["trackfitting_params"] = params_to_json(track_fitter.get_parameters());

    // Precompute boundary steiner indices (needs anode, which is available now)
    try {
        auto boundary = cluster.get_two_boundary_steiner_graph_idx(
            "steiner_graph", "steiner_pc");
        root["boundary_steiner_indices"] = Json::Value(Json::arrayValue);
        root["boundary_steiner_indices"].append(static_cast<Json::UInt>(boundary.first));
        root["boundary_steiner_indices"].append(static_cast<Json::UInt>(boundary.second));
    } catch (const std::exception& e) {
        logger().warn("Could not dump boundary steiner indices: {}", e.what());
    }

    // Main cluster (only if it differs from cluster)
    if (main_cluster && main_cluster != &cluster) {
        root["main_cluster"] = dump_cluster_data(*main_cluster);
    }

    Persist::dump(output_path, root, true);
    logger().info("Dumped init_first_segment inputs to {}", output_path);
}


// Helper: reconstruct a Dataset from JSON coordinate arrays
static Dataset json_to_dataset(const Json::Value& jpc, const std::vector<std::string>& coord_names)
{
    std::map<std::string, Array> arrays;
    for (const auto& name : coord_names) {
        if (!jpc.isMember(name)) continue;
        const auto& jarr = jpc[name];
        std::vector<double> vals;
        vals.reserve(jarr.size());
        for (const auto& v : jarr) {
            vals.push_back(v.asDouble());
        }
        arrays.emplace(name, Array(vals));
    }

    // Also restore flag_steiner_terminal if present
    if (jpc.isMember("flag_steiner_terminal")) {
        const auto& jterm = jpc["flag_steiner_terminal"];
        std::vector<int> vals;
        vals.reserve(jterm.size());
        for (const auto& v : jterm) {
            vals.push_back(v.asInt());
        }
        arrays.emplace("flag_steiner_terminal", Array(vals));
    }

    // Restore wpid array (WirePlaneId is same size as int)
    if (jpc.isMember("wpid")) {
        const auto& jwpid = jpc["wpid"];
        std::vector<int> vals;
        vals.reserve(jwpid.size());
        for (const auto& v : jwpid) {
            vals.push_back(v.asInt());
        }
        arrays.emplace("wpid", Array(vals));
    }

    return Dataset(arrays);
}

// Helper: reconstruct a graph from JSON
static Graphs::Weighted::Graph json_to_graph(const Json::Value& jg)
{
    size_t nv = jg["num_vertices"].asUInt();
    Graphs::Weighted::Graph graph(nv);

    const auto& jedges = jg["edges"];
    for (const auto& edge : jedges) {
        auto src = edge[0].asUInt();
        auto tgt = edge[1].asUInt();
        double weight = edge[2].asDouble();
        boost::add_edge(src, tgt, weight, graph);
    }
    return graph;
}

// Helper: create a minimal cluster node with blob child containing a "3d" PC
// and inject steiner data. Returns the cluster node (raw ptr, owned by parent).
static PointCloud::Tree::Points::node_t* make_cluster_node(
    PointCloud::Tree::Points::node_t& grouping_node,
    const Json::Value& jcluster)
{
    // Read scope coords
    std::vector<std::string> coord_names;
    if (jcluster.isMember("scope_coords")) {
        for (const auto& c : jcluster["scope_coords"]) {
            coord_names.push_back(c.asString());
        }
    }
    if (coord_names.empty()) {
        coord_names = {"x", "y", "z"};
    }

    // Build the cluster subtree: cluster_node -> blob_node(s)
    // Each blob needs a "3d" PC with ALL arrays (including corrected coords like x_t0cor)
    // and a "scalar" PC with blob metadata.
    auto cluster_ptr = std::make_unique<PointCloud::Tree::Points::node_t>();

    if (jcluster.isMember("blobs") && jcluster["blobs"].size() > 0) {
        // Restore real blob data from dump
        for (const auto& jblob_3d : jcluster["blobs"]) {
            // Reconstruct "3d" Dataset with all arrays
            std::map<std::string, Array> blob_arrays;
            Json::Value jtypes = jblob_3d.get("_array_types", Json::objectValue);
            for (const auto& key : jblob_3d.getMemberNames()) {
                if (key == "_array_types" || key == "_scalar") continue;
                const auto& jarr = jblob_3d[key];
                if (!jarr.isArray()) continue;

                std::string dtype = jtypes.get(key, "f8").asString();
                // dtypes are NumPy-style: "f4"=float, "f8"=double, "i4"=int, etc.
                if (dtype == "f4") {
                    std::vector<float> vals;
                    vals.reserve(jarr.size());
                    for (const auto& v : jarr) vals.push_back(v.asFloat());
                    blob_arrays.emplace(key, Array(vals));
                } else if (dtype[0] == 'i') {
                    std::vector<int> vals;
                    vals.reserve(jarr.size());
                    for (const auto& v : jarr) vals.push_back(v.asInt());
                    blob_arrays.emplace(key, Array(vals));
                } else {
                    // "f8", "d", or anything else -> double
                    std::vector<double> vals;
                    vals.reserve(jarr.size());
                    for (const auto& v : jarr) vals.push_back(v.asDouble());
                    blob_arrays.emplace(key, Array(vals));
                }
            }
            Dataset blob_3d_ds(blob_arrays);

            // Reconstruct "scalar" Dataset
            std::map<std::string, Array> scalar_arrays;
            if (jblob_3d.isMember("_scalar")) {
                const auto& jscalar = jblob_3d["_scalar"];
                Json::Value jsc_types = jscalar.get("_array_types", Json::objectValue);
                for (const auto& key : jscalar.getMemberNames()) {
                    if (key == "_array_types") continue;
                    const auto& jarr = jscalar[key];
                    if (!jarr.isArray()) continue;

                    std::string dtype = jsc_types.get(key, "f8").asString();
                    // dtypes are NumPy-style: "f4"=float, "f8"=double, "i4"=int, etc.
                    if (dtype == "f4") {
                        std::vector<float> vals;
                        vals.reserve(jarr.size());
                        for (const auto& v : jarr) vals.push_back(v.asFloat());
                        scalar_arrays.emplace(key, Array(vals));
                    } else if (dtype[0] == 'i') {
                        std::vector<int> vals;
                        vals.reserve(jarr.size());
                        for (const auto& v : jarr) vals.push_back(v.asInt());
                        scalar_arrays.emplace(key, Array(vals));
                    } else {
                        // "f8", "d", or anything else -> double
                        std::vector<double> vals;
                        vals.reserve(jarr.size());
                        for (const auto& v : jarr) vals.push_back(v.asDouble());
                        scalar_arrays.emplace(key, Array(vals));
                    }
                }
            }
            Dataset scalar_ds(scalar_arrays);

            cluster_ptr->insert(
                PointCloud::Tree::Points({
                    {"scalar", std::move(scalar_ds)},
                    {"3d", std::move(blob_3d_ds)},
                }));
        }
    } else {
        // Fallback: create a minimal dummy blob with one point
        using fa_float_t = Facade::float_t;
        using fa_int_t = Facade::int_t;
        std::map<std::string, Array> blob_arrays;
        for (const auto& name : coord_names) {
            blob_arrays.emplace(name, Array(std::vector<double>{0.0}));
        }
        Dataset scalar_ds({
            {"charge", Array(std::vector<fa_float_t>{1.0})},
            {"center_x", Array(std::vector<fa_float_t>{0.0})},
            {"center_y", Array(std::vector<fa_float_t>{0.0})},
            {"center_z", Array(std::vector<fa_float_t>{0.0})},
            {"wpid", Array(std::vector<fa_int_t>{0})},
            {"npoints", Array(std::vector<fa_int_t>{1})},
            {"slice_index_min", Array(std::vector<fa_int_t>{0})},
            {"slice_index_max", Array(std::vector<fa_int_t>{1})},
            {"u_wire_index_min", Array(std::vector<fa_int_t>{0})},
            {"u_wire_index_max", Array(std::vector<fa_int_t>{1})},
            {"v_wire_index_min", Array(std::vector<fa_int_t>{0})},
            {"v_wire_index_max", Array(std::vector<fa_int_t>{1})},
            {"w_wire_index_min", Array(std::vector<fa_int_t>{0})},
            {"w_wire_index_max", Array(std::vector<fa_int_t>{1})},
            {"max_wire_interval", Array(std::vector<fa_int_t>{1})},
            {"min_wire_interval", Array(std::vector<fa_int_t>{1})},
            {"max_wire_type", Array(std::vector<fa_int_t>{0})},
            {"min_wire_type", Array(std::vector<fa_int_t>{0})},
        });
        Dataset blob_3d_ds(blob_arrays);
        cluster_ptr->insert(
            PointCloud::Tree::Points({
                {"scalar", std::move(scalar_ds)},
                {"3d", std::move(blob_3d_ds)},
            }));
    }

    // Insert cluster into grouping
    auto* cluster_node = grouping_node.insert(std::move(cluster_ptr));

    // Get the Cluster facade
    auto* cluster = cluster_node->value.facade<Facade::Cluster>();

    // Inject steiner_pc into cluster's local_pcs
    if (jcluster.isMember("steiner_pc")) {
        auto steiner_ds = json_to_dataset(jcluster["steiner_pc"], coord_names);
        cluster->local_pcs()["steiner_pc"] = std::move(steiner_ds);
    }

    // Inject steiner_graph
    if (jcluster.isMember("steiner_graph")) {
        auto graph = json_to_graph(jcluster["steiner_graph"]);
        cluster->give_graph("steiner_graph", std::move(graph));
    }

    // Set flags
    if (jcluster.isMember("flags")) {
        const auto& jflags = jcluster["flags"];
        if (jflags.isMember("main_cluster") && jflags["main_cluster"].asInt()) {
            cluster->set_flag(Facade::Flags::main_cluster, 1);
        }
    }

    // Directly set the default scope with the dumped coordinate names.
    PointCloud::Tree::Scope loaded_scope{"3d", coord_names};
    cluster->set_default_scope(loaded_scope);

    // Restore the scope transform mapping so get_scope_transform() returns
    // the correct name (e.g. "T0Correction") instead of defaulting to "Unity".
    if (jcluster.isMember("scope_transform")) {
        std::string transform_name = jcluster["scope_transform"].asString();
        cluster->set_scope_transform(loaded_scope, transform_name);
    }

    return cluster_node;
}


PR::DebugIO::LoadedTestData
PR::DebugIO::load_init_first_segment_inputs(const std::string& input_path)
{
    auto root = Persist::load(input_path);

    LoadedTestData data;
    data.flag_back_search = root.get("flag_back_search", true).asBool();

    if (root.isMember("trackfitting_params")) {
        data.trackfitting_params = params_from_json(root["trackfitting_params"]);
    }

    // Precomputed boundary indices
    if (root.isMember("boundary_steiner_indices")) {
        const auto& jbi = root["boundary_steiner_indices"];
        data.boundary_steiner_indices = {jbi[0].asUInt(), jbi[1].asUInt()};
    }

    // Create grouping node (root of the tree)
    data.grouping_node = std::make_unique<PointCloud::Tree::Points::node_t>();
    auto* grouping = data.grouping_node->value.facade<Facade::Grouping>();
    (void)grouping;

    // Create the main cluster
    auto* cluster_node = make_cluster_node(*data.grouping_node, root["cluster"]);
    data.cluster = cluster_node->value.facade<Facade::Cluster>();

    // Create main_cluster (separate or same as cluster)
    if (root.isMember("main_cluster")) {
        auto* main_cluster_node = make_cluster_node(*data.grouping_node, root["main_cluster"]);
        data.main_cluster = main_cluster_node->value.facade<Facade::Cluster>();
    } else {
        data.main_cluster = data.cluster;
    }

    logger().info("Loaded init_first_segment inputs from {}", input_path);
    return data;
}
