#include "WireCellTruth/TruthPointLabeler.h"
#include "WireCellTruth/BlobDepoLabel.h"
#include "WireCellTruth/BlobDepoLabelTensor.h"
#include "WireCellTruth/Schema.h"
#include "WireCellTruth/TruthTable.h"

#include "WireCellAux/SamplingHelpers.h"
#include "WireCellAux/SimpleTensorSet.h"
#include "WireCellAux/TensorDMcommon.h"
#include "WireCellAux/TensorDMpointtree.h"
#include "WireCellIface/WirePlaneId.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/String.h"

#include <boost/filesystem.hpp>

WIRECELL_FACTORY(TruthPointLabeler, WireCell::Truth::TruthPointLabeler,
                 WireCell::INamed,
                 WireCell::ITruthPointLabeler, WireCell::IConfigurable)

using namespace WireCell;

Truth::TruthPointLabeler::TruthPointLabeler()
    : Aux::Logger("TruthPointLabeler", "truth")
{
}
Truth::TruthPointLabeler::~TruthPointLabeler()
{
}

Configuration Truth::TruthPointLabeler::default_configuration() const
{
    auto cfg = assoc_default_configuration(m_assoc);
    cfg["inpath"] = m_inpath;
    cfg["anode"] = "";
    cfg["face"] = m_face;
    cfg["detector_volumes"] = "";
    cfg["reco_time_offset"] = m_reco_time_offset;
    cfg["drift_speed"] = m_drift_speed;
    cfg["tick"] = m_label.tick;
    cfg["truth_table"] = m_truth_table;
    cfg["truth_prefix"] = m_truth_prefix;
    cfg["blobtruth_path"] = m_blobtruth_path;
    cfg["id_mode"] = m_id_mode;
    cfg["pfmerger"] = Json::objectValue;
    cfg["pfmerger"]["pf_ke_min"] = m_pfparams.pf_ke_min;
    cfg["pfmerger"]["q_min"] = m_pfparams.q_min;
    cfg["pfmerger"]["delta_len"] = m_pfparams.delta_len;
    cfg["pfmerger"]["compton_q"] = m_pfparams.compton_q;
    cfg["pfmerger"]["track_len_min"] = m_pfparams.track_len_min;
    cfg["pfmerger"]["michel_own"] = m_pfparams.michel_own;
    cfg["pfmerger"]["photon_prop"] = m_pfparams.photon_prop;
    cfg["pfmerger"]["void_radius"] = m_pfparams.void_radius;
    cfg["refine_by_nearest_depo"] = m_label.refine;
    cfg["d_cap"] = m_label.d_cap;
    cfg["end_window"] = m_label.end_window;
    cfg["kink_window"] = m_label.kink_window;
    cfg["kink_angle"] = m_label.kink_angle;
    cfg["heat_sigma"] = m_label.heat_sigma;
    cfg["shared_frac"] = m_label.shared_frac;
    cfg["labels_meta"] = m_labels_meta;
    cfg["label_dead"] = m_label_dead;
    return cfg;
}

void Truth::TruthPointLabeler::configure(const Configuration& cfg)
{
    assoc_configure(m_assoc, cfg);
    m_inpath = get(cfg, "inpath", m_inpath);
    m_face = get(cfg, "face", m_face);
    const std::string anode_tn = get<std::string>(cfg, "anode", "");
    if (anode_tn.empty()) {
        THROW(ValueError() << errmsg{"TruthPointLabeler: anode is required"});
    }
    m_anode = Factory::find_tn<IAnodePlane>(anode_tn);
    const std::string dv_tn = get<std::string>(cfg, "detector_volumes", "");
    m_reco_time_offset = get(cfg, "reco_time_offset", m_reco_time_offset);
    m_drift_speed = get(cfg, "drift_speed", m_drift_speed);
    m_label.tick = get(cfg, "tick", m_label.tick);
    if (!dv_tn.empty()) {
        m_dv = Factory::find_tn<IDetectorVolumes>(dv_tn);
        WirePlaneId wpid(kAllLayers, m_face, m_anode->ident());
        const auto md = m_dv->metadata(wpid);
        m_reco_time_offset = md["time_offset"].asDouble();
        m_drift_speed = md["drift_speed"].asDouble();
        m_label.tick = md["tick"].asDouble();
    }
    if (m_drift_speed <= 0.0) {
        THROW(ValueError() << errmsg{"TruthPointLabeler: drift_speed (or detector_volumes) is required"});
    }
    m_truth_table = get(cfg, "truth_table", m_truth_table);
    m_truth_prefix = get(cfg, "truth_prefix", m_truth_prefix);
    m_blobtruth_path = get(cfg, "blobtruth_path", m_blobtruth_path);
    m_id_mode = get(cfg, "id_mode", m_id_mode);
    if (cfg.isMember("pfmerger")) m_pfparams.configure(cfg["pfmerger"]);
    m_label.refine = get(cfg, "refine_by_nearest_depo", m_label.refine);
    m_label.d_cap = get(cfg, "d_cap", m_label.d_cap);
    m_label.end_window = get(cfg, "end_window", m_label.end_window);
    m_label.kink_window = get(cfg, "kink_window", m_label.kink_window);
    m_label.kink_angle = get(cfg, "kink_angle", m_label.kink_angle);
    m_label.heat_sigma = get(cfg, "heat_sigma", m_label.heat_sigma);
    m_label.shared_frac = get(cfg, "shared_frac", m_label.shared_frac);
    m_labels_meta = get(cfg, "labels_meta", m_labels_meta);
    m_label_dead = get(cfg, "label_dead", m_label_dead);
    if (m_label_dead) {
        log->warn("label_dead is reserved: dead nodes are not labelled in this version");
    }
    m_assoc.keep_depos = m_label.refine;
}

Configuration Truth::TruthPointLabeler::labels_meta(const TruthTable& tt, const PFResult& pf, const PointLabelSummary& sum) const
{
    Configuration m;
    m["schema_version"] = schema_version;
    m["schema"] = schema_version_string();
    m["pf_source"] = tt.source;
    m["pfmerger"] = m_pfparams.meta();
    m["point_label"] = m_label.meta();
    m["assoc"] = assoc_default_configuration(m_assoc);
    m["id_mode"] = m_id_mode;
    m["reco_time_offset"] = m_reco_time_offset;
    m["drift_speed"] = m_drift_speed;
    m["anode"] = m_anode->ident();
    m["face"] = m_face;
    if (!tt.meta.isNull()) m["truth_file_meta"] = tt.meta;
    m["process_names"] = Json::arrayValue;
    for (const auto& p : process_names()) m["process_names"].append(p);
    m["cls_names"] = Json::arrayValue;
    for (int c = 0; c <= (int)Cls::le; ++c) m["cls_names"].append(cls_name((Cls)c));
    m["truth_pc_keys"] = Json::arrayValue;
    for (const auto& k : truth_pc_keys()) m["truth_pc_keys"].append(k);
    m["summary"]["nclusters"] = (Json::UInt64)sum.nclusters;
    m["summary"]["nnodes"] = (Json::UInt64)sum.nnodes;
    m["summary"]["npoints"] = (Json::UInt64)sum.npoints;
    m["summary"]["nblobs_skipped"] = (Json::UInt64)sum.nblobs_skipped;
    m["summary"]["nghost"] = (Json::UInt64)sum.nghost;
    m["summary"]["nshared"] = (Json::UInt64)sum.nshared;
    m["summary"]["nrefined_differ"] = (Json::UInt64)sum.nrefined_differ;
    m["summary"]["npolylines"] = (Json::UInt64)sum.npolylines;
    m["summary"]["pf_nodes"] = (Json::UInt64)pf.nodes.size();
    m["summary"]["pf_merged_delta"] = (Json::UInt64)pf.nmerged_delta;
    m["summary"]["pf_merged_shower"] = (Json::UInt64)pf.nmerged_shower;
    m["summary"]["pf_merged_weak"] = (Json::UInt64)pf.nmerged_weak;
    m["summary"]["pf_merged_invisible"] = (Json::UInt64)pf.nmerged_invisible;
    m["summary"]["void_roots"] = Json::arrayValue;
    for (int r : pf.void_roots) m["summary"]["void_roots"].append(r);
    return m;
}

bool Truth::TruthPointLabeler::operator()(const input_tuple_type& intup, output_pointer& out)
{
    out = nullptr;
    auto ints = std::get<0>(intup);
    auto icluster = std::get<1>(intup);
    auto ideposet = std::get<2>(intup);
    if (!ints || !icluster || !ideposet) {
        if (ints) log->warn("join node got non-EOS ITensorSet {} on EOS", ints->ident());
        if (icluster) log->warn("join node got non-EOS ICluster {} on EOS", icluster->ident());
        if (ideposet) log->warn("join node got non-EOS IDepoSet {} on EOS", ideposet->ident());
        log->debug("EOS at call={}", m_count);
        ++m_count;
        return true;
    }
    const int ident = ints->ident();
    const std::string inpath = m_inpath.find('%') == std::string::npos ? m_inpath : String::format(m_inpath, ident);

    // 1. the truth table
    TruthTable tt;
    const auto& depos = *ideposet->depos();
    if (!m_truth_table.empty()) {
        const std::string fname = m_truth_table.find('%') == std::string::npos ? m_truth_table
            : String::format(m_truth_table, icluster->ident());
        if (boost::filesystem::exists(fname)) {
            tt = TruthTable::load(fname, log);
        }
        else {
            log->warn("no truth file {} for cluster {}: particle table from the depos", fname, icluster->ident());
        }
    }
    if (tt.empty()) {
        tt = TruthTable::from_depos(depos, m_assoc.id_from_prior);
    }
    // id map for Stage A index mode
    DepoBlobAssocParams assoc = m_assoc;
    bool index_mode = (m_id_mode == "index") || (m_id_mode == "auto" && !tt.sed_ids.empty());
    if (index_mode && tt.sed_ids.empty()) {
        log->warn("id_mode index but the truth table has no sed_ids: ids taken as track ids");
        index_mode = false;
    }
    if (index_mode) {
        auto map = std::make_shared<std::vector<int>>();
        map->reserve(tt.sed_ids.size());
        for (const auto& [tid, orig] : tt.sed_ids) map->push_back(tid);
        assoc.id_map = map;
        if (!tt.sed_energy.empty()) {
            assoc.energy_map = std::make_shared<std::vector<float>>(tt.sed_energy);
        }
    }

    // 2. the association (blob labels identical to BlobDepoLabel)
    auto res = depo_blob_assoc(icluster, depos, assoc, log);

    // 3. the particle flow with the measured visible charge per track
    std::map<int, double> q_vis;
    for (const auto& [vdesc, e] : res.blobs) {
        for (const auto& [tid, qe] : e.by_track) q_vis[tid] += qe.first;
    }
    PFMerger merger(m_pfparams);
    auto pf = merger.run(tt, q_vis);
    {
        // record the PF results in the table rows
        const auto& pfc = truth_pf_columns();
        const int c_tid = column_index(pfc, "trackid"), c_node = column_index(pfc, "pf_node"),
            c_cls = column_index(pfc, "cls"), c_merged = column_index(pfc, "merged_into"), c_nvis = column_index(pfc, "nvis_q");
        for (auto& row : tt.pf) {
            const int tid = (int)row[c_tid];
            const auto& gi = pf.info(tid);
            row[c_node] = gi.pf_node;
            row[c_cls] = (double)(int)gi.cls;
            row[c_merged] = gi.merged_into;
            auto qit = q_vis.find(tid);
            row[c_nvis] = qit == q_vis.end() ? 0.0 : qit->second;
        }
    }

    // 4. the tree
    Aux::TensorDM::TensorIndex ti(*ints->tensors());
    std::unique_ptr<PointCloud::Tree::Points::node_t> root;
    if (ti.get(inpath)) {
        try {
            root = Aux::TensorDM::as_pctree(ti, inpath);
        }
        catch (const KeyError& err) {
            // a tree tensor whose named point clouds are not in the set
            root = nullptr;
        }
    }
    if (!root) {
        // A face without a point tree: PointTreeBuilding emits no tree for
        // a face it never sampled (PDHD wall faces, pr_ml doc 20 R2), so
        // there is nothing to label -- the set passes through unchanged.
        log->debug("call={} ident {} cluster {}: no point tree at {}, passing the set through",
                   m_count, ident, icluster->ident(), inpath);
        out = ints;
        ++m_count;
        return true;
    }
    IAnodeFace::pointer face = m_anode->face(m_face);
    if (!face) {
        // SBND's per-APA clustering runs every APA as "face 0" while APA 1
        // only has face 1 (doc 20 sec 4.5): take the anode's first face.
        for (const auto& f : m_anode->faces()) {
            if (f) { face = f; break; }
        }
        if (face) {
            log->debug("no face {} on anode {}, using face {}", m_face, m_anode->ident(), face->which());
        }
    }
    if (!face) {
        THROW(ValueError() << errmsg{"TruthPointLabeler: no face " + std::to_string(m_face)});
    }
    const double toff = m_reco_time_offset, speed = m_drift_speed, depo_toff = m_assoc.time_offset;
    depo_position_f position = [face, toff, speed, depo_toff](const IDepo& d) {
        const double x = Aux::time2drift(face, toff, speed, d.time() + depo_toff);
        return Point(x, d.pos().y(), d.pos().z());
    };
    auto polylines = build_polylines(depos, assoc, position);
    auto sum = label_live_tree(*root, icluster, depos, res, tt, pf, polylines, m_label, position, log);

    // 5. the output set: pass-through tensors, the relabelled tree, the truth tensors
    auto tv = std::make_shared<ITensor::vector>();
    for (const auto& ten : *ints->tensors()) {
        const std::string dp = ten->metadata()["datapath"].asString();
        if (dp == inpath || dp.rfind(inpath + "/", 0) == 0) continue;
        tv->push_back(ten);
    }
    auto tree_tens = Aux::TensorDM::as_tensors(*root, inpath);
    tv->insert(tv->end(), tree_tens.begin(), tree_tens.end());
    const std::string bpath = m_blobtruth_path.find('%') == std::string::npos ? m_blobtruth_path
        : String::format(m_blobtruth_path, icluster->ident());
    auto btens = blobtruth_tensors(icluster, res, bpath);
    tv->insert(tv->end(), btens.begin(), btens.end());
    const std::string tprefix = m_truth_prefix.find('%') == std::string::npos ? m_truth_prefix
        : String::format(m_truth_prefix, ident);
    auto ttens = tt.as_tensors(ident, tprefix);
    tv->insert(tv->end(), ttens.begin(), ttens.end());

    Configuration setmd = ints->metadata();
    const auto meta = labels_meta(tt, pf, sum);
    setmd["labels_meta"] = meta;
    if (!m_labels_meta.empty()) {
        const std::string fname = m_labels_meta.find('%') == std::string::npos ? m_labels_meta : String::format(m_labels_meta, ident);
        Persist::dump(fname, meta, true);
    }
    out = std::make_shared<Aux::SimpleTensorSet>(ident, setmd, tv);
    log->debug("call={} ident {} cluster {} depos {}: {} nodes {} points labelled, {} pf nodes ({} source), {} ghost {} shared blobs, {} tensors out",
               m_count, ident, icluster->ident(), depos.size(), sum.nnodes, sum.npoints, pf.nodes.size(), tt.source,
               sum.nghost, sum.nshared, tv->size());
    ++m_count;
    return true;
}
