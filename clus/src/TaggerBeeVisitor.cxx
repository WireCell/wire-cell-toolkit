/// TaggerBeeVisitor -- the cosmic/containment tagger verdict Bee sets
/// (tagger_stm, tagger_tgm, tagger_fc, tagger_lm) written from INSIDE the PR
/// MultiAlgBlobClustering, so a standalone (no-art) PR job produces them
/// (ai-helper issue 33: the step-2 job wct-pr.jsonnet).
///
/// A port of the tagger branch of larwirecell's wclsTensorSetLabeler (the
/// "labeler_tagger" instance of the LArSoft 1-step chain), which reads the same
/// per-cluster quantities from the serialized post-PR tree.  Same content:
///   * only the "beam-window candidates" are dumped: clusters with
///     flag_main_cluster set and cluster_t0 in [beam_window[0], beam_window[1]);
///   * one point per blob "3d" point, in the "coords" arrays (the clustering
///     Bee coordinates, e.g. x_t0cor/y_cor/z_cor; raw x/y/z when unset or
///     absent), charge = max(blob charge / npoints, 1);
///   * cluster_id = 1 (tagged: flag_STM / flag_TGM / flag_FC set, lm_flag == 2)
///     or 0, real_cluster_id = 0 (Bee colours by the verdict).
/// Written to the shared "bee_sink" (the PR MABC's) under the Bee event index the
/// MABC publishes on the Ensemble (Ensemble::bee_index()) -- the MABC owns the
/// index, so the two writers cannot drift; an Ensemble without one is an error.
/// The MABC holds the sink open for the whole job, so this visitor neither
/// acquires nor releases it.
/// A blob whose "3d" PC lacks x/y/z (or has them of unequal length) or whose
/// "scalar" PC lacks "charge" is skipped, with one warning per event.
/// Must run after the taggers (TaggerCheckSTM/TGM/FC) and QLMatching's lm_flag.

#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/IBeeSink.h"
#include "WireCellClus/Facade_Ensemble.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"
#include "WireCellUtil/Bee.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Point.h"

#include <algorithm>

class TaggerBeeVisitor;
WIRECELL_FACTORY(TaggerBeeVisitor, TaggerBeeVisitor,
                 WireCell::IConfigurable, WireCell::Clus::IEnsembleVisitor)
using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

class TaggerBeeVisitor : public Aux::Logger, public IConfigurable, public Clus::IEnsembleVisitor {
public:
    TaggerBeeVisitor() : Aux::Logger("TaggerBeeVisitor", "clus") {}
    virtual ~TaggerBeeVisitor() {}

    virtual Configuration default_configuration() const {
        Configuration cfg;
        cfg["grouping"] = m_grouping;
        cfg["bee_sink"] = "";
        cfg["detector"] = m_detector;
        cfg["beam_window"][0] = m_bw_low;
        cfg["beam_window"][1] = m_bw_high;
        cfg["coords"] = Json::arrayValue;
        return cfg;
    }

    virtual void configure(const Configuration& cfg) {
        m_grouping = get(cfg, "grouping", m_grouping);
        m_detector = get(cfg, "detector", m_detector);
        const std::string sink = get<std::string>(cfg, "bee_sink", "");
        if (sink.empty()) raise<ValueError>("TaggerBeeVisitor: bee_sink is required");
        m_sink = Factory::find_tn<IBeeSink>(sink);
        if (cfg["beam_window"].isArray() && cfg["beam_window"].size() == 2) {
            m_bw_low = cfg["beam_window"][0].asDouble();
            m_bw_high = cfg["beam_window"][1].asDouble();
        }
        m_coords.clear();
        if (cfg["coords"].isArray() && cfg["coords"].size() == 3) {
            for (const auto& c : cfg["coords"]) m_coords.push_back(c.asString());
        }
    }

    virtual void visit(Ensemble& ensemble) const {
        const int index = ensemble.bee_index();
        if (index < 0) {
            raise<ValueError>("TaggerBeeVisitor: the Ensemble carries no Bee event index "
                              "(run it inside a MultiAlgBlobClustering pipeline)");
        }
        const int run = ensemble.rse_valid() ? ensemble.runNo() : ensemble.get_scalar<int>("runNo", 0);
        const int sub = ensemble.rse_valid() ? ensemble.subRunNo() : ensemble.get_scalar<int>("subRunNo", 0);
        const int evt = ensemble.rse_valid() ? ensemble.eventNo() : ensemble.get_scalar<int>("eventNo", 0);
        Bee::Points stm(m_detector, "tagger_stm"); stm.rse(run, sub, evt);
        Bee::Points tgm(m_detector, "tagger_tgm"); tgm.rse(run, sub, evt);
        Bee::Points fc (m_detector, "tagger_fc");  fc.rse(run, sub, evt);
        Bee::Points lm (m_detector, "tagger_lm");  lm.rse(run, sub, evt);

        auto groupings = ensemble.with_name(m_grouping);
        size_t ncand = 0, nskip = 0;
        if (!groupings.empty() && groupings.at(0)) {
            auto* gnode = groupings.at(0)->node();
            for (auto* cnode : gnode->children()) {
                bool candidate = false;
                int tag_stm = 0, tag_tgm = 0, tag_fc = 0, tag_lm = 0;
                {
                    auto cit = cnode->value.local_pcs().find("cluster_scalar");
                    if (cit == cnode->value.local_pcs().end()) continue;
                    auto& cs = cit->second;
                    auto rf = [&](const char* k) -> int {
                        auto a = cs.get(k); return (a && a->size_major() > 0 && a->elements<int>()[0] != 0) ? 1 : 0; };
                    const int is_main = rf("flag_main_cluster");
                    double t0 = 0;
                    { auto a = cs.get("cluster_t0"); if (a && a->size_major() > 0) t0 = a->elements<double>()[0]; }
                    candidate = is_main && t0 >= m_bw_low && t0 < m_bw_high;
                    tag_stm = rf("flag_STM");
                    tag_tgm = rf("flag_TGM");
                    tag_fc = rf("flag_FC");
                    { auto a = cs.get("lm_flag");
                      tag_lm = (a && a->size_major() > 0 && a->elements<int>()[0] == 2) ? 1 : 0; }
                }
                if (!candidate) continue;
                ++ncand;
                for (auto* bnode : cnode->children()) {
                    auto& lpcs = bnode->value.local_pcs();
                    auto sit = lpcs.find("scalar");
                    auto dit = lpcs.find("3d");
                    if (sit == lpcs.end() || dit == lpcs.end() || dit->second.size_major() == 0) continue;
                    auto& d3 = dit->second;
                    auto ax0 = d3.get("x"), ay0 = d3.get("y"), az0 = d3.get("z");
                    auto aq = sit->second.get("charge");
                    if (!ax0 || !ay0 || !az0 || !aq || aq->size_major() == 0
                        || ay0->size_major() != ax0->size_major() || az0->size_major() != ax0->size_major()) {
                        ++nskip;
                        continue;
                    }
                    const auto x = ax0->elements<double>();
                    const auto y = ay0->elements<double>();
                    const auto z = az0->elements<double>();
                    std::vector<double> txv, tyv, tzv;
                    if (m_coords.size() == 3) {
                        auto ax = d3.get(m_coords[0]), ay = d3.get(m_coords[1]), az = d3.get(m_coords[2]);
                        if (ax && ay && az && ax->size_major() == x.size() && ay->size_major() == x.size()
                            && az->size_major() == x.size()) {
                            auto sx = ax->elements<double>(); txv.assign(sx.begin(), sx.end());
                            auto sy = ay->elements<double>(); tyv.assign(sy.begin(), sy.end());
                            auto sz = az->elements<double>(); tzv.assign(sz.begin(), sz.end());
                        }
                    }
                    const bool have_tc = !txv.empty();
                    const double q = aq->elements<double>()[0];
                    const double qpp = x.size() ? std::max(q / x.size(), 1.0) : 1.0;
                    for (size_t i = 0; i < x.size(); ++i) {
                        const Point pt = have_tc ? Point(txv[i], tyv[i], tzv[i]) : Point(x[i], y[i], z[i]);
                        stm.append(pt, qpp, tag_stm, 0);
                        tgm.append(pt, qpp, tag_tgm, 0);
                        fc.append(pt, qpp, tag_fc, 0);
                        lm.append(pt, qpp, tag_lm, 0);
                    }
                }
            }
        }
        if (nskip) {
            log->warn("event ({},{},{}): skipped {} blob(s) lacking 3d x/y/z or scalar charge", run, sub, evt, nskip);
        }
        for (const auto* obj : {&stm, &tgm, &fc, &lm}) m_sink->write(*obj, index, run, sub, evt);
        log->debug("event ({},{},{}) bee index {}: {} beam-window candidate cluster(s)", run, sub, evt, index, ncand);
    }

private:
    std::string m_grouping{"live"};
    std::string m_detector{"sbnd"};
    double m_bw_low{0}, m_bw_high{0};
    std::vector<std::string> m_coords;
    IBeeSink::pointer m_sink;
};
