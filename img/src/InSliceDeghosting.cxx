
#include "WireCellImg/InSliceDeghosting.h"
#include "WireCellImg/CSGraph.h"
#include "WireCellImg/GeomClusteringUtil.h"
#include "WireCellAux/SimpleCluster.h"
#include "WireCellAux/ClusterHelpers.h"
#include "WireCellUtil/GraphTools.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/TimeKeeper.h"

#include <iterator>
#include <chrono>
#include <fstream>

WIRECELL_FACTORY(InSliceDeghosting, WireCell::Img::InSliceDeghosting, WireCell::INamed, WireCell::IClusterFilter,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Img;
using namespace WireCell::Aux;

Img::InSliceDeghosting::InSliceDeghosting()
  : Aux::Logger("InSliceDeghosting", "img")
{
}

Img::InSliceDeghosting::~InSliceDeghosting() {}

WireCell::Configuration Img::InSliceDeghosting::default_configuration() const
{
    WireCell::Configuration cfg;
    cfg["dryrun"] = m_dryrun;
    cfg["good_blob_charge_th"] = m_good_blob_charge_th;
    cfg["clustering_policy"] = m_clustering_policy;
    cfg["config_round"] = m_config_round;
    cfg["deghost_th"] = m_deghost_th;
    cfg["deghost_th1"] = m_deghost_th1;
    return cfg;
}

void Img::InSliceDeghosting::configure(const WireCell::Configuration& cfg)
{
    m_dryrun = get<bool>(cfg, "dryrun", m_dryrun);
    m_good_blob_charge_th = get<double>(cfg, "good_blob_charge_th", m_good_blob_charge_th);
    m_clustering_policy = get<std::string>(cfg, "clustering_policy", m_clustering_policy);
    m_config_round = get<int>(cfg, "config_round", m_config_round);
    m_deghost_th = get<float>(cfg, "deghost_th", m_deghost_th);
    m_deghost_th1 = get<float>(cfg, "deghost_th1", m_deghost_th1);

    log->debug("{}", cfg);
}

namespace {
    template <class Map, class Key, class Pos>
    void tag(Map& m, const Key& k, const Pos& p, const bool tag)
    {
        auto iter = m.find(k);
        if (iter == m.end()) m[k] = 0;
        int& pack = m[k];
        if (tag)
            pack |= (1 << p);
        else
            pack &= (0 << p);
    }

    template <class Map, class Key, class Pos>
    bool tag(const Map& m, const Key& k, const Pos& p)
    {
        auto iter = m.find(k);
        if (iter == m.end()) return false;
        int pack = iter->second >> p;
        return (pack & 1);
    }

    std::string dump_blob_tags(const InSliceDeghosting::vertex_tags_t& blob_tags)
    {
        std::stringstream ss;
        ss << "blob_tags ";
        std::unordered_map<int, size_t> counters;
        for (const auto& [vtx, pack] : blob_tags) {
            for (int pos = static_cast<int>(InSliceDeghosting::GOOD);
                 pos != static_cast<int>(InSliceDeghosting::TO_BE_REMOVED + 1); ++pos) {
                if (tag(blob_tags, vtx, pos)) counters[pos] += 1;
            }
        }
        for (const auto& [tag, counter] : counters) {
            ss << " " << tag << " " << counter;
        }
        return ss.str();
    }

    /// b -> connected channel idents
    std::unordered_map<WireCell::WirePlaneLayer_t, std::set<int> > connected_channels(const cluster_graph_t& cg,
                                                                                      const cluster_vertex_t& bdesc)
    {
        std::unordered_map<WireCell::WirePlaneLayer_t, std::set<int> > cidents;
        for (const auto& wdesc : neighbors_oftype<cluster_node_t::wire_t>(cg, bdesc)) {
            const auto iwire = get<cluster_node_t::wire_t>(cg[wdesc].ptr);
            cidents[iwire->planeid().layer()].insert(iwire->channel());
        }
        return cidents;
    }

    bool adjacent(std::unordered_map<WireCell::WirePlaneLayer_t, std::set<int> >& cid1,
                  std::unordered_map<WireCell::WirePlaneLayer_t, std::set<int> >& cid2)
    {
        std::map<WireCell::WirePlaneLayer_t, std::pair<int, int> > map1_plane_chs, map2_plane_chs;
        std::map<WireCell::WirePlaneLayer_t, int> map_plane_score;

        for (auto it = cid1.begin(); it != cid1.end(); it++) {
            map_plane_score[it->first] = 0;
            map1_plane_chs[it->first] = std::make_pair(*it->second.begin(), *it->second.rbegin());
        }

        for (auto it = cid2.begin(); it != cid2.end(); it++) {
            map2_plane_chs[it->first] = std::make_pair(*it->second.begin(), *it->second.rbegin());
        }

        int sum_score = 0;
        for (auto it = map_plane_score.begin(); it != map_plane_score.end(); it++) {
            if (map1_plane_chs[it->first].first == map2_plane_chs[it->first].second + 1 ||
                map2_plane_chs[it->first].first == map1_plane_chs[it->first].second + 1) {
                map_plane_score[it->first] = 1;
            }
            else if (map1_plane_chs[it->first].first <= map2_plane_chs[it->first].second &&
                     map2_plane_chs[it->first].first <= map1_plane_chs[it->first].second) {
                map_plane_score[it->first] = 2;
            }

            if (map_plane_score[it->first] == 0) return false;
            sum_score += map_plane_score[it->first];
        }

        if (sum_score >= 5) {
            return true;
        }
        else {
            return false;
        }
    }

}  // namespace

void InSliceDeghosting::blob_quality_ident(const cluster_graph_t& cg, vertex_tags_t& blob_tags)
{
    // loop over blob
    for (const auto& vtx : GraphTools::mir(boost::vertices(cg))) {
        const auto& node = cg[vtx];
        if (node.code() != 'b') {
            continue;
        }
        const auto iblob = get<cluster_node_t::blob_t>(node.ptr);
        //        double current_t = iblob->slice()->start();
        double current_q = iblob->value();  // charge ...
        if (current_q > m_good_blob_charge_th) {
            tag(blob_tags, vtx, GOOD, true);
            tag(blob_tags, vtx, POTENTIAL_GOOD, true);
        }
        //        else {
        //    blob_tags.insert({vtx, BAD});
        //}

        // loop over the b-b connections ...
        /// XIN: add between slice b-b connection info.
        /// follow b-b edges -> time, charge of connected blobs
        // std::vector<double> fronts, backs;
        for (auto v : neighbors_oftype<cluster_node_t::blob_t>(cg, vtx)) {
            const auto ib = get<cluster_node_t::blob_t>(cg[v].ptr);
            if (ib->value() > m_good_blob_charge_th) {
                tag(blob_tags, vtx, POTENTIAL_GOOD, true);
            }
            //  if (ib->slice()->start() > current_t) {
            //    backs.push_back(ib->value());
            // }
            // if (ib->slice()->start() < current_t) {
            //    fronts.push_back(ib->value());
            //}
        }

        /// DEBUGONLY:
        // log->debug("{} fronts {} backs {}", vtx, fronts.size(), backs.size());
    }
}
void InSliceDeghosting::local_deghosting1(const cluster_graph_t& cg, vertex_tags_t& blob_tags)
{
    // loop over time slices
    for (const auto& svtx : GraphTools::mir(boost::vertices(cg))) {
        if (cg[svtx].code() != 's') {
            continue;
        }

        std::unordered_map<int, float> wire_score_map;
        // Loop over the neighbor blobs of current slice. Keep the largest one.
        std::unordered_map<size_t, std::unordered_set<cluster_vertex_t> > view_groups;
        // std::unordered_map<RayGrid::grid_index_t, size_t> grid_occupancy;
        std::unordered_map<cluster_vertex_t, std::unordered_set<WireCell::WirePlaneLayer_t> > blob_planes;
        std::unordered_map<cluster_vertex_t, std::unordered_map<WireCell::WirePlaneLayer_t, std::set<int> > >
            blob_channels;

        for (auto bedge : GraphTools::mir(boost::out_edges(svtx, cg))) {
            auto bvtx = boost::target(bedge, cg);
            if (cg[bvtx].code() != 'b') {
                THROW(ValueError() << errmsg{
                          String::format("'s' node connnected to '%s' not implemented!", cg[bvtx].code())});
            }

            /// group blobs into #live view group
            /// ASSUMPTION: connected meas node == # live views
            auto meas = neighbors_oftype<cluster_node_t::meas_t>(cg, bvtx);
            //	    std::unordered_set<WireCell::WirePlaneLayer_t> live_planes;
            view_groups[meas.size()].insert(bvtx);
            // how do we know which plane meas_t is bad ... ???
            for (auto mea : meas) {
                const auto imeasure = get<cluster_node_t::meas_t>(cg[mea].ptr);
                // log->debug("planeid {}",imeasure->planeid());
                blob_planes[bvtx].insert(imeasure->planeid().layer());
                // live_planes.insert(imeasure->planeid().layer());
            }

            blob_channels[bvtx] = connected_channels(cg, bvtx);
            if (meas.size() == 3) {  // all three plane live ...
                auto& cidents = blob_channels[bvtx];
                for (auto it = cidents.begin(); it != cidents.end(); it++) {
                    for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++) {
                        if (wire_score_map.find(*it1) == wire_score_map.end()) {
                            wire_score_map[*it1] = 1;
                        }
                        else {
                            wire_score_map[*it1]++;
                        }
                    }
                }
            }
        }  // loop over blob ...

        // loop over two-wire blobs ...
        for (auto it = view_groups[2].begin(); it != view_groups[2].end(); it++) {
            auto& two_view_chs = blob_channels[*it];  // connected_channels(cg, *it);
            std::unordered_map<WireCell::WirePlaneLayer_t, bool> flag_plane;
            auto& live_planes = blob_planes[*it];
            for (auto it1 = two_view_chs.begin(); it1 != two_view_chs.end(); it1++) {
                if (live_planes.find(it1->first) != live_planes.end()) {
                    flag_plane[it1->first] = true;
                }
                else {
                    flag_plane[it1->first] = false;
                }
                if (flag_plane[it1->first]) {
                    for (auto it2 = it1->second.begin(); it2 != it1->second.end(); it2++) {
                        if (wire_score_map.find(*it2) == wire_score_map.end()) {
                            wire_score_map[*it2] = 1;
                        }
                        else {
                            wire_score_map[*it2]++;
                        }
                    }
                }
            }
        }

        std::unordered_set<cluster_vertex_t> cannot_remove;
        for (auto it = view_groups[2].begin(); it != view_groups[2].end(); it++) {
            auto& two_view_chs = blob_channels[*it];  // connected_channels(cg, *it);
            int count = 0;
            for (auto it1 = view_groups[3].begin(); it1 != view_groups[3].end(); it1++) {
                auto& three_view_chs = blob_channels[*it1];  // connected_channels(cg, *it1);
                if (tag(blob_tags, *it1, POTENTIAL_GOOD)) {
                    if (adjacent(two_view_chs, three_view_chs)) count++;
                    if (count == 2) {
                        cannot_remove.insert(*it);
                        break;
                    }
                }
            }
        }

        std::unordered_map<cluster_vertex_t, float> blob_high_score_map;

        for (auto it = view_groups[2].begin(); it != view_groups[2].end(); it++) {
            auto& two_view_chs = blob_channels[*it];  // connected_channels(cg, *it);
            std::unordered_map<WireCell::WirePlaneLayer_t, bool> flag_plane;
            auto& live_planes = blob_planes[*it];

            std::map<WireCell::WirePlaneLayer_t, float> score_plane;
            for (auto it1 = two_view_chs.begin(); it1 != two_view_chs.end(); it1++) {
                if (live_planes.find(it1->first) != live_planes.end()) {
                    flag_plane[it1->first] = true;
                }
                else {
                    flag_plane[it1->first] = false;
                }
                score_plane[it1->first] = -1;
                if (flag_plane[it1->first]) {
                    float sum1 = 0;
                    float sum2 = 0;
                    for (auto it2 = it1->second.begin(); it2 != it1->second.end(); it2++) {
                        sum2++;
                        sum1 += 1. / wire_score_map[*it2];
                    }
                    score_plane[it1->first] = sum1 / sum2;
                }
            }

            blob_high_score_map[*it] =
                std::max(std::max(score_plane[kUlayer], score_plane[kVlayer]), score_plane[kWlayer]);
        }

        for (auto it = view_groups[3].begin(); it != view_groups[3].end(); it++) {
            auto& two_view_chs = blob_channels[*it];  // connected_channels(cg, *it);
            std::unordered_map<WireCell::WirePlaneLayer_t, bool> flag_plane;
            auto& live_planes = blob_planes[*it];

            std::map<WireCell::WirePlaneLayer_t, float> score_plane;
            for (auto it1 = two_view_chs.begin(); it1 != two_view_chs.end(); it1++) {
                if (live_planes.find(it1->first) != live_planes.end()) {
                    flag_plane[it1->first] = true;
                }
                else {
                    flag_plane[it1->first] = false;
                }
                score_plane[it1->first] = -1;
                if (flag_plane[it1->first]) {
                    float sum1 = 0;
                    float sum2 = 0;
                    for (auto it2 = it1->second.begin(); it2 != it1->second.end(); it2++) {
                        sum2++;
                        sum1 += 1. / wire_score_map[*it2];
                    }
                    score_plane[it1->first] = sum1 / sum2;
                }
            }

            blob_high_score_map[*it] =
                std::max(std::max(score_plane[kUlayer], score_plane[kVlayer]), score_plane[kWlayer]);
            if (tag(blob_tags, *it, POTENTIAL_GOOD)) blob_high_score_map[*it] = 1;
        }

        for (auto it = view_groups[2].begin(); it != view_groups[2].end(); it++) {
            auto& two_view_chs = blob_channels[*it];  // connected_channels(cg, *it);
            std::unordered_map<WireCell::WirePlaneLayer_t, bool> flag_plane;
            const auto iblob1 = get<cluster_node_t::blob_t>(cg[*it].ptr);
            double current_q1 = iblob1->value();  // charge ...

            auto& live_planes = blob_planes[*it];
            for (auto it1 = two_view_chs.begin(); it1 != two_view_chs.end(); it1++) {
                if (live_planes.find(it1->first) != live_planes.end()) {
                    flag_plane[it1->first] = true;
                }
                else {
                    flag_plane[it1->first] = false;
                }
            }

            int count = 0;
            // find the merged wire associated with this blob ...
            auto meas = neighbors_oftype<cluster_node_t::meas_t>(cg, *it);
            for (auto mea : meas) {
                const auto imeasure = get<cluster_node_t::meas_t>(cg[mea].ptr);
                // blob_planes[bvtx].insert(imeasure->planeid().layer());
                if (flag_plane[imeasure->planeid().layer()]) {
                    int mcell_lwire = *two_view_chs[imeasure->planeid().layer()].begin();
                    int mcell_hwire = *two_view_chs[imeasure->planeid().layer()].rbegin();

                    // find the merged blob associated with the merged wire ...
                    auto blobs = neighbors_oftype<cluster_node_t::blob_t>(cg, mea);
                    for (auto blob : blobs) {
                        if (blob == *it) continue;
                        auto& blob_chs = blob_channels[blob];  // connected_channels(cg, blob);
                        const auto iblob2 = get<cluster_node_t::blob_t>(cg[blob].ptr);
                        double current_q2 = iblob2->value();  // charge ...

                        if (blob_high_score_map[blob] > m_deghost_th) {
                            int mcell1_lwire = *blob_chs[imeasure->planeid().layer()].begin();
                            int mcell1_hwire = *blob_chs[imeasure->planeid().layer()].rbegin();
                            int min_wire = std::max(mcell_lwire, mcell1_lwire);
                            int max_wire = std::min(mcell_hwire, mcell1_hwire);

                            // /// DEBUGONLY:
                            // const auto islice = get<cluster_node_t::slice_t>(cg[svtx].ptr);
                            // log->debug("start: {} three {} two {} cnrm {} score {} {} {} {} {} {} ", islice->start()
                            // /   islice->span(), 	       view_groups[3].size(), view_groups[2].size(), cannot_remove.size(),
                            // 	       current_q1, imeasure->planeid().layer(),
                            // 	       current_q2, min_wire, max_wire, (max_wire - min_wire + 1.0) / (mcell_hwire -
                            // mcell_lwire + 1.0) );

                            if ((max_wire - min_wire + 1.0) / (mcell_hwire - mcell_lwire + 1.0) >= m_deghost_th &&
                                current_q2 > current_q1 * m_deghost_th) {
                                count++;
                                break;
                            }
                        }
                    }
                }
            }

            if (count == 2 && cannot_remove.find(*it) == cannot_remove.end()) tag(blob_tags, *it, TO_BE_REMOVED, true);
        }

    }  // loop over slice
}

void InSliceDeghosting::local_deghosting(const cluster_graph_t& cg, vertex_tags_t& blob_tags)
{
    //    int count_keeper = 0;
    // loop over slices ...
    for (const auto& svtx : GraphTools::mir(boost::vertices(cg))) {
        if (cg[svtx].code() != 's') {
            continue;
        }

        // saved the used channels ...
        std::unordered_map<WireCell::WirePlaneLayer_t, std::set<int> > used_plane_channels;
        // Loop over the neighbor blobs of current slice. Keep the largest one.
        std::unordered_map<size_t, std::unordered_set<cluster_vertex_t> > view_groups;
        // std::unordered_map<RayGrid::grid_index_t, size_t> grid_occupancy;
        std::unordered_map<cluster_vertex_t, std::unordered_set<WireCell::WirePlaneLayer_t> > blob_planes;
        std::unordered_map<cluster_vertex_t, std::unordered_map<WireCell::WirePlaneLayer_t, std::set<int> > >
            blob_channels;

        for (auto bedge : GraphTools::mir(boost::out_edges(svtx, cg))) {
            /// XIN: bvtx -> RayGrid::Strip::grid_range_t, {2views}, {3views}
            /// CHECK: per slice map {RayGrid::grid_index_t -> float}
            auto bvtx = boost::target(bedge, cg);
            if (cg[bvtx].code() != 'b') {
                THROW(ValueError() << errmsg{
                          String::format("'s' node connnected to '%s' not implemented!", cg[bvtx].code())});
            }

            /// group blobs into #live view group
            /// ASSUMPTION: connected meas node == # live views
            auto meas = neighbors_oftype<cluster_node_t::meas_t>(cg, bvtx);
            //	    std::unordered_set<WireCell::WirePlaneLayer_t> live_planes;
            view_groups[meas.size()].insert(bvtx);
            // how do we know which plane meas_t is bad ... ???
            for (auto mea : meas) {
                const auto imeasure = get<cluster_node_t::meas_t>(cg[mea].ptr);
                // log->debug("planeid {}",imeasure->planeid());
                blob_planes[bvtx].insert(imeasure->planeid().layer());
                // live_planes.insert(imeasure->planeid().layer());
            }

            /// EXAMPLE: add up #blobs using each RayGrid::grid_index_t
            // RayGrid ... layer_index_t ...
            // const IBlob::pointer iblob = get<cluster_node_t::blob_t>(cg[bvtx].ptr);
            // const auto& strips = iblob->shape().strips();
            // for (const RayGrid::Strip& strip : strips) {
            //     for (RayGrid::grid_index_t gridx = strip.bounds.first; gridx < strip.bounds.second; ++gridx) {
            //         if (grid_occupancy.find(gridx) != grid_occupancy.end()) {
            //             grid_occupancy[gridx]++;
            //         }
            //         else {
            //             grid_occupancy[gridx] = 1;
            //         }
            //     }
            // }

            /// EXAMPLE: b -> connected channel idents
            // how do we know which plane is channel is?
            blob_channels[bvtx] = connected_channels(cg, bvtx);
            auto& cidents = blob_channels[bvtx];

            if (meas.size() == 3 && tag(blob_tags, bvtx, GOOD)) {  // all three plane live ...
                for (auto it = cidents.begin(); it != cidents.end(); it++) {
                    // int ich = it->first;
                    //		used_plane_channels[it->second].insert(it->first);
                    used_plane_channels[it->first].insert(it->second.begin(), it->second.end());
                }
            }
            //            log->debug("cidents {}", cidents.size());

            /// EXAMPLEONLY: temp tag all as TO_BE_REMOVED
            // blob_tags.insert({bvtx, TO_BE_REMOVED});
        }

        std::unordered_set<cluster_vertex_t> cannot_remove;
        for (auto it = view_groups[2].begin(); it != view_groups[2].end(); it++) {
            auto& two_view_chs = blob_channels[*it];  // connected_channels(cg, *it);
            int count = 0;
            for (auto it1 = view_groups[3].begin(); it1 != view_groups[3].end(); it1++) {
                auto& three_view_chs = blob_channels[*it1];  // connected_channels(cg, *it1);
                if (tag(blob_tags, *it1, GOOD)) {
                    if (adjacent(two_view_chs, three_view_chs)) count++;
                    if (count == 2) {
                        cannot_remove.insert(*it);
                    }
                }
            }
        }

        // loop over two-good wire cells
        for (auto it = view_groups[2].begin(); it != view_groups[2].end(); it++) {
            if (!tag(blob_tags, *it, POTENTIAL_GOOD)) {
                auto& two_view_chs = blob_channels[*it];  // connected_channels(cg, *it);
                std::unordered_map<WireCell::WirePlaneLayer_t, bool> flag_plane;
                auto& live_planes = blob_planes[*it];
                bool save_flag = false;
                for (auto it1 = two_view_chs.begin(); it1 != two_view_chs.end(); it1++) {
                    if (live_planes.find(it1->first) != live_planes.end()) {
                        flag_plane[it1->first] = true;
                    }
                    else {
                        flag_plane[it1->first] = false;
                    }
                    if (flag_plane[it1->first] && !save_flag) {
                        int nwire = 0;
                        int nwire_common = 0;
                        auto& used_chs = used_plane_channels[it1->first];
                        for (auto it2 = it1->second.begin(); it2 != it1->second.end(); it2++) {
                            nwire++;
                            if (used_chs.find(*it2) != used_chs.end()) nwire_common++;
                        }
                        if (nwire_common != nwire) save_flag = true;
                    }
                }

                if ((!save_flag) && (cannot_remove.find(*it) == cannot_remove.end())) {
                    tag(blob_tags, *it, TO_BE_REMOVED, true);
                }
            }
        }

        std::unordered_map<int, float> wire_score_map;
        for (auto it = view_groups[2].begin(); it != view_groups[2].end(); it++) {
            if (tag(blob_tags, *it, TO_BE_REMOVED)) continue;
            auto& two_view_chs = blob_channels[*it];  // connected_channels(cg, *it);
            std::unordered_map<WireCell::WirePlaneLayer_t, bool> flag_plane;
            auto& live_planes = blob_planes[*it];
            // std::unordered_map<WireCell::WirePlaneLayer_t, float> score_plane;
            for (auto it1 = two_view_chs.begin(); it1 != two_view_chs.end(); it1++) {
                if (live_planes.find(it1->first) != live_planes.end()) {
                    flag_plane[it1->first] = true;
                }
                else {
                    flag_plane[it1->first] = false;
                }
                // score_plane[it1->first] = -1;
                if (flag_plane[it1->first]) {
                    for (auto it2 = it1->second.begin(); it2 != it1->second.end(); it2++) {
                        if (wire_score_map.find(*it2) == wire_score_map.end()) {
                            wire_score_map[*it2] = 1;
                        }
                        else {
                            wire_score_map[*it2]++;
                        }
                    }
                }
            }
        }

        for (auto it = view_groups[2].begin(); it != view_groups[2].end(); it++) {
            if (tag(blob_tags, *it, TO_BE_REMOVED)) continue;
            auto& two_view_chs = blob_channels[*it];  // connected_channels(cg, *it);
            std::unordered_map<WireCell::WirePlaneLayer_t, bool> flag_plane;
            auto& live_planes = blob_planes[*it];
            std::unordered_map<WireCell::WirePlaneLayer_t, float> score_plane;
            for (auto it1 = two_view_chs.begin(); it1 != two_view_chs.end(); it1++) {
                if (live_planes.find(it1->first) != live_planes.end()) {
                    flag_plane[it1->first] = true;
                }
                else {
                    flag_plane[it1->first] = false;
                }
                score_plane[it1->first] = -1;
                if (flag_plane[it1->first]) {
                    float sum1 = 0;
                    float sum2 = 0;
                    auto& used_chs = used_plane_channels[it1->first];
                    for (auto it2 = it1->second.begin(); it2 != it1->second.end(); it2++) {
                        sum2++;
                        if (used_chs.find(*it2) == used_chs.end()) {
                            sum1 += 1. / wire_score_map[*it2];
                        }
                    }
                    score_plane[it1->first] = sum1 / sum2;
                }
            }
            if (score_plane[WireCell::kUlayer] <= m_deghost_th1 && score_plane[WireCell::kVlayer] <= m_deghost_th1 &&
                score_plane[WireCell::kWlayer] <= m_deghost_th1 && (cannot_remove.find(*it) == cannot_remove.end()) &&
                (!tag(blob_tags, *it, POTENTIAL_GOOD))) {
                tag(blob_tags, *it, POTENTIAL_BAD, true);

                /// DEBUGONLY:
                // const auto islice = get<cluster_node_t::slice_t>(cg[svtx].ptr);
                // log->debug("start: {} three {} two {} used {} {} {} cnrm {} score {} {} {}", islice->start() /
                // islice->span(),
                //            view_groups[3].size(), view_groups[2].size(), used_plane_channels[kUlayer].size(),
                //            used_plane_channels[kVlayer].size(), used_plane_channels[kWlayer].size(),
                //            cannot_remove.size(), score_plane[WireCell::kUlayer], score_plane[WireCell::kVlayer],
                //            score_plane[WireCell::kWlayer]);
            }
        }

        /// EXAMPLEONLY: rm 3views from TO_BE_REMOVED
        /// this is to prune all 3views from blob_tags, for single key deletion, use equal_range.
        //  for (auto i = blob_tags.begin(), last = blob_tags.end(); i != last;) {
        //    auto const& [key, value] = *i;
        // if (view_groups[3].find(key) != view_groups[3].end() && value == TO_BE_REMOVED) {
        //     i = blob_tags.erase(i);
        //     ++count_keeper;
        // }
        // else {
        //     ++i;
        // }
        // }
    }
    //    log->debug("count_keeper {}", count_keeper);
}
bool InSliceDeghosting::operator()(const input_pointer& in, output_pointer& out)
{
    out = nullptr;
    if (!in) {
        log->debug("EOS");
        return true;
    }

    TimeKeeper tk(fmt::format("InSliceDeghosting"));

    const auto in_graph = in->graph();
    log->debug("in_graph: {}", dumps(in_graph));

    // blob desc -> quality tag
    vertex_tags_t blob_tags;

    WireCell::cluster_graph_t cg_new_bb;

    if (m_config_round == 1) {
        /// 1, Ident good/bad blob groups. in: ICluster out: blob_quality_tags {blob_desc -> quality_tag}
        log->debug(tk(fmt::format("start blob tagging")));
        blob_quality_ident(in_graph, blob_tags);
        /// DEBUGONLY:
        log->debug("conf 1 step 1 {}", dump_blob_tags(blob_tags));

        /// 2, Local (in-slice) de-ghosting. in: ICluster + blob_quality_tags out: updated blob_quality_tags (remove or
        /// not)
        log->debug(tk(fmt::format("start local (in-slice) de-ghosting")));
        local_deghosting(in_graph, blob_tags);
        /// DEBUGONLY:
        log->debug("conf 1 step 2 {}", dump_blob_tags(blob_tags));

        /// 3, delete some blobs. in: ICluster + blob_quality_tags out: filtered ICluster
        log->debug(tk(fmt::format("start delete some blobs")));
        using VFiltered =
            typename boost::filtered_graph<cluster_graph_t, boost::keep_all, std::function<bool(cluster_vertex_t)> >;
        VFiltered fg_rm_bad_blobs(in_graph, {}, [&](auto vtx) { return !tag(blob_tags, vtx, TO_BE_REMOVED); });
        /// DEBUGONLY:
        log->debug("conf 1 step 3 {}", dump_blob_tags(blob_tags));

        /// 4, in-group clustering. in: ICluster + blob_quality_tags out: filtered ICluster
        /// TODO: do we need to call copy_graph?
        log->debug(tk(fmt::format("start per group clustering")));
        using desc_map_t = std::unordered_map<cluster_vertex_t, cluster_vertex_t>;
        using pm_desc_map_t = boost::associative_property_map<desc_map_t>;
        desc_map_t o2c1;
        WireCell::cluster_graph_t cg_old_bb;
        boost::copy_graph(fg_rm_bad_blobs, cg_old_bb, boost::orig_to_copy(pm_desc_map_t(o2c1)));
        log->debug("cg_old_bb: {}", dumps(cg_old_bb));
        /// rm old b-b edges
        using EFiltered =
            typename boost::filtered_graph<cluster_graph_t, std::function<bool(cluster_edge_t)>, boost::keep_all>;
        EFiltered fg_no_bb(cg_old_bb,
                           [&](auto edge) {
                               auto source = boost::source(edge, cg_old_bb);
                               auto target = boost::target(edge, cg_old_bb);
                               if (cg_old_bb[source].code() == 'b' and cg_old_bb[target].code() == 'b') {
                                   return false;
                               }
                               return true;
                           },
                           {});

        desc_map_t o2c2;
        boost::copy_graph(fg_no_bb, cg_new_bb, boost::orig_to_copy(pm_desc_map_t(o2c2)));

        /// reverse the mapping
        /// TODO: need protection?
        desc_map_t c2o;
        for (const auto [o, c1] : o2c1) {
            const auto c2 = o2c2[c1];
            c2o[c2] = o;
        }
        log->debug("cg_new_bb: {}", dumps(cg_new_bb));
        /// make new b-b edges within groups
        /// four separated groups, good&bad, good&good, bad&bad  (two of these might not be useful ...)
        std::vector<gc_filter_t> filters = {
            [&](const cluster_vertex_t& vtx) {
                return tag(blob_tags, vtx, POTENTIAL_GOOD) && !tag(blob_tags, vtx, POTENTIAL_BAD);
            },
            [&](const cluster_vertex_t& vtx) {
                return tag(blob_tags, vtx, POTENTIAL_GOOD) && tag(blob_tags, vtx, POTENTIAL_BAD);
            },
            [&](const cluster_vertex_t& vtx) {
                return !tag(blob_tags, vtx, POTENTIAL_GOOD) && !tag(blob_tags, vtx, POTENTIAL_BAD);
            },
            [&](const cluster_vertex_t& vtx) {
                return !tag(blob_tags, vtx, POTENTIAL_GOOD) && tag(blob_tags, vtx, POTENTIAL_BAD);
            },
        };
        /// each filter adds new b-b edges on cg_new_bb
        /// this one is too slow and stuck for 0 - 9592 full tick range
        // for (auto filter : filters) {
        //     geom_clustering(cg_new_bb, m_clustering_policy, filter);
        // }

        /// condense multi-map to map
        std::unordered_map<cluster_vertex_t, int> groups;
        for (size_t groupid = 0; groupid < filters.size(); ++groupid) {
            for (const auto& vtx : oftype(cg_new_bb, 'b')) {
                if (filters[groupid](c2o[vtx])) groups.insert({vtx, groupid});
            }
        }
        /// DEBUGONLY:
        log->debug("conf 1 step 4 {}", dump_blob_tags(blob_tags));
        /// DEBUGONLY:
        {
            std::unordered_map<int, int> group_sizes;
            for (const auto& [desc, id] : groups) {
                group_sizes[id] += 1;
            }
            for (const auto& [id, count] : group_sizes) {
                log->debug("group_sizes: {} {}", id, count);
            }
        }

        grouped_geom_clustering(cg_new_bb, m_clustering_policy, groups);
        log->debug(tk(fmt::format("end")));
        log->debug("cg_new_bb: {}", dumps(cg_new_bb));
    }
    else if (m_config_round == 2) {
        // after charge solving ...
        blob_quality_ident(in_graph, blob_tags);
        /// DEBUGONLY:
        log->debug("conf 2 {}", dump_blob_tags(blob_tags));

        // second round of local deghosting ...
        local_deghosting1(in_graph, blob_tags);

        using VFiltered =
            typename boost::filtered_graph<cluster_graph_t, boost::keep_all, std::function<bool(cluster_vertex_t)> >;
        VFiltered fg_rm_bad_blobs(in_graph, {}, [&](auto vtx) { return !tag(blob_tags, vtx, TO_BE_REMOVED); });
        // do we have to copy this every time ???
        boost::copy_graph(fg_rm_bad_blobs, cg_new_bb);
        log->debug("cg_new_bb 2nd: {}", dumps(cg_new_bb));
    }
    else if (m_config_round == 3) {
        // after charge solving ...
        blob_quality_ident(in_graph, blob_tags);
        /// DEBUGONLY:
        log->debug("conf 3 {}", dump_blob_tags(blob_tags));

        // at some points remove the Non-potential good ... (last round ...)
        using VFiltered =
            typename boost::filtered_graph<cluster_graph_t, boost::keep_all, std::function<bool(cluster_vertex_t)> >;
        VFiltered fg_rm_bad_blobs(in_graph, {}, [&](auto vtx) {
            if (in_graph[vtx].code() != 'b') return true;
            return tag(blob_tags, vtx, POTENTIAL_GOOD);
        });
        // do we have to copy this every time ???
        boost::copy_graph(fg_rm_bad_blobs, cg_new_bb);
        log->debug("cg_new_bb 3rd: {}", dumps(cg_new_bb));
    }

    /// output
    auto& out_graph = cg_new_bb;
    /// TODO: keep this?
    log->debug("in_graph: {}", dumps(in_graph));
    log->debug("out_graph: {}", dumps(out_graph));

    out = std::make_shared<Aux::SimpleCluster>(out_graph, in->ident());
    if (m_dryrun) {
        out = std::make_shared<Aux::SimpleCluster>(in_graph, in->ident());
    }
    return true;
}
