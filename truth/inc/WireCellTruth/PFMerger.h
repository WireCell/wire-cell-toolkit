/** PFMerger: G4 particles -> particle-flow (PF) nodes and classes.

    Header-only (owner decision O1).  Versioned: PFMergerParams is
    serialised into labels_meta.json by meta() so that a label set
    records the rules it was made with.

    Input: a TruthTable (every G4 particle kept) and, optionally, the
    visible charge per G4 track id measured in the reconstruction
    (sum of blobtruth q_bk; else the table's nvis_q from the
    SimEnergyDeposits).  Output: PF nodes and a per-G4-id map
    g4 -> (pf_node, cls, merged_into, instance_valid).

    Rules, in order, with their source (wcp-porting-img pr_ml/docs/19
    sec 3.6 and sec 4):

    1. Visible particles: e, mu, pi+-, K+-, Sigma+-, p, gamma, n
       (Pandora LArMCParticleHelper::IsVisible).  Every other particle
       (nuclei, pi0, K0, Lambda, neutrinos, ...) never gets a node;
       its ids map to the node of its nearest visible ancestor.
    2. Primary: the parent chain reaches the generator (mother 0)
       passing only through SPINE's invisible-parent list
       111/113/221/331/3212.
    3. Ionisation electrons (eIoni/muIoni/hIoni/ionIoni) from a
       track shorter than delta_len merge into the parent track;
       longer ones are their own node, class delta (Supera DeltaSize,
       pynuml muIoni -> muon).
    4. conv/compt/phot/eBrem/annihil/Rayl daughters of a shower
       particle merge into the shower root; so do ionisation electrons
       inside a shower (Supera merge order).
    5. A Decay or muMinusCaptureAtRest electron of a muon is a Michel,
       its own node (pynuml), when michel_own.
       With em_fold (default off) rule 4 also takes the EM daughters of
       a delta or Michel node: an electron's radiation stays with it
       (pr_ml doc 25 sec 2.8; off = a photon radiated by a delta or
       Michel electron is a shower root of its own).
       With michel_capture_ke_min > 0 (default 0 = off) a
       muMinusCaptureAtRest electron below that KE is not a Michel
       (the capture's atomic electrons); with michel_stopped_only
       (default off) a Decay electron is a Michel only when its muon
       stopped: never for a mu- (at rest a mu- decays under
       muMinusCaptureAtRest), and for a mu+ only when the decay waited
       at least michel_wait_min beyond the fastest possible flight
       (electron start_t - muon start_t - |muon end - start| / c).
       Both: pr_ml doc 25 sec 6 (owner 2026-10-24).  A rejected
       electron falls through to the rules below (shower root).
    6. A shower node with visible charge below compton_q is class LE
       with instance_valid 0 (SPINE: semantics kept, instance -1).
    7. A node survives when its visible charge >= q_min and its
       KE >= pf_ke_min and it is a primary, a track or delta longer
       than track_len_min, a shower root, or a Michel; otherwise it
       merges into its nearest surviving ancestor node.  A node
       without any ancestor node always survives.  An interaction with
       no surviving primary node carrying charge and starting within
       void_radius of its vertex is flagged void_root.
    8. Classes track / shower / michel / delta / le; shape precedence
       for shared points is Schema's (SPINE).

    Units: the TruthTable is in cm / ns / GeV (SBND attacher); the
    parameters are in WCT units.  Charges are electrons.
 */

#ifndef WIRECELLTRUTH_PFMERGER
#define WIRECELLTRUTH_PFMERGER

#include "WireCellTruth/Schema.h"
#include "WireCellTruth/TruthTable.h"
#include "WireCellUtil/Configuration.h"
#include "WireCellUtil/Units.h"

#include <cmath>
#include <cstdlib>
#include <functional>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace WireCell::Truth {

    struct PFMergerParams {
        int version{1};
        double pf_ke_min{10 * units::MeV};
        double q_min{10000.0};          // electrons; ~1 mm of a MIP
        double delta_len{1 * units::cm};
        double compton_q{200000.0};     // electrons; ~5 MeV
        double track_len_min{1 * units::cm};
        bool michel_own{true};
        double photon_prop{2.5 * units::cm};  // recorded; used by the matching of showers (WP0.2)
        double void_radius{1 * units::cm};
        /// Rule 4 also merges EM daughters into a delta or Michel node.
        /// Recorded in meta() only when on (labels_meta.json unchanged when off).
        bool em_fold{false};
        /// Rule 5: KE floor for muMinusCaptureAtRest electrons (0 = off).
        /// Recorded in meta() only when on.
        double michel_capture_ke_min{0.0};
        /// Rule 5: Decay electrons are Michels only for a stopped muon
        /// (mu- never; mu+ when the decay waited >= michel_wait_min).
        /// Recorded in meta() only when on.
        bool michel_stopped_only{false};
        double michel_wait_min{20 * units::ns};

        Configuration meta() const
        {
            Configuration m;
            m["version"] = version;
            m["pf_ke_min_MeV"] = pf_ke_min / units::MeV;
            m["q_min_e"] = q_min;
            m["delta_len_cm"] = delta_len / units::cm;
            m["compton_q_e"] = compton_q;
            m["track_len_min_cm"] = track_len_min / units::cm;
            m["michel_own"] = michel_own;
            m["photon_prop_cm"] = photon_prop / units::cm;
            m["void_radius_cm"] = void_radius / units::cm;
            m["visible_pdgs"] = "11 13 211 321 3222 3112 2212 22 2112 (Pandora IsVisible)";
            m["invisible_parents"] = "111 113 221 331 3212 (SPINE)";
            m["shape_precedence"] = "track > michel > shower > delta > le (SPINE)";
            if (em_fold) m["em_fold_delta_michel"] = true;
            if (michel_capture_ke_min > 0) m["michel_capture_ke_min_MeV"] = michel_capture_ke_min / units::MeV;
            if (michel_stopped_only) m["michel_stopped_only_wait_min_ns"] = michel_wait_min / units::ns;
            return m;
        }
        void configure(const Configuration& cfg)
        {
            pf_ke_min = get(cfg, "pf_ke_min", pf_ke_min);
            q_min = get(cfg, "q_min", q_min);
            delta_len = get(cfg, "delta_len", delta_len);
            compton_q = get(cfg, "compton_q", compton_q);
            track_len_min = get(cfg, "track_len_min", track_len_min);
            michel_own = get(cfg, "michel_own", michel_own);
            photon_prop = get(cfg, "photon_prop", photon_prop);
            void_radius = get(cfg, "void_radius", void_radius);
            em_fold = get(cfg, "em_fold", em_fold);
            michel_capture_ke_min = get(cfg, "michel_capture_ke_min", michel_capture_ke_min);
            michel_stopped_only = get(cfg, "michel_stopped_only", michel_stopped_only);
            michel_wait_min = get(cfg, "michel_wait_min", michel_wait_min);
        }
    };

    struct PFNode {
        int index{-1};          // position in PFResult::nodes
        int trackid{0};         // the root G4 id
        int pdg{0};
        Cls cls{Cls::unknown};
        int parent{-1};         // index of the parent node, -1 at the interaction root
        bool primary{false};
        int nu_row{-1};         // truth_nu row (interaction) of the root particle
        double ke{0};           // WCT units
        double length{0};       // chord start-end of the root particle, WCT units
        double q_vis{0};        // electrons, this node and every id merged into it
        bool survived{true};
        bool instance_valid{true};
        std::vector<int> ids;   // every G4 id mapped to this node, sorted
        std::string process;    // creation process of the root
    };

    struct G4Info {
        int pf_node{-1};
        Cls cls{Cls::unknown};
        int merged_into{0};     // trackid of the node root (self when own node)
        bool instance_valid{false};
    };

    struct PFResult {
        std::vector<PFNode> nodes;
        std::map<int, G4Info> g4;       // every trackid of the table
        std::set<int> void_roots;       // nu rows with no primary node at the vertex
        size_t nmerged_delta{0}, nmerged_shower{0}, nmerged_weak{0}, nmerged_invisible{0};

        const G4Info& info(int trackid) const
        {
            static const G4Info none;
            auto it = g4.find(trackid);
            return it == g4.end() ? none : it->second;
        }
    };

    class PFMerger {
      public:
        explicit PFMerger(const PFMergerParams& params = PFMergerParams()) : m_p(params) {}

        static bool is_visible(int pdg)
        {
            switch (std::abs(pdg)) {
            case 11: case 13: case 211: case 321: case 3222: case 3112: case 2212: case 22: case 2112:
                return true;
            default:
                return false;
            }
        }
        static bool is_skipped_parent(int pdg)
        {
            switch (std::abs(pdg)) {
            case 111: case 113: case 221: case 331: case 3212: return true;
            default: return false;
            }
        }
        static bool is_em(int pdg) { return std::abs(pdg) == 11 || pdg == 22; }
        static bool is_ionisation(const std::string& proc)
        {
            return proc == "eIoni" || proc == "muIoni" || proc == "hIoni" || proc == "ionIoni";
        }
        static bool is_shower_process(const std::string& proc)
        {
            return proc == "conv" || proc == "compt" || proc == "phot" || proc == "eBrem" ||
                proc == "annihil" || proc == "Rayl";
        }

        /// q_vis: measured visible charge per G4 id (electrons); when
        /// empty the table's nvis_q column is used.
        PFResult run(const TruthTable& tt, const std::map<int, double>& q_vis = {}) const
        {
            PFResult res;
            const auto& pfc = truth_pf_columns();
            const int c_tid = column_index(pfc, "trackid"), c_mom = column_index(pfc, "mother_trackid"),
                c_pdg = column_index(pfc, "pdg"), c_proc = column_index(pfc, "process"),
                c_ke = column_index(pfc, "KE"), c_nvis = column_index(pfc, "nvis_q"), c_nu = column_index(pfc, "nu_row"),
                c_sx = column_index(pfc, "start_x"), c_sy = column_index(pfc, "start_y"), c_sz = column_index(pfc, "start_z"),
                c_ex = column_index(pfc, "end_x"), c_ey = column_index(pfc, "end_y"), c_ez = column_index(pfc, "end_z"),
                c_st = column_index(pfc, "start_t");

            struct Row {
                int tid, mother, pdg, nu_row;
                std::string proc;
                double ke, len, q;
                double t0;          // start time (michel_stopped_only)
                Point start;
                // results
                int node{-1};       // own node index or -1
                int target{0};      // trackid whose node this id maps to (self when own node)
                bool primary{false};
            };
            std::map<int, Row> rows;    // trackid order
            for (const auto& r : tt.pf) {
                Row row;
                row.tid = (int)r[c_tid];
                row.mother = (int)r[c_mom];
                row.pdg = (int)r[c_pdg];
                row.nu_row = (int)r[c_nu];
                const int pcode = (int)r[c_proc];
                row.proc = (pcode >= 0 && pcode < (int)process_names().size()) ? process_names()[pcode] : "unknown";
                row.ke = r[c_ke] * units::GeV;
                const Point s(r[c_sx] * units::cm, r[c_sy] * units::cm, r[c_sz] * units::cm);
                const Point e(r[c_ex] * units::cm, r[c_ey] * units::cm, r[c_ez] * units::cm);
                row.len = (e - s).magnitude();
                row.start = s;
                row.t0 = r[c_st] * units::ns;
                auto qit = q_vis.find(row.tid);
                row.q = q_vis.empty() ? r[c_nvis] : (qit == q_vis.end() ? 0.0 : qit->second);
                row.target = row.tid;
                rows.emplace(row.tid, row);
            }
            auto parent_of = [&](const Row& r) -> const Row* {
                auto it = rows.find(r.mother);
                return it == rows.end() ? nullptr : &it->second;
            };

            // primaries (rule 2)
            for (auto& [tid, r] : rows) {
                const Row* p = parent_of(r);
                while (p && is_skipped_parent(p->pdg)) p = parent_of(*p);
                r.primary = (p == nullptr) || r.proc == "primary";
            }

            // classification and own-node decision, parents before
            // children (trackid order is creation order in Geant4; a
            // child with a smaller id than its parent is handled by the
            // recursive resolve below)
            std::function<void(Row&)> classify;
            std::set<int> done;
            auto make_node = [&](Row& r, Cls cls) {
                PFNode n;
                n.index = (int)res.nodes.size();
                n.trackid = r.tid;
                n.pdg = r.pdg;
                n.cls = cls;
                n.primary = r.primary;
                n.nu_row = r.nu_row;
                n.ke = r.ke;
                n.length = r.len;
                n.process = r.proc;
                res.nodes.push_back(n);
                r.node = n.index;
                r.target = r.tid;
            };
            classify = [&](Row& r) {
                if (done.count(r.tid)) return;
                done.insert(r.tid);
                Row* p = const_cast<Row*>(parent_of(r));
                if (p) classify(*p);
                const int ppdg = p ? p->pdg : 0;
                const bool parent_is_shower = p && p->node >= 0 && res.nodes[p->node].cls == Cls::shower;
                const bool parent_in_shower = p && (parent_is_shower ||
                                                    (p->node < 0 && p->target != p->tid &&
                                                     rows.at(p->target).node >= 0 &&
                                                     res.nodes[rows.at(p->target).node].cls == Cls::shower));
                // em_fold: a delta or Michel node (or an id merged into one) also takes EM daughters
                auto parent_in_cls = [&](Cls c) {
                    if (!p) return false;
                    if (p->node >= 0) return res.nodes[p->node].cls == c;
                    return p->target != p->tid && rows.at(p->target).node >= 0 &&
                        res.nodes[rows.at(p->target).node].cls == c;
                };
                const bool parent_in_em = parent_in_shower ||
                    (m_p.em_fold && (parent_in_cls(Cls::delta) || parent_in_cls(Cls::michel)));
                if (!is_visible(r.pdg)) {
                    // rule 1: map to the parent's target (or to nothing)
                    r.target = p ? p->target : r.tid;
                    ++res.nmerged_invisible;
                    return;
                }
                if (std::abs(r.pdg) == 11 && m_p.michel_own && std::abs(ppdg) == 13 &&
                    (r.proc == "Decay" || r.proc == "muMinusCaptureAtRest")) {
                    bool michel = true;
                    if (m_p.michel_capture_ke_min > 0 && r.proc == "muMinusCaptureAtRest" &&
                        r.ke < m_p.michel_capture_ke_min) {
                        michel = false;             // a capture's atomic electron
                    }
                    if (m_p.michel_stopped_only && r.proc == "Decay") {
                        if (ppdg == 13) {
                            michel = false;         // mu- Decay: in flight
                        }
                        else {
                            const double wait = r.t0 - p->t0 - p->len / units::clight;
                            if (wait < m_p.michel_wait_min) michel = false;   // mu+ decay in flight
                        }
                    }
                    if (michel) {
                        make_node(r, Cls::michel);  // rule 5
                        return;
                    }
                }
                if (std::abs(r.pdg) == 11 && is_ionisation(r.proc) && p) {
                    if (parent_in_em) {
                        r.target = p->target;        // rule 4: ionisation inside a shower
                        ++res.nmerged_shower;
                        return;
                    }
                    if (r.len < m_p.delta_len) {
                        r.target = p->target;        // rule 3: short delta into the track
                        ++res.nmerged_delta;
                        return;
                    }
                    make_node(r, Cls::delta);
                    return;
                }
                if (is_em(r.pdg) && p && parent_in_em && (is_shower_process(r.proc) || is_ionisation(r.proc) || r.proc == "Decay")) {
                    r.target = p->target;            // rule 4
                    ++res.nmerged_shower;
                    return;
                }
                if (is_em(r.pdg)) {
                    make_node(r, Cls::shower);       // shower root
                    return;
                }
                make_node(r, Cls::track);            // charged hadron, muon, neutron
            };
            for (auto& [tid, r] : rows) classify(r);

            // resolve targets to node roots, accumulate charge and ids
            auto root_of = [&](int tid) {
                int cur = tid;
                for (int guard = 0; guard < 100000; ++guard) {
                    const auto& r = rows.at(cur);
                    if (r.target == cur) return cur;
                    cur = r.target;
                }
                return cur;
            };
            for (auto& [tid, r] : rows) {
                r.target = root_of(tid);
            }
            auto accumulate = [&]() {
                for (auto& n : res.nodes) { n.q_vis = 0; n.ids.clear(); }
                for (auto& [tid, r] : rows) {
                    const auto& root = rows.at(r.target);
                    if (root.node < 0) continue;     // an invisible primary with no node
                    auto& n = res.nodes[root.node];
                    n.q_vis += r.q;
                    n.ids.push_back(tid);
                }
            };
            accumulate();

            // parent node of each node = node of the nearest ancestor with a node
            auto parent_node = [&](const PFNode& n) {
                const Row* p = parent_of(rows.at(n.trackid));
                while (p) {
                    const auto& root = rows.at(p->target);
                    if (root.node >= 0 && root.node != n.index) return root.node;
                    p = parent_of(*p);
                }
                return -1;
            };
            for (auto& n : res.nodes) n.parent = parent_node(n);

            // rule 7: survival, parents before children (node order is
            // creation order of the root ids after the recursion)
            for (auto& n : res.nodes) {
                if (n.parent < 0) { n.survived = true; continue; }
                const bool tracklike = n.cls == Cls::track || n.cls == Cls::delta;
                const bool strong = n.primary || (tracklike && n.length >= m_p.track_len_min) ||
                    n.cls == Cls::shower || n.cls == Cls::michel;
                n.survived = (n.q_vis >= m_p.q_min) && (n.ke >= m_p.pf_ke_min || n.cls == Cls::michel) && strong;
                if (!n.survived) {
                    int target = n.parent;
                    while (target >= 0 && !res.nodes[target].survived) target = res.nodes[target].parent;
                    if (target < 0) { n.survived = true; continue; }   // nothing to merge into
                    rows.at(n.trackid).target = res.nodes[target].trackid;
                    ++res.nmerged_weak;
                }
            }
            for (auto& [tid, r] : rows) r.target = root_of(tid);
            accumulate();
            for (auto& n : res.nodes) {
                n.parent = n.survived ? parent_node(n) : -1;
                if (n.survived && n.cls == Cls::shower && n.q_vis < m_p.compton_q) {
                    n.cls = Cls::le;                 // rule 6
                }
                n.instance_valid = n.survived && n.cls != Cls::le;
                std::sort(n.ids.begin(), n.ids.end());
            }

            // the g4 map
            for (const auto& [tid, r] : rows) {
                G4Info gi;
                gi.merged_into = r.target;
                const auto& root = rows.at(r.target);
                if (root.node >= 0) {
                    const auto& n = res.nodes[root.node];
                    gi.pf_node = n.index;
                    gi.cls = n.cls;
                    gi.instance_valid = n.instance_valid;
                }
                res.g4[tid] = gi;
            }

            // void roots: interactions without a surviving primary node at the vertex
            const auto& nuc = truth_nu_columns();
            const int n_idx = column_index(nuc, "nu_idx"), n_x = column_index(nuc, "vtx_x"),
                n_y = column_index(nuc, "vtx_y"), n_z = column_index(nuc, "vtx_z");
            for (size_t row = 0; row < tt.nu.size(); ++row) {
                if (tt.nu[row][n_idx] < 0) continue;
                const Point vtx(tt.nu[row][n_x] * units::cm, tt.nu[row][n_y] * units::cm, tt.nu[row][n_z] * units::cm);
                bool found = false;
                for (const auto& n : res.nodes) {
                    if (!n.survived || !n.primary || n.nu_row != (int)row || n.q_vis <= 0.0) continue;
                    if ((rows.at(n.trackid).start - vtx).magnitude() <= m_p.void_radius) { found = true; break; }
                }
                if (!found) res.void_roots.insert((int)row);
            }
            return res;
        }

        const PFMergerParams& params() const { return m_p; }

      private:
        PFMergerParams m_p;
    };
}

#endif
