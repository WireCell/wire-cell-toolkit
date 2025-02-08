/**
   Class and functions to help work with uboone TTrees, TC, TDC, Trun, T_light, T_match1.
 */

#ifndef WIRECELL_ROOT_UBOONETTREES
#define WIRECELL_ROOT_UBOONETTREES

#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Waveform.h"
#include "WireCellUtil/Logging.h"

#include <memory>
#include <vector>
#include <string>

#include "TFile.h"
#include "TTree.h"

namespace WireCell::Root {

    /**
       Sole interface to data in ROOT TTrees.
    */
    class UbooneTTrees {

        std::unique_ptr<TFile> m_tfile;
        TTree* m_activity{nullptr}; // always load
        TTree* m_live{nullptr};     // always load
        TTree* m_dead{nullptr};     // load only if "dead" in kinds
        TTree* m_bad{nullptr};      // optional
        TTree* m_flash{nullptr};    // load only if "light" in kinds
        TTree* m_match{nullptr};    // load only if "light" in kinds

        Long64_t m_nentries{0}; // number of entries in the trees.
        Long64_t m_entry{-1};   // not yet loaded.  Only use next() to advance

    public:

        // The data is available after construction and refreshed on each call
        // to next().  The data is invalid when/if next() throws IndexError and
        // "this" instance should be dropped shortly after.

        struct Header {
            // double eventTime{0};    // 1.45767e9. Unix time?
            double triggerTime{0};  // 1.48822e10. Units?
            int runNo{0};
            int eventNo{0};
            int nrebin{0};
            // float unit_dis{0};

            void set_addresses(TTree& tree) {
                tree.SetBranchAddress("triggerTime", &triggerTime);
                tree.SetBranchAddress("runNo", &runNo);
                tree.SetBranchAddress("eventNo", &eventNo);
                tree.SetBranchAddress("nrebin", &nrebin);
                // tree.SetBranchAddress("unit_dis",&unit_dis);
            }
        };
        Header header;

        struct Activity {

            // "activity"
            std::vector<int> *timesliceId{nullptr}; // [0,2399]
            std::vector<std::vector<int>> *timesliceChannel{nullptr}; // [0,8255]
            // Units are number of signal electrons.
            std::vector<std::vector<int>> *raw_charge{nullptr};
            std::vector<std::vector<int>> *raw_charge_err{nullptr};

            void set_addresses(TTree& tree) {
                tree.SetBranchAddress("timesliceId", &timesliceId);
                tree.SetBranchAddress("timesliceChannel", &timesliceChannel);
                tree.SetBranchAddress("raw_charge", &raw_charge);
                tree.SetBranchAddress("raw_charge_err", &raw_charge_err);
            }
        };
        Activity activity;

        // TC and TDC are nearly the same.
        struct Blob {
            std::vector<int> *cluster_id_vec{nullptr};

            // 0 value indicates the view for a blob is "dead"
            std::vector<int> *flag_u_vec{nullptr};
            std::vector<int> *flag_v_vec{nullptr};
            std::vector<int> *flag_w_vec{nullptr};

            std::vector<std::vector<int>> *wire_index_u_vec{nullptr};
            std::vector<std::vector<int>> *wire_index_v_vec{nullptr};
            std::vector<std::vector<int>> *wire_index_w_vec{nullptr};

            // Lord, forgive us our sins.
            std::vector<std::vector<int>*> flag_uvw;

            void set_addresses(TTree& tree) {
                tree.SetBranchAddress("cluster_id", &cluster_id_vec); 
                tree.SetBranchAddress("flag_u", &flag_u_vec);
                tree.SetBranchAddress("flag_v", &flag_v_vec);
                tree.SetBranchAddress("flag_w", &flag_w_vec);
                tree.SetBranchAddress("wire_index_u", &wire_index_u_vec);
                tree.SetBranchAddress("wire_index_v", &wire_index_v_vec);
                tree.SetBranchAddress("wire_index_w", &wire_index_w_vec);

                // derived
                flag_uvw = {flag_u_vec, flag_v_vec, flag_w_vec};
            }
        };

        struct TC : Blob { 
            std::vector<int> *time_slice_vec{nullptr};               // TC: singular
            std::vector<double> *q_vec{nullptr};

            void set_addresses(TTree& tree) {
                Blob::set_addresses(tree);
                tree.SetBranchAddress("q", &q_vec);
                tree.SetBranchAddress("time_slice", &time_slice_vec);
            }
        };
        TC live;

        struct TDC : Blob {
            std::vector<std::vector<int>> *time_slices_vec{nullptr}; // TDC: plural

            void set_addresses(TTree& tree) {
                Blob::set_addresses(tree);
                tree.SetBranchAddress("time_slice", &time_slices_vec);
            };
        };
        TDC dead;

        // Bad channels may override an activity map entry or delete it.
        class Bad {
            int chid{0}, plane{0};
            int tick_beg{0}, tick_end{0}; // one past?

            Waveform::ChannelMasks _masks;

        public:
            void set_addresses(TTree& tree) {
                tree.SetBranchAddress("chid", &chid);
                tree.SetBranchAddress("plane", &plane);
                tree.SetBranchAddress("start_time", &tick_beg);
                tree.SetBranchAddress("end_time", &tick_end);

                // derived
                _masks.clear();  // this is a one-shot fill
                const int nentries = tree.GetEntries();
                for (int entry = 0; entry<nentries; ++entry) {
                    tree.GetEntry(entry);
                    auto& brl = _masks[chid];
                    brl.emplace_back(tick_beg, tick_end);
                }
            }

            const Waveform::ChannelMasks& masks() const { return _masks; }

        };
        Bad bad;

        struct Flash {
            int type, flash_id;
            double time, tmin, tmax, qtot;
            double light[32], dlight[32];
            std::vector<int>* channels = nullptr;
            
            void set_addresses(TTree& tree) {
                tree.SetBranchAddress("time",&time);
                tree.SetBranchAddress("type",&type);
                tree.SetBranchAddress("flash_id",&flash_id);
                tree.SetBranchAddress("low_time",&tmin);
                tree.SetBranchAddress("high_time",&tmax);
                tree.SetBranchAddress("total_PE",&qtot);
                tree.SetBranchAddress("PE",light);
                tree.SetBranchAddress("PE_err",dlight);
                tree.SetBranchAddress("fired_channels",&channels);
            }
        };
        Flash flash;

        struct Match {
            int cluster_id, flash_id;

            void set_addresses(TTree& tree) {
                tree.SetBranchAddress("tpc_cluster_id", &cluster_id);
                tree.SetBranchAddress("flash_id", &flash_id);
            }            
        };
        Match match;

    public:

        // Construct the interface to the trees.
        //
        // This does NOT load any entry.  Must call next() to load first entry, etc.
        UbooneTTrees(const std::string& tfile_name, const std::vector<std::string>& kinds) {
            
            m_tfile.reset(TFile::Open(tfile_name.c_str(), "READ"));
            if (!m_tfile) {
                raise<IOError>("failed to open %s", tfile_name);
            }

            m_activity = open_ttree("Trun");
            m_nentries = m_activity->GetEntries();
            header.set_addresses(*m_activity);
            activity.set_addresses(*m_activity);

            m_live = open_ttree("TC");
            live.set_addresses(*m_live);

            // Needed for live or dead 2-view but optional
            m_bad = open_ttree("T_bad_ch", false);
            if (m_bad) {
                bad.set_addresses(*m_bad);
            }

            auto has = [&](const std::string& what) {
                for (const auto& k : kinds) {
                    if (what == k) return true;
                }
                return false;
            };

            if (has("dead")) {
                // only needed for kind=="dead"
                m_dead = open_ttree("TDC");
                dead.set_addresses(*m_dead);

            }
            if (has("light")) {
                m_flash = open_ttree("T_flash");
                flash.set_addresses(*m_flash);
                m_match = open_ttree("T_match");
                match.set_addresses(*m_match);
            }
        }

        void next() {
            ++m_entry;
            if (m_entry < m_nentries) { 
                m_activity->GetEntry(m_entry);
                m_live->GetEntry(m_entry);
                if (m_dead) {
                    m_dead->GetEntry(m_entry);
                }
                return;
            }
            raise<IndexError>("attempt to get entry %d past end of TTree with %d",
                              m_entry, m_nentries);
        }

        // Get live or dead blob data
        const Blob& blob(const std::string& which) {
            if (which == "live") return live;
            return dead;
        }

        int nentries() const { return m_nentries; }
        int entry() const { return m_entry; }
        // Number of slices represented in the data.
        int nslices_data() const { return activity.timesliceId->size(); }
        // Number of slices spanned from slice ID=0 to slice ID = max
        int nslices_span() const { return 1 + *std::max_element(activity.timesliceId->begin(),
                                                                activity.timesliceId->end()); }

    private:

        TTree* open_ttree(const std::string& tree_name, bool required = true) {
            auto ttree = reinterpret_cast<TTree*>(m_tfile->Get(tree_name.c_str()));
            if (!ttree && required) {
                raise<IOError>("failed to get TTree %s from %s", tree_name, m_tfile->GetName());
            }
            return ttree;
        }
    };

    // A class to iterate on TTree entries in an ordered list of files.
    class UbooneTFiles {
        std::vector<std::string> m_input;
        std::vector<std::string> m_kinds;
        WireCell::Log::logptr_t log;
        int m_calls{-1};
    public:
        UbooneTFiles(const std::vector<std::string>& fnames,
                     const std::vector<std::string>& kinds,
                    WireCell::Log::logptr_t log)
            : m_input(fnames)
            , m_kinds(kinds)
            , log(log) {

            // so we can use pop_back.
            std::reverse(m_input.begin(), m_input.end());
        }

        std::unique_ptr<UbooneTTrees> trees;

        bool next() {
            ++m_calls;
            if (!trees) {
                // We are starting a new file.

                if (m_input.empty()) {
                    // We have exhausted our input files.
                    return false;
                }

                auto fname = m_input.back();
                m_input.pop_back();

                try {
                    trees = std::make_unique<UbooneTTrees>(fname, m_kinds);
                }
                catch (IOError& err) {
                    log->warn("failed to open {}, skipping, call={}", fname, m_calls);
                    trees = nullptr;
                    return next();
                }
                catch (IndexError& err) {
                    trees = nullptr;
                    return next();
                }

                // sanity check file
                const int nentries = trees->nentries();
                if (! nentries) {
                    log->warn("no entries {}, skipping, call={}", fname, m_calls);
                    trees = nullptr;
                    return next();
                }
                log->debug("starting {} with {} entries, call={}", fname, nentries, m_calls);
            }

            // We have an open file with at least some entries, try to advance
            try {
                trees->next();
            }
            catch (IndexError& err) {
                trees = nullptr;
                return next();
            }

            // sanity check entry
            const int n_slices_data = trees->nslices_data();
            if (!n_slices_data) {
                log->warn("no slices in entry {}, skipping, call={}", trees->entry(), m_calls);
                return next();
            }

            log->debug("read {} slices in entry {}, call={}", n_slices_data, trees->entry(), m_calls);
            return true;
        }
    };
}

#endif 
