#ifndef WIRECELL_MATCH_QLMATCHING
#define WIRECELL_MATCH_QLMATCHING

#include "WireCellAux/Logger.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellIface/ITensorSetFilter.h"

#include "WireCellMatch/SemiAnalyticalModel.h"
#include "WireCellMatch/TimingTPCBundle.h"

#include <memory>
#include <string>
#include <vector>

namespace WireCell::Match {

    /// Charge-light matcher.
    ///
    /// Port of larwirecell/qlmatch/QLMatching to wire-cell-toolkit, without
    /// any larsoft/fhicl/art dependencies. The SBND semi-analytical optical
    /// model is loaded from a JSON file at configure time (see
    /// cfg/sbnd/semi-analytical-sbnd.json or equivalent).
    ///
    /// Input (ITensorSetFilter): pctree tensor set holding clusters, with the
    /// per-event optical flashes attached as the canonical "flash"/"light"/
    /// "flashlight" point clouds on the live root node (see
    /// Match::OpflashToFlashPCs, the same schema as root/UbooneClusterSource).
    /// Output: cluster tensor set with per-cluster t0 AND a per-cluster
    /// "flash" scalar (matched flash row index) set from the matched flash, so
    /// Clus::Facade::Cluster::get_flash() reflects the match.
    class QLMatching : public Aux::Logger,
                       public ITensorSetFilter,
                       public IConfigurable {
    public:
        QLMatching();
        virtual ~QLMatching();

        bool operator()(const input_pointer& in, output_pointer& out) override;

        void configure(const WireCell::Configuration& cfg) override;
        WireCell::Configuration default_configuration() const override;

    private:
        std::size_t m_count{0};
        int m_bee_index{0};

        // Number of optical-detector channels (per-flash PE vector length).
        // The only bit-safe source of nchan now that flashes arrive via the
        // canonical PCs (do NOT derive from max channel ident).
        int m_nchan{312};

        // PMTs on vs off (default true => use SBND PMTs); see ch_mask for the
        // disabled channels.
        bool m_pmts{true};
        bool m_data{true};
        std::vector<int> m_ch_mask;
        bool m_beamonly{false};
        double m_flash_minPE{500};
        double m_flash_mintime{-1500e3};
        double m_flash_maxtime{1500e3};
        double m_beam_mintime{-5e3};
        double m_beam_maxtime{5e3};
        double m_QtoL{0.5};
        // Drift speed used for the per-flash X correction (flash_x_offset =
        // sign * flash_time * drift_speed). In WCT internal units (length/time,
        // i.e. mm/ns numerically). Default is the historical hard-coded value;
        // configs should pass the common params.lar.drift_speed via jsonnet.
        double m_drift_speed{1.563e-3};
        // LASSO solution threshold below which a (flash, cluster) bundle is
        // dropped after each matching round. Was hard-coded 0.05 inline; pulled
        // out so it can be widened/narrowed from the jsonnet without rebuild.
        double m_strength_cutoff{0.05};

        // Default SBND VUV/VIS efficiency arrays, indexed by OpDet (312 entries).
        // Configuration "VUVEfficiency"/"VISEfficiency" arrays override these.
        std::vector<double> m_VUVEfficiency{
            0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.03570, 0.03570, 0.03570, 0.03570, 0.03570, 0.03570, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.03570, 0.03570, 0.03570, 0.03570, 0.03570, 0.03570, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.03570, 0.03570, 0.03570, 0.03570, 0.03570, 0.03570, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.03570, 0.03570, 0.03570, 0.03570, 0.03570, 0.03570, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.03920, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752, 0.01752
        };
        std::vector<double> m_VISEfficiency{
            0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.03570, 0.03570, 0.03570, 0.03570, 0.03570, 0.03570, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.03570, 0.03570, 0.03570, 0.03570, 0.03570, 0.03570, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.03570, 0.03570, 0.03570, 0.03570, 0.03570, 0.03570, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.03570, 0.03570, 0.03570, 0.03570, 0.03570, 0.03570, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.01264, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.02600, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271, 0.00271
        };

        IAnodePlane::pointer m_anode{nullptr};
        IDetectorVolumes::pointer m_dv;

        std::string m_inpath{"pointtrees/%d"};
        std::string m_outpath{"pointtrees/%d"};
        std::string m_bee_dir{"data"};
        float m_cluster_t0{-1e12};

        // Path to the JSON file holding VUVHits, VISHits, geometry and the
        // SBND OpDet array.
        std::string m_semimodel_file{"sbnd/photodet/semi-analytical-sbnd.json"};

        std::unique_ptr<SemiAnalyticalModel> m_semi_model;

        void remove_bundle_selection(TimingTPCBundleSelection to_be_removed,
                                     TimingTPCBundleSet& bundle_set);
        void remove_bundle_selection(
            TimingTPCBundleSelection to_be_removed,
            FlashBundlesMap& flash_bundles_map,
            ClusterBundlesMap& cluster_bundles_map,
            std::map<std::pair<Opflash*, WireCell::Clus::Facade::Cluster*>,
                     TimingTPCBundle::pointer>& flash_cluster_bundles_map);
        void organize_bundles(
            TimingTPCBundleSelection& results_bundles,
            std::map<std::pair<Opflash*, WireCell::Clus::Facade::Cluster*>,
                     TimingTPCBundle::pointer>& flash_cluster_bundles_map);
    };

} // namespace WireCell::Match

#endif
