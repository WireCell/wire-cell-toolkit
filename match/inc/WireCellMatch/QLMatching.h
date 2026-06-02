#ifndef WIRECELL_MATCH_QLMATCHING
#define WIRECELL_MATCH_QLMATCHING

#include "WireCellAux/Logger.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellIface/IFiducial.h"
#include "WireCellIface/ITensorSetFilter.h"

#include "WireCellMatch/SemiAnalyticalModel.h"
#include "WireCellMatch/TimingTPCBundle.h"

#include "WireCellUtil/Units.h"

#include <memory>
#include <string>
#include <unordered_map>
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
    /// Aux::FlashTensorToOpticalPCs, the same schema as root/UbooneClusterSource).
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

        // Number of optical-detector channels (per-flash PE vector length).
        // The only bit-safe source of nchan now that flashes arrive via the
        // canonical PCs (do NOT derive from max channel ident).
        int m_nchan{312};

        // PMTs on vs off (default true => apply the OpDet type mask); see
        // ch_mask for further per-channel disables.
        bool m_pmts{true};
        // OpDet types kept by the matching mask, matched against
        // SemiAnalyticalModel::OpticalDetector::type (1 = dome PMT,
        // 0 = (X)Arapuca). Default {1} = PMTs only, reproducing the historical
        // SBND selection. The mask is derived per-channel from the injected
        // OpDet table rather than a hard-coded SBND layout, so it generalizes to
        // other detectors via config.
        std::vector<int> m_active_opdet_types{1};
        bool m_data{true};
        std::vector<int> m_ch_mask;
        bool m_beamonly{false};
        double m_flash_minPE{500};
        double m_flash_mintime{-1.5 * units::ms};
        double m_flash_maxtime{1.5 * units::ms};
        double m_beam_mintime{-5 * units::us};
        double m_beam_maxtime{5 * units::us};
        double m_QtoL{0.5};
        // Drift speed used for the per-flash X correction (flash_x_offset =
        // sign * flash_time * drift_speed). In WCT internal units (length/time,
        // i.e. mm/ns numerically). Default is the historical hard-coded value;
        // configs should pass the common params.lar.drift_speed via jsonnet.
        double m_drift_speed{1.563 * units::mm / units::us};
        // LASSO solution threshold below which a (flash, cluster) bundle is
        // dropped after each matching round. Was hard-coded 0.05 inline; pulled
        // out so it can be widened/narrowed from the jsonnet without rebuild.
        double m_strength_cutoff{0.05};

        // ---- Tuning constants pulled out of the inline code (see
        // match/docs/improve_progress.md). The §B-§G defaults equal the historical
        // hard-coded literals (bit-identical if omitted); the §A cushions below
        // default to the MicroBooNE convention and intentionally differ. ----

        // §A active-volume cushions. The raw active-volume bounds (anode/cathode
        // X, Y span, Z span) now come from the IDetectorVolumes service
        // (m_dv->inner_bounds) so the code is detector-agnostic. These signed
        // cushions adjust the effective PE-inclusion and boundary-flag windows
        // relative to the true geometry, in the per-TPC anode->cathode drift
        // coordinate u (u=0 at the anode/PMT plane, u=u_cathode at the cathode
        // seam), so a single set serves both reversed-drift SBND APAs.
        //
        // Defaults follow the MicroBooNE prototype convention
        // (ToyMatching.cxx::calculate_pred_pe ~158-164): outward-positive in u.
        // They shift the inclusion window slightly vs the old SBND literals
        // (intended); overriding them from jsonnet can reproduce the old bounds
        // bit-identically. See match/docs/improve_progress.md §A.
        double m_anode_ext1{-2.0 * units::cm};    // PE-inclusion window edge, below the anode (low_x_cut_ext1)
        double m_anode_ext2{4.0 * units::cm};     // anode flag-window outer edge (low_x_cut_ext2)
        double m_cathode_ext1{1.2 * units::cm};   // PE-inclusion window edge, beyond the cathode (high_x_cut_ext1)
        double m_cathode_ext2{-2.0 * units::cm};  // cathode flag-window inner edge (high_x_cut_ext2)
        double m_y_cushion{0.0 * units::cm};       // signed inward(+)/outward(-) shift of each |y| edge
        double m_z_cushion{0.0 * units::cm};       // signed inward(+)/outward(-) shift of each z edge

        // §D pre-selection / bad-match gates.
        double m_mc_saturation_pe{5000};      // MC saturated-PMT mask trigger (total flash PE)

        // §E out-of-beam QA cuts.
        double m_outbeam_ks_max{0.2};
        double m_outbeam_chi2ndf_max{20};
        double m_outbeam_pe_frac{0.5};        // out-of-beam total-PE mismatch fraction

        // §C LASSO weights.
        double m_lasso_lambda{0.1};
        double m_delta_charge{0.01};
        double m_delta_light{0.025};
        double m_delta_shape{0.01};
        double m_bkg_weight{0.5};              // background-column weight (round 1)
        double m_pe_mismatch_knee{0.3};        // PE-mismatch weight knee (fraction of meas)
        double m_pe_mismatch_floor{0.3};       // PE-mismatch weight floor

        // §G flash PE-error model (forwarded to Opflash).
        double m_pe_err_floor{0.3};
        double m_pe_err_frac{0.3};
        double m_pe_err_knee{1.0};
        double m_flash_pe_threshold{0.0};      // Opflash "fired" threshold (PE)

        // §F bundle-quality thresholds (forwarded to TimingTPCBundle).
        double m_bundle_ks_merge_max{0.2};
        double m_bundle_chi2ndf_merge_max{20};
        double m_bundle_addmerge_exponent{0.8};
        double m_highconsist_ks_max{0.06};
        int    m_highconsist_min_ndf{3};
        double m_bundle_pe_ndf_knee{1.0};

        // §H raw readout-window truncation flag (T0-independent, APA-agnostic),
        // always computed. Flags a bundle whose cluster's leading/trailing time
        // slice sits within m_window_edge_ticks of the raw readout window
        // [0, m_readout_window_ticks]. Inert (no consumer yet); since nothing
        // reads it, always filling it leaves production output unchanged.
        int m_readout_window_ticks{3427};  // window end in raw ticks (SBND daq.nticks)
        int m_window_edge_ticks{4};        // edge-proximity threshold (~one slice = nticks_live_slice)

        // When true, discard any (flash, cluster) bundle whose cluster is not
        // contained in the TPC drift box once the flash-derived T0 x-offset is
        // applied (the prototype flag_good_bundle gate, ToyMatching.cxx 272-275).
        // Default OFF so existing production configs stay bit-identical; SBND
        // jsonnet sets require_containment: true.
        bool m_require_containment{false};

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

        // Optional SBND CPA structure-exclusion fiducial volume. When set (the
        // "cathode_fiducial" config tn), the cathode-end flag_at_x_boundary is
        // decided by a 3-D point-in-volume test against it; when null (default,
        // non-SBND) the cathode-end flag falls back to the original flat-cathode
        // 1-D drift-coordinate window, leaving those configs bit-identical.
        IFiducial::pointer m_cathode_fv{nullptr};

        // Per-cluster cache of the cluster's significant extreme points
        // (get_extreme_wcps, flattened). Used by the cathode-end fiducial test; the
        // geometric extremes are intrinsic to the cluster (offset-independent for
        // ranking) so they are computed once and reused across candidate flashes.
        // Cleared each event in operator().
        mutable std::unordered_map<const WireCell::Clus::Facade::Cluster*,
                                   std::vector<WireCell::Point>> m_extreme_cache;

        std::string m_inpath{"pointtrees/%d"};
        std::string m_outpath{"pointtrees/%d"};
        float m_cluster_t0{-1e12};

        // Path to the JSON file holding VUVHits, VISHits, geometry and the
        // SBND OpDet array.
        std::string m_semimodel_file{"sbnd/photodet/semi-analytical-sbnd.json"};

        std::unique_ptr<SemiAnalyticalModel> m_semi_model;

        // Per-OpDet table (type + position) loaded from semimodel_file. Used to
        // build the SemiAnalyticalModel and, at execute() time, to derive the
        // per-channel on/off mask and the per-TPC OpDet split.
        std::vector<SemiAnalyticalModel::OpticalDetector> m_opdets;
        double m_cathode_x{0.0};  // cathode-plane x; OpDets on its low/high side belong to TPC 0/1

        // Fill the per-bundle boundary flags (flag_close_to_PMT,
        // flag_at_x_boundary, flag_spec_end) for one (flash, cluster) bundle.
        // Ports the MicroBooNE prototype end-trimming walk (ToyMatching.cxx
        // ~176-290) into the per-TPC drift coordinate u. s/anode_x/u_cathode are
        // the geometry scalars computed once per anode in operator().
        //
        // Returns the prototype's flag_good_bundle / TPC-containment verdict: true
        // iff the (post-trim) cluster endpoints fall inside the TPC drift box (the
        // 4-part in-window guard, == ToyMatching.cxx 272-275). Returns false when
        // the cluster has no 3d points under this anode (it cannot be contained).
        // The caller discards uncontained bundles when m_require_containment.
        bool compute_endpoint_flags(TimingTPCBundle* bundle,
                                    WireCell::Clus::Facade::Cluster* cluster,
                                    double flash_x_offset,
                                    double s, double anode_x, double u_cathode) const;

        // Significant extreme points of a cluster (get_extreme_wcps, flattened),
        // memoized in m_extreme_cache. Used to locate the cathode endpoint.
        const std::vector<WireCell::Point>&
            cluster_extreme_points(WireCell::Clus::Facade::Cluster* cluster) const;

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
