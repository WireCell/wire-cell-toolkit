/** TruthPointLabeler: add the per-point "truth" point cloud to a live
    point tree.

    Join of (ITensorSet from PointTreeBuilding, ICluster, IDepoSet) ->
    ITensorSet.  The output set carries every input tensor not under
    the live tree's datapath unchanged, the live tree re-serialised
    with the added local point cloud on every blob node (its "3d" and
    "scalar" clouds untouched), the blobtruth / blobtruth_sum tensors
    of the association and the truth_nu / truth_pf tensors of the
    truth table, so that a downstream MultiAlgBlobClustering with
    aux_datatypes [truth_nu, truth_pf, blobtruth] forwards them.

    Configuration (every default is the legacy-neutral choice; the node
    is itself opt-in):

      the DepoBlobAssoc knobs of BlobDepoLabel: time_offset, nsigma,
        speed (required), pindex, id_from_prior, extra_smear_long,
        extra_smear_tran, extra_smear_tran_collection / _induction
      real_assoc      {extra_smear_tran_collection, extra_smear_tran_induction}: a
                      second association; only blobs the primary leaves empty
                      take its labels (default absent = off)
      inpath          datapath pattern of the live tree ("pointtrees/%d/live")
      anode, face     the AnodePlane and face of the tree (as PointTreeBuilding)
      detector_volumes  for the time_offset / drift_speed / tick metadata
                      (as PointTreeBuilding); when empty the knobs
                      reco_time_offset, drift_speed and tick are used
      truth_table     file pattern of the truth file ("%d" = cluster
                      ident); empty = build the table from the depos
      truth_prefix    datapath prefix of the truth tensors ("truth/%d")
      blobtruth_path  datapath prefix of blobtruth ("clusters/%d")
      id_mode         "auto" (sed_ids present => index), "track", "index"
      pfmerger        PFMergerParams keys (pf_ke_min, q_min, delta_len, ...)
      refine_by_nearest_depo, d_cap, end_window, kink_window, kink_angle,
        heat_sigma, shared_frac   PointLabelParams
      labels_meta     file to write the labels_meta JSON to ("" = none);
                      the same JSON is in the output set metadata
      label_dead      reserved (dead nodes are not labelled yet; true warns)
      polyline_drift_volume_only  build the trajectories from the depos of the
                      face's drift volume only (x between its anode and
                      cathode planes, any y / z; not behind the response
                      plane); false = every depo of the event, so a track
                      crossing the cathode gets a mirrored copy (pr_ml doc
                      23 sec 7.4)
      polyline_per_particle  one trajectory per Geant4 particle of a key
                      (truth table sed_ids orig_trackid in index mode, else
                      the depo id before folding); false = one time-ordered
                      polyline per key, a zigzag for folded keys (doc 23)

    wcp-porting-img pr_ml/docs/19 sec 3.6 (T6).
 */

#ifndef WIRECELLTRUTH_TRUTHPOINTLABELER
#define WIRECELLTRUTH_TRUTHPOINTLABELER

#include "WireCellTruth/DepoBlobAssoc.h"
#include "WireCellTruth/PFMerger.h"
#include "WireCellTruth/TruthPointLabel.h"
#include "WireCellAux/Logger.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellIface/ITruthPointLabeler.h"

#include <string>

namespace WireCell::Truth {

    class TruthPointLabeler : public Aux::Logger,
                              public ITruthPointLabeler,
                              public IConfigurable
    {
      public:
        TruthPointLabeler();
        virtual ~TruthPointLabeler();

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        virtual bool operator()(const input_tuple_type& intup,
                                output_pointer& out);

      private:
        Configuration labels_meta(const TruthTable& tt, const PFResult& pf, const PointLabelSummary& sum) const;

        DepoBlobAssocParams m_assoc;
        /// real_assoc (default off): a second association with extra
        /// transverse widths (extra_smear_tran_collection / _induction,
        /// the per-plane SP wire-filter widths; every other setting as
        /// the primary's).  A blob the primary leaves without charge but
        /// the second run fills takes the second run's entry; every other
        /// blob, the particle flow (q_vis) and the truth table stay the
        /// primary's (pr_ml doc 25 sec 7, owner 2026-10-25: the hybrid of
        /// the deghoster's W9 labels for the point labels).
        bool m_real_assoc{false};
        double m_real_tran_collection{0.0};
        double m_real_tran_induction{0.0};
        PFMergerParams m_pfparams;
        PointLabelParams m_label;
        std::string m_inpath{"pointtrees/%d/live"};
        std::string m_truth_table{""};
        std::string m_truth_prefix{"truth/%d"};
        std::string m_blobtruth_path{"clusters/%d"};
        std::string m_id_mode{"auto"};
        std::string m_labels_meta{""};
        bool m_label_dead{false};
        bool m_polyline_drift_volume_only{false};
        bool m_polyline_per_particle{false};
        IAnodePlane::pointer m_anode;
        int m_face{0};
        IDetectorVolumes::pointer m_dv;
        double m_reco_time_offset{0.0};
        double m_drift_speed{0.0};
        size_t m_count{0};
    };
}

#endif
