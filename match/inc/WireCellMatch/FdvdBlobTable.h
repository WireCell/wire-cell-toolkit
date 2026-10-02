/** FdvdBlobTable: the imaging blobs of all CRMs of an FD-VD event as one table.
 *
 * Fan-in of the per-CRM imaging ICluster streams (the same inputs the
 * clustering's PointTreeBuilding samples), emitting one ITensorSet with a
 * "blobs" tensor (f8 [nblob, NCOL]).  Each row is one blob with the light
 * prediction's blob centre computed exactly as the python reference reads it
 * from the cluster archive (fdvd_sim/stageB/img_eval_all.py blobs, lines
 * 45-67): y, z = mean of the blob corners as stored in the archive (the
 * ClusterArrays bnodes corner columns; needs ClusterFileSource
 * restore_corners: true, else the reloaded shape's corners are used), x from
 * the slice middle,
 *   x = (xw_mm - v_mm_per_us * ((start + span/2)/1e3 - toff_us)) / 10  [cm],
 * charge = IBlob::value().  The PC-tree blob scalars carry no corners
 * (PointTreeBuilding samples only "3d"), which is why the table is built from
 * the imaging stream.
 *
 * Rows follow the python blob order: inputs sorted by their "labels" string
 * (the python globs the archive names "clusters-apa-crm<N>-<pipe>" sorted),
 * blobs in graph vertex order within an input.  Columns:
 *   0 order  1 anode  2 face  3 slice_index_min  4 slice_index_max
 *   5 u_min 6 u_max 7 v_min 8 v_max 9 w_min 10 w_max   (strip bounds, layers 2-4)
 *   11 x_cm 12 y_cm 13 z_cm 14 charge 15 start_ns 16 span_ns  17 input port
 * Columns 1-10 are the join key with the PC-tree blob "scalar" PC
 * (Aux::fill_scalar_blob: slice_index = int(start / tick)).
 */
#ifndef WIRECELL_MATCH_FDVDBLOBTABLE
#define WIRECELL_MATCH_FDVDBLOBTABLE

#include "WireCellIface/IClusterFaninTensorSet.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"

#include <string>
#include <vector>

namespace WireCell::Match {

    class FdvdBlobTable : public Aux::Logger, public IClusterFaninTensorSet, public IConfigurable {
      public:
        static constexpr int NCOL = 18;

        FdvdBlobTable();
        virtual ~FdvdBlobTable();

        virtual std::vector<std::string> input_types();
        virtual bool operator()(const input_vector& invec, output_pointer& out);

        virtual WireCell::Configuration default_configuration() const;
        virtual void configure(const WireCell::Configuration& cfg);

      private:
        size_t m_multiplicity{1};
        std::vector<std::string> m_labels;   // one per input; rows ordered by sorted label
        double m_tick{500.0};                // ns, the slice_index unit of Aux::fill_scalar_blob
        double m_xw_mm{3250.7};              // img_pilot_eval.py line 29: XW, V, TOFF
        double m_v_mm_per_us{1.60563};
        double m_toff_us{8.7};
        // Use the imaging-time corners carried by ClusterFileSource(restore_corners: true) when present;
        // else the reloaded shape's re-derived corners.
        bool m_stored_corners{true};
        size_t m_nstored{0};
        size_t m_count{0};
    };

}  // namespace WireCell::Match

#endif
