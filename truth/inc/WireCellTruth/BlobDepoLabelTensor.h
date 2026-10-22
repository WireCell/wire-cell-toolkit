/**
   BlobDepoLabelTensor: (ICluster, IDepoSet) -> ITensorSet holding the
   per-(blob, G4 track) charge table of DepoBlobAssoc.

   Output tensors (all with md "datapath" under <datapath>):
   - optionally the filled cluster as TensorDM "cluster" tensors
     (as_tensors(filled_cluster, datapath));
   - "blobtruth": double (M, 6) with columns
     blob_ident, face, slice_ident, trackid, q_bk, e_bk
     one row per (blob, track) with q_bk > 0 (and q_bk/q_total >= min_frac);
   - "blobtruth_sum": double (B, 5) with columns
     blob_ident, q_true, q_top, tid_top, ntracks
     one row per blob of the cluster, in graph vertex order.
   Column names are in Truth::Schema and in each tensor's md["columns"].

   wcp-porting-img pr_ml/docs/19 sec 3.5.
 */

#ifndef WIRECELLTRUTH_BLOBDEPOLABELTENSOR
#define WIRECELLTRUTH_BLOBDEPOLABELTENSOR

#include "WireCellTruth/DepoBlobAssoc.h"
#include "WireCellAux/Logger.h"
#include "WireCellIface/IClusterDepoTensor.h"
#include "WireCellIface/IConfigurable.h"

namespace WireCell::Truth {

    class BlobDepoLabelTensor : public Aux::Logger,
                                public IClusterDepoTensor,
                                public IConfigurable
    {
      public:
        BlobDepoLabelTensor();
        virtual ~BlobDepoLabelTensor();

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        virtual bool operator()(const input_tuple_type& intup,
                                output_pointer& out);

      private:
        DepoBlobAssocParams m_params;
        /// Datapath pattern, "%d" replaced by the cluster ident.
        std::string m_datapath{"clusters/%d"};
        /// Include the filled cluster's TensorDM tensors.
        bool m_with_cluster{true};
        /// Drop (blob, track) rows below this fraction of the blob's total.
        double m_min_frac{0.0};
        size_t m_count{0};
    };
}

#endif
