/** DepoBlobAssoc: the BlobDepoFill integration, keyed by G4 track.

    This is img/src/BlobDepoFill.cxx duplicated (fork-by-duplication:
    the production component stays byte-identical) with one change:
    the charge each depo deposits in each blob is accumulated per
    track key as well as in total, so that a blob knows which G4
    particles made it and with which fractions (wcp-porting-img
    pr_ml/docs/19 sec 3.5; doc 15 sec 2.1 and 3.1).

    The total per blob is accumulated in the same (slice, depo, wire)
    order from the same doubles as BlobDepoFill, so filled_cluster()
    is bit-identical to BlobDepoFill's output when the extra smearing
    is off.

    The track key of a depo is prior()->id() when id_from_prior is
    set and a prior exists (the undrifted G4 depo), else id().  The
    energy share is the depo's charge share times the key depo's
    energy(); it is only meaningful when the depo file carried the
    truth sidecar (DepoFileSink truth_sidecar), else energy() is 1.
 */

#ifndef WIRECELLTRUTH_DEPOBLOBASSOC
#define WIRECELLTRUTH_DEPOBLOBASSOC

#include "WireCellIface/ICluster.h"
#include "WireCellIface/IDepo.h"
#include "WireCellUtil/Logging.h"

#include <map>
#include <unordered_map>

namespace WireCell::Truth {

    struct DepoBlobAssocParams {
        /// As BlobDepoFill: added to depo times before slicing.
        double time_offset{0.0};
        /// As BlobDepoFill: Gaussian extent considered, in sigma.
        double nsigma{3.0};
        /// As BlobDepoFill: drift speed converting extent_long to time.  Required.
        double speed{0.0};
        /// As BlobDepoFill: the primary plane the depos are diced onto.
        int pindex{2};
        /// Track key and energy from prior() when it exists.
        bool id_from_prior{true};
        /// Added in quadrature to the longitudinal extent (in time)
        /// and the transverse extent (in distance).  0 = off, in
        /// which case the integration is exactly BlobDepoFill's.
        double extra_smear_long{0.0};
        double extra_smear_tran{0.0};
    };

    struct DepoBlobAssocResult {
        struct BlobEntry {
            /// Sum over all depos, accumulated as BlobDepoFill does.
            double q_total{0.0};
            /// track key -> (charge, energy share), sorted by key.
            std::map<int, std::pair<double, double>> by_track;
        };
        /// Keyed by the blob's vertex descriptor in the input cluster graph.
        std::unordered_map<cluster_vertex_t, BlobEntry> blobs;
        size_t nblobs{0};
        size_t ndepos{0};       // depos in a sensitive volume of some face
    };

    /// Run the association.  Every b-node of the cluster gets an entry
    /// (possibly empty).  Throws ValueError if speed <= 0.
    DepoBlobAssocResult depo_blob_assoc(const ICluster::pointer& cluster,
                                        const IDepo::vector& depos,
                                        const DepoBlobAssocParams& params,
                                        Log::logptr_t log = nullptr);

    /// A copy of the cluster with each blob's value replaced by its
    /// q_total and uncertainty 0 -- what BlobDepoFill emits.
    ICluster::pointer filled_cluster(const ICluster::pointer& cluster,
                                     const DepoBlobAssocResult& res);
}

#endif
