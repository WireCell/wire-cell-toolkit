/** Produce point cloud representation of a set of blobs.

    See blob-sampling.org document for details.
 */

#ifndef WIRECELL_CLUS_BLOBSAMPLER
#define WIRECELL_CLUS_BLOBSAMPLER

#include "WireCellIface/IBlobSampler.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"
#include "WireCellUtil/Binning.h"

#include <regex>

namespace WireCell::Clus {

    class BlobSampler : public Aux::Logger, public IBlobSampler, public IConfigurable
    {
      public:

        BlobSampler();
        virtual ~BlobSampler() {};

        // IBlobSampler
        // virtual PointCloud::Dataset sample_blobs(const IBlob::vector& blobs);
        virtual std::tuple<PointCloud::Dataset, PointCloud::Dataset> sample_blob(const IBlob::pointer& blob,
                                                int bob_index = 0);
        
        // IConfigurable
        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        // Runtime configuration override interface
        virtual std::tuple<PointCloud::Dataset, PointCloud::Dataset> 
        sample_blob_with_config(const IBlob::pointer& blob,
                              int blob_index = 0,
                              const Configuration& runtime_config = Configuration());
        
        struct CommonConfig {

            /** Config: "time_offset".  A time value ADDED to the blob
                times prior to processing.  This may be required to remove
                any artificial "start time" from simulation or otherwise
                properly convert from drift time to drift location.
            */
            double time_offset{0*units::us};

            /** Config: "drift_speed".  The drift speed to assume in order
                to translate a blob time to a spacial position in the
                drift direction.
            */
            double drift_speed{1.6*units::mm/units::us};

            /** Config: "prefix".  Sets the string prepended to names of
             * the arrays in the produced dataset. */
            std::string prefix{""};

            /** Config: "tbins", "tmin" and "tmax" */
            Binning tbinning = Binning(1,0.0,1.0);

            /** Config: "extra".  Match potential extra arrays to
             * include in point cloud. */
            std::vector<std::string> extra = {};
            std::vector<std::regex> extra_re = {};

        };
        CommonConfig m_cc;

        /** Config "strategy".

            A full sampling is performed by applying a pipeline of
            functions each enacting a particular sampling strategy.
            The full result is the union of results from each applied
            strategy.

            The "strategy" configuration parameter value may be
            provided in one of these forms:

            - a string providing the strategy name and indicating the
              default configuration for that strategy be used.

            - an object providing a .name attribute naming the
              strategy.  Any other object attributes will be
              considered for providing per-strategy options.

            - an array of "strategy" configuration paramter values.

            The following strategy names are recognized.  They are
            distinquished in how points are sampled in the transverse
            plane.  They all have common options to determine how
            sampled points are projected along the time/drift
            dimension.
        
            - center :: the single point which is the average of the
              blob corner points.  This is the default.

            - corner :: the blob corner points.
            
            - edge :: the center of blob boundary edges.

            - grid :: uniformly spaced points on a ray grid.

              This accepts the following options:

              - step :: grid spacing given as relative to a grid step.
                (def=1.0).

              - planes :: an array of two plane indices, each in
                {0,1,2}, which determine the ray grid.

              Note, ray crossing points may be sampled with two
              instances of "grid" each with step=1.0 and with mutually
              unique pairs of plane indices.

            - stepped :: sample points on a stepped ray grid of two views.

              This accepts the following options:

              - min_step_size :: The minimium number of wires over
                which a step will be made.  default=3.

              - max_step_fraction :: The maximum fraction of a blob a
                step may take.  If non-positive, then all steps are
                min_step_size.  default=1/12.

              - offset :: How far along the diagonal from a crossing point to
                the crossing point of the next neighbor rays to place the point.
                A value of 0 places the point at the ray crossing.  A value of
                0.5 (default) places the point at the crossing of the ray's
                wires and makes this sampling equivalent to what Wire-Cell
                prototype uses.

            - bounds :: points sampled along the blob boundry edges.

              This accepts the following options:

              - step :: grid spacing given as relative to a grid step.
                (def=1.0).

              Note, first point on each boundary edge is one step away
              from the corner.  If "step" is larger than a boundary
              edge that edge will have no samples from this strategy.
              Include "corner" and/or "edge" strategy to include these
              two classes of points missed by "bounds".

            The points from any of the strategies above are projected
            along the time/drift dimension of the blob according to a
            linear grid given by a binning: (tbins, tmin, tmax).  The
            tmin/tmax are measured relative to the time span of the
            blob.  The default binning is (1, 0.0, 1.0) so that a
            single set of samples are made at the start of the blob.

            Each of the strategies may have a different time binning
            by setting:

            - tbins :: the number of time bins (def=1)

            - tmin :: the minimum relative time (def=0.0)

            - tmax :: the maximum relative time (def=1.0)

            These options may be set on the BlobSampler to change the
            default values for all strategies and may be set on each
            individual strategy.



              - charge_stepped :: sample points on a stepped ray grid with charge-based filtering.

              This is an enhanced version of the "stepped" strategy that adds
              charge-based filtering to improve sampling quality by considering
              wire signal strength. This strategy is based on the WCPPID sampling
              method from the Wire-Cell prototype. Bad/dead planes are automatically
              detected based on charge uncertainty values.

              This accepts the following options:

              - min_step_size :: The minimum number of wires over which a step 
                will be made. default=3.

              - max_step_fraction :: The maximum fraction of a blob a step may 
                take. If non-positive, then all steps are min_step_size. default=1/12.

              - offset :: How far along the diagonal from a crossing point to
                the crossing point of the next neighbor rays to place the point.
                A value of 0 places the point at the ray crossing. A value of
                0.5 (default) places the point at the crossing of the ray's
                wires. default=0.5.

              - tolerance :: Tolerance for pitch bounds checking. default=0.03.

              - charge_threshold_max :: Minimum charge threshold for wires in the
                plane with maximum wire coverage. default=4000.

              - charge_threshold_min :: Minimum charge threshold for wires in the
                plane with minimum wire coverage. default=4000.

              - charge_threshold_other :: Minimum charge threshold for wires in the
                third plane. default=4000.

              - max_wire_product_threshold :: When the product of max_wires × min_wires
                is less than or equal to this value, all wires are used instead of
                stepped sets. default=2500.

              - disable_mix_dead_cell :: Boolean flag controlling how zero-charge
                wires are handled. When true, zero-charge wires are treated as
                failing charge thresholds. default=false.

              - dead_threshold :: Charge uncertainty threshold for detecting dead/bad
                planes. If more than 50% of sampled channels in a plane have uncertainty
                above this threshold, the plane is considered bad and its charge
                threshold is set to 0. default=1e10.

              The charge-based filtering logic works as follows:
              1. Bad planes are automatically detected by analyzing charge uncertainty
              2. For bad planes, charge thresholds are automatically set to 0
              3. Wires from the stepped sets (mandatory wires) bypass some charge filtering
              4. Non-mandatory wires must meet charge thresholds
              5. If both crossed wires are mandatory, no additional charge filtering is applied
              6. If either crossed wire is non-mandatory, all three wire charges are checked
              7. Points where all three wires have zero charge are excluded

              Note: This strategy requires charge data with uncertainty information to be
              available in the blob's activity map through the slice interface.
        */
        struct Sampler;
        std::vector<std::unique_ptr<Sampler>> m_samplers;

        void add_strategy(Configuration cfg);
    };

}

#endif

