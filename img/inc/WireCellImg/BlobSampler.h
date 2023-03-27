/** Produce point cloud representation of blobs through application of
    various sampling strategies.
 
    This produces a dataset with the following array names:
 
    - <prefix>x :: point X-position in the global coordinate system
    - <prefix>y :: point Y-position in the global coordinate system
    - <prefix>z :: point Z-position in the global coordinate system
    - <prefix>strategy :: an integer count for the sampling strategy
    - <prefix>blobid :: the blob ident

    The prefix "<prefix>" may be changed by setting the "prefix"
    configuration parameter.  Default is empty string.
 */

#ifndef WIRECELL_IMG_BLOBSAMPLER
#define WIRECELL_IMG_BLOBSAMPLER

#include "WireCellIface/IBlobSampler.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"
#include "WireCellUtil/Binning.h"

namespace WireCell::Img {

    class BlobSampler : public Aux::Logger, public IBlobSampler, public IConfigurable
    {
      public:

        BlobSampler();
        virtual ~BlobSampler() = default;

        // IBlobSampler
        virtual PointCloud::Dataset sample_blobs(const IBlob::vector& blobs);
        
        // IConfigurable
        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;
        
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
                {0,1}, which determine the ray grid.

              Note, ray crossing points may be sampled with two
              instances of "grid" each with step=1.0 and with mutually
              unique pairs of plane indices.

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

        */
        struct Sampler;
        std::vector<std::unique_ptr<Sampler>> m_samplers;


        void add_strategy(Configuration cfg);

    };

};

#endif

