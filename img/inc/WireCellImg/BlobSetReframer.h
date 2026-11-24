/* Project the blobs in a blobset to a frame.
 *
 * The projection in each view of a blob will produce a trace with
 * total charge equal to the blob's charge value (or uncertainty).
 *
 * Time: ISlice::start() is relative to the time of the slice's frame.
 * The output frame takes its time from the slices' frame and a trace
 * starts at tbin = round(slice start / tick).  This is correct when the
 * blobs come straight from a slicer in the same job.  Clusters read
 * from file do not yet carry the original frame time (their slice
 * frames have time 0), so the output frame time is then 0.  See
 * docs/components/BlobSetReframer.md for the problem and suggested
 * remedies.
 */

#ifndef WIRECELLIMG_BLOBSETREFRAMER
#define WIRECELLIMG_BLOBSETREFRAMER

#include "WireCellAux/Logger.h"

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IBlobSetFramer.h"

#include "WireCellUtil/Units.h"

namespace WireCell::Img {

    class BlobSetReframer : public Aux::Logger, public IBlobSetFramer, public IConfigurable {
    public:
        BlobSetReframer();
        virtual ~BlobSetReframer();

        virtual void configure(const WireCell::Configuration& cfg);
        virtual WireCell::Configuration default_configuration() const;

        virtual bool operator()(const input_pointer& in, output_pointer& out);

    private:


        /// Config: "frame_tag" - if nonempty, apply this tag to
        /// output frame.  Note, the frame may span more than one
        /// anode.  Trace tags are applied with names like
        /// "anode<ident>" where "<ident>" gives the anode ident.
        std::string m_frame_tag{""};

        /// Config: "tick" - the sampling time for the output frame.
        /// The charge of a blob over a blob's time span will be
        /// resampled to "tick" so as to preserve current.  Eg, one
        /// unit of charge "q" per blob time slice of 2us will result
        /// in trace charge elements of 0.25q when resampled to a tick
        /// of 0.5us.
        double m_tick{0.5*units::us};

        /// Config: "source" - string to determine where to get "charge".
        ///
        /// May be one of: "blob" or "activity".  Default if unset is "blob".
        ///
        /// Config: "measure" - string to determine what kind of value to use as "charge".
        ///
        /// May be one of: "value" or "error".  Default if unset is "value".
        using trace_maker_f = std::function<ITrace::vector(const IBlob::pointer& iblob)>;
        trace_maker_f m_trace_maker;

        size_t m_count{0};
    };
}  // namespace WireCell::Img

#endif
