/**
   A reframer takes makes a "rectangular" frame filled with samples
   from the tagged traces of its input.  This new frame has exactly
   one trace for each channel and each trace is padded to span a
   uniform duration.  These configuration paramters control how the
   new frame is shaped:

   - tags :: list of trace tags to consider.  an empty list means to
     use all traces.  Default: [].

   - tbin :: the number of ticks measured from the *input* frame
     reference time to start each output trace.  Converted to time via
     the input frame's tick, this value will be applied to the output
     frame's reference time.  As such, the "tbin()" method of the
     output traces will return 0.  Default: 0.

   - nticks :: the total number of ticks in the output traces.
     Default: 0.

   - toffset :: an addition time offset arbitrarily added to the
     output frame's reference time.  This time is *asserts some
     fiction* and does not contribute to calculating the output tbin.
     Default: 0.0.

   - fill :: if a needed sample does not exist it will be set to this
     value.  Default: 0.0.

 */

#ifndef WIRECELL_GEN_REFRAMER
#define WIRECELL_GEN_REFRAMER

#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellAux/Logger.h"

#include <vector>
#include <string>

namespace WireCell {
    namespace Gen {

        class Reframer : public Aux::Logger,
                         public IFrameFilter, public IConfigurable {
           public:
            Reframer();
            virtual ~Reframer();

            virtual void configure(const WireCell::Configuration& config);
            virtual WireCell::Configuration default_configuration() const;

            virtual bool operator()(const input_pointer& inframe, output_pointer& outframe);

           private:
            /// ASSUMPTION: input traces and summary are in the same order
            std::pair<ITrace::vector, IFrame::trace_summary_t> process_one(const ITrace::vector& itraces, const IFrame::trace_summary_t& isummary);
            ITrace::vector process_one(const ITrace::vector& itraces); /// no summary version
            IAnodePlane::pointer m_anode;

            /// Configure: tags
            ///
            /// Will reframe each set of tagged traces separately, carrying
            /// forward the tag to the output.
            ///
            /// If empty, then all traces in the frame are reframed regardless
            /// of any trace tagging.
            std::vector<std::string> m_input_tags;

            /// Configure: ignore_tags
            ///
            /// If input_tags is empty and yet tagged traces are found in the
            /// frame it signifies the user likely made an error and a warning
            /// is printed.  Set this to true to supress the warning.  A true
            /// value is inconsistent with a non-empty "tags" array.
            bool m_ignore_tags{false};

            /// Configure: frame_tag
            ///
            /// If set, apply as a frame tag on the output frame.
            std::string m_frame_tag{""};

            /// Configure: toffset
            ///
            /// The output frame reference time is nominally set to the the
            /// input frame reference time and extended by "tbin" worth of
            /// "tick" in the output traces.  The "toffset" time will be
            /// arbitrarily added to that nominal time.
            double m_toffset{0};

            /// Configure: fill
            ///
            /// The value to assign to output trace samples that are not
            /// represented in the input frame.
            double m_fill{0};

            /// Configure: tbin
            ///
            /// The starting time bin (tick) of the traces.  This effectively
            /// truncates the output frame array.
            int m_tbin{0};

            /// Configure: nticks
            ///
            /// The size of the output frame array along the time dimension.
            int m_nticks{0};

            // count calls for more useful log messages.
            size_t m_count{0};
        };
    }  // namespace Gen
}  // namespace WireCell
#endif
