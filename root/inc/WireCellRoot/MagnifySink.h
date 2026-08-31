/** Sink data to a file format used by the "Magnify" GUI .
 *
 * This is technically a "filter" as it passes on its input.  This
 * allows an instance of the sink to sit in the middle of some longer
 * chain.
 *
 * FIXME: currently this class TOTALLY violates the encapsulation of
 * DFP by requiring the input file in order to transfer input data out
 * of band of the flow.
 */

#ifndef WIRECELLROOT_MAGNIFYFILESINK
#define WIRECELLROOT_MAGNIFYFILESINK

#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellUtil/Logging.h"

class TFile;

namespace WireCell {
    namespace Root {

        class MagnifySink : public IFrameFilter, public IConfigurable {
           public:
            MagnifySink();
            virtual ~MagnifySink();

            /// IFrameSink
            virtual bool operator()(const IFrame::pointer& in, IFrame::pointer& out);

            /// IConfigurable
            virtual WireCell::Configuration default_configuration() const;
            virtual void configure(const WireCell::Configuration& config);

           private:
            Configuration m_cfg;
            IAnodePlane::pointer m_anode;

            int m_nrebin;
            bool m_geo_written{false};
            void create_file();
            void do_shunt(TFile* output_tf);
            void write_geo_tree(TFile* output_tf);

            Log::logptr_t log;
        };
    }  // namespace Root
}  // namespace WireCell

#endif
