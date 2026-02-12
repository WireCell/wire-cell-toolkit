/** Sink tracking data to a ROOT file for MicroBooNE.
 *
 * This is a filter that takes grouping results and saves tracking information
 * to a ROOT file. It passes through its input, allowing it to sit in the
 * middle of a processing chain.
 */

#ifndef WIRECELLROOT_UBOONEMAGNIFYTRACKINGSINK
#define WIRECELLROOT_UBOONEMAGNIFYTRACKINGSINK

#include "WireCellIface/ITensorSetFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellUtil/Logging.h"

class TFile;
class TTree;

namespace WireCell {
    namespace Root {

        class UbooneMagnifyTrackingSink : public ITensorSetFilter, public IConfigurable {
           public:
            UbooneMagnifyTrackingSink();
            virtual ~UbooneMagnifyTrackingSink();

            /// ITensorSetFilter
            virtual bool operator()(const ITensorSet::pointer& in, ITensorSet::pointer& out);

            /// IConfigurable
            virtual WireCell::Configuration default_configuration() const;
            virtual void configure(const WireCell::Configuration& config);

           private:
            Log::logptr_t log;
            Configuration m_cfg;
            int m_runNo;
            int m_subRunNo;
            int m_eventNo;
            std::vector<IAnodePlane::pointer> m_anodes;
            IDetectorVolumes::pointer m_dv;

            void create_file();
            void write_bad_channels(TFile* output_tf, const ITensorSet::pointer& ts);
        };
    }  // namespace Root
}  // namespace WireCell

#endif