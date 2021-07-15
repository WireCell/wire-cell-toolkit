#ifndef WIRECELL_TRACKDEPOS
#define WIRECELL_TRACKDEPOS

#include "WireCellIface/IDepoSource.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Logging.h"

#include "WireCellGen/Cfg/TrackDepos/Structs.hpp"

#include <tuple>
#include <deque>

namespace WireCell {

    namespace Gen {

        /// A producer of depositions created from some number of simple, linear tracks.
        class TrackDepos : public IDepoSource, public IConfigurable {
           public:
            /// Create tracks with depositions every stepsize and assumed
            /// to be traveling at clight.
            TrackDepos(double stepsize = 1.0 * units::millimeter, double clight = 1.0);
            virtual ~TrackDepos();

            virtual void configure(const WireCell::Configuration& config);
            virtual WireCell::Configuration default_configuration() const;


            /// ISourceNode
            virtual bool operator()(IDepo::pointer& out);

          public:
            // Expose for use by unit tests
            WireCell::IDepo::vector depos();
            void add_track(double time, WireCell::Ray ray, double dedx = -1.0);
            using track_t = std::tuple<double, Ray, double>;
            std::vector<track_t> tracks() const { return m_tracks; }

          private:

            void apply_grouping();

            using config_t = WireCellGen::Cfg::TrackDepos::Config;
            config_t m_cfg;

            int m_count;
            Log::logptr_t l;

            std::vector<track_t> m_tracks; // keep to enable testing
            std::deque<WireCell::IDepo::pointer> m_depos;

        };

    }  // namespace Gen
}  // namespace WireCell

#endif
