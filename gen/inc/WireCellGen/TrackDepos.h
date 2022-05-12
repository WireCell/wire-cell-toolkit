#ifndef WIRECELL_TRACKDEPOS
#define WIRECELL_TRACKDEPOS

#include "WireCellAux/Logger.h"
#include "WireCellAux/Configurable.h"

#include "WireCellIface/IDepoSource.h"

#include "WireCellUtil/Units.h"

#include "WireCellGen/Cfg/TrackDepos.hpp"

#include <tuple>
#include <deque>

namespace WireCell {

    namespace Gen {

        using WireCell::Gen::Cfg::TrackDepos::Config;

        /// A producer of depositions created from some number of simple, linear tracks.
        class TrackDepos : public Aux::Logger,
                           public Aux::Configurable<Config>, 
                           public IDepoSource {
           public:
            /// Create tracks with depositions every stepsize and assumed
            /// to be traveling at clight.
            TrackDepos(double stepsize = 1.0 * units::millimeter, double clight = 1.0);
            virtual ~TrackDepos();

            // virtual void configure(const WireCell::Configuration& config);
            // virtual WireCell::Configuration default_configuration() const;


            /// ISourceNode
            virtual bool operator()(IDepo::pointer& out);

          public:
            // Expose for use by unit tests
            WireCell::IDepo::vector depos();
            void add_track(double time, WireCell::Ray ray, double dedx = -1.0);
            using track_t = std::tuple<double, Ray, double>;
            std::vector<track_t> tracks() const { return m_tracks; }

          private:

            virtual void configured();

            // using config_t = WireCell::Gen::Cfg::TrackDepos::Config;
            // config_t m_cfg;

            int m_count;
            // Log::logptr_t l;

            std::vector<track_t> m_tracks; // keep to enable testing
            std::deque<WireCell::IDepo::pointer> m_depos;

        };

    }  // namespace Gen
}  // namespace WireCell

#endif
