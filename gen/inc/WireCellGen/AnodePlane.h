/** An AnodePlane provides wire-related and some volumetric
 * geometrical information as well serving as the top of
 * anode/face/plane organizational hierarchy for accessing this info.
 * See also IAnodeFace and IWirePlane.
 */

#ifndef WIRECELLGEN_ANODEPLANE
#define WIRECELLGEN_ANODEPLANE

#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Logging.h"
#include <unordered_map>

#include "WireCellGen/Cfg/AnodePlane/Structs.hpp"

namespace WireCell {
    namespace Gen {

        class AnodePlane : public IAnodePlane, public IConfigurable {
        public:
            AnodePlane();
            virtual ~AnodePlane() {}

            // IConfigurable interface
            virtual void configure(const WireCell::Configuration& config);
            virtual WireCell::Configuration default_configuration() const;

            /// IAnodePlane interface
            virtual int ident() const { return m_cfg.ident; }
            virtual int nfaces() const { return m_faces.size(); }
            virtual IAnodeFace::pointer face(int ident) const;
            virtual IAnodeFace::vector faces() const { return m_faces; }
            virtual WirePlaneId resolve(int channel) const;
            virtual std::vector<int> channels() const;
            virtual IChannel::pointer channel(int chident) const;
            virtual IWire::vector wires(int channel) const;

        private:
            using config_t = WireCellGen::Cfg::AnodePlane::Config;
            config_t m_cfg;

            IAnodeFace::vector m_faces;

            std::unordered_map<int, int> m_c2wpid;
            std::unordered_map<int, IWire::vector> m_c2wires;
            std::vector<int> m_channels;
            std::unordered_map<int, IChannel::pointer> m_ichannels;
            Log::logptr_t l;
        };
    }  // namespace Gen

}  // namespace WireCell

#endif
