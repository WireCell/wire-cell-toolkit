#ifndef WIRECELLGEN_DUCTOR
#define WIRECELLGEN_DUCTOR

#include "WireCellUtil/Pimpos.h"
#include "WireCellUtil/Response.h"

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IDuctor.h"

#include "WireCellIface/IAnodeFace.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IPlaneImpactResponse.h"
#include "WireCellIface/IRandom.h"
#include "WireCellUtil/Logging.h"

#include "WireCellGen/Cfg/Ductor/Structs.hpp"

#include <vector>

namespace WireCell {
    namespace Gen {

        /** This IDuctor needs a Garfield2D field calculation data
         * file in compressed JSON format as produced by Python module
         * wirecell.sigproc.garfield.
         */
        class Ductor : public IDuctor, public IConfigurable {
           public:
            Ductor();
            virtual ~Ductor(){};

            // virtual void reset();
            virtual bool operator()(const input_pointer& depo, output_queue& frames);

            virtual void configure(const WireCell::Configuration& config);
            virtual WireCell::Configuration default_configuration() const;

           protected:

            // Note, at least DepoSplat inherits from Ductor and needs
            // access to the config.
            using config_t = WireCellGen::Cfg::Ductor::Config;
            config_t m_cfg{};
            
            
            IAnodePlane::pointer m_anode{nullptr};
            IRandom::pointer m_rng{nullptr};
            std::vector<IPlaneImpactResponse::pointer> m_pirs;
            IDepo::vector m_depos;
            std::string m_mode{"continuous"};

            // These two are initialized from config but we update
            // them so hold them outside of m_cfg.
            size_t m_frame_count{0};
            double m_start_time{0.0};

            virtual void process(output_queue& frames);
            virtual ITrace::vector process_face(IAnodeFace::pointer face, const IDepo::vector& face_depos);
            bool start_processing(const input_pointer& depo);
            Log::logptr_t l;
        };
    }  // namespace Gen
}  // namespace WireCell

#endif
