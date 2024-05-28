/** This class "splats" depos directly into a frame.
 *
 * It is approximately equivalent to combined simulation and sigproc
 * using the same response.
 *
 * See DepoSetSplat for a version that input depo sets.
 */

#ifndef WIRECELLGEN_DEPOSPLAT
#define WIRECELLGEN_DEPOSPLAT


#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IDuctor.h"

#include "WireCellIface/IAnodeFace.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IRandom.h"
#include "WireCellUtil/Logging.h"

#include <vector>

namespace WireCell {
    namespace Gen {
        class DepoSplat : public IDuctor, public IConfigurable {
          public:
            DepoSplat();
            virtual ~DepoSplat();

            // virtual void reset();
            virtual bool operator()(const input_pointer& depo, output_queue& frames);

            virtual void configure(const WireCell::Configuration& config);
            virtual WireCell::Configuration default_configuration() const;

          private:

            std::string m_anode_tn;

            IAnodePlane::pointer m_anode;

            IDepo::vector m_depos;

            double m_start_time;
            double m_readout_time;
            double m_tick;
            double m_drift_speed;
            double m_nsigma;
            std::string m_mode;

            // If set (config fluctuate=true), apply charge-preserving
            // fluctuation to each Gaussian sampling.  Default is no
            // fluctuation.
            IRandom::pointer m_rng;

            int m_frame_count;
            // if non-empty, set as tag on output frame
            std::string m_frame_tag{""};

            void process(output_queue& frames);
            ITrace::vector process_face(IAnodeFace::pointer face, const IDepo::vector& face_depos);
            bool start_processing(const input_pointer& depo);
            Log::logptr_t l;
        };
    }  // namespace Gen
}  // namespace WireCell



// #include "WireCellGen/Ductor.h"

// namespace WireCell {
//     namespace Gen {

//         // DepoSplat inherits from Ductor, replacing the heavy lifting
//         // with some lightweight laziness.
//         class DepoSplat : public Ductor {
//           public:
//             DepoSplat();
//             virtual ~DepoSplat();

//           protected:
//             virtual ITrace::vector process_face(IAnodeFace::pointer face, const IDepo::vector& depos);

//             /// SPD logger
//             Log::logptr_t l;
//         };
//     }  // namespace Gen
// }  // namespace WireCell

#endif
