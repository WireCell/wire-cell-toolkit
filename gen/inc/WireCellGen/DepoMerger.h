#ifndef WIRECELL_GEN_DEPOMERGER
#define WIRECELL_GEN_DEPOMERGER

#include "WireCellIface/IDepoMerger.h"
#include "WireCellIface/IConfigurable.h"

namespace WireCell {
    namespace Gen {

        class DepoMerger : public IDepoMerger, public IConfigurable {
           public:
            DepoMerger();
            virtual ~DepoMerger();

            // IDepoMerger
            virtual bool operator()(input_queues& inqs, output_queues& outqs);

            // IConfigurable
            virtual void configure(const WireCell::Configuration& config);
            virtual WireCell::Configuration default_configuration() const;

           private:
            // std::tuple<bool, bool> m_eos;
            int m_nin0, m_nin1, m_nout;
            bool m_eos;
        };
    }  // namespace Gen
}  // namespace WireCell

#endif
