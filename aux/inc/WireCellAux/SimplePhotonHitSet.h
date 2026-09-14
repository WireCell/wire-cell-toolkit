#ifndef WIRECELLAUX_SIMPLEPHOTONHITSET
#define WIRECELLAUX_SIMPLEPHOTONHITSET

#include "WireCellIface/IPhotonHitSet.h"

namespace WireCell::Aux {

    class SimplePhotonHitSet : public WireCell::IPhotonHitSet {
        int m_ident;
        WireCell::IPhotonHit::shared_vector m_hits;
        WireCell::Configuration m_metadata;

       public:
        SimplePhotonHitSet(int ident, const WireCell::IPhotonHit::vector& hits,
                           const WireCell::Configuration& metadata = WireCell::Configuration())
          : m_ident(ident)
          , m_hits(std::make_shared<WireCell::IPhotonHit::vector>(hits.begin(), hits.end()))
          , m_metadata(metadata)
        {
        }
        virtual ~SimplePhotonHitSet();
        virtual int ident() const { return m_ident; }
        virtual WireCell::IPhotonHit::shared_vector hits() const { return m_hits; }
        virtual WireCell::Configuration metadata() const { return m_metadata; }
    };

}  // namespace WireCell::Aux

#endif
