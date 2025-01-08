//  This adapts an IScaler to the IDepoSetScaler interface.
//
// This is a very dumb implementation as a per depo scaler is doing
// extra work to keep its output in time order where as we could do
// better by ignoring order during drifting and do a final sort().
//
// Further improvements would drift depos as a block and utilize SIMD.
// See gen-kokkos for smarter smarts.
//
// The only practical reason to use this is it will speed up Pgrapher
// (less so, TbbFlow) compared to using a bare per depo scaler.

#ifndef WIRECELLGEN_DEPOSETSCALER
#define WIRECELLGEN_DEPOSETSCALER

#include "WireCellIface/IDrifter.h"
#include "WireCellIface/IDepoSetFilter.h"
#include "WireCellIface/INamed.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"

namespace WireCell::Gen {

  class DepoSetScaler : public Aux::Logger,
    public IDepoSetFilter, public IConfigurable {     
  public:

    DepoSetScaler();
    virtual ~DepoSetScaler();

    // IDepoSetFilter
    virtual bool operator()(const input_pointer& in, output_pointer& out);

    /// WireCell::IConfigurable interface.
    virtual void configure(const WireCell::Configuration& config);
    virtual WireCell::Configuration default_configuration() const;

  private:

    IDrifter::pointer m_scaler{nullptr};
    size_t m_count{0};

  };

}



#endif
