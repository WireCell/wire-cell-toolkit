/** Implement a semaphore component interace. */

#ifndef WIRECELLAUX_SEMAPHORE
#define WIRECELLAUX_SEMAPHORE

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/ISemaphore.h"
#include "WireCellUtil/Semaphore.h"


namespace WireCell::Aux {
    class Semaphore : public ISemaphore,
                      public IConfigurable
    {
      public:
        Semaphore();
        virtual ~Semaphore();

        // IConfigurable interface
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;

        // ISemaphore
        virtual void acquire() const;
        virtual void release() const;
        virtual int get_count() const;

      private:

        mutable FastSemaphore m_sem;

    };
}  // namespace WireCell::Pytorch

#endif  // WIRECELLPYTORCH_TORCHSERVICE
