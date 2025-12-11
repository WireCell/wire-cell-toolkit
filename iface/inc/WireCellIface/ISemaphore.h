/** An interface to the semaphore pattern */

#ifndef WIRECELL_ISEMAPHORE
#define WIRECELL_ISEMAPHORE

#include "WireCellUtil/IComponent.h"

namespace WireCell {
    class ISemaphore : public IComponent<ISemaphore> {
      public:
        virtual ~ISemaphore();

        /// Block until available spot to hold the semaphore is
        /// available.
        virtual void acquire() const = 0;

        /// Release hold on the semaphore
        virtual void release() const = 0;

        /// Use Construct a Context on a semaphore in a local scope to
        /// automate release
        struct Context {
            ISemaphore::pointer sem;
            Context(ISemaphore::pointer sem) : sem(sem) { sem->acquire(); }
            ~Context() { sem->release(); }
        };

        /// Useful only for debugging.
        virtual int get_count() const { return -1; }

    };
}  // namespace WireCell

#endif  // WIRECELL_ITENSORFORWARD
