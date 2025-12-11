/** A thread semaphore.

    At least until we get C++20.

    Use this to limit the number of active threads in a code context,
    eg in a service method.

*/
#ifndef WIRECELL_SEMAPHORE
#define WIRECELL_SEMAPHORE

#include <mutex>
#include <atomic>
#include <condition_variable>
#include <cassert>

namespace WireCell {

    // https://vorbrodt.blog/2019/02/05/fast-semaphore/
    //
    // The "fast" version below is about 10x faster than the basic
    // version if count is more than 1.  The base (or "fast" with
    // count=1) runs at about a MHz.

    class Semaphore
    {
      public:

        /// Create a semaphore with "count" initial spaces available
        /// to grab.  If set to zero, no clients will be able to
        /// acquire().
        Semaphore(int count=0) noexcept
            : m_count(count) { assert(count > -1); }

        /// Set the semaphore count.  It is only safe to call this
        /// prior to any acquire() or release().  It is safe to
        /// increase the count via calls to release() not balanced by
        /// calls to acquire().
        void set_count(int c)
        {
            m_count = c;
        }

        /// Grab the semaphore, blocks thread until space is available.
        void acquire() noexcept
        {
            std::unique_lock<std::mutex> lock(m_mutex);
            m_cv.wait(lock, [&]() { return m_count != 0; });
            --m_count;
        }

        /// Make one more space available on the semaphore.  This must
        /// be called by a client that previously called acquire() in
        /// oder to release the semephore.  It may be called in
        /// isolation to increase space on the semaphore.  
        void release() noexcept
        {
            {
                std::unique_lock<std::mutex> lock(m_mutex);
                ++m_count;
            }
            m_cv.notify_one();
        }

        /// Get the number of semaphores left.  Note, this is subject to race
        /// but can be useful for debugging.
        virtual int get_count() const { return m_count; }

      private:
        int m_count;
        std::mutex m_mutex;
        std::condition_variable m_cv;
    };

    // See comments on semaphore above

    class FastSemaphore
    {
      public:
        FastSemaphore(int count) noexcept
            : m_count(count), m_semaphore(0) {}
 
        void set_count(int c)
        {
            m_count = c;
        }

        void acquire()
        {
            int count = m_count.fetch_sub(1, std::memory_order_relaxed);
            if (count < 1) {
                m_semaphore.acquire();
            }
            std::atomic_thread_fence(std::memory_order_acquire);
        }

        void release()
        {
            std::atomic_thread_fence(std::memory_order_release);
            int count = m_count.fetch_add(1, std::memory_order_relaxed);
            if (count < 0) {
                m_semaphore.release();
            }
        }
 
        /// Get the number of semaphores left.  Note, this is subject to race
        /// but can be useful for debugging.
        virtual int get_count() const { return m_count; }

      private:
        std::atomic<int> m_count;
        Semaphore m_semaphore;
    };
}


#endif
