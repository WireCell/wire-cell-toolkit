#ifndef WIRECELL_SHAREDFILESINK
#define WIRECELL_SHAREDFILESINK

#include "WireCellUtil/Spdlog.h"

#include "spdlog/sinks/base_sink.h"
#include "spdlog/details/file_helper.h"

#include <memory>
#include <mutex>
#include <string>

namespace WireCell {

    namespace Log {

        /** A file sink that shares ONE open file handle across every sink
            instance naming the same file.

            Why this exists.  Log::logger(name, share_sinks=false) mints a fresh
            sink object per logger so that each can carry its own pattern, and
            Aux::Logger::set_name() uses exactly that for every configured
            component.  With spdlog's basic_file_sink each of those objects
            opens its OWN FILE*, so a job with N components has N independent
            4096-byte stdio buffers appending to one file.  They fill and drain
            at different times, so one component's buffer is flushed into the
            middle of another component's line -- a log message cut mid-word
            with an unrelated (and often many seconds older) message spliced
            over it.  Measured on SBND production logs: every splice sat on an
            exact 4096-byte file offset, and the splice count tracked the
            component count (98 components -> 70 splices in 5240 lines).

            This sink keeps the per-instance formatter -- so the emitted line
            format is byte-for-byte what it was -- but routes every write
            through a single shared handle under a single shared mutex.  One
            buffer means writes can only ever land whole and in order.

            The handle is opened once per path and kept for the life of the
            process, which is how the per-logger sinks behaved anyway (they are
            never destroyed).  Append semantics match basic_file_sink's default.
         */
        class SharedFileSink
            : public spdlog::sinks::base_sink<std::mutex>
        {
          public:
            explicit SharedFileSink(const spdlog::filename_t& filename);
            ~SharedFileSink() = default;

            /// The file this sink writes to.
            const spdlog::filename_t& filename() const;

          protected:
            // Two locks, and both are needed.  base_sink<std::mutex> takes a
            // PER-INSTANCE lock around sink_it_(), exactly as basic_file_sink_mt
            // did -- spdlog's pattern_formatter caches cached_tm_/last_log_secs_
            // and mutates them inside format(), so dropping that lock would race
            // two threads sharing one sink.  The shared per-FILE mutex is then
            // taken inside it, around the write only.  The ordering is always
            // instance-then-shared and the shared one is never held while taking
            // an instance lock, so there is no deadlock.
            void sink_it_(const spdlog::details::log_msg& msg) override;
            void flush_() override;

          private:
            struct SharedFile;
            // The one-per-path table lives inside this member so the nested
            // type stays private to the sink.
            static std::shared_ptr<SharedFile> acquire(const spdlog::filename_t& filename);
            std::shared_ptr<SharedFile> m_file;
        };

    }  // namespace Log

}  // namespace WireCell

#endif
