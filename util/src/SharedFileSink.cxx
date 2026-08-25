#include "WireCellUtil/SharedFileSink.h"

#include <map>

using namespace WireCell;

/// One open handle and one mutex per file path.  Kept for the life of the
/// process, which is how the per-logger file sinks behaved anyway -- nothing
/// destroys them -- and it means a sink made late in a job appends to the same
/// buffer as one made at configure time instead of opening a second.
struct Log::SharedFileSink::SharedFile {
    std::mutex mtx;
    spdlog::details::file_helper fh;
    spdlog::filename_t name;
};

std::shared_ptr<Log::SharedFileSink::SharedFile>
Log::SharedFileSink::acquire(const spdlog::filename_t& filename)
{
    // Guards the table.  Sink construction happens while components are being
    // configured, which is not necessarily single-threaded.
    static std::mutex mtx;
    static std::map<spdlog::filename_t, std::shared_ptr<SharedFile>> table;

    std::lock_guard<std::mutex> lk(mtx);
    auto it = table.find(filename);
    if (it != table.end()) {
        return it->second;
    }
    auto sf = std::make_shared<SharedFile>();
    sf->name = filename;
    // truncate=false, matching basic_file_sink's default: a sink opens the
    // file for append, it does not clobber it.
    sf->fh.open(filename, false);
    table[filename] = sf;
    return sf;
}

Log::SharedFileSink::SharedFileSink(const spdlog::filename_t& filename)
    : m_file(acquire(filename))
{
}

const spdlog::filename_t& Log::SharedFileSink::filename() const
{
    return m_file->name;
}

void Log::SharedFileSink::sink_it_(const spdlog::details::log_msg& msg)
{
    // Format with THIS sink's own formatter -- that is the whole point of the
    // per-logger sink objects, and it keeps the emitted line byte-identical to
    // what basic_file_sink produced.  base_sink<std::mutex> already holds this
    // instance's lock here, which is what makes the formatter's mutable cache
    // safe; the shared lock below covers only the write.
    spdlog::memory_buf_t formatted;
    formatter_->format(msg, formatted);

    std::lock_guard<std::mutex> lk(m_file->mtx);
    m_file->fh.write(formatted);
}

void Log::SharedFileSink::flush_()
{
    std::lock_guard<std::mutex> lk(m_file->mtx);
    m_file->fh.flush();
}
