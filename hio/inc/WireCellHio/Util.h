/** A wrapper for h5cpp calls
 */

#ifndef WIRECELLHIO_UTIL
#define WIRECELLHIO_UTIL

#include <mutex>
#include <thread>
#include <chrono>

// wow, h5cpp is kind of crappy
#pragma GCC diagnostic push
#pragma GCC diagnostic warning "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic warning "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic warning "-Warray-bounds"
#pragma GCC diagnostic ignored "-Warray-bounds"
#pragma GCC diagnostic warning "-Wreturn-type"
#pragma GCC diagnostic ignored "-Wreturn-type"
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic warning "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic warning "-Wreorder"
#pragma GCC diagnostic ignored "-Wreorder"
#pragma GCC diagnostic warning "-Wswitch"
#pragma GCC diagnostic ignored "-Wswitch"
#pragma GCC diagnostic warning "-Wpedantic"
#pragma GCC diagnostic ignored "-Wpedantic"

#include <h5cpp/all>

#pragma GCC diagnostic pop

namespace WireCell {
    namespace Hio {

        extern std::mutex g_h5cpp_mutex;

    };  // namespace Hio
}  // namespace WireCell

#endif  // WIRECELLHIO_UTIL
