#include <iostream>
#include <memory>
#include <thread>
#include <chrono>
#include <mutex>

#include "WireCellUtil/String.h"
#include "WireCellUtil/Eigen.h"

// h5cpp NEED to be placed after Eigen to use h5::read<Eigen::Matrix>
#include <h5cpp/all>

std::mutex g_h5cpp_mutex;

template <class T>
using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

void hello(int& threadid)
{
    std::cout << "Hello World! Thread ID, " << threadid << std::endl;
    return;
}

namespace {
    h5::fd_t h5open(const int& threadid, const std::string& fname)
    {
        // comment the lock off to check thread safety
        std::lock_guard<std::mutex> guard(g_h5cpp_mutex);

        h5::fd_t fd;
        try {
            fd = h5::open(fname, H5F_ACC_RDONLY);  // H5F_ACC_RDONLY
            std::cout << "h5::open " << fname << std::endl;
        }
        catch (...) {
            std::cout << "Can't open " << fname << std::endl;
            throw;
        }
        return fd;
    }
    void h5read(const int& threadid, const h5::fd_t& fd)
    {
        const int sequence_max = 4;
        const int apa = threadid;
        const std::string tag = "orig";

        // comment the lock off to check thread safety
        std::lock_guard<std::mutex> guard(g_h5cpp_mutex);
        for (int sequence = 0; sequence < sequence_max; ++sequence) {
            auto key = WireCell::String::format("/%d/frame_%s%d", sequence, tag, apa);
            std::cout << "loading: " << key << std::endl;

            Eigen::MatrixXf d;
            try {
                d = h5::read<Eigen::MatrixXf>(fd, key);
            }
            catch (...) {
                std::cout << "Can't load " << key << std::endl;
                throw;
            }
        }
    }
}  // namespace

int main(int argc, const char* argv[])
{
    // if (argc != 3) {
    //   std::cout << "usage: threads <path-to-exported-script-module> nthreads\n";
    //   return -1;
    // }

    std::vector<std::string> fnames = {"orig-0.h5", "orig-1.h5", "orig-2.h5", "orig-3.h5"};

    std::vector<h5::fd_t> fds;
    for (auto fn : fnames) {
        fds.push_back(h5open(0, fn));
    }

    const int nthreads = fnames.size();

    std::cout << "test_h5cpp_threading\n";

    std::vector<std::thread> threads;
    for (int i = 0; i < nthreads; i++) {
        std::cout << "main() : creating thread, " << i << std::endl;
        // threads.push_back(std::thread(h5open, i, std::ref(fnames[i])));
        threads.push_back(std::thread(h5read, i, std::ref(fds[i])));
    }
    for (int i = 0; i < nthreads; i++) {
        threads[i].join();
    }

    return 0;
}