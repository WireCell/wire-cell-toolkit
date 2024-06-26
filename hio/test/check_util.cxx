#include "WireCellHio/Util.h"

using namespace WireCell;

namespace test {

    template <class T, class... args_t>
    h5::ds_t write(const h5::fd_t &fd, const std::string &dataset_path, const T *ptr, args_t &&... args)
    {
        return h5::write(fd, dataset_path, ptr, args...);
    }
}  // namespace test

int main()
{
    const std::string fn = "tmp.h5";
    h5::create(fn, H5F_ACC_TRUNC);
    h5::fd_t fd = h5::open(fn, H5F_ACC_RDWR);
    int d[2] = {42, 45};
    h5::write<int>(fd, "/d", d, h5::count{2});
    return 0;
}