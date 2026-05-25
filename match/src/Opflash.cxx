#include "WireCellMatch/Opflash.h"

#include "WireCellUtil/Exceptions.h"

// Boost 1.82's multi_array transitively includes boost/functional.hpp which
// uses std::unary_function/binary_function; gcc 12+ flags these as
// -Werror=deprecated-declarations. Silence around the boost include only.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <boost/multi_array.hpp>
#pragma GCC diagnostic pop

#include <algorithm>

using namespace WireCell;
using namespace WireCell::Match;

Opflash::Opflash(const ITensor::pointer ten, int idx, double threshold, int nchan)
{
    if (ten->shape().size() != 2) {
        raise<ValueError>("Opflash: input tensor dim %d != 2", ten->shape().size());
    }
    const int nrow = ten->shape()[0];
    const int ncol = ten->shape()[1];
    if (nrow < idx + 1) raise<ValueError>("Opflash: input tensor nrow %d < idx+1", nrow);
    if (ncol < nchan + 1) raise<ValueError>("Opflash: input tensor ncol %d < nchan+1", ncol);

    flash_id    = idx;
    m_nchan     = nchan;
    m_threshold = threshold;
    type        = 0;
    low_time    = 0;
    high_time   = 0;

    using MultiArray = boost::multi_array<double, 2>;
    boost::array<MultiArray::index, 2> shape = {nrow, ncol};
    boost::multi_array_ref<double, 2> ten_data((double*)ten->data(), shape);

    time     = ten_data[idx][0];
    total_PE = 0;
    PE.resize(nchan, 0);
    PE_err.resize(nchan, 1);

    for (int i = 0; i < nchan; ++i) {
        PE[i] = ten_data[idx][i + 1];
        PE_err[i] = (PE[i] < 1) ? 0.3 : 0.3 * PE[i];
        total_PE += PE[i];
        if (PE[i] > threshold) fired_channels.push_back(i);
    }
}

Opflash::~Opflash() = default;

bool Opflash::get_fired(int ch) const
{
    return std::find(fired_channels.begin(), fired_channels.end(), ch) != fired_channels.end();
}
