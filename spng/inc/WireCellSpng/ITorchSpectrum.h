/** FILL ME OUT
*/

#ifndef WIRECELLSPNG_ITORCHSPECTRUM
#define WIRECELLSPNG_ITORCHSPECTRUM

#include "WireCellUtil/IComponent.h"
#include <torch/torch.h>
// #include "WireCellSpng/hash_util.h"
#include <boost/compute/detail/lru_cache.hpp>

namespace WireCell {

class ITorchSpectrum : public IComponent<ITorchSpectrum> {
public:
    ITorchSpectrum() {};
    virtual ~ITorchSpectrum();

    /// Return the coldelec response data
    virtual torch::Tensor spectrum() const = 0;
    // virtual torch::Tensor spectrum(const std::vector<int64_t> & shape) = 0;

    virtual torch::Tensor spectrum(const std::vector<int64_t> & shape) = 0;

    /// Get the base shape of the response
    virtual const std::vector<int64_t> & shape() const {return m_shape;};

    /// Get any shifts of the response
    virtual std::vector<int64_t> shifts() const = 0;

protected:
    std::vector<int64_t> m_shape;
};

}  // namespace WireCell

#endif