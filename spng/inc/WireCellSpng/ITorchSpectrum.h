#ifndef WIRECELLSPNG_ITORCHSPECTRUM
#define WIRECELLSPNG_ITORCHSPECTRUM

#include "WireCellSpng/Torch.h"

#include "WireCellUtil/IComponent.h"

namespace WireCell {

/** ITorchSpectrum provides "spectra" as torch tensors of requested shape.
 *
 * On the surface, this interface merely returns tensors.  However,
 * implementations are expected to satisfy the following contract:
 * 
 * - Tensors provide data in Fourier-space representation.  Aka "frequency" or
 * "periodicity" and NOT "interval" aka "channel" aka "wire" aka "time"
 * representation.  When the tensor has real value, the spectrum is interpreted
 * as a spectral amplitude (eg the abs() of complex-valued samples).
 *
 * - The spectrum tensor spans the full Fourier domain and not the "half" below
 * the Nyquist frequency.  The caller is free to truncate half in order to use
 * irfft().
 *
 * - The sampling period, if meaningful, must be communicated out of band.
 *
 * - The spectrum(shape) method produces a tensor that represents an
 * interpolation (and not a padding) in Fourier space.  
 * 
 * FIXME: the shifts() method is deprecated and will be removed.  All
 * implementations SHOULD provide spectra that do not induce artificial shifts
 * when used as a convolution kernel.  To assure that, the spectra should be
 * prepared following these rules:
 *
 * - A kernel tensor dimension along which the content has mirror symmetry must
 * place the reflection point in sample zero.
 *
 * - A kernel tensor dimension along which the content has an origin must place
 * that origin point at sample zero.
 *
 * In the example of an FR, the channel dimension has mirror symmetry and the
 * time dimension has an origin.
 *
 * After existing implementations have their artificial shifts removed, shifts()
 * method will be deleted.
 *
 * Note: shape() provides required information related handling to "logical"
 * shifts incurred by interpolation.
 *
 * See the spng/docs/decon.org document for more details on proper padding and
 * shifts and how to deal with any wrap-around when a spectrum is used.
 */
class ITorchSpectrum : public IComponent<ITorchSpectrum> {
public:
    ITorchSpectrum() {};
    virtual ~ITorchSpectrum();

    /// Return the spectrum interpolated into the given shape.
    ///
    using shape_t = std::vector<int64_t>;
    virtual torch::Tensor spectrum(const shape_t& shape) const = 0;

    /// Return the "natural shape" of the spectrum.
    ///
    /// A spectrum based on an initial sampling should provide the shape of that
    /// sampled tensor.  A spectrum that is produced by sampling a continuous
    /// function should return zeros for each dimension.
    virtual shape_t shape() const = 0;


    /// deprecated.
    virtual shape_t shifts() const = 0;

};

}  // namespace WireCell

#endif
