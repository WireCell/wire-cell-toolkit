#ifndef WIRECELL_SPNG_RATIOKERNEL
#define WIRECELL_SPNG_RATIOKERNEL

#include "WireCellSpng/ContextBase.h"
#include "WireCellSpng/Logger.h"
#include "WireCellSpng/ITorchSpectrum.h"

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IFilterWaveform.h"


namespace WireCell::SPNG {
    /// A kernel that is the simple Fourier-space ratio of two kernels.
    /// No convolution is done.
    struct RatioKernel : public ContextBase, 
                         public Logger,
                         public ITorchSpectrum, virtual public IConfigurable {
        RatioKernel();
        virtual ~RatioKernel() = default;

        /// Return the ratio spectrum.
        virtual torch::Tensor spectrum(const shape_t & shape) const;

        /// Ratios are sampled from analytic functions which have no native
        /// size/shape.  
        virtual shape_t shape() const;

        // Deprecated.  Returns zeros.
        virtual shape_t shifts() const { return shape_t(2, 0); }

        // IConfigurable
        ///
        /// Configuration:
        ///
        /// - numerator :: the WCT type/name of an ITorchSpectrum supplying numerator.
        /// - denominator :: the WCT type/name of an ITorchSpectrum supplying denominator.
        virtual void configure(const WireCell::Configuration& config);
        // Rely on ContextBase
        virtual WireCell::Configuration default_configuration() const;

    private:

        ITorchSpectrum::pointer m_num, m_den;

    };
}

#endif
