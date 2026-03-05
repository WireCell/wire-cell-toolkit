#pragma once

#include "WireCellSpng/TensorFilter.h"

#include "WireCellUtil/HanaJsonCPP.h"
#include "WireCellIface/IConfigurable.h"

namespace WireCell::SPNG {
    struct RebaselinerConfig {

        /// The dimension over which to rebaseline.
        int dim = -1;

        /// The number of consecutive zeros required to delimit an ROI.  A
        /// zero-crossing inside an ROI could very rarely hit exactly zero and
        /// even more rarely hit exactly zero twice in a row, etc.  Making this
        /// number too big will combine ROIs that are truly separated.  
        int consequtive_zeros = 2;

        /// Define a "small" ROI.
        int min_roi_size = 1;

        /// Reduce the ROI by this much on either size.
        int shrink_size = 0;

        /// If true, remove "small" ROIs by setting their values to zero.
        bool remove_small = false;

        /// If true, after rebaseline, clamp all values to be non-negative.
        bool remove_negative = false;


    };
}

BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::RebaselinerConfig,
                        dim, consequtive_zeros,
                        min_roi_size, shrink_size,
                        remove_small, remove_negative);

namespace WireCell::SPNG {

    /// Apply ROI baseline.  See rebaseline_zero() in Rebaseline.h for details of algorithm.
    struct Rebaseliner : public TensorFilter,
                         virtual public IConfigurable {

        Rebaseliner();
        virtual ~Rebaseliner() = default ;

        /// TensorFilter
        virtual ITorchTensor::pointer filter_tensor(const ITorchTensor::pointer& in);

        // IConfigurable - see KernelConvolveConfig for configuration documentation.
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;
        
    private:

        RebaselinerConfig m_config;

    };
}
