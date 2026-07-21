#include "WireCellSpng/TensorPacker.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"
#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(SPNGTensorPacker,
                 WireCell::SPNG::TensorPacker,
                 WireCell::SPNG::TensorPacker::fan_type,
                 WireCell::IConfigurable);

namespace WireCell::SPNG {
    
    TensorPacker::TensorPacker()
        : FaninBase<ITorchTensor, ITorchTensorSet>("TensorPacker", "spng")
    {
    }

    void TensorPacker::fanin_combine(const input_vector& inv, output_pointer& out)
    {
        Configuration md;
        out = std::make_shared<SimpleTorchTensorSet>(m_count, md, inv);
        logit(out, "tensor packer");        
    };

}
