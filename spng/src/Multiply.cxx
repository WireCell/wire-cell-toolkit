#include "WireCellSpng/Multiply.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(SPNGMultiply,
                 WireCell::SPNG::Multiply,
                 WireCell::SPNG::Multiply::fan_type,
                 WireCell::IConfigurable);

namespace WireCell::SPNG {

    Multiply::Multiply() 
        : FaninBase<ITorchTensor>("Multiply", "spng")
    {
    }

    void Multiply::fanin_combine(const input_vector& inv, output_pointer& out)
    {
        Configuration md = inv[0]->metadata();
        torch::Tensor product = inv[0]->tensor().clone();

        for (size_t ind=1; ind < inv.size(); ++ind) {
            update(md, inv[ind]->metadata());
            product = torch::mul(product, inv[ind]->tensor());
        }
        out = std::make_shared<SimpleTorchTensor>(product, md);
        logit(out, "tensor product");        
    }



}
