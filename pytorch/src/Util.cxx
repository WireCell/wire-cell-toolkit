
#include "WireCellPytorch/Util.h"

#include "WireCellUtil/Exceptions.h"
#include "WireCellAux/SimpleTensorSet.h"
#include "WireCellAux/SimpleTensor.h"

using namespace WireCell;

std::string Pytorch::dump(const torch::Tensor &ten)
{
    std::stringstream ss;
    ss << "torch::Tensor: ";
    ss << " shape {";
    for (int idim = 0; idim < ten.dim(); ++idim) {
        ss << ten.size(idim) << " ";
    }
    ss << "} ";
    return ss.str();
}

ITensorSet::pointer Pytorch::to_itensor(const std::vector<torch::IValue> &inputs)
{
    ITensor::vector *itv = new ITensor::vector;

    int ind = 0;
    for (auto ival : inputs) {
        auto ten = ival.toTensor().cpu();
        if (ten.dim() != 4) {
            THROW(ValueError() << errmsg{"ten.dim() != 4"});
        }
        std::vector<size_t> shape = {(size_t) ten.size(0), (size_t) ten.size(1), (size_t) ten.size(2),
                                     (size_t) ten.size(3)};
        // TODO need to figure out type from dtyp
        Aux::SimpleTensor *st = new Aux::SimpleTensor(shape, (float*)nullptr);
        size_t nbyte = 4;
        for (auto n : shape) nbyte *= n;
        auto data = (float *) st->data();
        memcpy(data, (float *) ten[0][0].data_ptr<float>(), nbyte);
        itv->push_back(ITensor::pointer(st));
        ++ind;
    }

    // FIXME use a correct seqno
    int seqno = 0;
    Configuration md;

    return std::make_shared<Aux::SimpleTensorSet>(seqno, md, ITensor::shared_vector(itv));
}

std::vector<torch::IValue> Pytorch::from_itensor(const ITensorSet::pointer &inputs, const bool gpu)
{
    std::vector<torch::IValue> ret;

    for (auto iten : *inputs->tensors()) {
        if (iten->shape().size() != 4) {
            THROW(ValueError() << errmsg{"iten->shape().size()!=4"});
        }
        // TODO determine data type from metadata
        auto ten = torch::from_blob((float *) iten->data(), {(long) iten->shape()[0], (long) iten->shape()[1],
                                                             (long) iten->shape()[2], (long) iten->shape()[3]});
        if (gpu) {
            ten = ten.to(torch::Device(torch::kCUDA, 0));
            assert(ten.device().type() == torch::kCUDA);
            ret.push_back(ten);
            //ret.push_back(ten.cuda());
        }
        else {
            ret.push_back(ten);
        }
    }

    return ret;
}
