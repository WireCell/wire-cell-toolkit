#ifdef HAVE_LIBTORCH
#include "WireCellUtil/String.h"

#include "WireCellAux/SimpleTensorSet.h"
#include "WireCellAux/SimpleTensor.h"
#include "WireCellIface/ITensorSet.h"
#include "WireCellZio/FlowConfigurable.h"
#include "WireCellUtil/Exceptions.h"

#include "zio/domo/client.hpp"
#include "zio/tens.hpp"

#include "WireCellPytorch/Torch.h"  // One-stop header.

#include <iostream>
#include <memory>
#include <thread>
#include <chrono>

using namespace WireCell;

namespace {
    ITensorSet::pointer to_itensor(const std::vector<torch::IValue> &inputs)
    {
        ITensor::vector *itv = new ITensor::vector;

        int ind = 0;
        for (auto ival : inputs) {
            auto ten = ival.toTensor().cpu();
            // std::cout << "ten.shape: {"
            //           << ten.size(0) << ", "
            //           << ten.size(1) << ", "
            //           << ten.size(2) << ", "
            //           << ten.size(3) << "} "
            //           << std::endl;
            std::vector<size_t> shape = {(size_t) ten.size(0), (size_t) ten.size(1), (size_t) ten.size(2),
                                         (size_t) ten.size(3)};
            // TODO need to figure out type from dtyp
            Aux::SimpleTensor *st = new Aux::SimpleTensor(shape, (float*)nullptr);
            // std::cout << "st.shape: {"
            //           << st->shape()[0] << ", "
            //           << st->shape()[1] << ", "
            //           << st->shape()[2] << "} \n";
            size_t nbyte = 4;
            for (auto n : shape) nbyte *= n;
            // std::cout << "nbyte: " << nbyte << std::endl;
            auto data = (float *) st->data();
            // std::cout
            //     << "ten[0][0].data<float>(): " << ten[0][0].data<float>() << "\n"
            //     << "ten[0].data<float>(): " << ten[0].data<float>() << "\n"
            //     << "(float*)ten[0].data<float>(): " << (float *)ten[0].data<float>() << "\n";
            // auto src = (float *)ten[0].data<float>();
            // std::cout
            //     << "src[0]: " << src[0] << "\n"
            //     << "src[nbyte/4-1]: " << src[nbyte / 4 - 1] << "\n";
            memcpy(data, (float *) ten[0][0].data<float>(), nbyte);
            itv->push_back(ITensor::pointer(st));
            ++ind;
        }

        int seqno = 0;
        Configuration md;
        // Json::Reader reader;
        // reader.parse(label[zio::tens::form]["metadata"].dump(), md);

        return std::make_shared<Aux::SimpleTensorSet>(seqno, md, ITensor::shared_vector(itv));
    }
    torch::IValue from_itensor(const ITensorSet::pointer &inputs)
    {
        if (inputs->tensors()->size() != 1) {
            THROW(ValueError() << errmsg{"inputs->tensors()->size()!=1"});
        }
        auto iten = inputs->tensors()->front();
        if (iten->shape().size() != 4) {
            THROW(ValueError() << errmsg{"iten->shape().size()!=4"});
        }
        // TODO determine data type from metadata
        auto ten = torch::from_blob((float *) iten->data(), {(long) iten->shape()[0], (long) iten->shape()[1],
                                                             (long) iten->shape()[2], (long) iten->shape()[3]});

        // std::cout << "ten.shape: {"
        //           << ten.size(0) << ", "
        //           << ten.size(1) << ", "
        //           << ten.size(2) << ", "
        //           << ten.size(3) << "} "
        //           << std::endl;
        return ten;
    }
    std::string dump(const zio::Message &msg)
    {
        std::stringstream ss;
        ss << "zio.Message: ";
        ss << "ZIO" << msg.level() << msg.form() << msg.label();
        ss << " + [0x" << msg.origin() << "," << msg.granule() << "," << msg.seqno() << "]";
        ss << " + [" << msg.payload().str() << "]";
        return ss.str();
    }
    std::string dump(const ITensorSet::pointer &itens)
    {
        std::stringstream ss;
        ss << "ITensorSet: ";
        Json::FastWriter jwriter;
        ss << itens->ident() << ", " << jwriter.write(itens->metadata());
        for (auto iten : *itens->tensors()) {
            ss << "shape: [";
            for (auto l : iten->shape()) {
                ss << l << " ";
            }
            ss << "]\n";
        }
        return ss.str();
    }
}  // namespace

int main()
{
    // client
    std::shared_ptr<zio::domo::Client> m_client;
    zmq::context_t ctx;
    zmq::socket_t sock(ctx, ZMQ_CLIENT);
    try {
        m_client = std::make_shared<zio::domo::Client>(sock, "tcp://localhost:5555");
    }
    catch (...) {
        THROW(RuntimeError() << errmsg{"Client init faileds!"});
    }

    //
    std::vector<torch::jit::IValue> inputs;
    inputs.push_back(torch::ones({1, 3, 800, 600}, torch::dtype(torch::kFloat32).device(torch::kCUDA, 0)));

    auto iitens = to_itensor(inputs);
    std::cout << "torch -> itens ... OK\n";

    auto msg = Zio::FlowConfigurable::pack(iitens);
    std::cout << "itens -> msg ... OK\n";
    // std::cout << dump(msg) << "\n";

    zmq::multipart_t mmsg(msg.toparts());
    // std::cout << mmsg.str() << "\n";
    m_client->send("torch:dnnroi", mmsg);
    std::cout << "Client::send ... OK\n";

    mmsg.clear();
    m_client->recv(mmsg);
    std::cout << "Client::recv ... OK\n";
    // std::cout << mmsg.str() << "\n";

    msg.fromparts(mmsg);
    // std::cout << dump(msg) << "\n";
    auto oitens = Zio::FlowConfigurable::unpack(msg);
    std::cout << "msg -> itens ... OK\n";
    // std::cout << dump(oitens);

    auto out_ival = from_itensor(oitens);
    std::cout << "itens -> torch ... OK\n";

    return 0;
}
#else
int main() { return 0; }
#endif
