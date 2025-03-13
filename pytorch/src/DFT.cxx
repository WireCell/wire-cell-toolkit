#include "WireCellPytorch/DFT.h"
#include "WireCellUtil/NamedFactory.h"

#include "WireCellPytorch/Torch.h"
#include <torch/csrc/api/include/torch/fft.h>


WIRECELL_FACTORY(TorchDFT, WireCell::Pytorch::DFT,
                 WireCell::IDFT,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Pytorch;

DFT::DFT()
{
}

DFT::~DFT()
{
}

Configuration DFT::default_configuration() const
{
    Configuration cfg;

    // one of: {cpu, gpu, gpuN} where "N" is a GPU number.  "gpu"
    // alone will use GPU 0.
    cfg["device"] = "cpu";
    return cfg;
}

void DFT::configure(const WireCell::Configuration& cfg)
{
    auto dev = get<std::string>(cfg, "device", "cpu");
    m_ctx.connect(dev);
}


using torch_transform = std::function<torch::Tensor(const torch::Tensor&)>;

static
void doit(const TorchContext& ctx,
          const IDFT::complex_t* in, IDFT::complex_t* out,
          int64_t nrows, int64_t ncols, // 1d vec should have nrows=1
          torch_transform func)
{
    TorchSemaphore sem(ctx);
    torch::NoGradGuard no_grad;

    int64_t size = nrows*ncols;

    auto dtype = torch::TensorOptions().dtype(torch::kComplexFloat);

    // 1) in->src
    if (in != out) {            // from_blob() doesn't like const data
        memcpy(out, in, sizeof(IDFT::complex_t)*size);
    }

    torch::Tensor src = torch::from_blob(out, {nrows, ncols}, dtype);

    // 2) dst = func(src)
    src = src.to(ctx.device());
    auto dst = func(src);

    // Making contiguous costs a copy but gets the data in row-major
    // so the (2nd) copy next actually gives correct results.  This
    // corrects optimizations that libtorch makes for transpose (and
    // others) eg when our func is 1b.  Alternatively, may avoid both
    // copies by iterating over indices but presumably (?) that is
    // slower.  Likewise we make contiguous on the device as that is
    // presumably (?) faster when the device is GPU.
    dst = dst.contiguous().cpu();

    // 3) dst->out
    memcpy(out, dst.data_ptr(), sizeof(IDFT::complex_t)*size);
}


void DFT::fwd1d(const IDFT::complex_t* in, IDFT::complex_t* out, int size) const
{
    doit(m_ctx, in, out, 1, size,
         [](const torch::Tensor& src) { return torch::fft::fft(src); });
}


void DFT::inv1d(const IDFT::complex_t* in, IDFT::complex_t* out, int size) const
{
    doit(m_ctx, in, out, 1, size, // fixme: check norm
         [](const torch::Tensor& src) { return torch::fft::ifft(src); });
}


void DFT::fwd1b(const IDFT::complex_t* in, IDFT::complex_t* out,
                int nrows, int ncols, int axis) const
{
    doit(m_ctx, in, out, nrows, ncols, [&](const torch::Tensor& src) {
        return torch::fft::fft2(src, torch::nullopt, {axis}); });
}


void DFT::inv1b(const IDFT::complex_t* in, IDFT::complex_t* out,
                int nrows, int ncols, int axis) const
{
    doit(m_ctx, in, out, nrows, ncols, [&](const torch::Tensor& src) {
        return torch::fft::ifft2(src, torch::nullopt, {axis}); });
}       

        
void DFT::fwd2d(const IDFT::complex_t* in, IDFT::complex_t* out,
                int nrows, int ncols) const
{
    doit(m_ctx, in, out, nrows, ncols,
         [](const torch::Tensor& src) { return torch::fft::fft2(src); });
}


void DFT::inv2d(const IDFT::complex_t* in, IDFT::complex_t* out,
                int nrows, int ncols) const
{
    doit(m_ctx, in, out, nrows, ncols,
         [](const torch::Tensor& src) { return torch::fft::ifft2(src); });
}

