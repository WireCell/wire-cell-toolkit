#include "WireCellPytorch/Torch.h"
#include <torch/csrc/api/include/torch/fft.h>

#include <vector>
#include <complex>
#include <iostream>

using complex_t = std::complex<float>;

void dump(const std::vector<complex_t>& v, int nrows, int ncols, const std::string& msg="")
{
    std::cerr << msg << ": ("<<nrows<<","<<ncols<<")\n";
    for (int irow = 0; irow<nrows; ++irow) {
        for (int icol = 0; icol<ncols; ++icol) {
            auto c =  v[irow*ncols + icol];
            std::cerr<< "(" << std::real(c) << "," << std::imag(c) << ")\t";
        }
        std::cerr << "\n";
    }
}

void dump(const at::Tensor& ten, const std::string& msg="")
{
    int nrows = ten.size(0);
    int ncols = ten.size(1);
    std::cerr << msg << ": ("<<nrows<<","<<ncols<<")\n"
              << ten << "\n";
}

int main()
{
    const int nrows=3, ncols=8;
    int size = nrows*ncols;
    std::vector<complex_t> v(size, 0);
    //std::iota(v.begin(), v.end(), 0);
    v[ncols+2] = 1.0;
    dump(v, nrows, ncols, "v");

    // Note: gpu is almost 10x SLOWER than CPU due to kernel load time!
    // auto device = at::Device(at::kCPU);
    auto device = at::Device(at::kCUDA);

    auto typ_options = at::TensorOptions().dtype(at::kComplexFloat);
    // auto dev_options = typ_options.device(device);
    
    for (int axis = 0; axis < 2; ++ axis) {
        at::Tensor src = at::from_blob(v.data(), {nrows, ncols}, typ_options);
        dump(src, "src");

        src = src.to(device);
        // In pytorch dim=(0,) transforms along columns, ie follows
        // numpy.fft convention.  Both directions work on both CPU and
        // CUDA.
        at::Tensor dst = torch::fft::fft2(src, {}, {axis,});

        // BUT, BEWARE that the underlying storage will NOT reflect
        // logical row-major ordering.  Indexing is as expected but
        // memory returned by data_ptr() will reflect transpose
        // optimizations.  At the expense of a copy the contiguous()
        // method provides expected row-major storage order.  
        dst = dst.contiguous().cpu();

        dump(dst, "dst");

        std::vector<complex_t> v2(size, 0);
        memcpy(v2.data(), dst.data_ptr(), sizeof(complex_t)*size);

        std::cerr << "axis=" << axis << " dim=" << axis
                  << " shape=(" << src.size(0) << "," << src.size(1) << ")\n";
        dump(v2, nrows, ncols, "dft(v)");
    }
    return 0;
}
