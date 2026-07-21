#include "WireCellSpng/Stream.h"
#include "WireCellUtil/Stream.h"

using namespace WireCell;

// Return the pigenc dtype string.
std::string Stream::dtype_of(torch::Tensor ten)
{
    auto dt = ten.dtype();

    if (dt == at::kByte) return "<u1";
    if (dt == at::kChar) return "<i1";
    if (dt == at::kShort) return "<i2";
    if (dt == at::kInt) return "<i4";
    if (dt == at::kLong) return "<i8";
    if (dt == at::kHalf) return "<f2";
    if (dt == at::kFloat) return "<f4";
    if (dt == at::kDouble) return "<f8";
    if (dt == at::kComplexHalf) return "<c2";
    if (dt == at::kComplexFloat) return "<c4";
    if (dt == at::kComplexDouble) return "<c8";
    return "";
}


torch::TensorOptions Stream::tensoroptions_from_dtype(const std::string& dtype)
{
    torch::TensorOptions to;
    if (dtype == "<u1") return to.dtype(torch::kByte);
    if (dtype == "<i1") return to.dtype(torch::kChar);
    if (dtype == "<i2") return to.dtype(torch::kShort);
    if (dtype == "<i4") return to.dtype(torch::kInt);
    if (dtype == "<i8") return to.dtype(torch::kLong);
    if (dtype == "<f2") return to.dtype(torch::kHalf);
    if (dtype == "<f4") return to.dtype(torch::kFloat);
    if (dtype == "<f8") return to.dtype(torch::kDouble);
    if (dtype == "<c2") return to.dtype(torch::kComplexHalf);
    if (dtype == "<c4") return to.dtype(torch::kComplexFloat);
    if (dtype == "<c8") return to.dtype(torch::kComplexDouble);
    return to;
}


Stream::shape_t Stream::shape_of(torch::Tensor ten)
{
    auto sizes = ten.sizes();
    Stream::shape_t shape(sizes.begin(), sizes.end());
    return shape;
}

std::ostream& Stream::write(std::ostream& so, const std::string& fname, torch::Tensor ten)
{
    // Stream protocol: each entry has:
    // 1) header: file name + file size
    // 2) payload
    //
    // npz payload has two parts:
    // 2.a) pig header
    // 2.b) pig data

    pigenc::Header head;
    auto shape = Stream::shape_of(ten);
    head.set(shape, dtype_of(ten));

    // Member header
    size_t fsize = head.file_size();
    custard::write(so, fname, fsize);

    // payload header
    head.write(so);

    // payload data
    // Avoid the copy needed to use pigenc::File and write more directly.
    auto here = ten.to(torch::kCPU).contiguous();
    so.write((const char*)here.data_ptr(), here.element_size() * here.numel());
    return so;
}

std::istream& Stream::read(std::istream& si, std::string& fname, torch::Tensor& ten)
{
    size_t fsize=0;
    custard::read(si, fname, fsize);
    if (!si) return si;
    pigenc::Header head;
    head.read(si);
    if (!si) return si;

    auto shape = head.shape();
    std::vector<long int> sizes(shape.begin(), shape.end());

    ten = torch::empty(sizes, tensoroptions_from_dtype(head.dtype()));
    si.read((char*)ten.data_ptr(), head.data_size());
    return si;
}

