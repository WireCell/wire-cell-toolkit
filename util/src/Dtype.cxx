#include "WireCellUtil/Dtype.h"

size_t WireCell::dtype_size(const std::string& dt)
{
    if (dt == "c") return sizeof(char);
    if (dt == "B") return sizeof(std::byte);
    if (dt == "i1" or dt == "u1") return sizeof(int8_t);
    if (dt == "i2" or dt == "u2") return sizeof(int16_t);
    if (dt == "i4" or dt == "u4") return sizeof(int32_t);
    if (dt == "i8" or dt == "u8") return sizeof(int64_t);
    if (dt == "f4") return sizeof(float);
    if (dt == "f8") return sizeof(double);
    if (dt == "c8") return sizeof(std::complex<float>);
    if (dt == "c16") return sizeof(std::complex<double>);
    return 0;
}
std::string WireCell::dtype(const std::type_info& ti)
{
    if (ti == typeid(char)) return "c";
    if (ti == typeid(std::byte)) return "B";

    if (ti == typeid(int8_t))  return "i1";
    if (ti == typeid(int16_t)) return "i2";
    if (ti == typeid(int32_t)) return "i4";
    if (ti == typeid(int64_t)) return "i8";

    if (ti == typeid(uint8_t))  return "u1";
    if (ti == typeid(uint16_t)) return "u2";
    if (ti == typeid(uint32_t)) return "u4";
    if (ti == typeid(uint64_t)) return "u8";

    if (ti == typeid(float)) return "f4";
    if (ti == typeid(double)) return "f8";
    if (ti == typeid(std::complex<float>)) return "c8";
    if (ti == typeid(std::complex<double>)) return "c16";
    return "";
}
const std::type_info& WireCell::dtype_info(const std::string& dt)
{
    if (dt == "c") return typeid(char);
    if (dt == "B") return typeid(std::byte);
    if (dt == "i1" or dt == "u1") return typeid(int8_t);
    if (dt == "i2" or dt == "u2") return typeid(int16_t);
    if (dt == "i4" or dt == "u4") return typeid(int32_t);
    if (dt == "i8" or dt == "u8") return typeid(int64_t);
    if (dt == "f4") return typeid(float);
    if (dt == "f8") return typeid(double);
    if (dt == "c8") return typeid(std::complex<float>);
    if (dt == "c16") return typeid(std::complex<double>);
    return typeid(nullptr);
}
