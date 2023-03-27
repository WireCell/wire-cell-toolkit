#include "WireCellIface/ITensor.h"
#include "WireCellUtil/Dtype.h"

std::string WireCell::ITensor::dtype() const
{
    return WireCell::dtype(element_type());
}
