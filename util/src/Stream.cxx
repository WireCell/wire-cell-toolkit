#include "WireCellUtil/Stream.h"


std::ostream& WireCell::Stream::write(
    std::ostream& so,
    const std::string& fname,
    const std::byte* bytes,
    const WireCell::Stream::shape_t& shape,
    const std::string& dtype_name)
{
    pigenc::File pig;
    pigenc::Header& phead = pig.header();
    const bool fortran_order = false;
    phead.set(shape, dtype_name, fortran_order);
    auto& dat = pig.data();
    const char* cbytes = (const char*)bytes; // match pigenc byte type
    dat.insert(dat.begin(), cbytes, cbytes + phead.data_size());
    size_t fsize = phead.file_size();
    custard::write(so, fname, fsize);
    if (!so) return so;
    return pig.write(so);
}


