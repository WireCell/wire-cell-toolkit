#include "WireCellUtil/PointCloudArray.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Dtype.h"

using namespace WireCell::PointCloud;

//// array

Array::~Array()
{
}

Array::Array()
{
}

Array::Array(const Array& rhs)
    : m_shape(rhs.m_shape)
    , m_ele_size(rhs.m_ele_size)
    , m_dtype(rhs.m_dtype)
    , m_store(rhs.m_store)
    , m_metadata(rhs.m_metadata)
{
    update_span();
}

Array::Array(Array&& rhs)
    : m_shape(std::move(rhs.m_shape))
    , m_ele_size(std::move(rhs.m_ele_size))
    , m_dtype(std::move(rhs.m_dtype))
    , m_store(std::move(rhs.m_store))
    , m_bytes(std::move(rhs.m_bytes))
    , m_metadata(std::move(rhs.m_metadata))
{
}

Array::Array(const std::byte* data, const std::string& dtype, const shape_t& shape)
{
    assign(data,dtype,shape);
}
Array::Array(std::byte* data, const std::string& dtype, const shape_t& shape, bool share)
{
    assign(data,dtype,shape,share);
}

void Array::assign(const std::byte* data, const std::string& dtype, const shape_t& shape)
{
    assign(const_cast<std::byte*>(data), dtype, shape, false);
}

void Array::assign(std::byte* data, const std::string& dtype, const shape_t& shape, bool share)
{
    m_store.clear();
    m_shape = shape;
    m_dtype = dtype;
    size_t nbytes = m_ele_size = dtype_size(dtype);

    for (const auto& n : shape) {
        nbytes *= n;
    }
    if (share) {
        m_bytes = span_t<std::byte>(data, nbytes);
    }
    else {
        m_store.assign(data, data+nbytes);
        update_span();
    }
}

// Assignment 
Array& Array::operator=(const Array& rhs)
{
    clear();
    m_shape = rhs.m_shape;
    m_ele_size = rhs.m_ele_size;
    m_dtype = rhs.m_dtype;
    if (rhs.m_store.empty()) { // shared
        m_bytes = span_t<std::byte>(rhs.m_bytes.data(), rhs.m_bytes.size());
    }
    else {
        m_store = rhs.m_store;
        update_span();
    }
    m_metadata = rhs.m_metadata;
    return *this;
}

Array& Array::operator=(Array&& rhs)
{
    clear();
    std::swap(m_shape, rhs.m_shape);
    std::swap(m_ele_size, rhs.m_ele_size);
    std::swap(m_dtype, rhs.m_dtype);
    std::swap(m_bytes, rhs.m_bytes);
    std::swap(m_metadata, rhs.m_metadata);
    return *this;
}

Array Array::slice(size_t position, size_t count, bool share) 
{
    const size_t ndims = m_shape.size();
    shape_t shape = m_shape;
    shape[0] = count;

    // The number of bytes into the flattened (possibly N-d) array where major
    // axis element position starts.
    size_t start_bytes = m_ele_size * position;
    for (size_t idim=1; idim < ndims; ++idim) {
        start_bytes *= m_shape[idim];
    }

    // Build array on the slice 
    return Array(m_store.data() + start_bytes, m_dtype, shape, share);
}
Array Array::slice(size_t position, size_t count) const
{
    const size_t ndims = m_shape.size();
    shape_t shape = m_shape;
    shape[0] = count;

    // The number of bytes into the flattened (possibly N-d) array where major
    // axis element position starts.
    size_t start_bytes = m_ele_size * position;
    for (size_t idim=1; idim < ndims; ++idim) {
        start_bytes *= m_shape[idim];
    }

    // Build array on the slice 
    return Array(m_store.data() + start_bytes, m_dtype, shape);
}

Array Array::zeros_like(size_t nmaj)
{
    Array ret;
    ret.m_metadata = m_metadata;
    ret.m_dtype = m_dtype;
    ret.m_ele_size = m_ele_size;
    if (m_shape.empty()) {
        ret.update_span();
        return ret;
    }
    auto shape = m_shape;
    shape[0] = nmaj;
    size_t size = m_ele_size;
    for (auto s : shape) {
        size *= s;
    }
    ret.m_store.resize(size, std::byte(0));
    ret.m_shape = shape;
    ret.update_span();
    return ret;
}

/// Return the total number of elements.
size_t Array::num_elements() const
{
    if (m_bytes.empty()) {
        return 0;
    }
    if (m_shape.empty()) {
        return 0;
    }

    size_t num = 1;
    for (const auto& s : m_shape) {
        num *= s;
    }
    return num;
}            

static bool append_compatible(const Array::shape_t& s1, const Array::shape_t& s2)
{
    size_t size = s1.size();

    if (size != s2.size()) return false; // same dimensions
    if (size == 1) return true;          // all 1D append okay
    for (size_t ind=1; ind<size; ++ind) {
        if (s1[ind] != s2[ind]) return false;
    }
    return true;
}

void Array::append(const Array& tail)
{
    const auto tshape = tail.shape();
    if (not append_compatible(m_shape, tshape)) {
        THROW(ValueError() << errmsg{"array append with incompatible shape"});
    }
    append(tail.bytes());
}

void Array::append(const std::byte* data, size_t nbytes)
{
    if (nbytes % m_ele_size) {
        THROW(ValueError() << errmsg{"byte append not compatible with existing type"});
    }
    const size_t nelem = nbytes / m_ele_size;

    size_t rows = 0, notrows=1;
    for (auto s : m_shape) {
        if (!rows) {
            rows = s;
            continue;
        }
        notrows *= s;
    }

    if (nelem % notrows) {
        THROW(ValueError() << errmsg{"byte append not compatible with existing shape"});
    }
    const size_t nrows = nelem/notrows;

    assure_mutable();
    m_store.insert(m_store.end(), data, data+nbytes);
    m_shape[0] += nrows;
    update_span();
}


bool Array::operator==(const Array& other) const
{
    // same size, shape, type, content and metadata.  Cheapest
    // test first.  bail early.
    if (this->size_major() != other.size_major()) return false;
        
    const auto s1 = this->shape();
    const auto s2 = other.shape();
    if (s1.size() != s2.size()) return false;
    if (! std::equal(s1.begin(), s1.end(), s2.begin())) return false;
    const auto b1 = this->bytes();
    const auto b2 = other.bytes();
    if (! std::equal(b1.begin(), b1.end(), b2.begin())) return false;
    if (this->metadata() != other.metadata()) return false;
    return true;
}
bool Array::operator!=(const Array& other) const 
{
    return !(*this == other);
}
