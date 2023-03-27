#include "WireCellUtil/PointCloud.h"

using namespace WireCell::PointCloud;

//// array

// Copy constructor
Array::Array(const Array& rhs)
{
    (*this) = rhs;
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



//// dataset

bool Dataset::add(const std::string& name, const Array& arr)
{
    if (m_store.empty()) {
        m_store.insert({name,arr});
        return true;
    }
    if (size_major() != arr.size_major()) {
        THROW(WireCell::ValueError() << WireCell::errmsg{"major axis size mismatch when adding \""+name+"\" to dataset"});
    }
    auto it = m_store.find(name);
    if (it == m_store.end()) {
        m_store.insert({name,arr});
        return true;
    }
    it->second = arr;
    return false;
}

selection_t Dataset::selection(const std::vector<std::string>& names) const
{
    selection_t ret;
    for (const auto& name : names) {
        auto it = m_store.find(name);
        if (it == m_store.end()) {
            return selection_t();
        }
        ret.push_back(std::cref(it->second));
    }
    return ret;
}


Dataset::keys_t store_keys(const Dataset::store_t& store)
{
    Dataset::keys_t ret;
    for (const auto& [k,v] : store) {
        ret.insert(k);
    }
    return ret;
}

Dataset::keys_t difference(const Dataset::keys_t& a, const Dataset::keys_t& b)
{
    Dataset::keys_t diff;
    std::set_difference(a.begin(), a.end(),
                        b.begin(), b.end(),
                        std::inserter(diff, diff.end()));
    return diff;
}

Dataset Dataset::zeros_like(size_t nmaj)
{
    Dataset ret;
    for (auto& [key, arr] : m_store) {
        ret.add(key, arr.zeros_like(nmaj));
    }
    return ret;
}

void Dataset::append(const std::map<std::string, Array>& tail)
{
    if (tail.empty()) return;

    const size_t nmaj_before = size_major();

    auto diff = difference(keys(), store_keys(tail));
    if (diff.size() > 0) {
        THROW(ValueError() << errmsg{"missing keys in append"});
    }

    size_t nmaj = 0;            // check rectangular
    for (auto& [key, arr] : tail) {
        if (nmaj == 0) {        // first
            nmaj = arr.size_major();
            continue;
        }
        if (arr.size_major() != nmaj) {
            THROW(ValueError() << errmsg{"variation in major axis size in append"});
        }
    }

    for (auto& [key, arr] : m_store) {
        const auto& [ok, oa] = *tail.find(key);
        arr.append(oa.bytes());
    }

    const size_t nmaj_after = size_major();

    for (auto cb : m_append_callbacks) {
        cb(nmaj_before, nmaj_after);
    }
}


