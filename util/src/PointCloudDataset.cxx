#include "WireCellUtil/PointCloudDataset.h"
#include "WireCellUtil/Logging.h"

#include <memory>               // make_shared()

using namespace WireCell::PointCloud;

/*
 *  Dataset
 */

static void debug(const Dataset* ds, const std::string& ctx = "")
{
    return;
    // spdlog::debug("Dataset: {}x{} @ {} {}",
    //               ds->size(), ds->size_major(), (void*)ds, ctx);
}

Dataset::Dataset()
{
    debug(this, "default constructor");
}
    
Dataset::Dataset(const std::map<std::string, Array>& arrays)
{
    for (const auto& [n,a] : arrays) {
        add(n,a);
    }
    debug(this, "arrays constructor");
}

Dataset::Dataset(const Dataset& other)
    : m_store(other.m_store)
    , m_append_callbacks(other.m_append_callbacks)
    , m_metadata(other.m_metadata)
{
    debug(this, "copy constructor");
}

Dataset::Dataset(Dataset&& other)
    : m_store(std::move(other.m_store))
    , m_append_callbacks(std::move(other.m_append_callbacks))
    , m_metadata(std::move(other.m_metadata))

{
   debug(this, "move constructor");
}

Dataset Dataset::slice(size_t position, size_t count) const
{
    Dataset dslc;
    for (const auto& name : keys()) {
        dslc.add(name, get(name)->slice(position, count));
    }
    return dslc;
}
Dataset Dataset::slice(size_t position, size_t count, bool share)
{
    Dataset dslc;
    for (const auto& name : keys()) {
        dslc.add(name, get(name)->slice(position, count, share));
    }
    return dslc;
}

Dataset& Dataset::operator=(const Dataset& other)
{
    m_store = other.m_store;
    m_append_callbacks = other.m_append_callbacks;
    m_metadata = other.m_metadata;
    debug(this, "copy assignment");
    return *this;
}

Dataset& Dataset::operator=(Dataset&& other)
{
    std::swap(m_store, other.m_store);
    std::swap(m_append_callbacks, other.m_append_callbacks);
    std::swap(m_metadata, other.m_metadata);
    return *this;
}

Dataset::~Dataset()
{
    debug(this, "destructor");
}

void Dataset::add(const std::string& name, const Array& arr)
{
    if (m_store.size() && size_major() != arr.size_major()) {
        raise<ValueError>("major axis size mismatch when adding \"%s\" to dataset with %d arrays: have=%d got=%d", name, m_store.size(), size_major(), arr.size_major());
    }
    auto [it,inserted] = m_store.insert({name, std::make_shared<Array>(arr)});
    if (inserted) {
        return;
    }
    raise<ValueError>("array named \"%s\" already exists", name);
}
void Dataset::add(const std::string& name, Array&& arr)
{
    if (m_store.size() && size_major() != arr.size_major()) {
        raise<ValueError>("major axis size mismatch when adding \"%s\" to dataset with %d arrays: have=%d got=%d", name, m_store.size(), size_major(), arr.size_major());
    }
    auto [it,inserted] = m_store.insert({name, std::make_shared<Array>(std::move(arr))});
    if (inserted) {
        return;
    }
    raise<ValueError>("array named \"%s\" already exists", name);
}

Dataset::selection_t Dataset::selection(const name_list_t& names) 
{
    selection_t ret;
    for (const auto& name : names) {
        auto it = m_store.find(name);
        if (it == m_store.end()) {
            return selection_t();
        }
        ret.push_back(it->second);
    }
    return ret;
}
Dataset::const_selection_t Dataset::selection(const name_list_t& names) const
{
    const_selection_t ret;
    for (const auto& name : names) {
        auto it = m_store.find(name);
        if (it == m_store.end()) {
            return selection_t();
        }
        ret.push_back(it->second);
    }
    return ret;
}


template<typename MapType>
auto map_keys(const MapType& m)
{
    std::vector<typename MapType::key_type> ret;
    for (const auto& [k,_] : m) {
        ret.push_back(k);
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
        ret.add(key, arr->zeros_like(nmaj));
    }
    return ret;
}

Dataset::keys_t Dataset::keys() const
{
    keys_t ret;
    for (const auto& [k,_] : m_store) {
        ret.push_back(k);
    }
    std::sort(ret.begin(), ret.end());
    return ret;
}
bool Dataset::has(const std::string& key) const
{
    for (const auto& [k,_] : m_store) {
        if (k == key) return true;
    }
    return false;
}
size_t Dataset::size_major() const
{
    if (m_store.empty()) {
        return 0;
    }
    return m_store.begin()->second->size_major();
}

Dataset::array_ptr Dataset::get(const std::string& name) const
{
    auto it = m_store.find(name);
    if (it == m_store.end()) {
        return nullptr;
    }
    return it->second;            
}

Dataset::mutable_array_ptr Dataset::get(const std::string& name) 
{
    auto it = m_store.find(name);
    if (it == m_store.end()) {
        return nullptr;
    }
    return it->second;            
}

void Dataset::append(const Dataset& tail)
{
    if (tail.empty()) {
        return;                 // nothing to append
    }

    const size_t nmaj_before = size_major();

    if (m_store.empty()) {      // we are nothing and become tail
        for (const auto& key : tail.keys()) {
            auto arr = tail.get(key);
            add(key, *arr);
        }
    }
    else {                      // actual append
        auto diff = difference(keys(), tail.keys());
        if (diff.size() > 0) {
            raise<ValueError>("missing keys in append: %d missing", diff.size());
        }

        size_t nmaj = 0;            // check rectangular
        for (const auto& key : tail.keys()) {
            auto arr = tail.get(key);
            if (nmaj == 0) {        // first
                nmaj = arr->size_major();
                continue;
            }
            if (arr->size_major() != nmaj) {
                raise<ValueError>("variation in major axis size in append: can not mix sizes %d and %d",
                                  nmaj, arr->size_major());
            }
        }

        for (auto& [key, arr] : m_store) {
            auto oa = tail.get(key);
            arr->append(oa->bytes());
        }
    }

    const size_t nmaj_after = size_major();

    for (auto cb : m_append_callbacks) {
        cb(nmaj_before, nmaj_after);
    }
}

// fixme: the append() above is a horrible copy-paste-hack from the
// append() below!  They are in danger of drifting apart.  Please fix
// so they share common code.  In haste, I leave this problem.

void Dataset::append(const std::map<std::string, Array>& tail)
{
    if (tail.empty()) {
        return;                 // nothing to append
    }

    const size_t nmaj_before = size_major();

    if (m_store.empty()) {      // we are nothing and become tail
        for (const auto& [key, arr] : tail) {
            add(key, arr);
        }
    }
    else {                      // actual append
        auto diff = difference(keys(), map_keys(tail));
        if (diff.size() > 0) {
            raise<ValueError>("missing keys in append: %d missing", diff.size());
        }

        size_t nmaj = 0;            // check rectangular
        for (auto key : keys()) {
            auto it = tail.find(key);
            if (nmaj == 0) {        // first
                nmaj = it->second.size_major();
                continue;
            }
            const auto& arr = it->second;
            if (arr.size_major() != nmaj) {
                raise<ValueError>("variation in major axis size in append: can not mix sizes %d and %d",
                                  nmaj, arr.size_major());
            }
        }
        for (auto& [key, arr] : m_store) {
            auto it = tail.find(key);
            arr->append(it->second.bytes());
        }
    }

    const size_t nmaj_after = size_major();

    for (auto cb : m_append_callbacks) {
        cb(nmaj_before, nmaj_after);
    }
}

bool Dataset::operator==(Dataset const& other) const
{
    if (m_store.size() != other.m_store.size()) return false;
    for (const auto& [key, arr] : m_store) {
        auto oa = other.get(key);
        if (!oa) return false;
        if (*arr != *oa) {
            return false;
        }
    }
    if (this->metadata() != other.metadata()) {
        return false;
    }
    return true;
}

bool Dataset::operator!=(Dataset const& other) const
{
    return !(*this == other);
}
