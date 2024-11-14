#if defined __clang__
 #pragma clang diagnostic push
 #pragma clang diagnostic ignored "-Wunused-local-typedefs"
#endif

#include "boost/icl/interval_map.hpp"
#include "boost/icl/split_interval_map.hpp"
#include "boost/icl/separate_interval_set.hpp"
#include "boost/icl/split_interval_set.hpp"

#if defined __clang__
 #pragma clang diagnostic pop
#endif

#include <iostream>
#include <vector>
#include <typeinfo>

template<typename ISet>
void doit_set()
{
    using key_t = typename ISet::domain_type;

    using interval_t = typename boost::icl::interval<key_t>::interval_type;
    // using interval_set_t = typename boost::icl::interval_set<Key>;

    const std::vector<std::pair<key_t, key_t>> intervals = {
        {1,4}, {2,5}, {5,7}
    };
    ISet is;
    for (const auto& [f,s] : intervals) {
        is += interval_t::right_open(f,s);
    }
    std::cerr << is << std::endl;
}

template<typename IMap>
void doit_map(bool common)
{
    using key_t = typename IMap::domain_type;
    using set_t = typename IMap::codomain_type;
    using interval_t = typename IMap::interval_type;
    // using interval_map_t = typename boost::icl::interval_map<key_t, set_t>;

    const std::vector<std::pair<key_t, key_t>> intervals = {
        {1,4}, {2,5}, {5,7}
    };
    IMap im;
    size_t count{1};
    for (const auto& [f,s] : intervals) {
        im += std::make_pair(interval_t::right_open(f,s), set_t{count});
        ++count;
        if (common) --count;
    }
    // this causes compilation errors with gcc 12 and clang 14.
    // std::cerr << im << std::endl;
}

template<typename Key>
void doit_key()
{
    std::cerr << "key type: " << typeid(Key).name() << "\n";

    std::cerr << "interval_set joining:\n";
    doit_set<boost::icl::interval_set<Key>>();
    std::cerr << "interval_set separating:\n";
    doit_set<boost::icl::separate_interval_set<Key>>();
    std::cerr << "interval_set splitting:\n";
    doit_set<boost::icl::split_interval_set<Key>>();

    std::cerr << "interval_map joining same values:\n";
    doit_map<boost::icl::interval_map<Key, std::set<size_t>>>(true);
    std::cerr << "interval_map splitting same values:\n";
    doit_map<boost::icl::split_interval_map<Key, std::set<size_t>>>(true);

    std::cerr << "interval_map joining different values:\n";
    doit_map<boost::icl::interval_map<Key, std::set<size_t>>>(false);
    std::cerr << "interval_map splitting different values:\n";
    doit_map<boost::icl::split_interval_map<Key, std::set<size_t>>>(false);

    std::cerr << "interval_map joining same values sum:\n";
    doit_map<boost::icl::interval_map<Key, size_t>>(true);
    std::cerr << "interval_map splitting same values sum:\n";
    doit_map<boost::icl::split_interval_map<Key, size_t>>(true);

    std::cerr << "interval_map joining different values sum:\n";
    doit_map<boost::icl::interval_map<Key, size_t>>(false);
    std::cerr << "interval_map splitting different values sum:\n";
    doit_map<boost::icl::split_interval_map<Key, size_t>>(false);

    

}

int main()
{
    doit_key<int>();
    // doit_key<float>();
}
