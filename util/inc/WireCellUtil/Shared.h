#ifndef WIRECELL_SHARED
#define WIRECELL_SHARED

#include <memory>
#include <vector>
namespace WireCell {

    // Provide base class defining things for objects that are expected to be
    // held by shared pointer.  Use like:
    //
    // class Thing : public Shared<Thing> {
    //     // ... etc ...
    // };
    template<typename T>
    struct Shared {

        // The wrapped sub class type.
        using shared_type = T;

        // For holding a facade by a shared pointer.
        using pointer = std::shared_ptr<T>;
        using const_pointer = std::shared_ptr<const T>;

        // Simple collection of shared facades.
        using vector = std::vector<pointer>;
        using const_vector = std::vector<const_pointer>;
    };
    
}
#endif
