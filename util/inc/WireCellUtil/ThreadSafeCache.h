#ifndef WIRECELL_THREADSAFECACHE
#define WIRECELL_THREADSAFECACHE

#include <string>
#include <mutex>         // For std::unique_lock and std::shared_lock
#include <shared_mutex>  // For std::shared_mutex
#include <optional>      // For std::optional in C++17
#include <list>          // For std::list (doubly linked list for LRU order)
#include <unordered_map> // For std::unordered_map (for key lookup)
#include <stdexcept>     // For std::invalid_argument

namespace WireCell {

    /// A thread-safe LRU cache with map-like semantics.
    ///
    /// This cache has a fixed capacity and evicts the Least Recently Used (LRU)
    /// item when capacity is reached.
    ///
    /// Key type must be Hashable and EqualityComparable.
    /// Value type must be CopyConstructible and CopyAssignable.
    template <typename Key, typename Value>
    class ThreadSafeCache {
    public:
        // Type aliases for convenience
        using shared_mutex = std::shared_mutex;
        using shared_lock = std::shared_lock<shared_mutex>;
        using unique_lock = std::unique_lock<shared_mutex>;

        // Internal type for list nodes
        using list_node_type = std::pair<Key, Value>;
        using list_iterator = typename std::list<list_node_type>::iterator;
        using map_type = std::unordered_map<Key, list_iterator>;
        using list_type = std::list<list_node_type>;

        /**
         * @brief Constructs a ThreadSafeCache with a specified capacity.
         * @param capacity The maximum number of items the cache can hold. Must be > 0.
         * @throws std::invalid_argument if capacity is 0.
         */
        explicit ThreadSafeCache(size_t capacity)
            : capacity_(capacity) {
            if (capacity == 0) {
                throw std::invalid_argument("ThreadSafeCache capacity must be greater than 0.");
            }
        }

        void reset_capacity(size_t capacity) {
            capacity_ = capacity;
            shrink(capacity);
        }

        // The size() method is const and only reads cache_map_ size, so it can use shared_lock.
        size_t size() const {
            shared_lock lock(mutex_);
            return cache_map_.size();
        }

        // Methods that modify the cache state require a unique_lock.

        void clear() {
            unique_lock lock(mutex_);
            cache_map_.clear();
            cache_list_.clear();
        }

        std::optional<Value> get(const Key& key) {
            unique_lock lock(mutex_); // Unique lock required as get modifies LRU order

            auto it_map = cache_map_.find(key);
            if (it_map == cache_map_.end()) {
                return std::nullopt; // Key not found
            }

            // Key found, update LRU order: move the item to the front of the list
            list_iterator it_list = it_map->second;
            cache_list_.splice(cache_list_.begin(), cache_list_, it_list); // Move to front

            return it_list->second; // Return the value
        }

        /// Shrink cache to given size.  If size is zero, shrink to capacity - 1 to allow room for one put().
        void shrink(size_t size = 0) {
            if (!size) size = capacity_ - 1;
            while (cache_map_.size() > size) {
                Key lru_key = cache_list_.back().first;
                cache_map_.erase(lru_key);
                cache_list_.pop_back();
            }
        }

        void put(const Key& key, const Value& value) {
            unique_lock lock(mutex_); // Exclusive lock for writing

            auto it_map = cache_map_.find(key);
            if (it_map != cache_map_.end()) {
                // Key exists: update value and move to front (MRU)
                list_iterator it_list = it_map->second;
                it_list->second = value; // Update value
                cache_list_.splice(cache_list_.begin(), cache_list_, it_list); // Move to front
            } else {
                // Key does not exist: add new item

                shrink();

                // Insert new item at the front (MRU)
                cache_list_.emplace_front(key, value);
                cache_map_[key] = cache_list_.begin(); // Map key to new list iterator
            }
        }

        void remove(const Key& key) {
            unique_lock lock(mutex_); // Exclusive lock for writing

            auto it_map = cache_map_.find(key);
            if (it_map != cache_map_.end()) {
                list_iterator it_list = it_map->second;
                cache_list_.erase(it_list); // Erase from list
                cache_map_.erase(it_map);   // Erase from map
            }
        }

        /**
         * @brief Attempts to retrieve a value for the given key.
         * If the key is found, the existing value is returned as an std::optional.
         * If the key is NOT found, an exclusive lock on the map is acquired
         * and transferred to =out_lock=, and an empty std::optional is returned.
         *
         * The caller is then responsible for computing the value and calling
         * =put_and_release= with the provided =out_lock=.
         *
         * This method ensures that only one thread will perform the computation
         * for a given missing key by holding an exclusive lock.
         *
         * @param key The key to look for.
         * @param out_lock An empty unique_lock reference to receive the exclusive lock
         *                 if the key is not found.
         * @return An std::optional<Value> containing the value if found, or empty if not found.
         */
        std::optional<Value> get_or_lock(const Key& key, unique_lock& out_lock) {
            // Acquire an exclusive lock immediately.
            // This ensures that if the key is not present, this thread
            // prevents any other thread from trying to compute/put the same key
            // or modify the map in a conflicting way.
            out_lock = unique_lock(mutex_);

            auto it_map = cache_map_.find(key);
            if (it_map != cache_map_.end()) {
                // Key found, update LRU order: move the item to the front of the list
                list_iterator it_list = it_map->second;
                cache_list_.splice(cache_list_.begin(), cache_list_, it_list); // Move to front

                Value result = it_list->second;
                out_lock.unlock(); // Key found, release the lock as we don't need to put.
                return result;
            } else {
                // Key not found. 'out_lock' now holds the exclusive lock.
                // The caller must compute the value and then call put_and_release.
                return std::nullopt;
            }
        }

        /**
         * @brief Puts a value into the cache using an already held exclusive lock, then releases the lock.
         * This method is intended to be used in conjunction with get_or_lock()
         * when the latter returns an empty optional (indicating the key was not found
         * and an exclusive lock was acquired).
         *
         * @param key The key to insert or update.
         * @param value The value to associate with the key.
         * @param held_lock An exclusive lock previously acquired by get_or_lock().
         *                  It will be released by this method.
         * @throws std::logic_error if =held_lock= does not own a lock on the correct mutex.
         */
        void put_and_release(const Key& key, const Value& value, unique_lock& held_lock) {
            // Sanity check: Ensure the caller actually holds the lock.
            if (!held_lock.owns_lock() || held_lock.mutex() != &mutex_) {
                throw std::logic_error("put_and_release called without owning the correct exclusive lock.");
            }

            // At this point, the exclusive lock is held. We can safely put the value,
            // managing LRU eviction if needed.
            auto it_map = cache_map_.find(key);
            if (it_map != cache_map_.end()) {
                // Key exists: update value and move to front (MRU)
                list_iterator it_list = it_map->second;
                it_list->second = value; // Update value
                cache_list_.splice(cache_list_.begin(), cache_list_, it_list); // Move to front
            } else {
                // Key does not exist: add new item

                shrink();

                // Insert new item at the front (MRU)
                cache_list_.emplace_front(key, value);
                cache_map_[key] = cache_list_.begin(); // Map key to new list iterator
            }

            held_lock.unlock(); // Release the exclusive lock.
        }

        /**
         * @brief Atomically gets a value for the given key, or computes and inserts it if not found.
         *
         * If the key is found, its value is returned, and its LRU status is updated.
         * If the key is not found, the provided computation function is called (while holding an
         * exclusive lock), its result is inserted into the cache (evicting LRU if needed),
         * and then the computed value is returned.
         *
         * This entire operation is atomic with respect to other cache operations.
         *
         * @param key The key to look for or compute.
         * @param compute_func A callable (e.g., lambda, std::function) that takes no arguments
         *                     and returns a Value. This function is only called if the key
         *                     is not found in the cache.
         * @return The value associated with the key (either retrieved or computed).
         */
        template <typename F> // F is a callable that returns Value
        Value get_or_compute(const Key& key, F compute_func) {
            unique_lock lock(mutex_); // Exclusive lock for the entire operation

            auto it_map = cache_map_.find(key);
            if (it_map != cache_map_.end()) {
                // Key found: update LRU order and return value
                list_iterator it_list = it_map->second;
                cache_list_.splice(cache_list_.begin(), cache_list_, it_list); // Move to front
                return it_list->second;
            } else {
                // Key not found. Compute while still holding the exclusive lock.
                Value computed_value = compute_func();

                shrink();

                // Insert new item at the front (MRU)
                cache_list_.emplace_front(key, computed_value);
                cache_map_[key] = cache_list_.begin(); // Map key to new list iterator
                return computed_value;
            }
        }

    private:
        map_type cache_map_;  // Maps Key to iterator in cache_list_
        list_type cache_list_; // Stores (Key, Value) pairs in LRU order (front=MRU, back=LRU)
        size_t capacity_;      // Maximum number of items in the cache
        mutable shared_mutex mutex_; // Must be mutable for const size() method
    };
} // namespace WireCell

#endif
