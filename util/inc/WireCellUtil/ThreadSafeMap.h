#ifndef WIRECELL_THREADSAFECACHE
#define WIRECELL_THREADSAFECACHE

#include <boost/unordered_map.hpp>

#include <string>
#include <mutex> // For std::unique_lock and std::shared_lock
#include <shared_mutex>
#include <optional>      // For std::optional in C++17


namespace WireCell {

    /// A thread-safe cache with map-like semantics.
    ///
    /// This can grow in size without bounds
    ///
    /// Key type must be Hashable and EqualityComparable.
    /// Value type must be CopyConstructible and CopyAssignable.
    template <typename Key, typename Value>
    class ThreadSafeMap {
    public:
        using shared_mutex = std::shared_mutex;
        using shared_lock = std::shared_lock<shared_mutex>;
        using unique_lock = std::unique_lock<shared_mutex>;

        ThreadSafeMap() = default;

        size_t size() const {
            shared_lock lock(mutex_);
            return cache_.size();
        }

        void clear() {
            unique_lock lock(mutex_);
            cache_.clear();
        }

        void put(const Key& key, const Value& value) {
            unique_lock lock(mutex_); // Exclusive lock for writing
            cache_[key] = value;
        }

        std::optional<Value> get(const Key& key) const {
            shared_lock lock(mutex_); // Shared lock for reading
            auto it = cache_.find(key);
            if (it != cache_.end()) {
                return it->second;
            }
            return std::nullopt;
        }

        void remove(const Key& key) {
            unique_lock lock(mutex_); // Exclusive lock for writing
            cache_.erase(key);
        }

        /**
         * @brief Attempts to retrieve a value for the given key.
         * If the key is found, the existing value is returned as an std::optional.
         * If the key is NOT found, an exclusive lock on the map is acquired
         * and transferred to =out_lock=, and an empty std::optional is returned.
         *
         * The caller is then responsible for computing the value and calling
         * =put_and_release_lock= with the provided =out_lock=.
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

            auto it = cache_.find(key);
            if (it != cache_.end()) {
                Value result = it->second;
                out_lock.unlock(); // Key found, release the lock as we don't need to put.
                return result;
            } else {
                // Key not found. 'out_lock' now holds the exclusive lock.
                // The caller must compute the value and then call put_and_release_lock.
                return std::nullopt;
            }
        }

        /**
         * @brief Puts a value into the map using an already held exclusive lock, then releases the lock.
         * This method is intended to be used in conjunction with get_or_lock().
         * when the latter returns an empty optional (indicating the key was not found
         * and an exclusive lock was acquired).
         *
         * @param key The key to insert or update.
         * @param value The value to associate with the key.
         * @param held_lock An exclusive lock previously acquired by =find_or_lock_for_compute=.
         *                  It will be released by this method.
         * @throws std::logic_error if =held_lock= does not own a lock on the correct mutex.
         */
        void put_and_release(const Key& key, const Value& value, unique_lock& held_lock) {
            // Sanity check: Ensure the caller actually holds the lock.
            // This helps prevent accidental misuse leading to undefined behavior or deadlocks.
            if (!held_lock.owns_lock() || held_lock.mutex() != &mutex_) {
                throw std::logic_error("put_and_release_lock called without owning the correct exclusive lock.");
            }

            // At this point, the exclusive lock is held. We can safely put the value.
            cache_[key] = value;
            held_lock.unlock(); // Release the exclusive lock.
        }

        // Alternative to the get_or_compute()/put_and_release_lock() pair that
        // will call compute_func if needed.
        template <typename F> // F is a callable that returns Value
        Value get_or_compute(const Key& key, F compute_func) {
            unique_lock lock(mutex_); // Exclusive lock for the entire operation

            auto it = cache_.find(key);
            if (it != cache_.end()) {
                return it->second;
            } else {
                // Key not found. Compute while still holding the exclusive lock.
                Value computed_value = compute_func();
                cache_[key] = computed_value;
                return computed_value;
            }
        }

    private:
        boost::unordered_map<Key, Value> cache_;
        mutable shared_mutex mutex_;
    };
}    

#endif
