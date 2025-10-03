#include "WireCellUtil/Logging.h"
#include "WireCellUtil/doctest.h"
#include "WireCellUtil/ThreadSafeCache.h"

// Include necessary Standard Library headers
#include <string>
#include <mutex>         // For std::unique_lock and std::shared_lock
#include <shared_mutex>  // For std::shared_mutex
#include <optional>      // For std::optional in C++17
#include <vector>
#include <thread>
#include <numeric> // For std::iota
#include <atomic>  // For checking counters in concurrent tests
#include <chrono>  // For std::this_thread::sleep_for in some tests
#include <functional> // For std::function
#include <stdexcept> // For std::invalid_argument

// Required headers for ThreadSafeCache internal implementation
#include <list>          // For std::list (doubly linked list for LRU order)
#include <unordered_map> // For std::unordered_map (for key lookup)

// --- Unit Tests ---

using namespace WireCell;

// --- Basic Functionality Tests ---

TEST_CASE("util thread safe cache Constructor and Basic State") {
    SUBCASE("Valid capacity") {
        ThreadSafeCache<int, std::string> cache(5);
        CHECK(cache.size() == 0);
    }

    // SUBCASE("Zero capacity throws exception") {
    //     REQUIRE_THROWS_AS(ThreadSafeCache<int, std::string>(0), std::invalid_argument);
    // }
}

TEST_CASE("util thread safe cache Basic Put, Get, and LRU eviction") {
    ThreadSafeCache<std::string, int> cache(3); // Capacity of 3
    std::optional<int> value_out;

    SUBCASE("Put items below capacity") {
        cache.put("k1", 10);
        cache.put("k2", 20);
        CHECK(cache.size() == 2);
        value_out = cache.get("k1");
        REQUIRE(value_out);
        CHECK(*value_out == 10);
    }

    SUBCASE("Put items exceeding capacity and check eviction") {
        cache.put("k1", 10); // MRU: k1
        cache.put("k2", 20); // MRU: k2, k1
        cache.put("k3", 30); // MRU: k3, k2, k1 (cache full)
        CHECK(cache.size() == 3);

        cache.put("k4", 40); // MRU: k4, k3, k2 (k1 should be evicted)
        CHECK(cache.size() == 3);
        CHECK_FALSE(cache.get("k1")); // k1 should be gone
        value_out = cache.get("k4");
        REQUIRE(value_out);
        CHECK(*value_out == 40);
        value_out = cache.get("k3");
        REQUIRE(value_out);
        CHECK(*value_out == 30);
    }

    SUBCASE("Get updates LRU order") {
        cache.put("k1", 10); // MRU: k1
        cache.put("k2", 20); // MRU: k2, k1
        cache.put("k3", 30); // MRU: k3, k2, k1 (cache full)

        cache.get("k1");     // MRU: k1, k3, k2 (k2 is now LRU)

        cache.put("k4", 40); // MRU: k4, k1, k3 (k2 should be evicted)
        CHECK(cache.size() == 3);
        CHECK_FALSE(cache.get("k2")); // k2 should be gone
        value_out = cache.get("k1");
        REQUIRE(value_out);
        CHECK(*value_out == 10);
    }

    SUBCASE("Update existing item moves to front") {
        cache.put("k1", 10); // MRU: k1
        cache.put("k2", 20); // MRU: k2, k1
        cache.put("k3", 30); // MRU: k3, k2, k1 (cache full)

        cache.put("k1", 100); // Update k1, it becomes MRU. MRU: k1, k3, k2
        CHECK(cache.size() == 3);

        cache.put("k4", 40); // MRU: k4, k1, k3 (k2 should be evicted)
        CHECK(cache.size() == 3);
        CHECK_FALSE(cache.get("k2"));
        value_out = cache.get("k1");
        REQUIRE(value_out);
        CHECK(*value_out == 100);
    }
}

TEST_CASE("util thread safe cache Remove and Clear functionality") {
    ThreadSafeCache<int, std::string> cache(5);
    cache.put(1, "one");
    cache.put(2, "two");
    cache.put(3, "three");
    CHECK(cache.size() == 3);

    SUBCASE("Remove an existing item") {
        cache.remove(2);
        CHECK(cache.size() == 2);
        CHECK_FALSE(cache.get(2));
        REQUIRE(cache.get(1));
        REQUIRE(cache.get(3));
    }

    SUBCASE("Remove a non-existent item") {
        cache.remove(99); // Should have no effect
        CHECK(cache.size() == 3);
    }

    SUBCASE("Clear all items") {
        cache.clear();
        CHECK(cache.size() == 0);
        CHECK_FALSE(cache.get(1));
        CHECK_FALSE(cache.get(2));
    }
}

// --- Concurrent Access Tests ---

TEST_CASE("util thread safe cache Concurrent Puts with LRU Eviction") {
    const int capacity = 100;
    ThreadSafeCache<int, int> cache(capacity);
    const int num_threads = 4;
    const int puts_per_thread = 500; // Total puts: 2000, capacity 100

    std::vector<std::thread> threads;
    for (int i = 0; i < num_threads; ++i) {
        threads.emplace_back([&cache, i, puts_per_thread]() {
            for (int j = 0; j < puts_per_thread; ++j) {
                // Keys will be from 0 to num_threads*puts_per_thread - 1
                int key = i * puts_per_thread + j;
                cache.put(key, key * 2);
            }
        });
    }

    for (auto& t : threads) {
        t.join();
    }

    // Final size must be equal to capacity
    CHECK(cache.size() == capacity);

    // It's non-deterministic which exact 100 items remain,
    // but we can at least check that they are valid items that *could* be in the cache.
    // This part is tricky. A more robust test would be to test with a small capacity
    // and specific patterns. For large N, we just verify size and no crashes.
    for (int k = 0; k < capacity; ++k) {
        // Try to get some items (non-deterministic which ones will be hit)
        std::optional<int> val_opt = cache.get(k);
        if (val_opt) {
            // If it's there, check value consistency
            CHECK(*val_opt == k * 2);
        }
    }
}

TEST_CASE("util thread safe cache Concurrent Gets and LRU Update (Unique_lock impact)") {
    const int capacity = 500;
    ThreadSafeCache<int, int> cache(capacity);
    for (int i = 0; i < capacity; ++i) {
        cache.put(i, i * 10);
    }
    CHECK(cache.size() == capacity);

    const int num_threads = 8;
    const int gets_per_thread = 2000;
    std::vector<std::thread> threads;
    std::atomic<int> successful_reads(0);

    for (int i = 0; i < num_threads; ++i) {
        threads.emplace_back([&cache, &successful_reads, capacity, gets_per_thread]() {
            for (int j = 0; j < gets_per_thread; ++j) {
                int key = j % capacity; // Cycle through existing keys
                std::optional<int> value_opt = cache.get(key);
                if (value_opt) {
                    if (*value_opt == key * 10) {
                        successful_reads++;
                    }
                }
            }
        });
    }

    for (auto& t : threads) {
        t.join();
    }

    // All gets should succeed and retrieve correct values, even with unique_lock contention
    CHECK(successful_reads == num_threads * gets_per_thread);
    CHECK(cache.size() == capacity); // Size should remain constant
}

TEST_CASE("util thread safe cache Concurrent get_or_lock / put_and_release") {
    const size_t capacity = 10;
    ThreadSafeCache<int, double> cache(capacity);
    const int num_keys = 20; // More keys than capacity to test eviction
    const int num_threads = 4;
    const int access_per_thread = 100;

    std::atomic<int> expensive_function_call_count(0);

    auto call_expensive_function = [&](int key_val) -> double {
        expensive_function_call_count++;
        std::this_thread::sleep_for(std::chrono::microseconds(100)); // Simulate work
        return static_cast<double>(key_val * 10.0 + 0.5);
    };

    std::vector<std::thread> threads;
    for (int i = 0; i < num_threads; ++i) {
        threads.emplace_back([&]() {
            for (int k_idx = 0; k_idx < access_per_thread; ++k_idx) {
                int key = k_idx % num_keys; // Cycle through keys, some will exceed capacity

                ThreadSafeCache<int, double>::unique_lock lock_handle;
                std::optional<double> found_value = cache.get_or_lock(key, lock_handle);

                if (found_value) {
                    // Value was already in the cache, no computation needed
                    CHECK(*found_value == static_cast<double>(key * 10.0 + 0.5));
                } else {
                    // Key not found, lock_handle holds the exclusive lock.
                    double computed_value = call_expensive_function(key);
                    cache.put_and_release(key, computed_value, lock_handle);
                }
            }
        });
    }

    for (auto& t : threads) {
        t.join();
    }

    // Final size must be equal to capacity
    CHECK(cache.size() == capacity);

    // This sub-test checks that expensive_function_call_count is exactly once per key
    // when keys accessed are within capacity, guaranteeing no redundant computations.
    const int test_num_keys_for_count = 10;
    ThreadSafeCache<int, double> cache_for_count(test_num_keys_for_count);
    expensive_function_call_count = 0; // Reset counter for new sub-test

    std::vector<std::thread> threads_for_count;
    for (int i = 0; i < num_threads; ++i) {
        threads_for_count.emplace_back([&]() {
            for (int k_idx = 0; k_idx < access_per_thread; ++k_idx) {
                int key = k_idx % test_num_keys_for_count; // Access keys within capacity
                ThreadSafeCache<int, double>::unique_lock lock_handle;
                std::optional<double> found_value = cache_for_count.get_or_lock(key, lock_handle);
                if (!found_value) {
                    double computed_value = call_expensive_function(key);
                    cache_for_count.put_and_release(key, computed_value, lock_handle);
                }
            }
        });
    }
    for (auto& t : threads_for_count) {
        t.join();
    }
    CHECK(cache_for_count.size() == test_num_keys_for_count);
    CHECK(expensive_function_call_count == test_num_keys_for_count); // Should be called exactly once per key

    for (int key = 0; key < test_num_keys_for_count; ++key) {
        std::optional<double> value_out = cache_for_count.get(key);
        REQUIRE(value_out);
        CHECK(*value_out == static_cast<double>(key * 10.0 + 0.5));
    }
}


TEST_CASE("util thread safe cache Concurrent get_or_compute") {
    const size_t capacity = 10;
    ThreadSafeCache<int, double> cache(capacity);
    const int num_keys = 20; // More keys than capacity to test eviction
    const int num_threads = 4;
    const int access_per_thread = 100;

    std::atomic<int> expensive_function_call_count(0);

    auto compute_value_func = [&](int key_val) -> double {
        expensive_function_call_count++;
        std::this_thread::sleep_for(std::chrono::microseconds(100)); // Simulate work
        return static_cast<double>(key_val * 100.0 + 0.123);
    };

    std::vector<std::thread> threads;
    for (int i = 0; i < num_threads; ++i) {
        threads.emplace_back([&]() {
            for (int k_idx = 0; k_idx < access_per_thread; ++k_idx) {
                int key = k_idx % num_keys; // Cycle through keys, some will exceed capacity

                double retrieved_value = cache.get_or_compute(key, [&]() {
                    return compute_value_func(key);
                });

                CHECK(retrieved_value == static_cast<double>(key * 100.0 + 0.123));
            }
        });
    }

    for (auto& t : threads) {
        t.join();
    }

    // Final size must be equal to capacity
    CHECK(cache.size() == capacity);

    // This sub-test checks that expensive_function_call_count is exactly once per key
    // when keys accessed are within capacity, guaranteeing no redundant computations.
    const int test_num_keys_for_count = 10;
    ThreadSafeCache<int, double> cache_for_count(test_num_keys_for_count);
    expensive_function_call_count = 0; // Reset counter for new sub-test

    std::vector<std::thread> threads_for_count;
    for (int i = 0; i < num_threads; ++i) {
        threads_for_count.emplace_back([&]() {
            for (int k_idx = 0; k_idx < access_per_thread; ++k_idx) {
                int key = k_idx % test_num_keys_for_count; // Access keys within capacity
                cache_for_count.get_or_compute(key, [&]() {
                    return compute_value_func(key);
                });
            }
        });
    }
    for (auto& t : threads_for_count) {
        t.join();
    }
    CHECK(cache_for_count.size() == test_num_keys_for_count);
    CHECK(expensive_function_call_count == test_num_keys_for_count); // Should be called exactly once per key

    for (int key = 0; key < test_num_keys_for_count; ++key) {
        std::optional<double> value_out = cache_for_count.get(key);
        REQUIRE(value_out);
        CHECK(*value_out == static_cast<double>(key * 100.0 + 0.123));
    }
}


TEST_CASE("util thread safe cache Mixed Concurrent Operations (Put, Get, Remove, Compute)") {
    const size_t capacity = 20;
    ThreadSafeCache<int, std::string> cache(capacity);
    const int num_threads = 10;
    const int operations_per_thread = 500;
    const int total_keys_range = 100; // Keys will be from 0 to 99

    std::atomic<int> expensive_compute_count(0);

    auto compute_string_func = [&](int key_val) -> std::string {
        expensive_compute_count++;
        // Simulate real work
        std::this_thread::sleep_for(std::chrono::microseconds(50));
        return "computed_val_" + std::to_string(key_val);
    };

    std::vector<std::thread> threads;
    for (int i = 0; i < num_threads; ++i) {
        threads.emplace_back([&]() {
            for (int j = 0; j < operations_per_thread; ++j) {
                int key = j % total_keys_range; // Cycle through keys

                int choice = (j + i) % 4; // Vary operations more

                if (choice == 0) { // Put
                    cache.put(key, "put_val_" + std::to_string(key));
                } else if (choice == 1) { // Get
                    std::optional<std::string> val = cache.get(key);
                    // Can't reliably check value as it might be overwritten concurrently
                    // But ensures no crash.
                } else if (choice == 2) { // Remove
                    cache.remove(key);
                } else { // get_or_compute
                    std::string val = cache.get_or_compute(key, [&]() {
                        return compute_string_func(key);
                    });
                    // Similarly, can't reliably check against a fixed expected value due to concurrent puts/removes
                }
                // A small sleep to increase context switching and contention
                if (j % 50 == 0) {
                    std::this_thread::sleep_for(std::chrono::microseconds(1));
                }
            }
        });
    }

    for (auto& t : threads) {
        t.join();
    }

    // After all operations, verify no crashes and a reasonable final state.
    // The exact final state (which items are present, which are LRU, exact size)
    // is highly non-deterministic with mixed operations.
    CHECK(cache.size() <= capacity); // Max size should not exceed capacity
    CHECK(cache.size() >= 0);        // Sanity check

    // Check one fixed key for stability after everything
    cache.put(999, "final_test_value");
    std::optional<std::string> final_val_opt = cache.get(999);
    REQUIRE(final_val_opt);
    CHECK(*final_val_opt == "final_test_value");
    CHECK(expensive_compute_count >= 0); // At least some computations might have happened
}
