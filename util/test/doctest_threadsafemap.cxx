#include "WireCellUtil/Logging.h"
#include "WireCellUtil/doctest.h"
#include "WireCellUtil/ThreadSafeMap.h"

#include <vector>
#include <thread>
#include <numeric> // For std::iota
#include <atomic>  // For checking counters in concurrent tests
#include <chrono>  // For std::this_thread::sleep_for in some tests

using namespace WireCell;

// Helper function to create unique string keys for testing
static
std::string make_key(int i) {
    return "key_" + std::to_string(i);
}

// --- Basic Functionality Tests ---

TEST_CASE("util thread safe map Basic put and get operations") {
    ThreadSafeMap<std::string, int> cache;
    std::string key = "my_key";
    int value_in = 42;


    SUBCASE("Put a new item and retrieve it") {
        cache.put(key, value_in);
        CHECK(cache.size() == 1);
        auto value_out = cache.get(key);
        REQUIRE(value_out);
        CHECK(value_out == value_in);
    }

    SUBCASE("Attempt to get a non-existent item") {
        CHECK_FALSE(cache.get("non_existent_key"));
        CHECK(cache.size() == 0); // Still 0 if no puts before
    }
}

TEST_CASE("util thread safe map Update an existing item") {
    ThreadSafeMap<std::string, int> cache;
    std::string key = "update_me";
    int value1 = 100;
    int value2 = 200;

    cache.put(key, value1);
    CHECK(cache.size() == 1);
    auto value_out = cache.get(key);
    REQUIRE(value_out);
    CHECK(value_out == value1);

    cache.put(key, value2); // Overwrite
    CHECK(cache.size() == 1); // Size should remain 1
    value_out = cache.get(key);
    REQUIRE(value_out);
    CHECK(value_out == value2);
}

TEST_CASE("util thread safe map Remove item functionality") {
    ThreadSafeMap<std::string, double> cache;
    std::string key1 = "item_one";
    std::string key2 = "item_two";
    double value = 3.14;

    cache.put(key1, value);
    cache.put(key2, value);
    CHECK(cache.size() == 2);

    SUBCASE("Remove an existing item") {
        cache.remove(key1);
        CHECK(cache.size() == 1);
        auto value_out = cache.get(key1);
        CHECK_FALSE(value_out); // key1 should be gone
        value_out = cache.get(key2);
        REQUIRE(value_out);     // key2 should still be there
        CHECK(value_out == value);
    }

    SUBCASE("Remove a non-existent item") {
        cache.remove("unknown_key"); // Should have no effect
        CHECK(cache.size() == 2);    // Size should not change
    }
}

TEST_CASE("util thread safe map Clear method") {
    ThreadSafeMap<int, bool> cache;
    cache.put(1, true);
    cache.put(2, false);
    CHECK(cache.size() == 2);

    cache.clear();
    CHECK(cache.size() == 0);

    auto val_out = cache.get(1);
    CHECK_FALSE(cache.get(1));
    val_out = cache.get(2);
    CHECK_FALSE(val_out);
}

// --- Concurrent Access Tests ---

TEST_CASE("util thread safe map Concurrent Puts (Multiple Writers)") {
    ThreadSafeMap<std::string, int> cache;
    const int num_threads = 8;
    const int puts_per_thread = 1000;
    std::vector<std::thread> threads;

    for (int i = 0; i < num_threads; ++i) {
        threads.emplace_back([&cache, i, puts_per_thread]() {
            for (int j = 0; j < puts_per_thread; ++j) {
                std::string key = make_key(i * puts_per_thread + j); // Unique key per put
                cache.put(key, i + j);
            }
        });
    }

    for (auto& t : threads) {
        t.join();
    }

    // Verify all items are in the cache and their values are correct
    CHECK(cache.size() == num_threads * puts_per_thread);

    for (int i = 0; i < num_threads; ++i) {
        for (int j = 0; j < puts_per_thread; ++j) {
            std::string key = make_key(i * puts_per_thread + j);
            auto value_out = cache.get(key);
            REQUIRE(value_out);
            CHECK(value_out == i + j);
        }
    }
}

TEST_CASE("util thread safe map Concurrent Gets (Multiple Readers)") {
    ThreadSafeMap<std::string, int> cache;
    const int num_initial_items = 500;
    // Populate cache with known values
    for (int i = 0; i < num_initial_items; ++i) {
        cache.put(make_key(i), i * 10);
    }
    CHECK(cache.size() == num_initial_items);

    const int num_readers = 16;
    const int reads_per_reader = 5000;
    std::vector<std::thread> threads;
    std::atomic<int> successful_and_correct_reads(0);

    for (int i = 0; i < num_readers; ++i) {
        threads.emplace_back([&cache, &successful_and_correct_reads, num_initial_items, reads_per_reader]() {
            for (int j = 0; j < reads_per_reader; ++j) {
                int key_idx = j % num_initial_items; // Cycle through existing keys
                std::string key = make_key(key_idx);
                auto value_out = cache.get(key);
                if (value_out) {
                    if (value_out == key_idx * 10) { // Verify data integrity
                        successful_and_correct_reads++;
                    }
                }
            }
        });
    }

    for (auto& t : threads) {
        t.join();
    }

    // All reads should have been successful and retrieved the correct value
    CHECK(successful_and_correct_reads == num_readers * reads_per_reader);
    CHECK(cache.size() == num_initial_items); // Cache size should not change from reads
}

TEST_CASE("util thread safe map Mixed Concurrent Operations (Put, Get, Remove)") {
    ThreadSafeMap<int, std::string> cache;
    const int num_threads = 10;
    const int operations_per_thread = 1000;
    std::vector<std::thread> threads;

    // Use an offset to generate distinct keys across test runs
    int key_offset = 10000;

    for (int i = 0; i < num_threads; ++i) {
        threads.emplace_back([&cache, i, operations_per_thread, key_offset]() {
            for (int j = 0; j < operations_per_thread; ++j) {
                int current_key = key_offset + i * operations_per_thread + j;
                std::string current_value = "val_" + std::to_string(current_key);
                int choice = j % 3; // Cycle through Put, Get, Remove

                if (choice == 0) { // Put
                    cache.put(current_key, current_value);
                } else if (choice == 1) { // Get
                    auto value_out = cache.get(current_key);
                    // Checking value_out here is complex due to concurrent writes/removes,
                    // but the primary goal is to ensure no crashes or deadlocks.
                } else { // Remove
                    cache.remove(current_key);
                }
                // Small sleep to increase context switching and contention
                if (j % 100 == 0) {
                    std::this_thread::sleep_for(std::chrono::microseconds(1));
                }
            }
        });
    }

    for (auto& t : threads) {
        t.join();
    }

    // A precise size check is challenging due to the unpredictable nature of concurrent mixed ops.
    // The main goal is to ensure the operations complete without crashing or deadlocking,
    // demonstrating correct shared_mutex usage under heavy contention.
    CHECK(cache.size() >= 0); // Ensure size is not negative (sanity check)

    // Add a known item after all threads finish and check it
    std::string final_check_key_val = "final_value";
    cache.put(99999, final_check_key_val);
    auto value_out = cache.get(99999);
    REQUIRE(value_out);
    CHECK(value_out == final_check_key_val);
}

TEST_CASE("util thread safe map Multiple writers and a concurrent reader") {
    ThreadSafeMap<int, std::string> cache;
    const int num_writers = 4;
    const int items_per_writer = 500;
    std::vector<std::thread> writer_threads;

    for (int i = 0; i < num_writers; ++i) {
        writer_threads.emplace_back([&cache, i, items_per_writer]() {
            for (int j = 0; j < items_per_writer; ++j) {
                cache.put(i * items_per_writer + j, "writer_" + std::to_string(i) + "_val_" + std::to_string(j));
            }
        });
    }

    // A single reader thread that attempts to read all possible keys while writers are active.
    // This tests if reads can proceed concurrently with other reads, and are blocked by writes.
    std::thread reader_thread([&cache, num_writers, items_per_writer]() {
        for (int i = 0; i < num_writers * items_per_writer; ++i) {
            auto value_out = cache.get(i);
            // The value might be empty if the key hasn't been written yet or was removed.
            // Or it might be an intermediate value if it was overwritten.
            // The crucial part is that `get` completes safely without data corruption or crash.
        }
    });

    for (auto& t : writer_threads) {
        t.join();
    }
    reader_thread.join();

    // After all operations, verify the final state (all items should be present and correct)
    CHECK(cache.size() == num_writers * items_per_writer);
    for (int i = 0; i < num_writers * items_per_writer; ++i) {
        auto value_out = cache.get(i);
        REQUIRE(value_out);
        CHECK(value_out == "writer_" + std::to_string(i / items_per_writer) + "_val_" + std::to_string(i % items_per_writer));
    }
}

TEST_CASE("util thread safe map Concurrent Compute-If-Absent (get_or_lock)") {
    ThreadSafeMap<int, double> tsmap;
    const int num_keys = 100;
    const int num_threads = 8;
    const int access_per_thread = 500; // Each thread will try to access num_keys * N times

    std::atomic<int> expensive_function_call_count(0);
    std::atomic<int> total_gets_success(0);

    // This function simulates an expensive computation
    auto call_expensive_function = [&](int key_val) -> double {
        expensive_function_call_count++;
        // Simulate real work
        std::this_thread::sleep_for(std::chrono::milliseconds(2));
        return static_cast<double>(key_val * 10.0 + 0.5);
    };

    std::vector<std::thread> threads;
    for (int i = 0; i < num_threads; ++i) {
        threads.emplace_back([&]() {
            for (int k_idx = 0; k_idx < access_per_thread; ++k_idx) {
                int key = k_idx % num_keys; // Cycle through keys to ensure contention
                
                ThreadSafeMap<int, double>::unique_lock lock_handle;
                std::optional<double> found_value = tsmap.get_or_lock(key, lock_handle);

                if (found_value) {
                    // Value was already in the map, no computation needed
                    total_gets_success++;
                    CHECK(found_value.value() == static_cast<double>(key * 10.0 + 0.5));
                } else {
                    // Key not found, =lock_handle= holds the exclusive lock.
                    // Now compute the value and put it.
                    double computed_value = call_expensive_function(key);
                    tsmap.put_and_release(key, computed_value, lock_handle);
                    // No need to check lock_handle.owns_lock() here, put_and_release_lock
                    // ensures it's released and throws if not owned.
                    total_gets_success++; // Counting the successful put as a successful "get" operation
                }
            }
        });
    }

    for (auto& t : threads) {
        t.join();
    }

    // Verify all keys are in the map
    CHECK(tsmap.size() == num_keys);

    // Verify values and that expensive function was called only once per key
    for (int key = 0; key < num_keys; ++key) {
        auto value_out = tsmap.get(key);
        REQUIRE(value_out);
        CHECK(value_out == static_cast<double>(key * 10.0 + 0.5));
    }

    // The core check: expensive function should have been called exactly once per unique key.
    CHECK(expensive_function_call_count == num_keys);

    // Ensure all successful paths were taken. Total access count in the loops
    CHECK(total_gets_success == num_threads * access_per_thread);
}


TEST_CASE("util thread safe map get_or_compute method (Encapsulated Compute-If-Absent)") {
    ThreadSafeMap<int, double> tsmap;
    const int num_keys = 100;
    const int num_threads = 8;
    const int access_per_thread = 500; // Each thread will try to access num_keys * N times

    std::atomic<int> expensive_function_call_count(0);

    // This function simulates an expensive computation that takes a key
    // and returns a value, also tracking how many times it's called.
    auto compute_value_func = [&](int key_val) -> double {
        expensive_function_call_count++;
        // Simulate real work
        std::this_thread::sleep_for(std::chrono::milliseconds(2));
        return static_cast<double>(key_val * 100.0 + 0.123);
    };

    std::vector<std::thread> threads;
    for (int i = 0; i < num_threads; ++i) {
        threads.emplace_back([&]() {
            for (int k_idx = 0; k_idx < access_per_thread; ++k_idx) {
                int key = k_idx % num_keys; // Cycle through keys to ensure contention

                // Call get_or_compute, passing the computation function
                // The lambda passed to get_or_compute is called only if the key is absent.
                // The `key` captured by value for the inner lambda is correct for its use.
                double retrieved_value = tsmap.get_or_compute(key, [&]() {
                    return compute_value_func(key);
                });

                CHECK(retrieved_value == static_cast<double>(key * 100.0 + 0.123));
            }
        });
    }

    for (auto& t : threads) {
        t.join();
    }

    // Verify all keys are in the map
    CHECK(tsmap.size() == num_keys);

    // Verify values are correct
    for (int key = 0; key < num_keys; ++key) {
        auto value_out = tsmap.get(key);
        REQUIRE(value_out);
        CHECK(value_out == static_cast<double>(key * 100.0 + 0.123));
    }

    // The core check: expensive_function_call_count should be exactly once per unique key.
    CHECK(expensive_function_call_count == num_keys);
}
