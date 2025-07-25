#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Flagged.h"

#include "WireCellUtil/doctest.h"

#include <iostream>

using namespace WireCell;

TEST_CASE("flagged with enums") {
    enum E {a=1<<1, b=1<<2, c=1<<3, a_and_b=a|b, b_and_c=b|c};
    Flagged f;

    f.set_flags(E::a);
    REQUIRE(f.flags() == E::a);
    REQUIRE(f.flags_any(E::a));
    REQUIRE(f.flags_all(E::a));
    f.set_flags(E::b);
    REQUIRE(f.flags() == E::a_and_b);
    f.set_flag(42);
    REQUIRE(f.flags_all(E::a_and_b));

    
}

TEST_CASE("flagged with scoped enums") {
    // Normally, use of scoped enums are annoying, err, I mean strongly typed
    // because we must sprinkle static casts.  But Flagged is templated and can
    // accept the Enum type so things work correctly.

    enum class E {a=1<<1, b=1<<2, c=1<<3, a_and_b=a|b, b_and_c=b|c};
    Flagged<E> f;
    f.set_flags(E::a);
    REQUIRE(f.flags() == E::a);
    REQUIRE(f.flags_any(E::a));
    REQUIRE(f.flags_all(E::a));
    f.set_flags(E::b);
    REQUIRE(f.flags() == E::a_and_b);
    f.set_flag(42);
    REQUIRE(f.flags_all(E::a_and_b));



    
}



// The rest is mostly gemini, with small cleanups/fixes.  It tests flagged at in
// terms of bit-level tests.  Normal users probably never want to see all these
// 0's and 1's and should stick with enums as above.

TEST_CASE("flagged initial state") {
    Flagged f;
    REQUIRE(f.flags() == 0);
    REQUIRE(!f.flags_any(1));
    REQUIRE(!f.flag_at(0));
}

TEST_CASE("flagged set_flag and flag_at") {
    Flagged f;
    f.set_flag(0);
    REQUIRE(f.flags() == (1ULL << 0));
    REQUIRE(f.flag_at(0));
    REQUIRE(!f.flag_at(1));

    f.set_flag(5);
    REQUIRE(f.flags() == ((1ULL << 0) | (1ULL << 5)));
    REQUIRE(f.flag_at(0));
    REQUIRE(f.flag_at(5));
    REQUIRE(!f.flag_at(4));

    const auto max_index = Flagged<>::max_index;

    f.set_flag(max_index); // Max index
    
    REQUIRE(f.flags() == ((1ULL << 0) | (1ULL << 5) | (1ULL << max_index)));
    REQUIRE(f.flag_at(max_index));

    // Test out of bounds index - should be a no-op
    f.set_flag(Flagged<>::max_flags);
    REQUIRE(f.flags() == ((1ULL << 0) | (1ULL << 5) | (1ULL << max_index))); // Should not change
    f.set_flag(255); // max index_t, but oobs.
    REQUIRE(f.flags() == ((1ULL << 0) | (1ULL << 5) | (1ULL << max_index))); // Should not change
}

TEST_CASE("flagged set_flags and flags_all") {
    Flagged f;
    f.set_flags(0b0101); // Binary for 5
    REQUIRE(f.flags() == 5);

    f.set_flags(0b1000); // Binary for 8
    REQUIRE(f.flags() == (5 | 8)); // 13 (0b1101)

    f.set_flags(0ULL); // Setting no flags should not change existing ones
    REQUIRE(f.flags() == 13);

    f.clear_flags();
    f.set_flags(~0ULL); // Set all bits
    REQUIRE(f.flags() == ~0ULL);
}

TEST_CASE("flagged clear_flags") {
    Flagged f;
    f.set_flags((1ULL << 10) | (1ULL << 20));
    REQUIRE(f.flags() != 0);

    f.clear_flags();
    REQUIRE(f.flags() == 0);

    f.clear_flags(); // Clearing already clear flags
    REQUIRE(f.flags() == 0);
}

TEST_CASE("flagged flags_set") {
    Flagged f;
    f.set_flags(0b10110); // 22 (bits 1, 2, 4 set)

    REQUIRE(f.flags_set(0b00110) == 0b00110); // bits 1, 2 are common
    REQUIRE(f.flags_set(0b10000) == 0b10000); // bit 4 is common
    REQUIRE(f.flags_set(0b01000) == 0b00000); // bit 3 not common
    REQUIRE(f.flags_set(0b11111) == 0b10110); // all common with self
    REQUIRE(f.flags_set(0b00000) == 0b00000); // no common flags
    REQUIRE(f.flags_set(~0ULL) == 0b10110); // common with all set
}

TEST_CASE("flagged flags_any") {
    Flagged f;
    f.set_flags(0b10100); // 20 (bits 2, 4 set)

    REQUIRE(f.flags_any(0b00100)); // bit 2 is set
    REQUIRE(f.flags_any(0b10000)); // bit 4 is set
    REQUIRE(f.flags_any(0b10100)); // bits 2 and 4 are set
    REQUIRE(f.flags_any(0b11100)); // bits 2 and 4 are set (and bit 0, 1 not)

    REQUIRE(!f.flags_any(0b00001)); // bit 0 not set
    REQUIRE(!f.flags_any(0b00010)); // bit 1 not set
    REQUIRE(!f.flags_any(0b01000)); // bit 3 not set
    REQUIRE(!f.flags_any(0b00000)); // no flags to check
}

TEST_CASE("flagged keep_flags") {
    Flagged f;
    f.set_flags(0b110110); // 54 (bits 1, 2, 4, 5 set)
    REQUIRE(f.flags() == 0b110110);

    f.keep_flags(0b10110); // Keep bits 1, 2, 4
    REQUIRE(f.flags() == 0b10110); // bit 5 not kept

    f.keep_flags(0b110001); // keep bits 0, 4, 5 but only 1,2,4 there
    REQUIRE(f.flags() == 0b10000); // bit 4 only left

    f.clear_flags();
    f.set_flags(0b1111); // Bits 0, 1, 2, 3 set

    f.keep_flags(0b0101); // Keep bits 0, 2
    REQUIRE(f.flags() == 0b0101);

    f.keep_flags(0b0000); // Keep no flags
    REQUIRE(f.flags() == 0);

    f.set_flags(0b101010);
    f.keep_flags(~0ULL); // Keep all flags
    REQUIRE(f.flags() == 0b101010);
}

TEST_CASE("flagged Chaining operations") {
    Flagged f;
    f.set_flag(0);
    f.set_flags(1ULL << 2);
    f.keep_flags((1ULL << 0) | (1ULL << 1) | (1ULL << 2)); // Should keep 0 and 2
    REQUIRE(f.flags() == ((1ULL << 0) | (1ULL << 2)));
    f.clear_flags();
    REQUIRE(f.flags() == 0);
}

TEST_CASE("flagged Edge cases for index_t") {
    const auto max_index = Flagged<>::max_index;
    const auto max_flags = Flagged<>::max_flags;
    using index_t = Flagged<>::index_t;

    Flagged f;
    f.set_flag(index_t(0));
    REQUIRE(f.flag_at(index_t(0)));

    f.set_flag(index_t(max_index));
    REQUIRE(f.flag_at(index_t(max_index)));

    // Out of bounds check for set_flag
    f.clear_flags();
    f.set_flag(index_t(max_flags));
    REQUIRE(f.flags() == 0);

    // Out of bounds check for flag_at
    REQUIRE(!f.flag_at(index_t(max_flags)));
    REQUIRE(!f.flag_at(index_t(255))); // Max value for std::byte
}
