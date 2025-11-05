#ifndef WIRECELL_UTIL_FLAGGED
#define WIRECELL_UTIL_FLAGGED

#include <cstddef>              // for size_t

namespace WireCell {

          
    /** Flagged is a simple class to manage "flags".
     *
     * It may be used as a mixin base class.
     *
     * Simple rules of flags:
     * - A flag is true ("1", "on", "set", etc) or false ("0", "off", "unset", etc).
     * - Up to 64 flags are supported.
     * - A flag state is represented by one bit in a bit mask of size_t.
     * - The state of all flags can be bit-wise queried with a bit map.
     * - An individual flag can also be addressed by its "bit index".
     *
     * Each friendly subclass is expected to provide its users an enum that
     * defines the possible user-level flags.  These definitions should
     * respect the bit-map nature of flags.
     *
     * enum {a=1<<1,b=1<<2,c=1<<3,b_and_c=(1<<2|1<<3)};
     * MyFlags f = b_and_c;
     * set_flags(f);
     * assert(flags_set(f));
     * assert(flags_set(b));
     * assert(flag_at(3));  // "c"
     *
     * When using "scoped enums" aka "strong enum" aka "enum class", the user
     * can make a Flagged with that type to avoid having to pepper user code
     * with static_cast<size_t>.
     *
     */
    template<typename FlagsType = size_t>
    class Flagged {
    public:
        
        /// The storage type for the flags bit map.
        using flags_t = size_t;

        /// A "bit index" must take values in [0,63], inclusive.
        using index_t = unsigned char;

        /// Maximum number of flags we can support.
        static constexpr index_t max_flags = 64u;

        /// The maximum index. 
        static constexpr index_t max_index = 63u;

        /// Flag setters

        /// Set a specific bit to true.  If index is out of bounds, this is a quiet no-op.
        void set_flag(index_t index) {
            if (index < max_flags) { m_flags |= ( 1ULL << index); };
        }

        /// Set any number of flags.  The bits set in "flags" are set in our
        /// internal flags.  Any existing set flags are kept.
        void set_flags(FlagsType flags) { m_flags |= static_cast<flags_t>(flags); }

        /// Keep only the set flags in "flags" that are also set in our flags.
        /// New flags are not set, this disable existing flags not in set in the
        /// input.
        void keep_flags(FlagsType flags) { m_flags &= static_cast<flags_t>(flags); }

        /// Set all flagged bits to zero.
        void clear_flags() { m_flags = 0; }

        /// Query flag state

        /// Get full flag bit map
        FlagsType flags() const { return static_cast<FlagsType>(m_flags); }

        /// Return which flags in the given "flags" are set in our flags.
        FlagsType flags_set(FlagsType flags) const {
            return static_cast<FlagsType>( m_flags & static_cast<flags_t>(flags) );
        }

        /// Return true if any of the given "flags" are set in our flags.
        bool flags_any(FlagsType flags) const { return (m_flags & static_cast<flags_t>(flags)) != 0; }

        /// Return true if all flags in "flags" are set in our flags.
        bool flags_all(FlagsType flags) const {
            return (m_flags & static_cast<flags_t>(flags)) == static_cast<flags_t>(flags);
        }

        /// Return true if flag at given bit index is set.
        bool flag_at(index_t index) const {
            if (index < max_flags) { return flags_any( 1ULL << index); };
            return false;
        }

    private:
        flags_t m_flags{0};         // bitmap of up to 64 flags
        

    };            
}

#endif
