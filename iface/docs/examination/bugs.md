# Potential Bugs in iface

This document records potential bugs found during code examination of
the `iface/` directory.  The iface package is predominantly a set of
abstract interface classes (134 headers) with minimal implementation
(12 source files).  Only files containing actual logic are examined.

---

## BUG-1: Copy-paste error in `ascending_time` comparator (IDepo.cxx:17-26) **[CONFIRMED BUG]**

```cpp
bool WireCell::ascending_time(const WireCell::IDepo::pointer& lhs, const WireCell::IDepo::pointer& rhs)
{
    if (lhs->time() == rhs->time()) {
        if (lhs->pos().x() == lhs->pos().x()) {   // BUG: compares lhs to itself
            return lhs.get() < rhs.get();
        }
        return lhs->pos().x() < lhs->pos().x();   // BUG: compares lhs to itself
    }
    return lhs->time() < rhs->time();
}
```

**Problem**: Both the inner `if` condition and the `return` statement compare
`lhs->pos().x()` to **itself** instead of to `rhs->pos().x()`.  This means:

- The inner `if` is always true (a value always equals itself, barring NaN),
  so the x-position tiebreaker is never actually used.
- Even if the `if` were false, the return line would always return `false`.

**Impact**: When two depositions have the same time, the comparator falls
through to pointer-address comparison instead of using x-position as
the intended tiebreaker.  Sorting by this comparator may produce a
different (non-physical) order than intended.  In practice, the
pointer-address fallback still provides a strict weak ordering, so it
won't crash, but the ordering is arbitrary rather than
physics-motivated.

**Fix**: Change both occurrences of `lhs->pos().x()` on the RHS of the
comparison to `rhs->pos().x()`.

---

## BUG-2: Same copy-paste error in `descending_time` (IDepo.cxx:29-38) **[CONFIRMED BUG]**

```cpp
bool WireCell::descending_time(const WireCell::IDepo::pointer& lhs, const WireCell::IDepo::pointer& rhs)
{
    if (lhs->time() == rhs->time()) {
        if (lhs->pos().x() == lhs->pos().x()) {   // BUG: same as BUG-1
            return lhs.get() > rhs.get();
        }
        return lhs->pos().x() > lhs->pos().x();   // BUG: same as BUG-1
    }
    return lhs->time() > rhs->time();
}
```

**Problem**: Identical to BUG-1.  Copy-paste propagated the same mistake.

---

## BUG-3: Copy-paste error in `ascending_index` comparator (IWire.cxx:11-17) **[CONFIRMED BUG]**

```cpp
bool WireCell::ascending_index(IWire::pointer lhs, IWire::pointer rhs)
{
    if (lhs->planeid() == rhs->planeid()) {
        return lhs->index() < rhs->index();
    }
    return lhs->planeid() < lhs->planeid();  // BUG: compares lhs to itself
}
```

**Problem**: The final line compares `lhs->planeid()` against **itself**
instead of `rhs->planeid()`.  This always returns `false`, which means
when two wires are on different planes, neither is considered "less
than" the other.  This could break strict weak ordering requirements
for `std::sort` if wires from multiple planes are mixed, potentially
leading to undefined behavior.

**Fix**: Change the last line to `return lhs->planeid() < rhs->planeid();`.

---

## BUG-4: `ISequence::end()` returns `cbegin()` instead of `cend()` (ISequence.h:46-48) **[CONFIRMED BUG]**

```cpp
virtual const_iterator begin() const { return cbegin(); }
virtual const_iterator end() const { return cbegin(); }   // BUG: should be cend()
virtual iterator begin() { return iterator(cbegin()); }
virtual iterator end() { return iterator(cbegin()); }     // BUG: should be cend()
```

**Problem**: Both `end()` overloads (const and non-const) return `cbegin()`
instead of `cend()`.  Any iteration over an `ISequence` via `begin()/end()`
will see an empty range.

**Impact**: `ISequence` appears to be a legacy abstraction (the `IData`
types use `shared_vector` instead, and `ITrace`/`IFrame` etc. don't
inherit from `ISequence`).  If nothing currently uses
`ISequence::begin()/end()` directly, the bug is dormant.  But it is
still a trap for any future code that tries to iterate an `ISequence`.

**Fix**: Change both `end()` implementations to return `cend()`.

---

## BUG-5: `cluster_node_t::code()` out-of-bounds access when `ptr.index() == 0` (ICluster.h:95-102)

```cpp
inline char code() const
{
    auto ind = ptr.index();
    if (ind == std::variant_npos) {
        return 0;
    }
    return known_codes[ind-1];  // if ind==0, accesses known_codes[-1]
}
```

**Problem**: The `ptr_t` variant has `size_t` at index 0 (the default/null
state).  When `ind == 0`, the expression `known_codes[ind-1]` wraps to
`known_codes[SIZE_MAX]` which is undefined behavior (out-of-bounds
access on a `std::string`).

**Mitigating factor**: Index 0 represents a default-constructed node with
no valid pointer, so it should be rare in practice.  However, the
comment at line 86 says "the index to a code is one less than the
ptr.index()", implying this offset is intentional for valid types
(1-5), but the null case (0) is not handled.

**Fix**: Add a check for `ind == 0` and return `\0` (or another sentinel).

---

## BUG-6: `cluster_node_t::ident()` builds vector of closures on every call (ICluster.h:113-129)

```cpp
int ident() const
{
    using ident_f = std::function<int()>;
    std::vector<ident_f> ofs {
        [](){return 0;},
        [&](){return std::get<channel_t>(ptr)->ident();},
        ...
    };
    const auto ind = ptr.index();
    if (ind == std::variant_npos) {
        return 0;
    }
    return ofs[ind]();
}
```

**Problem (correctness)**: If `ind == 0`, `ofs[0]()` returns 0 which is
fine.  But there is no check that `ind` is within range of `ofs`.
In theory, `std::variant` index is always within `[0, N)`, but the
`std::variant_npos` check suggests defensive coding; the same
defense should check `ind < ofs.size()`.

**Problem (efficiency)**: See EFFICIENCY section.

---

## BUG-7: Floating-point equality comparison in `IDepoDriftCompare` (IDepo.h:69)

```cpp
if (t1 == t2) {
    return lhs.get() < rhs.get();
}
```

**Problem**: Comparing floating-point times with `==` is fragile.  Two
depositions that are extremely close in drift-adjusted time but not
bit-identical will not be caught by this tiebreaker, and instead be
ordered by their computed `t1 < t2`, which for nearly-equal values is
sensitive to floating-point rounding.

**Mitigating factor**: The comparator is used in a `std::set` which
requires strict weak ordering.  Exact equality is sufficient for
the ordering contract.  The concern is practical: the pointer-address
tiebreaker may not trigger when physics would suggest it should.

**Severity**: Low.  The same pattern exists in `ascending_time` and
`descending_time`.

---

## BUG-8: `WirePlaneId::face()` extracts only 1 bit (WirePlaneId.cxx:36)

```cpp
int WireCell::WirePlaneId::face() const { return (m_pack & (1 << face_shift)) >> 3; }
```

**Problem**: `face_shift = 3` and the mask is `(1 << 3) = 0x8`, extracting
only a single bit.  This limits face values to 0 or 1.  Meanwhile,
the constructor packs `face << face_shift` which allows face values
larger than 1.  If a face value of 2 or higher were passed, bit 4
would be set but `face()` would return 0 (since bit 3 is 0).

**Mitigating factor**: In current WCT usage, faces are always 0 or 1
(MicroBooNE has 1 face, DUNE has 2), so this is not triggered in
practice.  But the encoding is fragile if future detectors have
more faces.

---

## BUG-9: `WirePlaneId::valid()` can never detect negative apa/face (WirePlaneId.cxx:39-46)

```cpp
bool WireCell::WirePlaneId::valid() const
{
    if (apa() < 0) return false;
    if (face() < 0) return false;
    ...
}
```

**Problem**: `face()` returns a non-negative value (it extracts bits via
mask-and-shift), and `apa()` returns `m_pack >> apa_shift` which for
negative packed values would produce a large positive value (since
`m_pack` is `int` and the shift is arithmetic).  These checks can
never actually trigger for the way the value is encoded.

**Severity**: Low.  Defensive but unreachable code.

---

## BUG-10: `WirePlaneId::convert` template assumes cfg has at least 3 elements (WirePlaneId.h:88-94)

```cpp
template <>
inline WireCell::WirePlaneId convert<WireCell::WirePlaneId>(const Configuration& cfg,
                                                            const WireCell::WirePlaneId& def)
{
    return WireCell::WirePlaneId(iplane2layer[convert<int>(cfg[0])], convert<int>(cfg[1], 0),
                                 convert<int>(cfg[2], 0));
}
```

**Problem**: No bounds check on `cfg`.  If `cfg` has fewer than 3 elements,
`cfg[0]`, `cfg[1]`, or `cfg[2]` may return null/invalid JSON values
leading to unexpected behavior.  Also, `iplane2layer[convert<int>(cfg[0])]`
has no bounds check on the index; if `cfg[0]` is not 0, 1, or 2 then
`iplane2layer` is accessed out of bounds.

**Severity**: Medium.  This depends on whether configuration is always
validated upstream.

---

## BUG-11: Duplicate include in IfaceDesctructors.cxx (line 81)

```cpp
#include "WireCellIface/IQueuedoutNode.h"
#include "WireCellIface/IQueuedoutNode.h"   // duplicate
```

**Impact**: Harmless due to include guards.  Cosmetic issue only.

---

## BUG-12: Duplicate include in WirePlaneId.h (lines 7 and 9)

```cpp
#include "WireCellUtil/Spdlog.h"
#include <ostream>
#include <functional>
#include "WireCellUtil/Spdlog.h"   // duplicate
```

**Impact**: Harmless due to include guards.  Cosmetic issue only.

---

## Summary

| ID     | File            | Severity | Category         |
|--------|-----------------|----------|------------------|
| BUG-1  | IDepo.cxx       | High     | Copy-paste typo  |
| BUG-2  | IDepo.cxx       | High     | Copy-paste typo  |
| BUG-3  | IWire.cxx       | High     | Copy-paste typo  |
| BUG-4  | ISequence.h     | Medium   | Logic error      |
| BUG-5  | ICluster.h      | Medium   | Off-by-one/UB    |
| BUG-6  | ICluster.h      | Low      | Missing bounds   |
| BUG-7  | IDepo.h         | Low      | Float equality   |
| BUG-8  | WirePlaneId.cxx | Low      | Bit-packing      |
| BUG-9  | WirePlaneId.cxx | Low      | Dead code        |
| BUG-10 | WirePlaneId.h   | Medium   | Missing bounds   |
| BUG-11 | IfaceDesctructors.cxx | Cosmetic | Duplicate include |
| BUG-12 | WirePlaneId.h   | Cosmetic | Duplicate include |
