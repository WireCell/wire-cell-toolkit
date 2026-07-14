# Consolidated Efficiency Summary

All efficiency observations across the apps module.

---

## Overview

The `apps` module is **not performance-critical**. It runs once at startup (Main class)
or as one-shot diagnostic tools (dumpers, CLI tools). None of the code is in a hot
loop or processes large data volumes. The efficiency observations below are therefore
minor -- they are noted for code quality and correctness of patterns, not because
they materially affect WCT runtime performance.

---

## Efficiency Observations

| ID | File:Line | Description | Impact |
|----|-----------|-------------|--------|
| EFF-MAIN-1 | Main.cxx:304-412 | **Four iterations over m_cfgmgr.all()**: Instantiation (304), naming (323), configuring (342), and finalize (399) each iterate the full configuration array. Each pass re-extracts type/name from JSON. Could be consolidated into fewer passes or cached. | Low -- O(100) components typical |
| EFF-MAIN-2 | Main.cxx:263-269 | **Persist::Parser created per config file**: A new parser instance is constructed inside the loop for each config file. If parser construction involves expensive initialization (Jsonnet VM), this repeats that cost unnecessarily. | Low to Medium -- depends on Parser internals |
| EFF-MAIN-3 | Main.cxx:118+ | **String copies in range-for loops**: `for (auto x : ...)` copies each element. Using `const auto& x` would avoid copying. Appears ~18 times across cmdline() and initialize(). | Negligible -- short strings at startup |
| EFF-MAIN-4 | Main.cxx:362-393 | **Two-pass app lookup and execution**: First loop builds `app_objs` vector, second loop executes. Could be a single loop. | Negligible -- 1-3 apps typical |

---

## Positive Patterns

Several aspects of the code are well-designed from an efficiency standpoint:

1. **Factory caching**: The `NamedFactory` system caches component instances by (type, name).
   Repeated lookups return the same pointer -- no redundant construction.

2. **Configuration merge**: The `update()` function performs in-place recursive JSON merge
   rather than rebuilding configurations from scratch.

3. **Plugin lazy loading**: Plugins are loaded only when explicitly requested, not eagerly
   discovered.

4. **TBB thread control**: Thread limiting via `tbb::global_control` is scoped to
   the `operator()()` method using `unique_ptr`, ensuring proper cleanup.

5. **Diagnostic dumpers**: Single-pass hierarchy traversal in AnodeDumper, efficient
   factory lookups, reasonable JSON construction patterns.

---

## Recommendations (if changes are considered)

1. **Highest value**: Add bounds checking to `String::split` results in Main.cxx
   (BUG-MAIN-1). This is a correctness fix, not efficiency, but it's the most
   impactful change in this module.

2. **Low effort**: Change `for (auto x : ...)` to `for (const auto& x : ...)` across
   all loops. A mechanical change that follows modern C++ best practices.

3. **Consider**: Move `Persist::Parser` construction out of the config file loop in
   `initialize()`. Requires understanding Parser internals to confirm it's safe to reuse.

4. **Not recommended**: Consolidating the four `m_cfgmgr.all()` passes. The current
   separation has clear semantic meaning (instantiate -> name -> configure -> finalize)
   and the ordering matters (naming must precede log config, which must precede
   configuration). Merging them would sacrifice clarity for negligible performance gain.
