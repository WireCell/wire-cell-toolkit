---
name: wct-testing
description: How to write and run tests for the Wire-Cell Toolkit. Use when adding tests, reproducing a bug before fixing it, or choosing between doctest, atomic, and BATS integration tests.
---

# Testing in the Wire-Cell Toolkit

Write tests with good coverage for new code. **Before fixing a bug, write a test that
reproduces it by failing, then fix code until the test passes.**

Pick the test type by what the code under test needs:

## Unit tests (function-level) — Doctest

- File: `<pkg>/test/doctest-<name>.cxx`, using the C++ Doctest framework.
- Run the entire suite: `build/wcdoctest`.
- Prefer these whenever the behavior can be checked by calling functions/classes directly.

## Atomic tests — a `main()` with no arguments

- Use when the test needs its own `main()` but takes **no** command-line arguments.
- File: `<pkg>/test/test-<name>.cxx` → compiles to `build/<pkg>/test-<name>`.

## Integration tests — BATS (bats-core)

Use to drive WCT CLIs such as `wire-cell`.

- File: `<pkg>/test/test-<name>.bats`.
- Start the file with:
  ```bash
  bats_load_library wct-bats.sh
  ```
  which pulls in the helper library `test/wct-bats.sh`.
- Run with the vendored runner: `test/bats/bin/bats <pkg>/test/test-<name>.bats`.
- Helper reference and logging levels: `test/docs/bats.org` and the header of
  `test/wct-bats.sh`.

## Checklist

- [ ] For a bug fix, the failing test was written and seen to fail *first*.
- [ ] New behavior is covered by the most direct test type that fits.
- [ ] Unit tests pass under `build/wcdoctest`; integration tests pass under the bats runner.
