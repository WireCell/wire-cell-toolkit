# Consolidated Bug Summary

All potential bugs found across the apps module, sorted by severity.

---

## Medium Severity

These bugs may cause incorrect behavior or unhelpful failures under plausible conditions.

| ID | File:Line | Description |
|----|-----------|-------------|
| BUG-MAIN-1 | Main.cxx:132-154 | **String split without bounds checking**: `String::split(vev, "=")` result accessed at `vv[0]` and `vv[1]` without size check. If user passes `-V foo` (no `=`), accessing `vv[1]` is out-of-bounds (crash). Affects `--ext-str`, `--ext-code`, `--tla-str`, `--tla-code`. |
| BUG-WC-1 | wire-cell.cxx:16-27 | **Non-WireCell exceptions not caught**: Only `WireCell::Exception` is caught. `boost::program_options::error` from malformed arguments, `std::runtime_error` from component failures, etc. all produce unhelpful `terminate()` messages instead of readable errors. |
| BUG-WCS-1 | wcsonnet.cxx:91-92 | **Output file write errors silently ignored**: `ofstream` not checked for open/write failure. Program returns 0 (success) even if output cannot be written (permission denied, disk full, etc.). |
| BUG-WCS-2 | wcsonnet.cxx:90 | **Jsonnet parse exceptions not caught**: `parser.load()` throws on syntax errors or missing files. Exception propagates to `terminate()` with no readable error message. |
| BUG-WCS-3 | wcsonnet.cxx:20-21 | **Values containing "=" are truncated**: `String::split(one, "=")` splits on ALL `=` chars. `--ext-str 'key=a=b'` is rejected as invalid instead of being parsed as key=`a=b`. Same issue in Main.cxx but crashes there instead of erroring. |

## Low Severity

Minor issues or edge cases unlikely to trigger in normal operation.

| ID | File:Line | Description |
|----|-----------|-------------|
| BUG-AD-1 | AnodeDumper.h:1 | **Header guard collision**: `AnodeDumper.h` uses `WIRECELLAPPS_CONFIGDUMPER` guard, same as `ConfigDumper.h`. Including both in one translation unit silently skips the second. Currently no file does this. |
| BUG-MAIN-2 | Main.cxx:173-182 | **Logsink with colons silently dropped**: Split on `:` produces 3+ elements for paths like `C:\logs\out.log`; the sink is silently ignored (no warning). |
| BUG-MAIN-3 | Main.cxx:373-384 | **Duplicate debug log message**: When compiled with TBB, the "executing N apps" message is logged twice (unconditional at line 373, then again inside `#if HAVE_TBB_LIB`). Without TBB, the unconditional message misleadingly includes "thread limit". |
| BUG-MAIN-4 | Main.cxx:46,409 | **finalize() uses logger during destruction**: Destructor calls `finalize()` which uses `log->debug()`. If spdlog is already torn down (static destruction order), this could crash. Unlikely with typical stack-allocated Main. |
| BUG-MAIN-5 | Main.cxx:198-200 | **No validation of --threads value**: Negative values accepted and passed to `tbb::global_control::max_allowed_parallelism`, which is undefined behavior for negative inputs. |
| BUG-WCW-1 | wcwires.cxx:69-81 | **Silent no-op without flags**: Running `wcwires file.json` without `-o` or `-v` silently does nothing and returns success. |
| BUG-WCW-2 | wcwires.cxx:75,114 | **catch(...) discards error details**: Validation failure details lost; user only sees "input invalid" with no indication of what failed. |
| BUG-WCW-4 | wcwires.cxx:53-57 | **type_size(1) on flag options**: `add_flag` with `type_size(1)` is unusual and may confuse CLI11's argument parsing for boolean flags. |

## Informational

Code quality issues that don't affect behavior.

| ID | File:Line | Description |
|----|-----------|-------------|
| BUG-CD-1 | ConfigDumper.cxx:49 | **Redundant convert<string>**: `c` is already a string; `convert<string>(c)` is unnecessary. |
| BUG-ND-1 | NodeDumper.cxx:68-72 | **Misleading loop variable**: `intype` used for both input and output type loops. Readability issue only. |
| BUG-WCS-4 | wcsonnet.cxx:53-66 | **Dead commented-out code**: 14 lines of commented-out alternative CLI parsing logic. |
| BUG-WCW-3 | wcwires.cxx:15-31 | **Dead commented-out code**: 17 lines of commented-out `parse_param` function. |
