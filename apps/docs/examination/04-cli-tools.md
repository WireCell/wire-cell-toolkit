# CLI Tools Examination

Files examined:
- `apps/wire-cell.cxx` (29 lines)
- `apps/wcsonnet.cxx` (96 lines)
- `apps/wcwires.cxx` (124 lines)

---

## wire-cell.cxx

### Algorithm

Minimal entry point wrapper around `Main`:

```
main(argc, argv):
    Create Main instance
    rc = m.cmdline(argc, argv)     // parse CLI, returns 1 if --help/--version
    if rc == 0:
        m.initialize()             // load config, plugins, components
        m()                        // execute all apps
    catch WireCell::Exception:
        print error, return 1
    return 0
```

### Potential Bugs

#### BUG-WC-1: Non-WireCell exceptions not caught (wire-cell.cxx:16-27) -- Medium

```cpp
try {
    rc = m.cmdline(argc, argv);
    // ...
}
catch (Exception& e) {
    cerr << errstr(e) << endl;
    return 1;
}
return 0;
```

Only `WireCell::Exception` is caught. Standard exceptions (`std::runtime_error`,
`std::bad_alloc`, `boost::program_options::error`, etc.) propagate uncaught,
producing an unhelpful `terminate()` message. The `cmdline()` method uses
`boost::program_options` which throws `boost::program_options::error` subtypes
on invalid arguments. A malformed command line would produce a cryptic abort
instead of a readable error message.

#### BUG-WC-2: Return code 0 even after --help/--version (wire-cell.cxx:27) -- Low

When `cmdline()` returns 1 (for `--help` or `--version`), the program falls through
to `return 0`. This is actually correct Unix behavior (help and version should
return 0), but the path is not obvious -- the `rc` variable is set to 1 but never
used as the return value.

### Efficiency Notes

No concerns. This is a trivial wrapper.

---

## wcsonnet.cxx

### Algorithm

Compiles Jsonnet configuration files to JSON, honoring WCT conventions:

```
main(argc, argv):
    Parse CLI options (output file, search paths, ext-str, ext-code, tla-str, tla-code)
    Parse each variable parameter as "name=value" pairs
    Validate filename is provided
    Create Persist::Parser with paths and variables
    Load and compile Jsonnet file
    Write JSON output to file
```

The `parse_param()` helper splits on `=` and validates that exactly 2 parts result.

### Potential Bugs

#### BUG-WCS-1: Output file not checked for write errors (wcsonnet.cxx:91-92) -- Medium

```cpp
std::ofstream out(output);
out << jdat << std::endl;
```

The `ofstream` constructor does not throw on failure by default. If the output file
cannot be opened (e.g., permission denied, directory doesn't exist), `out` is in
a fail state and the write silently does nothing. The program returns 0 (success)
despite producing no output.

Should check `out.is_open()` or `out.good()` after construction and after writing.

#### BUG-WCS-2: Persist::Parser::load() exceptions not caught (wcsonnet.cxx:90) -- Medium

```cpp
auto jdat = parser.load(filename);
```

If the Jsonnet file has syntax errors or the file doesn't exist, `parser.load()`
throws. This exception is not caught anywhere, producing an unhelpful `terminate()`
message with no indication of what went wrong.

#### BUG-WCS-3: Values containing "=" are truncated (wcsonnet.cxx:20-21) -- Medium

```cpp
auto two = String::split(one, "=");
if (two.size() != 2) {
```

`String::split` splits on ALL occurrences of `=`. A value like `--ext-str 'key=a=b'`
produces 3 elements and is rejected as invalid. The correct behavior would be to
split only on the first `=`, allowing values to contain `=` characters.

This same issue affects Main.cxx (BUG-MAIN-1), but there it crashes instead of
producing an error message.

#### BUG-WCS-4: Dead commented-out code (wcsonnet.cxx:53-66) -- Low

Lines 53-66 contain commented-out alternative CLI parsing logic. This is dead code
that adds noise but doesn't affect behavior.

### Efficiency Notes

No concerns. This is a simple compile-and-dump tool.

---

## wcwires.cxx

### Algorithm

Loads, validates, and optionally converts wire geometry files:

```
main(argc, argv):
    Parse CLI options (output, correction level, validate, fail-fast, epsilon)
    if no output specified:
        if validate requested:
            Load file with Correction::load
            Run validate(raw, epsilon, fail_fast)
            Print "input valid" or "input invalid"
    else (output specified):
        Map correction int to Correction enum
        Load file with specified correction level
        if validate requested:
            Run validate(store, epsilon, fail_fast)
        Dump corrected store to output file
```

Correction levels:
- `empty` (0): treated as `load`
- `load` (1): raw load, no corrections
- `order` (2): fix wire ordering
- `direction` (3): fix wire direction vectors
- `pitch` (4, default): fix wire pitch values

### Potential Bugs

#### BUG-WCW-1: No action when no output and no validation (wcwires.cxx:69-81) -- Low

```cpp
if (output.empty()) {
    if (do_validate) {
        // ... validate
    }
    // else: nothing happens
}
```

If the user runs `wcwires myfile.json` without `-o` or `-v`, the program silently
does nothing and returns 0. It would be more helpful to either print a warning or
require at least one action flag.

#### BUG-WCW-2: catch(...) discards error details (wcwires.cxx:75, 114) -- Low

```cpp
try {
    validate(raw, repsilon, do_fail_fast);
}
catch (...) {
    std::cerr << "input invalid\n";
    return 1;
}
```

The catch-all discards the exception message. The `validate()` function likely
throws with details about what validation failed. These details are lost, making
it harder for users to fix their wire files.

#### BUG-WCW-3: Commented-out parse_param function (wcwires.cxx:15-31) -- Low

Dead commented-out code (copy of `parse_param` from wcsonnet.cxx). Adds clutter.

#### BUG-WCW-4: type_size(1) on flag options (wcwires.cxx:53-57) -- Low

```cpp
app.add_flag("-v,--validate", do_validate, ...)->type_size(1)->allow_extra_args(false);
app.add_flag("-f,--fail-fast", do_fail_fast, ...)->type_size(1)->allow_extra_args(false);
```

`add_flag` with `type_size(1)` is unusual -- flags normally don't take arguments.
CLI11 may handle this correctly by treating it as a boolean option, but it's
inconsistent with typical flag usage. Users might be confused about whether
`-v true` or just `-v` is expected.

### Efficiency Notes

#### EFF-WCW-1: Redundant file load when validating with output (wcwires.cxx:69-121)

When both validation and output are requested, the file is loaded once with the
specified correction level (line 109). The corrected version is validated and then
dumped. This is actually correct -- you want to validate the corrected output.

However, the validate-only path (line 71) loads with `Correction::load` (raw),
while the output path (line 109) loads with the user-specified correction. These
are intentionally different behaviors: validate-only checks the raw input, while
output mode validates the corrected result. This design is sound.

---

## Common Observations

### CLI Library Split
- `wire-cell.cxx` uses `boost::program_options` (via Main class)
- `wcsonnet.cxx` and `wcwires.cxx` use CLI11 (header-only, vendored as `CLI11.hpp`)

This split means two CLI parsing libraries are dependencies. The CLI11 tools are
standalone and don't go through the Main class, so the split is reasonable.

### Error Handling Gap
All three CLI tools have incomplete error handling for non-WireCell exceptions.
The `wire-cell` wrapper catches only `WireCell::Exception`, and the other two
tools don't wrap their main logic in try/catch at all. Standard exceptions
(from JSON parsing, file I/O, etc.) would produce unhelpful abort messages.
