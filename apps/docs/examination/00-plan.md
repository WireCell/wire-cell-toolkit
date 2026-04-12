# apps/ Code Examination Plan

## Scope

Systematic examination of the Wire Cell Toolkit `apps` module (~800 lines of C++
across 12 source files) for:

1. **Potential bugs** -- logic errors, missing error handling, null dereferences, resource leaks
2. **Running efficiency** -- unnecessary copies, redundant iterations, missing checks
3. **Algorithm documentation** -- what each component does, how it fits into the WCT pipeline

## File Inventory

### Group 1: Main Class (~520 lines) -- highest priority
| File | Lines | Role |
|------|-------|------|
| `Main.h` | 118 | Main class interface and member declarations |
| `Main.cxx` | 414 | Core orchestrator: cmdline parsing, initialization, execution |

### Group 2: Application Components (~220 lines)
| File | Lines | Role |
|------|-------|------|
| `AnodeDumper.h` | 25 | AnodeDumper interface |
| `AnodeDumper.cxx` | 77 | Anode geometry dumper |
| `ConfigDumper.h` | 25 | ConfigDumper interface |
| `ConfigDumper.cxx` | 65 | Component default configuration dumper |
| `NodeDumper.h` | 25 | NodeDumper interface |
| `NodeDumper.cxx` | 82 | Node type metadata dumper |

### Group 3: CLI Executables (~250 lines)
| File | Lines | Role |
|------|-------|------|
| `wire-cell.cxx` | 29 | Main WCT entry point |
| `wcsonnet.cxx` | 96 | Jsonnet compiler CLI |
| `wcwires.cxx` | 124 | Wire geometry validator/converter CLI |

### Group 4: Tests
| File | Lines | Role |
|------|-------|------|
| `test_dlopen.cxx` | 54 | Plugin/factory singleton tests |
| `test_apps.bats` | 23 | CLI smoke tests |
| `anode-dumper.jsonnet` | 28 | AnodeDumper configuration example |

## Output Documents

| Document | Contents |
|----------|----------|
| `00-plan.md` | This plan |
| `01-overview.md` | Architecture overview, design patterns, pipeline flow |
| `02-main-class.md` | Main.h/Main.cxx deep examination |
| `03-app-components.md` | AnodeDumper, ConfigDumper, NodeDumper examination |
| `04-cli-tools.md` | wire-cell.cxx, wcsonnet.cxx, wcwires.cxx examination |
| `05-bug-summary.md` | Consolidated bug table |
| `06-efficiency-summary.md` | Consolidated efficiency recommendations |
