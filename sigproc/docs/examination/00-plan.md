# Signal Processing Code Examination Plan

## Scope

Systematic examination of the Wire Cell Toolkit `sigproc` module (~16,600 lines of C++
across 30 source files) for:

1. **Potential bugs** -- logic errors, off-by-one, uninitialized values, numeric pitfalls
2. **Running efficiency** -- unnecessary copies, redundant computation, memory allocation patterns
3. **Algorithm documentation** -- what each component does, how it fits into the pipeline

## File Grouping and Examination Order

### Group 1: Core Signal Processing Pipeline (~6,000 lines) -- highest priority
| File | Lines | Role |
|------|-------|------|
| `OmnibusSigProc.cxx` | 1,925 | Main signal processing orchestrator |
| `ROI_refinement.cxx` | 3,155 | Region-of-interest refinement (largest file) |
| `ROI_formation.cxx` | 937 | Initial ROI identification |

### Group 2: L1SP Filter (~700 lines)
| File | Lines | Role |
|------|-------|------|
| `L1SPFilter.cxx` | 692 | L1-norm sparse deconvolution filter |

### Group 3: Noise Filtering Subsystem (~2,000 lines)
| File | Lines | Role |
|------|-------|------|
| `OmnibusNoiseFilter.cxx` | 294 | Top-level noise filter orchestrator |
| `OmnibusPMTNoiseFilter.cxx` | 405 | PMT-specific noise removal |
| `OmniChannelNoiseDB.cxx` | 702 | Configurable channel noise database |
| `SimpleChannelNoiseDB.cxx` | 522 | Simple channel noise database |
| `NoiseModeler.cxx` | 257 | Noise spectrum modeling |
| `NoiseRanker.cxx` | (in src/) | Noise ranking utility |

### Group 4: Detector-Specific Implementations (~5,100 lines)
| File | Lines | Role |
|------|-------|------|
| `Protodune.cxx` | 1,020 | ProtoDUNE-SP configuration |
| `ProtoduneHD.cxx` | 1,043 | ProtoDUNE-HD configuration |
| `ProtoduneVD.cxx` | 1,354 | ProtoDUNE-VD configuration |
| `Microboone.cxx` | 1,352 | MicroBooNE configuration |
| `Icarus.cxx` | (in src/) | ICARUS configuration |
| `DuneCrp.cxx` | 369 | DUNE CRP configuration |

### Group 5: Utilities and Support (~1,800 lines)
| File | Lines | Role |
|------|-------|------|
| `PeakFinding.cxx` | 474 | Peak detection algorithm |
| `ChannelSelector.cxx` | 194 | Channel selection/routing |
| `FrameMerger.cxx` | 147 | Frame merging logic |
| `Diagnostics.cxx` | 144 | Diagnostic output |
| `Derivations.cxx` | ~100 | Derived quantities |
| Filters (Hf/Lf), Responses, etc. | ~300 | Filter/response implementations |

## Methodology

For each file/group:

1. **Read the header** to understand the class interface and configuration
2. **Read the implementation** looking for:
   - Uninitialized variables, use-after-move, dangling references
   - Off-by-one errors in loop bounds and array indexing
   - Numeric issues: division by zero, unclamped acos/asin, float comparison
   - Resource leaks, missing error handling at boundaries
   - Logic errors in conditional branches
3. **Assess efficiency**:
   - Unnecessary copies (large containers passed by value, repeated string ops)
   - Redundant computation (repeated lookups, recomputable values in loops)
   - Memory allocation patterns (repeated allocations in hot loops)
   - FFT usage patterns (plan reuse, unnecessary transforms)
4. **Document the algorithm** at a level useful for a physicist/developer working on the code

## Output Documents

All documents go in `sigproc/docs/examination/`:

| Document | Contents |
|----------|----------|
| `00-plan.md` | This plan |
| `01-overview.md` | High-level algorithm overview and pipeline architecture |
| `02-core-sigproc.md` | Core signal processing: OmnibusSigProc, ROI formation/refinement |
| `03-l1sp-filter.md` | L1SP filter algorithm and examination |
| `04-noise-filtering.md` | Noise filtering subsystem examination |
| `05-detector-specific.md` | Detector-specific code examination |
| `06-utilities.md` | Utility components examination |
| `07-bug-summary.md` | Consolidated bug/issue list across all files |
| `08-efficiency-summary.md` | Consolidated efficiency recommendations |
