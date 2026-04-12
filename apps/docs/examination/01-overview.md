# apps/ Architecture Overview and Algorithm Documentation

## Role in the Wire-Cell Toolkit

The `apps` module is the **top-level application layer** of the Wire-Cell Toolkit (WCT).
It sits above the interface (`iface`), utility (`util`), and all component modules
(e.g., `gen`, `sigproc`, `img`). Its responsibilities are:

1. **Bootstrap the toolkit** -- parse CLI arguments, load plugins, read configuration,
   instantiate and configure components
2. **Execute application components** -- run registered `IApplication` objects in sequence
3. **Provide CLI programs** -- `wire-cell` (general entry point), `wcsonnet` (Jsonnet compiler),
   `wcwires` (wire geometry validator)
4. **Provide diagnostic dumpers** -- `AnodeDumper`, `ConfigDumper`, `NodeDumper` for
   inspecting component state

## Architecture

```
CLI Programs          Diagnostic Apps          External Embedders
(wire-cell,           (AnodeDumper,            (LArSoft art/FHiCL)
 wcsonnet, wcwires)    ConfigDumper,
      |                NodeDumper)
      |                    |
      v                    v
  +-------------------------------+
  |          Main class           |  <-- Central orchestrator
  |  cmdline() -> initialize()    |
  |  -> operator()() -> finalize()|
  +-------------------------------+
      |         |         |
      v         v         v
  Persist    Plugin     Factory
  (Jsonnet   Manager    (Named
   config)   (dlopen)   Factory)
      |         |         |
      v         v         v
  +-------------------------------+
  |   WCT Component Ecosystem     |
  |  (IConfigurable, IApplication,|
  |   INode, ITerminal, INamed)   |
  +-------------------------------+
```

## Main Class Lifecycle

The `Main` class orchestrates the entire WCT boot-up sequence. It can be driven
either through `cmdline()` (for CLI use) or through fine-grained setup methods
(for embedding in frameworks like LArSoft).

### Phase 1: Setup (cmdline or explicit calls)

Configuration is accumulated into member variables:
- `m_cfgfiles`: configuration file paths (Jsonnet/JSON)
- `m_plugins`: plugin library names to load
- `m_apps`: application components to execute
- `m_load_path`: search paths for configuration file resolution
- `m_extvars`, `m_extcode`, `m_tlavars`, `m_tlacode`: Jsonnet external variables

### Phase 2: Initialization (initialize)

Sequential steps:
1. **Set up logging** -- configure spdlog pattern and create "main" logger
2. **Load configuration** -- parse each Jsonnet/JSON file using `Persist::Parser`
   with the accumulated load paths and external variables, merge into `ConfigManager`
3. **Extract "wire-cell" entry** -- special config entry carries plugin and app lists
4. **Load plugins** -- use `PluginManager` to dlopen shared libraries
5. **Instantiate components** -- iterate all config entries, call `Factory::lookup<Interface>`
   to create component instances
6. **Set component names** -- find components implementing `INamed`, call `set_name()`
7. **Apply log configuration** -- set log levels for named loggers
8. **Configure components** -- for each `IConfigurable`, get default config, merge with
   user config via `update()`, and call `configure()`

### Phase 3: Execution (operator())

1. Look up all `IApplication` objects by their type:name strings
2. If TBB is available and thread limit is set, create `tbb::global_control`
3. Execute each application sequentially via `aobj->execute()`

### Phase 4: Finalization (finalize, also called by destructor)

1. Iterate all config entries, find components implementing `ITerminal`
2. Call `finalize()` on each terminal component for cleanup

## Design Patterns

### Factory Pattern
All WCT components are created through `NamedFactory`. Components register themselves
via `WIRECELL_FACTORY` macros. The factory system supports:
- Type-based lookup: `Factory::lookup<IFace>(type, name)`
- Named instances: same type can have multiple named instances
- Singleton caching: same (type, name) always returns the same instance

### Plugin Pattern
Component implementations live in shared libraries (plugins). The `PluginManager`
uses `dlopen()` to load them at runtime. This allows the apps layer to be generic --
it doesn't link against any specific component implementations.

### Configuration Pattern
Configuration flows through the system as `Json::Value` objects:
1. Jsonnet files are compiled to JSON
2. JSON is organized as an array of `{type, name, data}` objects
3. Each component's `default_configuration()` provides defaults
4. User config is merged on top via `update()` (recursive merge)
5. The merged config is passed to `configure()`

The special `wire-cell` config entry (type="wire-cell") carries top-level
orchestration data: which plugins to load and which apps to run.

### Interface Segregation
Components implement only the interfaces they need:
- `IApplication`: has `execute()` -- top-level runnable
- `IConfigurable`: has `configure()`, `default_configuration()` -- accepts configuration
- `ITerminal`: has `finalize()` -- needs cleanup
- `INamed`: has `set_name()` -- accepts a string name
- `INode`: has `input_types()`, `output_types()`, `concurrency()`, `category()` -- data flow node

## CLI Programs

### wire-cell
The main WCT entry point. Ultra-thin wrapper: creates `Main`, calls `cmdline()`,
`initialize()`, and `operator()()`. Catches `WireCell::Exception` and prints error details.
This is the program that end-users typically invoke to run WCT processing pipelines.

### wcsonnet
A Jsonnet compiler that is aware of WCT conventions. Uses `Persist::Parser` which
internally uses the Go-based Jsonnet library (faster than the C++ reference implementation).
Honors `WIRECELL_PATH` environment variable for import resolution. Useful for debugging
Jsonnet configuration files independently of WCT execution.

### wcwires
A wire geometry utility that loads, validates, and optionally converts wire description
files. Supports correction levels (load, order, direction, pitch) that progressively
fix common issues in wire geometry definitions. The validator checks geometric properties
with a configurable epsilon tolerance.

## Diagnostic App Components

### AnodeDumper
Dumps anode plane geometry to JSON. Traverses the anode -> face -> plane hierarchy
and extracts identifiers, channel counts, and wire counts. Useful for verifying
detector geometry is correctly configured.

### ConfigDumper
Dumps default configurations for components. Can dump all known `IConfigurable`
components or a specified subset. Useful for discovering what configuration options
a component accepts.

### NodeDumper
Dumps node type metadata including input/output port types, concurrency model,
and category. Useful for understanding the data flow graph structure.

## Test Infrastructure

The module has minimal testing:
- `test_dlopen.cxx`: Verifies plugin loading, factory singleton behavior, and named
  instance distinctness. This is a critical infrastructure test.
- `test_apps.bats`: BATS shell tests verifying that CLI executables exist and produce
  expected output (smoke tests).
- `anode-dumper.jsonnet`: Example Jsonnet configuration for AnodeDumper, parameterized
  by detector (default: "pdsp"). Demonstrates the configuration pattern.
