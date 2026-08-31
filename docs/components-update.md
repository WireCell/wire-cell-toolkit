# Component Documentation Update Prompt

This file is an LLM prompt for updating `docs/components.md` and
`docs/components/NAME.md` when Wire-Cell Toolkit components are added or
modified.  Read it in full before acting.

---

## Your task

You will be given one or more component source files to document (or re-document).
For each component you must produce a **block** in the format shown below and then
call the update script to write the files.

---

## Codebase layout

```
toolkit/
  iface/inc/WireCellIface/   interface headers (INode hierarchy)
  <pkg>/inc/WireCell<Pkg>/   component class headers
  <pkg>/src/                 component class sources (.cxx)
  docs/components.md         summary table (one row per component)
  docs/components/NAME.md    per-component detail page
  docs/components-update.py  maintenance script
```

A component is any class that:
1. Is registered with `WIRECELL_FACTORY(NAME, ConcreteClass, Interface1, ...)` in a `.cxx` file.
2. Includes at least one `INode`-derived interface (see table below) in the macro arguments.

---

## How to find the files for a component

Given a source file path `<pkg>/src/Foo.cxx` containing
`WIRECELL_FACTORY(Foo, WireCell::Bar::Foo, WireCell::IFrameFilter, ...)`:

- **Source**: `<pkg>/src/Foo.cxx`
- **Header**: `<pkg>/inc/WireCell<Pkg>/Foo.h`  (capitalise the package name, e.g. `gen` → `WireCellGen`)
- **Interface header**: `iface/inc/WireCellIface/IFrameFilter.h`

Read both the header and the source.  Multiple components can live in one `.cxx` file;
read all `WIRECELL_FACTORY` lines to identify them.

---

## What to extract

### Description
One sentence describing what the component *does* operationally.  Focus on the
data transformation, not the class hierarchy.  Do not start with "This component".

### Node category and types
Derive these from the primary `INode`-derived interface using the table below.
Do not guess; if an interface is not in the table, report it as `unknown`.

| Interface | Category | Input | Output |
| --- | --- | --- | --- |
| IFrameFilter | function | IFrame | IFrame |
| IFrameSource | source | (none) | IFrame |
| IFrameSink | sink | IFrame | (none) |
| IFrameFanin | fanin | vector\<IFrame\> | IFrame |
| IFrameFanout | fanout | IFrame | vector\<IFrame\> |
| IFrameJoiner | join | (IFrame, IFrame) | IFrame |
| IFrameMerge | hydra | vector\<IFrame\> | IFrame |
| IFrameSlicer | queuedout | IFrame | ISlice |
| IFrameSlices | function | IFrame | ISliceSet |
| IFrameSplitter | split | IFrame | (IFrame, IFrame) |
| IFrameTensorSet | function | IFrame | ITensorSet |
| IDepoFilter | function | IDepo | IDepo |
| IDepoSource | source | (none) | IDepo |
| IDepoSink | sink | IDepo | (none) |
| IDepoFanout | fanout | IDepo | vector\<IDepo\> |
| IDepoMerger | hydra | (IDepo, IDepo) | IDepo |
| IDepoCollector | queuedout | IDepo | IDepoSet |
| IDepoFramer | function | IDepoSet | IFrame |
| IDepoSetFilter | function | IDepoSet | IDepoSet |
| IDepoSetSource | source | (none) | IDepoSet |
| IDepoSetSink | sink | IDepoSet | (none) |
| IDepoSetFanin | fanin | vector\<IDepoSet\> | IDepoSet |
| IDepoSetFanout | fanout | IDepoSet | vector\<IDepoSet\> |
| IDepos2DeposOrFrame | hydra | IDepoSet | IDepoSet or IFrame |
| IClusterFilter | function | ICluster | ICluster |
| IClusterSource | source | (none) | ICluster |
| IClusterSink | sink | ICluster | (none) |
| IClusterFanin | fanin | vector\<ICluster\> | ICluster |
| IClusterFanout | fanout | ICluster | vector\<ICluster\> |
| IClusterFramer | function | ICluster | IFrame |
| IClusterTensorSet | function | ICluster | ITensorSet |
| IClusterFaninTensorSet | fanin | vector\<ICluster\> | ITensorSet |
| IBlobSetProcessor | function | IBlobSet | IBlobSet |
| IBlobSetSink | sink | IBlobSet | (none) |
| IBlobSetSource | source | (none) | IBlobSet |
| IBlobSetFanin | fanin | vector\<IBlobSet\> | IBlobSet |
| IBlobSetFanout | fanout | IBlobSet | vector\<IBlobSet\> |
| IBlobSetFramer | function | IBlobSet | IFrame |
| IBlobDeclustering | function | ICluster | IBlobSet |
| IBlobDepoFill | join | (ICluster, IDepoSet) | ICluster |
| IBlobSampling | function | IBlobSet | ITensorSet |
| IBlobTensoring | queuedout | IBlobSet | ITensorSet |
| ITensorSetFilter | function | ITensorSet | ITensorSet |
| ITensorSetSource | source | (none) | ITensorSet |
| ITensorSetSink | sink | ITensorSet | (none) |
| ITensorSetFanin | fanin | vector\<ITensorSet\> | ITensorSet |
| ITensorSetFrame | function | ITensorSet | IFrame |
| ITensorSetCluster | function | ITensorSet | ICluster |
| ITensorSetUnpacker | fanout | ITensorSet | vector\<ITensor\> |
| ITensorPacker | fanin | vector\<ITensor\> | ITensorSet |
| ISliceFanout | fanout | ISlice | vector\<ISlice\> |
| ISliceFrameSink | sink | (ISlice, IFrame) | (none) |
| ISliceStriper | function | ISlice | IStripeSet |
| ITiling | function | ISlice | IBlobSet |
| IWireGenerator | function | IWireParameters | IWire::vector |
| IWireSource | source | (none) | IWire::vector |
| IWireSummarizer | function | IWire::vector | IWireSummary |
| IClustering | queuedout | IBlobSet | ICluster |
| IDuctor | queuedout | IDepo | IFrame |
| IDiffuser | queuedout | IDepo | IDiffusion |
| IDrifter | queuedout | IDepo | IDepo |

### Configuration parameters
A component is configurable if `IConfigurable` appears in its `WIRECELL_FACTORY` arguments.
If so, read the `default_configuration()` method body — every key set there is a parameter.
Also read the `configure()` method to understand what each key controls.

For each parameter write: `param_name: one-line description (default: value)`.
If there is no obvious default, write `(required)` instead.

---

## Output block format

Produce exactly one block per component.  Separate multiple blocks with `---`.
Do not add any text outside the blocks.

```
NAME: Foo
CONCRETE: WireCell::Bar::Foo
CATEGORY: function
INPUT: IFrame
OUTPUT: IFrame
IFACE: IFrameFilter
DESCRIPTION: Applies a spatial smoothing filter to every trace in a frame.
CONFIGURABLE: yes
CONFIG_PARAMS:
sigma: Gaussian smoothing width in ticks (default: 3.0)
tag: input trace tag to process; empty means all traces (default: "")
---
NAME: FooSink
CONCRETE: WireCell::Bar::FooSink
CATEGORY: sink
INPUT: IFrame
OUTPUT: (none)
IFACE: IFrameSink
DESCRIPTION: Writes each frame to a binary file for offline inspection.
CONFIGURABLE: yes
CONFIG_PARAMS:
filename: output file path (required)
---
```

Rules:
- `NAME` is the first argument to `WIRECELL_FACTORY` exactly as written.
- `CONCRETE` is the second argument exactly as written (may include namespace).
- `CATEGORY`, `INPUT`, `OUTPUT` come from the interface table above.
- `IFACE` is the bare interface class name (no namespace prefix).
- `DESCRIPTION` is a single sentence, no trailing period.
- `CONFIGURABLE` is `yes` or `no`.
- `CONFIG_PARAMS` lists one parameter per line as `name: description`.
  Omit the section entirely if `CONFIGURABLE: no`.
- If a parameter description would use a pipe character `|`, replace it with `or`.

---

## Workflow

### Adding a new component

1. Run the inventory to confirm the component is new:
   ```
   python3 docs/components-update.py inventory
   ```

2. Read the component's `.cxx` and `.h` files and produce a block (see format above).

3. Save the block to a temporary file, then write the docs:
   ```
   python3 docs/components-update.py write /tmp/block.txt
   ```
   This writes `docs/components/NAME.md` and regenerates `docs/components.md`.

### Updating an existing component

Same as above.  The `write` subcommand overwrites the existing `NAME.md` and
regenerates the table row.

### Removing a component

1. Delete the file:
   ```
   rm docs/components/NAME.md
   ```
2. Rebuild the table:
   ```
   python3 docs/components-update.py regen
   ```

### Rebuilding the table only (no content changes)

If only `docs/components.md` is stale (e.g. after a merge):
```
python3 docs/components-update.py regen
```

### Detecting what changed after a git pull or merge

```
python3 docs/components-update.py inventory
```

This reports:
- **NEW** — components present in source but missing from docs
- **REMOVED** — docs entries with no corresponding `WIRECELL_FACTORY`
- **STRUCTURALLY CHANGED** — interface or class name differs between source and docs
- **SOURCE MODIFIED SINCE HEAD~1** — source files touched in the last commit

Process NEW and STRUCTURALLY CHANGED components by reading their source and
producing blocks, then run `write` for each.  Process REMOVED by deleting the
`.md` file and running `regen`.

---

## Example session

A developer adds `toolkit/gen/src/MyFilter.cxx` containing:
```cpp
WIRECELL_FACTORY(MyFilter, WireCell::Gen::MyFilter,
                 WireCell::IFrameFilter, WireCell::IConfigurable)
```

You would:
1. Read `gen/src/MyFilter.cxx` and `gen/inc/WireCellGen/MyFilter.h`.
2. Identify: NAME=MyFilter, CONCRETE=WireCell::Gen::MyFilter,
   IFACE=IFrameFilter → category=function, input=IFrame, output=IFrame,
   CONFIGURABLE=yes.
3. Read `default_configuration()` and `configure()` for parameters.
4. Produce and save the block, then run:
   ```
   python3 docs/components-update.py write /tmp/block.txt
   ```
