# Pgraph: Algorithm and Architecture Explanation

Examined: 2026-04-12
Scope: All headers in `inc/WireCellPgraph/` and sources in `src/`

---

## 1. Overview

Pgraph implements a **single-threaded Data Flow Programming (DFP) engine** for
the Wire-Cell Toolkit (WCT). It orchestrates the execution of processing
components (nodes) connected by data-carrying edges, following a directed
acyclic graph (DAG) topology. The design prioritizes low memory usage and
simplicity over parallelism.

The key design principle: data "flows" from sources (producers) through
processing nodes to sinks (consumers). The engine determines execution order
and repeatedly fires nodes that have data to process, until no more work
remains.

---

## 2. Core Abstractions

### 2.1 Data and Type Erasure

All data in the graph is type-erased using `boost::any`:

```
Data = boost::any      (a single data item, can hold any type)
Queue = deque<Data>    (a FIFO buffer of data items)
Edge = shared_ptr<Queue>  (a shared, reference-counted queue)
```

This allows heterogeneous data types (depositions, frames, traces, etc.) to
flow through the same graph infrastructure. Type safety is maintained at the
**port signature** level — each port declares a string signature (a mangled C++
type name), and connections are only allowed between ports with matching
signatures.

### 2.2 Port

A Port mediates between a Node and an Edge. Each Port has:
- **Type**: input (head) or output (tail)
- **Signature**: a C++ type string for connection validation
- **Name**: optional human-readable identifier
- **Edge**: the shared queue this port is connected to

Key operations:
- `put(data)`: push data onto the output port's edge queue
- `get(pop=true)`: pop (or peek at) data from the input port's edge queue
- `plug(edge)`: connect an edge to this port

Ports enforce directionality: `put()` throws if called on an input port;
`get()` throws if called on an output port.

### 2.3 Node

A Node is the unit of computation. It has:
- Zero or more **input ports** and **output ports** (stored in `m_ports[2]`)
- A unique **instance** number (from a static counter)
- A pure virtual `operator()()` — the main execution method, returns true if
  it produced/consumed data, false if it had nothing to do
- A pure virtual `ident()` — returns an identifying string

The `connected()` method checks that all ports have edges plugged in.

### 2.4 Edge

An Edge is simply a `shared_ptr<deque<boost::any>>`. Two ports (one output,
one input) share the same Edge. The output port pushes data onto the back;
the input port pops data from the front. This is a classic producer-consumer
FIFO channel.

The use of `deque` (rather than `queue`) allows iteration, which the Hydra
wrapper exploits to bulk-copy queue contents.

### 2.5 Graph

The Graph class manages the collection of nodes and their connections, and
implements the execution engine. It maintains:
- `m_nodes`: all registered nodes (keyed by instance number)
- `m_edges`: list of (tail, head) node pairs
- `m_edges_forward[node]`: list of downstream nodes
- `m_edges_backward[node]`: list of upstream nodes
- `m_nodes_timer[node]`: cumulative CPU time per node

---

## 3. Graph Construction

### 3.1 Node Registration

Nodes are added via `add_node(Node*)`, which stores the node keyed by its
instance number. Duplicate adds are silently absorbed (same key → overwrites
with same pointer).

### 3.2 Connection

`connect(tail, head, tpind, hpind)` performs:

1. **Signature validation**: checks that the tail's output port signature
   matches the head's input port signature
2. **Edge creation**: creates a new shared Queue
3. **Plugging**: plugs the same Edge into both the tail's output port and the
   head's input port
4. **Registration**: adds both nodes to the graph, records the edge in
   `m_edges`, `m_edges_forward`, and `m_edges_backward`

After all connections, `connected()` verifies every port on every node has an
edge.

---

## 4. Topological Sort (Kahn's Algorithm)

**File:** `src/Graph.cxx:59-90`

Before execution, the graph is topologically sorted using Kahn's algorithm.
This establishes an execution priority order.

### Algorithm Steps:

1. **Compute in-degrees**: For each node, count the number of incoming edges.
   Source nodes (no incoming edges) have in-degree 0.

2. **Initialize seed set**: Collect all nodes with in-degree 0. These are the
   "roots" of the DAG — nodes that can execute without waiting for input from
   other graph nodes. The seeds are stored in a `std::map<size_t, Node*>`
   keyed by instance number, providing **deterministic ordering** among
   equal-priority nodes.

3. **Iterative processing**: Repeatedly:
   - Remove a node from the seed set (smallest instance number first)
   - Append it to the sorted output
   - For each downstream neighbor, decrement its in-degree
   - If any neighbor's in-degree reaches 0, add it to the seed set

4. **Return**: The sorted vector, ordered from sources to sinks.

**Determinism note:** Using `std::map` for seeds (rather than a set or queue)
ensures that among nodes with the same topological priority, the one with the
smallest instance number is processed first. This makes the sort deterministic
across runs, which is important for reproducible scientific computations.

---

## 5. Execution Model

**File:** `src/Graph.cxx:111-148`

### 5.1 Main Execution Loop

The execution follows a "demand-driven, bottom-up" strategy:

```
sorted = sort_kahn()  // source → ... → sink order

loop forever:
    for each node in REVERSE sorted order (sinks first):
        if call_node(node) succeeds:
            break and restart from sinks
    if no node fired:
        done
```

Key design choices:

1. **Reverse iteration (sinks first)**: The loop iterates from sinks toward
   sources. This is demand-driven: it tries to drain data from the bottom of
   the pipeline first, which keeps memory usage low by consuming buffered data
   before producing more.

2. **Single-fire-then-restart**: When any node successfully fires (returns
   true), the loop immediately restarts from the sinks. This ensures that
   downstream consumers get a chance to drain data as soon as it's produced,
   preventing queue buildup.

3. **Termination**: When a full scan from sinks to sources finds no node that
   can fire, the graph has quiesced — all sources are exhausted and all
   queues are drained.

### 5.2 Node Firing Semantics

A node's `operator()()` returns:
- **true**: the node did useful work (consumed input and/or produced output)
- **false**: the node could not do anything (no input available, output buffer
  full, or source exhausted)

The wrappers implement specific readiness logic:
- **Source**: fires if output port is empty (backpressure)
- **Sink**: fires if input port has data
- **Function**: fires if input has data AND output is empty
- **Queuedout**: fires if input has data (may produce multiple outputs)
- **JoinFanin**: fires if ALL input ports have data AND output is empty
- **SplitFanout**: fires if input has data AND not all outputs are full
- **Hydra**: fires if ALL input ports have data

### 5.3 execute_upstream (Alternative Execution Path)

`execute_upstream(node)` is a recursive, pull-based execution mode:
- For each parent of the target node, try to call it
- If a parent can't fire, recursively execute ITS parents first
- After processing parents, try to call the target node

This provides a way to "pull" data through the graph on demand for a specific
node, rather than running the full graph. It is not used in the main
`Pgrapher::execute()` path but is available for specialized use cases.

---

## 6. Wrapper and Factory Pattern

### 6.1 The Problem

WCT components implement typed interfaces (e.g., `IDepoSource`,
`IFrameFilter`), but the Pgraph engine operates at the type-erased
`boost::any` level. The wrapper pattern bridges this gap.

### 6.2 Wrapper Hierarchy

```
Pgraph::Node (abstract base)
  |
  +-- PortedNode (handles port construction and ident() from INode)
        |
        +-- Source      (wraps ISourceNodeBase)
        +-- Sink        (wraps ISinkNodeBase)
        +-- Function    (wraps IFunctionNodeBase)
        +-- Queuedout   (wraps IQueuedoutNodeBase)
        +-- JoinFanin<T> (template, wraps IJoinNodeBase or IFaninNodeBase)
        |     +-- Join  = JoinFanin<IJoinNodeBase>
        |     +-- Fanin = JoinFanin<IFaninNodeBase>
        +-- SplitFanout<T> (template, wraps ISplitNodeBase or IFanoutNodeBase)
        |     +-- Split  = SplitFanout<ISplitNodeBase>
        |     +-- Fanout = SplitFanout<IFanoutNodeBase>
        +-- Hydra       (wraps IHydraNodeBase)
```

**PortedNode** is the key base class. Its constructor:
1. Takes an `INode::pointer`
2. Queries `input_types()` and `output_types()` from the INode interface
3. Creates the corresponding Pgraph::Port objects with matching type signatures

Each concrete wrapper:
1. Downcasts the INode to the specific base type (e.g., `ISourceNodeBase`)
2. Implements `operator()()` with type-erasure bridging:
   - Gets `boost::any` data from input ports
   - Calls the underlying INode's operator with the type-erased data
   - Puts the `boost::any` result onto output ports

### 6.3 Factory

The Factory maps INode categories to Maker objects:

```
Factory constructor:
    bind_maker<Source>(INode::sourceNode)
    bind_maker<Sink>(INode::sinkNode)
    bind_maker<Function>(INode::functionNode)
    ... (9 categories total)

Factory::operator()(INode::pointer):
    if already wrapped → return cached wrapper
    look up Maker by node category
    create wrapper via Maker
    cache and return
```

The caching ensures each INode gets exactly one Pgraph::Node wrapper, even if
the same INode appears in multiple edge definitions.

### 6.4 Node Category Summary

| Category    | Inputs | Outputs | Semantics                              |
|-------------|--------|---------|----------------------------------------|
| Source      | 0      | 1       | Generates data items                   |
| Sink        | 1      | 0       | Consumes data items                    |
| Function    | 1      | 1       | Transforms one input to one output     |
| Queuedout   | 1      | 1       | One input may produce multiple outputs |
| Join        | N      | 1       | Merges N typed inputs into one output  |
| Fanin       | N      | 1       | Merges N same-type inputs into one     |
| Split       | 1      | N       | Splits one input into N typed outputs  |
| Fanout      | 1      | N       | Copies one input to N same-type outputs|
| Hydra       | N      | M       | General N-to-M with queue access       |

---

## 7. Configuration via Pgrapher

### 7.1 Jsonnet Configuration

The Pgrapher is configured via Jsonnet with this structure:

```jsonnet
{
    type: "Pgrapher",
    data: {
        edges: [
            { tail: {node: "Type:name", port: 0},
              head: {node: "Type:name", port: 0} },
            ...
        ],
        verbosity: 1  // 0=silent, 1=timer summary, 2=ExecMon tracing
    }
}
```

Each edge specifies a tail (output side) and head (input side), each with a
node typename and optional port index (default 0).

### 7.2 Configuration Flow

`Pgrapher::configure()`:
1. Creates a local `Pgraph::Factory`
2. For each edge in the configuration:
   a. Looks up the tail and head INode components by typename (they must
      already exist in the WCT named factory registry)
   b. Wraps them via the Pgraph::Factory (creating or reusing wrappers)
   c. Connects them in the Graph with the specified port indices
3. Verifies the graph is fully connected (all ports have edges)

### 7.3 Execution Flow

`Pgrapher::execute()`:
1. Calls `m_graph.execute()` which:
   a. Topologically sorts the graph (Kahn's algorithm)
   b. Runs the main execution loop (demand-driven, bottom-up)
2. Prints timing information based on verbosity level

---

## 8. Data Flow Example

Consider this simple graph from `test_pgraph.jsonnet`:

```
TrackDepos:cosmics ──(port 0)──> DepoMerger ──> DumpDepos
TrackDepos:beam ────(port 1)──/
```

Execution proceeds as follows:

1. **Sort**: Kahn's algorithm produces: [cosmics, beam, merger, dump]
2. **Execute** (reverse order iteration):
   - Try DumpDepos (sink) — input empty → false
   - Try DepoMerger (join) — both inputs empty → false
   - Try beam (source) — output empty → produces depo → **true** → restart
   - Try DumpDepos — input empty → false
   - Try DepoMerger — only one input has data → false (needs both)
   - Try cosmics (source) — output empty → produces depo → **true** → restart
   - Try DumpDepos — input empty → false
   - Try DepoMerger — both inputs have data → merges, produces output → **true**
   - Try DumpDepos — input has data → consumes → **true**
   - ... (continues until sources are exhausted)
3. **Terminate**: When no node can fire, execution ends.

This demand-driven approach naturally handles backpressure: sources only produce
when their output queue is empty, preventing unbounded memory growth.

---

## 9. Timing and Profiling

The graph provides two levels of performance monitoring:

### 9.1 Per-Node CPU Timing (verbosity >= 1)

Every `call_node()` invocation is bracketed with `std::clock()` calls. The
cumulative CPU time per node is stored in `m_nodes_timer`. After execution,
`print_timers()` reports nodes sorted by descending CPU time.

### 9.2 ExecMon Tracing (verbosity == 2)

When enabled, each successful node call records an `ExecMon` event with
detailed memory and timing information. The `ExecMon::summary()` output
provides a chronological trace of the execution for debugging.

---

## 10. Comparison with TbbFlow

Per the README, pgraph is one of two DFP engines in WCT:
- **Pgraph**: Single-threaded, low-memory, simple scheduling
- **TbbFlow** (in `tbb/` sub-package): Multi-threaded via Intel TBB

Pgraph's single-threaded design means:
- No synchronization overhead on queues or node execution
- Deterministic execution order (important for reproducibility)
- Lower peak memory (demand-driven consumption)
- No parallelism (cannot exploit multi-core hardware)

The tradeoff is appropriate for workloads where I/O or single-node computation
dominates, or where deterministic reproducibility is required.
