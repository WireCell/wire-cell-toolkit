# Pgraph: Efficiency Examination

Examined: 2026-04-12
Scope: All headers in `inc/WireCellPgraph/` and sources in `src/`

---

## E1. Graph::execute — Restart-from-Bottom Strategy (Impact: Medium-High)

**File:** `src/Graph.cxx:111-148`

The main execution loop in `Graph::execute()` uses a restart-from-bottom
approach:

```cpp
while (true) {
    int count = 0;
    bool did_something = false;

    for (auto nit = nodes.rbegin(); nit != nodes.rend(); ++nit, ++count) {
        Node* node = *nit;
        bool ok = call_node(node);
        // ...
        if (ok) {
            did_something = true;
            break;  // start again from bottom of graph
        }
    }
    if (!did_something) {
        return true;
    }
}
```

Every time a node successfully executes, the inner loop `break`s and restarts
from the bottom (sinks) of the topological order. For a graph with N nodes and
M total data items to process, this results in O(N * M) node-call attempts in
the worst case. Each restart scans from the bottom, calling nodes that may have
nothing to do (returning false), before reaching one that can fire.

**Specific concern:** For large graphs (dozens of nodes) processing many data
items (thousands of events), the overhead of repeatedly scanning non-ready
nodes adds up. Each `call_node()` invocation that returns false still costs
a virtual function call plus whatever readiness checking the wrapper does.

**Potential improvement:** A ready-queue approach could avoid re-scanning: when
a node produces output, enqueue its downstream consumers. This would reduce
the scheduling overhead from O(N) per fired node to O(1) amortized.

---

## E2. boost::any Type Erasure Overhead (Impact: Medium)

**File:** `inc/WireCellPgraph/Port.h:17`, used throughout all wrappers

All data flowing through the graph is wrapped in `boost::any`:

```cpp
typedef boost::any Data;
```

Every `put()` and `get()` operation involves `boost::any` copy semantics.
For small types, `boost::any` may heap-allocate (depending on the
implementation's small-object optimization threshold). For large objects like
`shared_ptr<IFrame>`, the `boost::any` copy copies the shared_ptr (incrementing
the reference count), which involves an atomic operation.

In the hot path (the main execution loop), every data item passes through:
1. `boost::any` construction in the wrapper's `operator()()`
2. `push_back` into the queue (may copy the `any`)
3. `front()` + `pop_front()` in the consumer (another copy + destroy)

**Note:** Modern alternatives like `std::any` (C++17) or `std::variant` could
reduce overhead, but `boost::any` is deeply embedded in the INode interface.
Move semantics on `boost::any` could help if the queue operations used
`std::move`, but they currently don't — `Port::put()` takes `Data&` (lvalue
reference), preventing move:

```cpp
void Port::put(Data& data);  // takes lvalue ref, not Data&&
```

---

## E3. Port::put Takes Lvalue Reference — Prevents Move Semantics (Impact: Medium)

**File:** `inc/WireCellPgraph/Port.h:53`, `src/Port.cxx:73-82`

```cpp
void put(Data& data);
```

Taking `Data&` (non-const lvalue reference) prevents callers from moving
temporaries into the queue. In multiple wrappers, the data is a local variable
that will not be used again:

```cpp
// Source wrapper (Wrappers.h:99-105)
boost::any obj;
m_ok = (*m_wcnode)(obj);
oport().put(obj);  // obj could be moved
```

```cpp
// Function wrapper (Wrappers.h:152-159)
boost::any out;
auto in = ip.get();
bool ok = (*m_wcnode)(in, out);
op.put(out);  // out could be moved
```

If `put()` accepted `Data&&` (or used a forwarding reference), the final copy
into the deque could be replaced with a move, avoiding potential heap allocation
for the `boost::any` copy.

---

## E4. PortedNode::ident() Creates Stringstream Per Call (Impact: Low-Medium)

**File:** `inc/WireCellPgraph/Wrappers.h:56-73`

```cpp
virtual std::string ident() {
    std::stringstream ss;
    ss << "<Node " << " type:" << WireCell::type(*(m_wcnode.get()))
       << " cat:" << m_wcnode->category()
       << " sig:" << demangle(m_wcnode->signature());
    // ... more formatting ...
    return ss.str();
}
```

`ident()` is called in `Graph::call_node()` (via SPDLOG_LOGGER_TRACE) and
in `Graph::print_timers()`. The TRACE calls are typically compiled out in
release builds, but if TRACE is enabled, `ident()` is called on every
`call_node()` invocation (potentially thousands of times per event). Each call
constructs a `std::stringstream`, performs multiple `demangle()` calls (which
may allocate), and returns a new string.

**Potential improvement:** Cache the ident string in the constructor or compute
it lazily once, since the node's identity doesn't change after construction.

---

## E5. Graph::print_timers — String Parsing Instead of Structured Data (Impact: Low)

**File:** `src/Graph.cxx:180-201`

```cpp
std::string iden = it->second->ident();
std::vector<std::string> tags;
boost::split(tags, iden, [](char c) { return c == ' '; });
l_timer->info("Timer: {} : {} sec", tags[2].substr(5), it->first);
```

The timer printing reconstructs the node's type name by parsing the `ident()`
string. This involves:
- Calling `ident()` (see E4 — creates stringstream, demangling, etc.)
- `boost::split` to tokenize by spaces (allocates a vector of strings)
- `substr(5)` to strip the "type:" prefix

This is only called once at the end of execution so the runtime impact is
minimal, but it's worth noting as a fragile and wasteful approach compared to
storing the node type name as structured data.

---

## E6. sort_kahn — Copies Edge Pairs (Impact: Low)

**File:** `src/Graph.cxx:62`

```cpp
for (auto th : m_edges) {
    nincoming[th.first] += 0;
    nincoming[th.second] += 1;
}
```

`m_edges` is `vector<pair<Node*, Node*>>`. The range-for copies each pair.
Using `const auto&` would avoid the copies. The pairs are small (two pointers),
so the impact is negligible, but it's a general pattern issue that appears in
multiple places.

---

## E7. Graph::execute — std::clock for Timing (Impact: Low)

**File:** `src/Graph.cxx:127-131`

```cpp
auto start = std::clock();
bool ok = call_node(node);
m_nodes_timer[node] += (std::clock() - start) / (double) CLOCKS_PER_SEC;
```

`std::clock()` measures CPU time, not wall time. This is actually correct for
profiling CPU-bound processing. However, `std::clock()` has low resolution on
many systems (often 1ms granularity on Linux with CLOCKS_PER_SEC=1000000 but
actual resolution ~1ms). For very fast nodes, the measured time may be
dominated by the timer overhead itself.

More importantly, `std::clock()` is called on **every** `call_node()` attempt,
including the many attempts that return false (node not ready). The timing
accumulation is only meaningful for `ok == true` calls, but the clock overhead
is paid for every attempt. Moving the clock calls inside the `if (ok)` block
would reduce overhead, though it would slightly misattribute time spent in
failed readiness checks.

---

## E8. ExecMon Always Constructed (Impact: Low)

**File:** `inc/WireCellPgraph/Graph.h:80`

```cpp
ExecMon m_em;
```

The `ExecMon` object is always constructed as a member of `Graph`, regardless
of whether `m_enable_em` is true. If `ExecMon` construction is expensive (e.g.,
it captures initial system state), this is wasted work for the common case of
verbosity < 2.

**Potential improvement:** Use `std::optional<ExecMon>` or a pointer, and only
construct it when `m_enable_em` is set to true.

---

## E9. Hydra Wrapper — Full Queue Copy for Input Processing (Impact: Medium)

**File:** `inc/WireCellPgraph/Wrappers.h:327-338`

```cpp
IHydraNodeBase::any_queue_vector inqv(nin);
for (size_t ind = 0; ind < nin; ++ind) {
    Edge edge = iports[ind].edge();
    // ...
    inqv[ind].insert(inqv[ind].begin(), edge->begin(), edge->end());
}
```

The Hydra copies the **entire contents** of each input edge queue into a local
queue vector. For nodes with large input buffers (e.g., buffering many frames),
this copies every `boost::any` item, each potentially involving a heap
allocation for the `boost::any` copy.

After calling the IHydraNode, it then pops items from the original queue to
match the consumed amount (line 358-362). The copy-then-diff approach is
inherently wasteful compared to passing the edge queues directly or using
move semantics.

---

## E10. Pgrapher::configure — Local Factory Object (Impact: Negligible)

**File:** `src/Pgrapher.cxx:55`

```cpp
Pgraph::Factory fac;
```

A new `Factory` is created every time `configure()` is called. This constructs
9 `Maker` objects via `new`. Since `configure()` is typically called once, this
is not a performance issue, but combined with B1 (memory leak), repeated
`configure()` calls would leak more memory each time.

---

## E11. m_nodes Uses std::map Instead of unordered_map (Impact: Negligible)

**File:** `inc/WireCellPgraph/Graph.h:74`

```cpp
std::map<size_t, Node*> m_nodes;
```

There is a commented-out `unordered_set<Node*>` above it. The `std::map`
provides O(log N) lookup vs O(1) for `unordered_map`. Given that `m_nodes` is
primarily iterated (in `connected()`) and looked up during `add_node()`, and
typical graph sizes are small (tens of nodes), the difference is negligible.
The `std::map` does provide deterministic iteration order, which may have been
the reason for the change.
