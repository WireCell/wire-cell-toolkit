# Pgraph: Potential Bug Examination

Examined: 2026-04-12
Scope: All headers in `inc/WireCellPgraph/` and sources in `src/`

---

## B1. Memory Leak in Factory — Raw `new` Without Ownership (Severity: Medium)

**File:** `inc/WireCellPgraph/Factory.h:34`, `src/Factory.cxx:41`

In `Factory::bind_maker()`:
```cpp
m_factory[cat] = new MakerT<Wrapper>;
```
And in the Maker's `operator()`:
```cpp
virtual Node* operator()(INode::pointer wcnode) { return new Wrapper(wcnode); }
```

Both allocate with raw `new` but the owning maps (`m_factory` of type
`map<NodeCategory, Maker*>` and `m_nodes` of type `map<INode::pointer, Node*>`)
store raw pointers. Neither `Factory`'s destructor nor any cleanup code calls
`delete` on these allocations.

**Impact:** Every `Maker` object (9 total, one per node category) and every
wrapped `Node*` created during the lifetime of a `Factory` instance is leaked
when the `Factory` goes out of scope. In `Pgrapher::configure()` (line 55),
a local `Pgraph::Factory fac` is created, populates `m_factory` with 9 Makers
and `m_nodes` with all wrapped nodes, then goes out of scope without cleanup.

The `Node*` objects survive because `Graph` also stores raw pointers to them,
but the `Maker*` objects in `m_factory` are unconditionally leaked.

**Note:** The `Node*` leak is partially mitigated by the fact that `Graph`
continues to use them, but there is no clear ownership transfer, and if the
`Graph` is destroyed before the nodes, there is no cleanup path at all.

---

## B2. Nfan/SplitFanout: Iterating Output Ports by Value — Data Not Actually Sent (Severity: High)

**File:** `test/test_pipegraph.cxx:191`, `test/doctest_graph.cxx:198`

In the test `Nfan::operator()()`:
```cpp
for (auto p : output_ports()) {
    p.put(obj);
}
```

`output_ports()` returns `PortList&` (a `vector<Port>`). The range-for uses
`auto p` (by value), so each iteration creates a **copy** of the `Port`. The
copy holds a copy of the `shared_ptr<Queue>` edge, so `put()` does push data
onto the shared queue, and this works correctly because `Edge` is a
`shared_ptr<Queue>`.

However, this pattern is still fragile: it creates unnecessary copies of Port
objects (including string members), and if Port's semantics ever changed to
use a non-shared edge or acquire move-only resources, the code would silently
break. The production wrapper `SplitFanout` (Wrappers.h:277) correctly uses
index-based access `oports[ind].put(outv[ind])`, which is the safer pattern.

**Reclassification:** This is not a data-loss bug due to shared_ptr semantics,
but it is a latent correctness risk and inefficiency in the test code.

---

## B3. SplitFanout: Input Consumed Even When Output Ports Are All Full (Severity: Medium)

**File:** `inc/WireCellPgraph/Wrappers.h:249-280`

In `SplitFanout::operator()()`:
```cpp
Port& ip = iport();
if (ip.empty()) {
    return false;  // don't call me if there is not any new input
}

auto& oports = output_ports();
size_t nout = oports.size();

bool full = true;
for (size_t ind = 0; ind < nout; ++ind) {
    if (oports[ind].empty()) {
        full = false;
    }
}
if (full) {
    return false;  // don't call me if all my output has something
}

auto in = ip.get();   // <--- input is CONSUMED here
```

The `full` check only returns false when **all** output ports have data. If
even one output port is empty, execution proceeds and the input is consumed
via `ip.get()`. However, the underlying `m_wcnode` is then called with the
consumed input.

The potential issue: the "full" check should arguably verify that **no** output
port has unconsumed data (i.e., all are empty), not just that at least one is
empty. As written, if some output ports still have data from a previous call but
one has been drained, the node will fire again. Whether this is a bug depends on
the expected semantics of the downstream Split/Fanout node — if it expects to
fire only when all outputs are drained, this is incorrect.

---

## B4. Hydra Wrapper: Fragile Input Queue Synchronization (Severity: Medium-High)

**File:** `inc/WireCellPgraph/Wrappers.h:312-372`

The Hydra wrapper has a documented "BIG FAT WARNING" (line 355-357):

```cpp
// 4) pop dfp input queues to match.  BIG FAT
// WARNING: this trimming assumes caller only
// pop_front's.  Really should hunt for which ones
// have been removed.
for (size_t ind = 0; ind < nin; ++ind) {
    size_t want = inqv[ind].size();
    while (iports[ind].size() > want) {
        iports[ind].get();
    }
}
```

The Hydra copies the entire input edge queue into a local `inqv`, passes it to
the concrete IHydraNode, then uses the **size difference** to determine how many
items were consumed. This assumes the IHydraNode only removes items from the
front (`pop_front`). If the node removes items from the middle or end, or
reorders them, the wrong items will be popped from the actual DFP queue.

Additionally, between steps 0 (checking all inputs non-empty) and step 1
(copying queue contents), and between step 1 (copying) and step 4 (trimming),
the code accesses the edge queues directly via `iports[ind].edge()` (line 329),
bypassing the Port API. While this is single-threaded, it creates a maintenance
hazard if the Port or Edge abstraction ever adds invariants.

---

## B5. Hydra Wrapper: Direct Edge Access Bypasses Port API (Severity: Low)

**File:** `inc/WireCellPgraph/Wrappers.h:329, 367`

```cpp
Edge edge = iports[ind].edge();
// ...
inqv[ind].insert(inqv[ind].begin(), edge->begin(), edge->end());
```

And for output:
```cpp
Edge edge = oports[ind].edge();
edge->insert(edge->end(), outqv[ind].begin(), outqv[ind].end());
```

The Hydra directly manipulates the underlying `deque` of both input and output
edges, bypassing the `Port::get()` and `Port::put()` methods. The `Port::get()`
method validates that the port is an input port and that the edge exists;
`Port::put()` validates the port is an output port. Bypassing these checks means
errors in the Hydra's port setup would not be caught.

---

## B6. Graph::connect — No Bounds Check on Port Indices (Severity: Medium)

**File:** `src/Graph.cxx:31-32`

```cpp
Port& tport = tail->output_ports()[tpind];
Port& hport = head->input_ports()[hpind];
```

This directly indexes into the PortList vectors without bounds checking. If
`tpind` or `hpind` is out of range, this is undefined behavior (vector
out-of-bounds access). The `Node::port()` method (Node.h:30-35) does have
bounds checking and throws `ValueError`, but `connect()` does not use it.

**Impact:** An incorrect jsonnet configuration specifying a port index beyond
the node's actual port count will cause undefined behavior rather than a clear
error message. In `Pgrapher::configure()`, port indices come directly from
user configuration: `int port = get(jone, "port", 0)`.

---

## B7. Graph::sort_kahn — Cycle Detection Missing (Severity: Low-Medium)

**File:** `src/Graph.cxx:59-90`

Kahn's algorithm can detect cycles: if the returned sorted list has fewer nodes
than the total node count, the graph contains a cycle. The implementation does
not perform this check:

```cpp
std::vector<Node*> ret;
// ... Kahn's algorithm ...
return ret;
```

If a user misconfigures a cyclic graph, `sort_kahn()` will silently return a
partial ordering missing the cycle participants. The `execute()` method would
then run only the non-cyclic nodes, producing silently incorrect/incomplete
results without any error or warning.

---

## B8. Graph::sort_kahn — Only Considers Edge-Connected Nodes (Severity: Low)

**File:** `src/Graph.cxx:59-66`

```cpp
std::unordered_map<Node*, int> nincoming;
for (auto th : m_edges) {
    nincoming[th.first] += 0;
    nincoming[th.second] += 1;
}
```

The `nincoming` map is populated solely from `m_edges`. Nodes added via
`add_node()` but not connected to any edge are not included in `nincoming` and
thus not in the topological sort output. While `Graph::connected()` would catch
fully disconnected nodes, there is a window between `add_node()` and full
connectivity checking where `sort_kahn()` could silently drop nodes.

---

## B9. Graph::execute — Unreachable Return Statement (Severity: Cosmetic)

**File:** `src/Graph.cxx:147`

```cpp
while (true) {
    // ...
    if (!did_something) {
        return true;
    }
}
return true;  // shouldn't reach
```

The `return true` after the infinite loop is unreachable. The comment
acknowledges this. Not a runtime bug but a code clarity issue.

---

## B10. Node Instance Counter Not Thread-Safe (Severity: Low)

**File:** `inc/WireCellPgraph/Node.h:12-13`, `src/Node.cxx:3`

```cpp
Node() { m_instance = m_instances++; }
// ...
static size_t m_instances;
```

The static counter `m_instances` is incremented without synchronization. While
pgraph is designed for single-threaded execution, if nodes are ever constructed
from multiple threads (e.g., during parallel configuration), this would be a
data race. The comment in README.org explicitly states "single-threaded", so
this is low severity but worth noting.

---

## B11. Graph::print_timers — Fragile String Parsing for Node Identity (Severity: Low)

**File:** `src/Graph.cxx:190-193`

```cpp
std::string iden = it->second->ident();
std::vector<std::string> tags;
boost::split(tags, iden, [](char c) { return c == ' '; });
l_timer->info("Timer: {} : {} sec", tags[2].substr(5), it->first);
```

This parses the `ident()` string by splitting on spaces and accessing `tags[2]`
with a hardcoded `.substr(5)` to extract the type name. The format of `ident()`
is defined in `Wrappers.h:57-72`:
```
<Node  type:SomeType cat:3 sig:...
```

If `ident()` format changes (e.g., extra or fewer spaces, reordered fields),
`tags[2]` or the `.substr(5)` offset could be wrong, causing an out-of-range
exception or garbled timer output. For non-PortedNode implementations (like test
nodes), the ident format is completely different (e.g., `"src[0]"`), which would
crash `print_timers()`.

---

## B12. execute_upstream — Potential Infinite Recursion on Cycles (Severity: Low)

**File:** `src/Graph.cxx:92-108`

```cpp
int Graph::execute_upstream(Node* node) {
    for (auto parent : m_edges_backward[node]) {
        bool ok = call_node(parent);
        if (ok) { ++count; continue; }
        count += execute_upstream(parent);
    }
    // ...
}
```

If a graph has cycles (which shouldn't happen in a DAG, but there's no
enforcement), `execute_upstream()` will recurse infinitely and overflow the
stack. Combined with B7 (no cycle detection), a misconfigured graph could cause
a crash here.

---

## B13. `dynamic_pointer_cast` Results Not Checked in Wrappers (Severity: Medium)

**File:** `inc/WireCellPgraph/Wrappers.h:88, 116, 139, 170, 200, 246, 308`

All wrapper constructors do:
```cpp
m_wcnode = std::dynamic_pointer_cast<ISourceNodeBase>(wcnode);
```
but never check if the result is null. If the wrong INode type is passed to a
wrapper (e.g., a Sink INode to a Source wrapper), `m_wcnode` will be null and
`(*m_wcnode)(obj)` will dereference a null pointer. The Factory selects
wrappers by category so this should not happen in practice, but there is no
defensive check.
