# Wire & Anode Geometry -- Code Examination

## Summary

This group of 11 files implements the wire-cell detector geometry hierarchy:
**AnodePlane > AnodeFace > WirePlane**, plus wire generation, parameterization,
summarization, schema loading, spatial filtering, and sourcing utilities.

### Key Findings

| ID | File | Severity | Description |
|----|------|----------|-------------|
| B1 | AnodePlane.cxx | Bug | `response_x`/`anode_x` read from null JSON when `sensitive_face=false` |
| B2 | AnodePlane.cxx | Bug | `channel_wire_collector_t` leaks `SimpleChannel` objects on error paths |
| B3 | AnodePlane.cxx | Minor | Pitch extension logic ignores oblique pitch directions |
| B4 | WireSummary.cxx | Bug | Memory leak: `m_cache` is never deleted in destructor |
| B5 | WireSummary.cxx | Bug | `closest()` truncates index (should round), off-by-one possible |
| B6 | WireGenerator.cxx | Bug | Raw `new GenWire` / `new IWire::vector` -- ownership unclear on error |
| B7 | WireBoundedDepos.cxx | Bug | Config reads `cfg["regions"]` but `default_configuration` names it `"wires"` |
| B8 | WirePlane.h | Bug | Raw `Pimpos*` pointer -- never deleted, memory leak |
| E1 | AnodePlane.cxx | Efficiency | Linear scan in `face(int ident)` -- trivially O(n), n small |
| E2 | MegaAnodePlane.cxx | Efficiency | `channels()` copies and concatenates all sub-anode channels each call |
| E3 | WireSummary.cxx | Efficiency | `by_chan()` uses `operator[]` which inserts empty entry on miss |
| E4 | AnodeFace.cxx | Minor | `get_raypairs` uses only first two wires -- fragile if ordering wrong |

---

## File-by-File Analysis

---

### AnodePlane.h / AnodePlane.cxx

**Algorithm overview:**
AnodePlane is the top-level configurable geometry component. It reads a
`WireSchema::Store` (via an `IWireSchema` component) and configuration for
face locations (response, anode, cathode X-planes). For each face it constructs
`WirePlane` objects (with `Pimpos` impact-position maps), collects channels,
and builds per-plane sensitive bounding boxes that are intersected to form
the overall face sensitive volume.

**Potential bugs:**

**(B1) Reading JSON values from null face (lines 169-184):**
When `sensitive_face = false` (line 162), the code still executes:
```cpp
const double response_x = jface["response"].asDouble();  // line 169
const double anode_x = get(jface, "anode", response_x);  // line 170
```
`jface` is `Json::nullValue` here. `Json::nullValue.operator[]("response")`
returns another null, and `asDouble()` on null returns 0.0, which is silently
used. While the values are not used for `bbvols` (guarded by `sensitive_face`),
`dirx` at line 186 and the `Pimpos` origin at line 234 still use `response_x`.
This means a non-sensitive face gets `response_x = 0` and `pimpos_origin.x() = 0`,
which may produce incorrect pitch distance calculations if that face's wire
planes are ever queried.

**(B2) Memory management in `channel_wire_collector_t` (lines 78-102):**
The struct allocates `SimpleChannel` objects with raw `new` (line 92). The
comment says "This struct explicitly does NOT free its new channels as it is
assumed they will all become wrapped in IChannel::pointers." However, if
`configure()` throws an exception between channel creation and wrapping into
`shared_ptr` at line 253, those `SimpleChannel` objects leak. Additionally,
channels for wires with `segment() > 0` (line 245) that have no segment-0
wire will have a `SimpleChannel` allocated via line 92 but never wrapped in
a `shared_ptr`.

**(B3) Pitch extension only handles axis-aligned pitch (lines 289-292):**
```cpp
if (std::fabs(pitch_dir.z()) > 0.999)
    pext = Point(0, 0, 0.5 * mean_pitch);
else if (std::fabs(pitch_dir.y()) > 0.999)
    pext = Point(0, 0.5 * mean_pitch, 0);
```
If the pitch direction is oblique (e.g., at 60 degrees for U/V planes), neither
condition is satisfied and `pext` remains `(0,0,0)` -- the bounding box gets no
half-pitch extension for those planes. The comment says "FIXME: probably need a
better comparison". This means the sensitive volume intersection may be slightly
too tight for the U and V wire planes.

**Efficiency concerns:**

**(E1) Linear scan in `face()` (lines 326-334):**
`face(int ident)` does a linear search over `m_faces`. With typically 2 faces
this is fine, but the pattern is worth noting.

**Key algorithmic details:**

- Wire plane ordering is assumed U/V/W via `iplane2layer[]` (line 197).
- The sensitive volume for a face is the *intersection* of per-plane bounding
  boxes (lines 304-312), not the union. This ensures only the region covered
  by all three planes is considered sensitive.
- `dirx` (line 186) encodes which direction is "into" the drift volume:
  `+1` if response is at larger X than anode, `-1` otherwise.
- Channels are deduplicated via sort+unique (lines 321-323).

---

### AnodeFace.h / AnodeFace.cxx

**Algorithm overview:**
AnodeFace stores the wire planes, sensitive bounding box, and constructs a
`RayGrid::Coordinates` object from ray pairs. The ray pairs comprise two
boundary pairs (horizontal Y-direction, vertical Z-direction) plus one pair
per wire plane derived from the first two wires of each plane.

**Potential bugs:**

**(E4) `get_raypairs` uses only first two wires (lines 36-37):**
```cpp
const auto wray0 = wires[0]->ray();
const auto wray1 = wires[1]->ray();
```
This assumes `wires[0]` and `wires[1]` exist and are properly ordered. If a
plane has only one wire (degenerate case), this will crash. No bounds check.

The pitch is computed from these two wires and used to offset `r1` and `r2`
by half a pitch (lines 40-41). This is correct for uniform-pitch planes but
assumes the first two wires are representative of the full plane's pitch.

**Key algorithmic details:**

- The `RayGrid::Coordinates` is initialized with 2 + N_planes ray pairs
  (typically 5 total: 2 bounds + 3 wire planes).
- The bounding box X-coordinate is zeroed (lines 18-21), projecting the
  sensitive volume to the Y-Z plane for the ray grid.
- The commented-out `dirx()` method (lines 72-83) shows a previous attempt
  to derive drift direction from wire geometry -- abandoned in favor of
  explicit configuration.

---

### MegaAnodePlane.h / MegaAnodePlane.cxx

**Algorithm overview:**
MegaAnodePlane aggregates multiple `IAnodePlane` instances to provide a
unified channel-resolution interface. It delegates `resolve()`, `channel()`,
`wires()`, and `channels()` by iterating over sub-anodes. Face-related
methods return dummies (empty vectors, nullptr).

**Potential bugs:**

- `ident()` always returns `-1` (line 39) and `faces()` returns empty (line 41).
  These are documented as dummies, but callers checking `ident() >= 0` will
  incorrectly conclude this is an invalid anode.
- `nfaces()` (line 36) indexes `m_anodes[0]` without checking if `m_anodes`
  is empty. If configured with no anodes, this is undefined behavior.

**Efficiency concerns:**

**(E2) `channels()` concatenates on every call (lines 49-57):**
```cpp
std::vector<int> ret;
for (auto& anode : m_anodes) {
    auto chans = anode->channels();
    ret.insert(ret.end(), chans.begin(), chans.end());
}
return ret;
```
This allocates and copies on every call. If called frequently (e.g., in a loop),
the result should be cached. Also, unlike `AnodePlane::channels()`, the result
is not deduplicated -- if anodes share channels, duplicates will appear.

---

### WirePlane.h / WirePlane.cxx

**Algorithm overview:**
Simple data holder for a wire plane: stores ident, a `Pimpos*` pointer, wire
vector, and channel vector.

**Potential bugs:**

**(B8) Raw `Pimpos*` is never freed (line 33 of .h):**
`WirePlane` stores `Pimpos* m_pimpos` as a raw pointer. The `Pimpos` is
allocated with `new` in `AnodePlane::configure()` (line 239). `WirePlane`'s
destructor (line 13 of .cxx) is empty -- it does not delete `m_pimpos`. Since
`AnodePlane` doesn't track or free it either, this is a memory leak every time
`configure()` is called (the old `Pimpos` objects from a previous configuration
are lost).

- `m_bbox` member (line 38 of .h) is declared but never initialized or used
  anywhere. Dead code.

---

### WireParams.h / WireParams.cxx

**Algorithm overview:**
Provides configurable wire parameters (bounding box, pitch vectors, angles) for
an idealized triple-plane wire geometry. Two `set()` overloads: one fully
explicit, one simplified assuming symmetric U/V angles and uniform pitch.

**Potential bugs:**

- No validation that pitch values are positive. Zero or negative pitch will
  produce degenerate cross products and zero-length pitch vectors.

**Key algorithmic details:**

- Pitch direction is computed as `xaxis.cross(wireDir)` (lines 74-76, 124-126).
  For W wires along Y, this gives pitch along Z. For U/V wires at +/-60 degrees,
  the pitch is perpendicular to the wire in the Y-Z plane.
- The simplified `set()` places planes at 0.375, 0.25, 0.125 of `dx` along X
  (lines 132-134), meaning U is farthest from the origin and W is closest.
- Wire offsets (`offset_mm`) shift the wire origin along the pitch direction
  from center (lines 83-85). This controls which wire lands on the Y-axis.

---

### WireGenerator.h / WireGenerator.cxx

**Algorithm overview:**
Generates wire geometry from `IWireParameters`. For each plane, it starts at
the parameter-defined origin and walks in both positive and negative pitch
directions, creating wires by intersecting infinite wire rays with the bounding
box. Wires that don't intersect the box are discarded (termination condition).

**Potential bugs:**

**(B6) Raw `new` usage (lines 54-69, 132):**
`make_wire()` returns a raw `GenWire*` allocated with `new` (line 69). If
`box_intersection` returns non-3, it returns `0` (null). The caller checks for
null but the successfully-created wires in `these_wires` vector are raw pointers.
If an exception occurs between creation and wrapping into `IWire::pointer`
(line 118), they leak.

`new IWire::vector` at line 132 is wrapped into a `shared_ptr` via
`output_pointer` at line 136, which is fine, but the pattern is unusual.

**Key algorithmic details:**

- Wire ident encoding: `(planeIndex+1) * 100000 + wireIndex` (lines 31-33).
  This limits each plane to < 100000 wires.
- The `proto` vector (line 81) is computed as `pitch.cross(xaxis).norm()` -- a
  unit vector along the wire direction, derived from the pitch direction crossed
  with the drift (X) axis.
- After generating wires in both directions, they are sorted by original index,
  then re-indexed contiguously starting from 0 (lines 111-119).
- Wire rays are oriented with `hits.first.y() < hits.second.y()` (lines 66-68),
  ensuring consistent Y-direction ordering.

---

### WireSummary.h / WireSummary.cxx

**Algorithm overview:**
Provides spatial queries over a wire collection: closest wire to a point,
bounding wire pair, pitch distance. Internally caches per-plane wire vectors
and pitch geometry. Pitch is determined from two wires near the center of each
plane.

**Potential bugs:**

**(B4) Memory leak in destructor (line 189):**
```cpp
WireSummary::~WireSummary() {}
```
The destructor does not `delete m_cache`. The cache is allocated with `new` at
line 187. This leaks the entire `WireSummaryCache` (including three
`WirePlaneCache` objects and all cached data) every time a `WireSummary`
is destroyed.

**(B5) Index truncation in `closest()` (lines 117-118):**
```cpp
double dist = wpc->pitch_distance(point);
return wpc->wire_by_index(dist / wpc->pitch_mag);
```
`dist / pitch_mag` is a floating-point value that gets implicitly truncated
(not rounded) when passed to `wire_by_index(int)`. For points near the boundary
between two wires, truncation toward zero gives the wrong result. The correct
approach would be `std::round(dist / pitch_mag)` or `std::lround(...)`. For
example, if the point is at pitch distance 2.9 * pitch_mag, the index will be 2
instead of the correct 3.

**Efficiency concerns:**

**(E3) `by_chan()` side effect (line 91):**
```cpp
IWireSegmentSet& got = chan2wire[chan];
```
Using `operator[]` on `std::map` inserts a default entry if the key doesn't
exist. If queried with unknown channel IDs, the map grows unboundedly. Should
use `find()` instead.

**Key algorithmic details:**

- Pitch is determined from two wires near the center of the collection
  (`wires.size()/2`), not from the first two. This is more robust for
  non-uniform wire spacing near edges.
- `pitch_distance()` uses `ray_dist()`, measuring signed distance from the
  pitch ray origin (first wire center) to the query point along the pitch
  direction.

---

### WireSummarizer.h / WireSummarizer.cxx

**Algorithm overview:**
Trivial function node that constructs a `WireSummary` from an `IWire::vector`.
Acts as a converter/adapter in a data flow graph.

**No bugs or efficiency concerns found.** This is a 22-line pass-through.

---

### WireSchemaFile.h / WireSchemaFile.cxx

**Algorithm overview:**
Loads wire geometry from a file into a `WireSchema::Store`. Supports a
configurable correction mode (default: pitch correction). The store is loaded
during `configure()` and served via `wire_schema_store()`.

**Potential bugs:**

- If `configure()` is never called but `wire_schema_store()` is accessed, the
  returned store is default-constructed (empty). No error or warning is given.
  The constructor accepts a filename but does not load it -- loading only happens
  in `configure()`.

**No efficiency concerns.** The store is loaded once and cached.

---

### WireBoundedDepos.h / WireBoundedDepos.cxx

**Algorithm overview:**
Spatial filter for depositions based on wire regions. Depos are accepted or
rejected based on whether their nearest wire (in pitch coordinates) falls within
configured wire-number ranges. Supports accept/reject modes and multiple
regions (each region is an AND of per-plane wire bounds; regions are OR'd).

**Potential bugs:**

**(B7) Config key mismatch (lines 39 vs 67):**
`default_configuration()` defines:
```cpp
cfg["wires"] = Json::arrayValue;  // line 39
```
But `configure()` reads:
```cpp
auto jregions = cfg["regions"];  // line 67
```
The key `"regions"` is never defined in the default config. If a user relies on
the default config structure, the regions array will be empty (null JSON) and no
filtering will occur. The documented parameter name in the header comment
(lines 22-26) says `"wires"`, matching the default config but not the actual
reader.

- `m_pimpos[iplane]` (line 99) is indexed by `iplane` from the wire bounds
  config. If the user specifies a plane index >= `m_pimpos.size()`, this is an
  out-of-bounds access (no bounds check).

- The `break` at line 63 means only the first non-empty face's planes are used
  for pimpos. The `// fixme:` comment acknowledges this limitation. For
  multi-face detectors where different faces have different plane geometries,
  this will silently use incorrect pimpos.

**Key algorithmic details:**

- Wire number is determined via `Pimpos::distance()` (pitch coordinate) and
  `region_binning().bin()` (discretization to wire index).
- The inclusive range `[imin, imax]` (line 105) matches the header documentation
  ("inclusive of the max").
- Lazy calculation of closest wire per plane (lines 88-101) avoids redundant
  pitch calculations when multiple regions reference the same plane.

---

### WireSource.h / WireSource.cxx

**Algorithm overview:**
Convenience source node combining `WireParams` and `WireGenerator` into a
single configurable component. Delegates configuration to `WireParams` and
wire generation to `WireGenerator`.

**No bugs or efficiency concerns found.** This is a thin facade.

---

## Cross-Cutting Concerns

### Memory Management Patterns
Several classes use raw `new` without corresponding `delete`:
- `WirePlane::m_pimpos` (B8) -- `Pimpos*` leaked on reconfiguration
- `WireSummary::m_cache` (B4) -- `WireSummaryCache*` leaked on destruction
- `WireGenerator::make_wire` (B6) -- raw `GenWire*` before wrapping
- `channel_wire_collector_t` (B2) -- raw `SimpleChannel*` on error paths

Modern C++ would use `std::unique_ptr` for owned resources. These are
legacy patterns but represent real leak risks, especially under exception paths.

### Coordinate System Assumptions
The entire group assumes:
- X axis = drift direction
- Wire planes are perpendicular to X
- Y axis = roughly vertical (gravity)
- Z axis = roughly horizontal (beam)
- Wire rays point in +Y direction by convention
- Pitch direction = X cross wireDirection

### Configuration Consistency
The `WireBoundedDepos` config key mismatch (B7) is a functional bug that could
cause silent no-op filtering. The `AnodePlane` null-face JSON access (B1) is
a latent issue that happens to work due to JSONCPP returning 0.0 for null
numeric access.
