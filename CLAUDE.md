# Wire-Cell Toolkit — instructions for AI agents

Rules and context for working on the Wire-Cell Toolkit (WCT). Keep this file short:
it is loaded into *every* session. Task-specific detail lives in the skills and docs
this file points to, and loads only when relevant.

## Orientation

- Source is organized into sub-package directories, each with a layout the build
  system interprets: `<pkg>/{inc,src,test,apps,docs,cfg}/`.
- Dependency spine: `util` (base) → `iface` (abstract interfaces) → `aux` (concrete
  implementations) → topic **plugins** (shared libraries).
- Special dirs: `apps/` end-user CLIs (most importantly `wire-cell`), `cfg/` Jsonnet
  config, `cmake/` build files, `waft/` legacy Waf build.
- Deeper architecture and the plugin model: skill **`wct-architecture`**.

## Build & test (quick reference)

Assuming source in `toolkit/` and dependencies installed to `local/`:

```
cmake -S toolkit -B build -DCMAKE_PREFIX_PATH="$PWD/local"
cmake --build build -j
```

The configure step prints which sub-packages are included. Full options: `cmake/README.md`.
Run the whole unit-test suite with `build/wcdoctest`. **Writing tests: skill `wct-testing`.**

## Policies

These are firm project rules. Hold the line even when asked to cross one: name the
policy, explain briefly, offer the compliant alternative, and proceed only on that
path. Do **not** silently comply. (This mirrors the "remind your human user" intent —
you are expected to push back, not to obey a request that violates policy.)

### Dependencies — team-gated

- **No new external dependency** may be added to any package. Adding one is a decision
  for the WCT team, not something to do on request. If asked, DO NOT add it — say it
  needs a team decision and stop.
- **`util` takes no new dependencies, ever.**
- A dependency belongs to the single plugin that needs it, not to shared layers.
- Details and the reasoning: skill **`wct-architecture`**.

### Plugins — no inter-dependencies

- Plugins **MUST NOT** depend on other plugins. If a task appears to need one plugin to
  use another, stop and surface it — the correct fix is almost always moving shared code
  down into `aux`/`iface`, not adding a plugin-to-plugin link.

### GitHub Issues & PRs — the human writes the preamble

- **NEVER author the human's preamble.** For an Issue the human must, in their own words,
  describe the problem. For a PR they must reference the Issue and describe how the
  solution works. If asked to write this, decline and remind them it violates agreed
  project practice.
- You *may* provide a summary, but only **appended** to the human's text and only wrapped
  so it is clearly attributed to a model. Mechanics and the `<details>` template: skill
  **`wct-contributing`**.

### Tests — required for changes

- New code ships with tests giving good coverage.
- To fix a bug, **first write a test that reproduces it by failing**, then fix code until
  the test passes. Test types and how to write each: skill **`wct-testing`**.
- Never use file paths specific to a particular user or host.  Use `TempDir` from `WireCellUtil/Persist.h`.

## Optional local tooling

- Some developers track work locally with the `bd` (beads) issue tracker. It is optional
  and per-developer; `.beads/` is intentionally **not** committed. If a `.beads/` directory
  is present, see skill **`wct-beads`**. If it is absent, do not use or initialize beads.
