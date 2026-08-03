---
name: wct-architecture
description: The Wire-Cell Toolkit dependency and plugin architecture — the util/iface/aux spine, how plugins aggregate dependencies at runtime, and the rules against new and inter-plugin dependencies. Use when adding a package, wiring dependencies, or deciding where code belongs.
---

# WCT architecture: layers, plugins, dependencies

WCT factors dependencies by package. Abstract interface classes let concrete
implementations live in leaf **plugin** shared libraries, so dependencies are
aggregated at *runtime* rather than baked into shared layers. Plugin shared
libraries also factor code by major application topic.

## The dependency spine

- **`util/`** — base library. **No new dependencies allowed, ever.**
- **`iface/`** — abstract interface base classes; almost no concrete code; depends only
  on `util`.
- **`aux/`** — general concrete implementations built on `iface`.
- **Plugins** — the remaining package directories produce plugin shared libraries that
  depend on `aux` or `iface`.

## Hard rules

- **Plugins MUST NOT depend on other plugins.** If code in one plugin seems to need
  another plugin, that is a signal to move the shared code *down* into `aux` or `iface`
  (behind an interface), not to add a plugin-to-plugin dependency. Stop and surface this
  rather than wiring the link.
- **Do not add a new external dependency** to any package — and never create or modify a
  plugin in a way that pulls one in. New dependencies are decided by the WCT team. If a
  task seems to require one, stop and raise it with the user; do not add it on request.
- A dependency a plugin legitimately already uses stays scoped to that single plugin.

## Placement guide

| You have… | Put it in… |
|---|---|
| A pure abstract interface | `iface/` |
| Concrete code usable across topics, no new deps | `aux/` |
| Topic-specific concrete code | the relevant plugin |
| Code two plugins both need | move to `aux`/`iface` — not a plugin cross-link |
| Anything needing a new external dependency | stop — team decision |
