---
name: wct-beads
description: Use the `bd` (beads) issue tracker for local task tracking in WCT. Use ONLY when a .beads/ directory already exists in the tree — beads is optional and per-developer. Covers day-to-day bd commands; does not init beads or impose any git workflow.
---

# Beads (`bd`) for local task tracking

Beads is **optional and per-developer** in WCT. Some developers track work locally with
`bd`; others do not.

## Precondition — check before doing anything

Engage this skill only if a **`.beads/` directory exists** at the top of the working
tree. If it does not, the developer is not using beads here: do **not** run `bd`, and do
**not** create the directory or run any init/setup command. Just track work however the
task otherwise calls for.

```bash
test -d .beads && echo "beads in use" || echo "no beads here — skip"
```

## Scope limits (do not cross these)

- **Never** run `bd init`, `bd migrate`, or any command that creates/initializes a beads
  database. If `.beads/` is absent, beads is simply not in use here.
- **Impose no git workflow.** Beads' own defaults suggest a mandatory commit/push routine;
  ignore that. `.beads/` is intentionally **not committed** to this repo — do not `git add`
  it, and do not push, pull, or gate work on beads state.
- Treat `bd` as a private, local scratch tracker only.

## Everyday commands

```bash
bd ready                 # issues ready to work (no open blockers)
bd list                  # list issues (add filters, e.g. --status open)
bd show <id>             # full detail for one issue
bd create "Title"        # create an issue (add -d "desc", -p <prio>, -t <type>)
bd update <id> --claim   # atomically claim an issue to work on it
bd update <id> ...       # change fields (status, priority, assignee, ...)
bd note <id> "text"      # append a working note
bd close <id>            # close a finished issue
bd dep add <id> <blocker-id>   # record that <id> is blocked by <blocker-id>
bd search "text"         # find issues by text
bd status                # database overview / stats
```

Run `bd <command> --help` for flags. Use beads to track multi-step or cross-session work
locally; keep it out of anything shared (commits, PRs, Issues) unless the developer asks.
