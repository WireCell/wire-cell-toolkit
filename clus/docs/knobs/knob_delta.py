#!/usr/bin/env python3
"""Report where a detector's compiled config departs from the C++ knob defaults.

WHY.  doctest_clus_knob_defaults.cxx pins the *C++* default of every knob the
SBND pattern-recognition port added, and it passes -- but it says nothing about
what any detector actually runs.  Reading the C++ (or watching that test go
green) therefore gives the wrong answer for the parameters a detector overrides.
This script prints that delta so the operating point is discoverable without
changing a single default.

The C++ defaults are read from the test itself, which is the authoritative
record precisely because it is enforced on every build: if a default moved, the
test fails rather than this table going quietly stale.

Usage:
    knob_delta.py <doctest_clus_knob_defaults.cxx> <compiled-config.json>

The compiled config must be a REAL compiled job (wcsonnet output or a pinned
ref/*/prod_prjob.json), not a bare module compile -- compiling an imported
module yields {} and every knob would look un-overridden.
"""
import json
import re
import sys
from collections import OrderedDict

NUM = re.compile(r'CHECK_KNOB_NUM\s*\(\s*\w+\s*,\s*"([^"]+)"\s*,\s*([-0-9.eE+]+)\s*\)')
BOOL = re.compile(r'CHECK_KNOB_BOOL\s*\(\s*\w+\s*,\s*"([^"]+)"\s*,\s*(true|false)\s*\)')


def cxx_defaults(path):
    """key -> (value, kind) as pinned by the defaults test."""
    text = open(path).read()
    out = {}
    for key, val in NUM.findall(text):
        out[key] = (float(val), "num")
    for key, val in BOOL.findall(text):
        out[key] = (val == "true", "bool")
    return out


def config_values(path):
    """key -> {value: [component names]} over every node's data block."""
    doc = json.load(open(path))
    nodes = doc if isinstance(doc, list) else [doc]
    seen = {}
    for node in nodes:
        if not isinstance(node, dict):
            continue
        data = node.get("data")
        if not isinstance(data, dict):
            continue
        who = "%s:%s" % (node.get("type", "?"), node.get("name", ""))
        for key, val in data.items():
            if isinstance(val, (bool, int, float)):
                seen.setdefault(key, OrderedDict()).setdefault(_norm(val), []).append(who)
    return seen


def _norm(v):
    return bool(v) if isinstance(v, bool) else float(v)


def _fmt(v):
    if isinstance(v, bool):
        return "true" if v else "false"
    return ("%g" % v)


def main(argv):
    if len(argv) != 3:
        sys.exit(__doc__)
    defaults = cxx_defaults(argv[1])
    values = config_values(argv[2])
    if not values:
        sys.exit("ERROR: no component data found in %s -- is it a compiled JOB, "
                 "or an imported module that compiled to {}?" % argv[2])

    overridden, agreeing, unpinned = [], [], []
    for key, (dflt, _kind) in sorted(defaults.items()):
        if key not in values:
            continue
        for val, who in values[key].items():
            if val != dflt:
                overridden.append((key, dflt, val, len(who)))
            else:
                agreeing.append(key)
    for key in sorted(values):
        if key not in defaults:
            unpinned.append(key)

    print("# Knob delta: C++ default vs compiled config")
    print()
    print("- C++ defaults pinned by the test: **%d**" % len(defaults))
    print("- Of those, present in this config: **%d**" % (len(overridden) + len(agreeing)))
    print("- **Overridden (config != C++ default): %d**" % len(overridden))
    print("- Set to the C++ default explicitly: %d" % len(agreeing))
    print("- Config keys with no pinned C++ default: %d" % len(unpinned))
    print()
    print("Every row below is a value you get WRONG by reading the C++ default.")
    print()
    print("| knob | C++ default | this config | nodes |")
    print("|---|---|---|---|")
    for key, dflt, val, n in sorted(overridden):
        print("| `%s` | %s | **%s** | %d |" % (key, _fmt(dflt), _fmt(val), n))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
