# What a detector actually runs, versus the C++ default

`clus/test/doctest_clus_knob_defaults.cxx` pins the **C++ default** of every
knob the SBND pattern-recognition port added, and it passes on every build.
That test is load-bearing: each knob shipped under the rule that with its key
absent the component reproduces the pre-knob behaviour byte-for-byte, and a
silently moved default would make every past A/B gate vacuous while still
showing green.

It is also easy to misread. It says nothing about what any detector runs.
**192 of the knobs present in the SBND production config differ from their C++
default** -- so reading the C++ (or watching that test go green) gives the
wrong answer for most of the operating point.

`sbnd-operating-point.md` in this directory is that delta, generated. Examples:

| knob | C++ default | SBND production |
|---|---|---|
| `sccc_max_gap` | 5 | 10 |
| `sccc_kink_max` | 15 | 18 |
| `kine_shower_fudge_factor` | 0.8 | 0.86 |

## Regenerate

```
python3 clus/docs/knobs/knob_delta.py \
    clus/test/doctest_clus_knob_defaults.cxx \
    <a compiled job or pinned ref/*/prod_prjob.json>
```

The second argument must be a **compiled job**, not an imported module: SBND's
`clus.jsonnet` and the per-detector `clus.jsonnet` files are modules whose
fields are all hidden (`::`), so compiling one directly yields `{}` and every
knob would look un-overridden. The script refuses an input with no component
`data` blocks for that reason.

## Why the defaults are not simply changed to match

Because they are not SBND's to change. `qlport/uboone-mabc.jsonnet` sets none
of these keys -- it inherits the C++ defaults -- so moving a default to SBND's
production value would silently change uBooNE, the frozen reference chain, and
would break the byte-identical-when-absent guarantee the whole knob discipline
rests on. The delta is documentation, not a migration plan.
