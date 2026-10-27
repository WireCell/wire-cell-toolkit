// Test multifan.jsonnet: the auto-factoring fan-out/fan-in tree builder.
// Run (via the WCT test harness) through wcsonnet; all checks are `assert`s so a
// failure aborts compilation.  Emits a small real fanout->sinks graph at the end.
local g = import 'pgraph.jsonnet';
local mf = import 'multifan.jsonnet';

// ---- factoring: leaves sum to exactly N, every fan node <= max --------------
local leaves_ok(n, max) = mf.leaves(mf.factor(n, max)) == n;

// deepest fan multiplicity in a factor tree (a leaf int, or array length + kids)
local maxfan(tree) =
  if std.type(tree) == 'number' then tree
  else std.foldl(function(a, t) std.max(a, maxfan(t)), tree, std.length(tree));

local check_factor(n, max) =
  assert leaves_ok(n, max) :
    'factor(%d,%d) leaves=%d != %d' % [n, max, mf.leaves(mf.factor(n, max)), n];
  assert maxfan(mf.factor(n, max)) <= max :
    'factor(%d,%d) has a fan > max' % [n, max];
  true;

// a spread of sizes incl. powers, primes, off-by-one, and non-default max
local factor_checks = [
  check_factor(n, 10)
  for n in [1, 2, 7, 10, 11, 12, 13, 20, 97, 100, 101, 320, 999, 1000]
] + [ check_factor(n, m) for n in [1, 5, 50, 320] for m in [2, 3, 7] ];

// ---- built graphs: correct port counts, every configured fan <= max ---------
local fan_mults(pnode) = [
  x.data.multiplicity
  for x in g.uses(pnode)
  if std.objectHas(x, 'data') && std.objectHas(x.data, 'multiplicity')
];

local check_fanout(n, max) =
  local fo = mf.fanout('to%d_%d' % [n, max], 'DepoSetFanout', n, max);
  assert std.length(fo.iports) == 1 : 'fanout(%d) iports != 1' % n;
  assert std.length(fo.oports) == n : 'fanout(%d) oports=%d != %d' % [n, std.length(fo.oports), n];
  assert std.all([m <= max for m in fan_mults(fo)]) : 'fanout(%d,%d) a node exceeds max' % [n, max];
  true;

local check_fanin(n, max) =
  local fi = mf.fanin('in%d_%d' % [n, max], 'FrameFanin', n, max);
  assert std.length(fi.iports) == n : 'fanin(%d) iports=%d != %d' % [n, std.length(fi.iports), n];
  assert std.length(fi.oports) == 1 : 'fanin(%d) oports != 1' % n;
  assert std.all([m <= max for m in fan_mults(fi)]) : 'fanin(%d,%d) a node exceeds max' % [n, max];
  true;

local graph_checks = [ check_fanout(n, 10) for n in [1, 4, 12, 100, 320] ]
                   + [ check_fanin(n, 10) for n in [1, 4, 12, 100, 320] ];

// fan size <= max is returned as a SINGLE atomic node (no intern wrapper)
local atomic_check =
  local fo = mf.fanout('atom', 'DepoSetFanout', 7, 10);
  assert std.length(g.uses(fo)) == 1 : 'fanout(7,10) should be a single atomic node';
  true;

// ---- a real (compilable) fanout -> N sinks graph ----------------------------
local N = 23;  // > max -> exercises a tree (23 = [[6,6],[6,5]] etc.)
local source = g.pnode({ type: 'TrackDepos' }, nin=0, nout=1);
local fo = mf.fanout('demo', 'DepoSetFanout', N, 10);
local sinks = [ g.pnode({ type: 'DumpDepos', name: 'sink%02d' % i }, nin=1, nout=0)
               for i in std.range(0, N - 1) ];
local capoff = g.intern(innodes=[fo], centernodes=sinks,
                        edges=[ g.edge(fo, sinks[i], i, 0) for i in std.range(0, N - 1) ],
                        name='capoff');
local graph = g.pipeline([source, capoff], 'main');

assert std.all(factor_checks) : 'factor checks failed';
assert std.all(graph_checks) : 'graph checks failed';
assert atomic_check : 'atomic check failed';

g.uses(graph) + [{ type: 'Pgrapher', data: { edges: g.edges(graph) } }]
