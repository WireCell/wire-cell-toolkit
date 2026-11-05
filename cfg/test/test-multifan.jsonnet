// Test multifan.jsonnet: the auto-factoring fan-out/fan-in tree builder.
// Run (via the WCT test harness) through wcsonnet; all checks are `assert`s so a
// failure aborts compilation.  Emits a small real fanout->sinks graph at the end.
local g = import 'pgraph.jsonnet';
local mf = import 'multifan.jsonnet';

// all() needs jsonnet >= 0.19; this works with 0.18.
local all(xs) = std.foldl(function(a, b) a && b, xs, true);

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
  assert all([m <= max for m in fan_mults(fo)]) : 'fanout(%d,%d) a node exceeds max' % [n, max];
  true;

local check_fanin(n, max) =
  local fi = mf.fanin('in%d_%d' % [n, max], 'FrameFanin', n, max);
  assert std.length(fi.iports) == n : 'fanin(%d) iports=%d != %d' % [n, std.length(fi.iports), n];
  assert std.length(fi.oports) == 1 : 'fanin(%d) oports != 1' % n;
  assert all([m <= max for m in fan_mults(fi)]) : 'fanin(%d,%d) a node exceeds max' % [n, max];
  true;

local graph_checks = [ check_fanout(n, 10) for n in [1, 4, 12, 100, 320] ]
                   + [ check_fanin(n, 10) for n in [1, 4, 12, 100, 320] ];

// fan size <= max is returned as a SINGLE atomic node (no intern wrapper)
local atomic_check =
  local fo = mf.fanout('atom', 'DepoSetFanout', 7, 10);
  assert std.length(g.uses(fo)) == 1 : 'fanout(7,10) should be a single atomic node';
  true;

// ---- per-node data from a function: each node sees its own tree ports -------
// Record node.ports in the data and check that it partitions 0..N-1 in order at
// every node, that leaves map their ports one-to-one and that each node's port
// count is its multiplicity.
local fan_nodes(pnode) = [
  x for x in g.uses(pnode)
  if std.objectHas(x, 'data') && std.objectHas(x.data, 'multiplicity')
];
local record = function(node) { ports: node.ports, leaf: node.leaf, offset: node.offset };

local check_node_ports(x) =
  local d = x.data;
  local m = d.multiplicity;
  assert std.length(d.ports) == m : '%s: %d port spans for multiplicity %d' % [x.name, std.length(d.ports), m];
  assert std.flattenArrays(d.ports) == std.range(d.offset, d.offset + std.length(std.flattenArrays(d.ports)) - 1) :
    '%s: port spans are not contiguous from offset %d' % [x.name, d.offset];
  assert (!d.leaf) || d.ports == [ [d.offset + i] for i in std.range(0, m - 1) ] :
    '%s: leaf ports are not one per tree port' % x.name;
  true;

local check_data_func(n, max, builder) =
  local nodes = fan_nodes(builder('df%d_%d' % [n, max], n, max));
  local leaf_ports = std.flattenArrays([ std.flattenArrays(x.data.ports) for x in nodes if x.data.leaf ]);
  local roots = [ x for x in nodes if x.data.offset == 0
                  && std.length(std.flattenArrays(x.data.ports)) == n ];
  assert all([ check_node_ports(x) for x in nodes ]) : 'data(%d,%d) node ports' % [n, max];
  assert std.sort(leaf_ports) == std.range(0, n - 1) : 'data(%d,%d) leaves do not cover 0..N-1 once' % [n, max];
  assert std.length(roots) >= 1 : 'data(%d,%d) no node covers all tree ports' % [n, max];
  true;

local data_checks = [
  check_data_func(n, max, function(name, n, max) mf.fanout(name, 'FrameFanout', n, max, record))
  for n in [1, 7, 23, 97] for max in [3, 10]
] + [
  check_data_func(n, max, function(name, n, max) mf.fanin(name, 'FrameFanin', n, max, record))
  for n in [1, 7, 23, 97] for max in [3, 10]
];

// per-port tag_rules from a function have one rule per port of each node, and
// the leaves carry the rule of their own tree port
local tagged_check =
  local n = 23;
  local rules = [ { frame: { '.*': 'p%d' % p } } for p in std.range(0, n - 1) ];
  local fo = mf.fanout('tagged', 'FrameFanout', n, 10,
                       function(node)
                         if node.leaf then { tag_rules: [ rules[p[0]] for p in node.ports ] }
                         else {});
  local nodes = fan_nodes(fo);
  local leaves = [ x for x in nodes if std.objectHas(x.data, 'tag_rules') ];
  assert all([ std.length(x.data.tag_rules) == x.data.multiplicity for x in leaves ]) :
    'tag_rules length differs from multiplicity';
  assert std.sort([ r.frame['.*'] for x in leaves for r in x.data.tag_rules ])
         == std.sort([ 'p%d' % p for p in std.range(0, n - 1) ]) :
    'leaf tag_rules do not cover each tree port once';
  // tree port i leaves a leaf whose rule for that port is rules[i]
  local by_tn = { [x.type + ':' + x.name]: x for x in leaves };
  assert all([ by_tn[fo.oports[i].node].data.tag_rules[fo.oports[i].port] == rules[i]
               for i in std.range(0, n - 1) ]) :
    'tree port does not get its own rule';
  true;

// object data is merged unchanged into every node
local object_check =
  local fo = mf.fanout('obj', 'FrameFanout', 23, 10, { foo: 42 });
  assert all([ x.data.foo == 42 for x in fan_nodes(fo) ]) : 'object data not on every node';
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

assert all(factor_checks) : 'factor checks failed';
assert all(graph_checks) : 'graph checks failed';
assert atomic_check : 'atomic check failed';
assert all(data_checks) : 'data function checks failed';
assert tagged_check : 'tagged check failed';
assert object_check : 'object data check failed';

g.uses(graph) + [{ type: 'Pgrapher', data: { edges: g.edges(graph) } }]
