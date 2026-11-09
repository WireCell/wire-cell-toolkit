// multifan.jsonnet -- fan-out / fan-in subgraphs that honor a maximum per-node
// fan multiplicity by AUTOMATICALLY factoring a large fan into a tree of smaller
// fans.
//
// Why: a single WCT fan node (DepoSetFanout, FrameFanout, FrameFanin, ...) and
// in particular the TbbFlow FanoutWrapper cap the number of ports of one node
// (TBB asserts nout <= 10).  To fan to/from more spines than that one must build
// a tree ("fans of fans").  Unlike pgrapher/common/funcs.jsonnet's multifan*,
// which require the caller to hand-specify the node count and multiplicity of
// every tree layer, this version derives the whole tree from only the requested
// total multiplicity and a single `max` -- nothing extra to thread through the
// FCL -> Jsonnet boundary.
//
//   local mf = import 'multifan.jsonnet';
//   local fo = mf.fanout('fan', 'DepoSetFanout', 320);  // pnode: 1 in, 320 out
//   local fi = mf.fanin ('fan', 'FrameFanin',    320);  // pnode: 320 in, 1 out
//
// If multiplicity <= max a single atomic fan node is returned; otherwise a
// subgraph (g.intern) whose exposed ports are exactly the 1 + N (fanout) or
// N + 1 (fanin) leaf ports, in order.  Every fan node in the tree has
// multiplicity <= max.
//
// `data` sets the rest of each fan node's `data` (the computed `multiplicity` is
// always set).  It is either:
//
// - an object, merged unchanged into EVERY fan node of the tree.  Use this only
//   for port-independent settings.  Per-port arrays such as FrameFanout or
//   FrameFanin `tag_rules` (or FrameFanin `tags`) have one entry per port of a
//   node, and the nodes of a tree have different multiplicities and cover
//   different ranges of the N tree ports, so an object can not express them.
//
// - a function(node) returning an object, called once per fan node with
//
//     node.multiplicity  the node's number of fan ports (m)
//     node.leaf          true if the node's fan ports are tree ports
//     node.offset        the tree port index of the node's first covered port
//     node.ports         array of m arrays: the tree port indices reached
//                        through each of the node's fan ports, in order
//
//   For a leaf, node.ports[i] == [node.offset + i].  Eg, for a FrameFanout tree
//   with one rule per tree port in `rules` (length N), apply them at the leaves
//   and give inner nodes no tag_rules (FrameFanout then copies frames as-is):
//
//     data = function(node)
//       if node.leaf then { tag_rules: [ rules[p[0]] for p in node.ports ] }
//       else {}
//
//   A FrameFanin tree is the mirror image, but its inner nodes need rules (or
//   `tags`) that carry forward the tags made by the nodes upstream of them, as a
//   FrameFanin without them drops the input tags.

local g = import 'pgraph.jsonnet';

{
  // integer power b^e (e >= 0)
  local ipow(b, e) = std.foldl(function(a, _) a * b, std.range(1, e), 1),

  // smallest depth d (>=1) such that max^d >= n
  local depth(n, max) = std.filter(function(d) ipow(max, d) >= n, std.range(1, 64))[0],

  // --- the shared factoring --------------------------------------------------
  // factor(n, max) -> a tree describing how to split an n-way fan into fans no
  // larger than `max`.  A LEAF is an integer m (1 <= m <= max): a single fan of
  // size m.  An INTERNAL node is an ARRAY of child trees: a fan whose size is the
  // array length, each of whose ports feeds one child subtree.  Reading leaves
  // left-to-right they sum to EXACTLY n.  When n does not factor evenly we round
  // the implied depth up and spread the remainder, so leaf fans differ by at
  // most one.  Pure data (no graph) so fanout() and fanin() share it.
  factor(n, max=10)::
    assert n >= 1 : 'multifan.factor: multiplicity must be >= 1';
    assert max >= 2 : 'multifan.factor: max must be >= 2';
    if n <= max then n
    else
      local d = depth(n, max);
      local cap = ipow(max, d - 1);      // most leaves one child subtree can hold
      local b = std.ceil(n / cap);       // minimal top branching (2 <= b <= max)
      local base = std.floor(n / b);
      local rem = n - base * b;          // the first `rem` children get one extra
      [ self.factor(base + (if i < rem then 1 else 0), max)
        for i in std.range(0, b - 1) ],

  // count the leaf slots (= total fan multiplicity) a factor tree provides
  leaves(tree):: if std.type(tree) == 'number' then tree
                 else std.foldl(function(a, t) a + self.leaves(t), tree, 0),

  // --- node helper -----------------------------------------------------------
  // the tree port indices reached through each top-level port of `tree`, whose
  // first tree port is `offset`
  local port_spans(tree, offset) =
    if std.type(tree) == 'number' then [ [offset + i] for i in std.range(0, tree - 1) ]
    else
      local sizes = [ $.leaves(t) for t in tree ];
      local starts = std.foldl(function(a, n) a + [a[std.length(a) - 1] + n],
                               sizes, [offset]);
      [ std.range(starts[i], starts[i + 1] - 1) for i in std.range(0, std.length(tree) - 1) ],

  local fan_data(tree, offset, data) =
    if std.type(data) == 'function'
    then data({
      multiplicity: if std.type(tree) == 'number' then tree else std.length(tree),
      leaf: std.type(tree) == 'number',
      offset: offset,
      ports: port_spans(tree, offset),
    })
    else data,

  local fan_node(fantype, name, tree, offset, nin, nout, data) =
    local mult = if std.type(tree) == 'number' then tree else std.length(tree);
    g.pnode({ type: fantype, name: name,
              data: { multiplicity: mult } + fan_data(tree, offset, data) },
            nin=nin, nout=nout),

  // --- fanout: 1 input port -> N output ports --------------------------------
  local build_out(tree, fantype, basename, path, data, offset=0) =
    if std.type(tree) == 'number'
    then fan_node(fantype, basename + path, tree, offset, 1, tree, data)
    else
      local b = std.length(tree);
      local spans = port_spans(tree, offset);
      local root = fan_node(fantype, basename + path, tree, offset, 1, b, data);
      local kids = [ build_out(tree[i], fantype, basename, path + '_%d' % i, data,
                               spans[i][0])
                     for i in std.range(0, b - 1) ];
      g.intern(innodes=[root], outnodes=kids,
               edges=[ g.edge(root, kids[i], i, 0) for i in std.range(0, b - 1) ],
               name=basename + path),

  fanout(name, fantype, multiplicity, max=10, data={})::
    build_out(self.factor(multiplicity, max), fantype, name, '', data),

  // As fanout() but driven by an EXPLICIT factor tree (see factor() for the
  // leaf-int / array-of-children format) instead of auto-factoring a count.
  // Useful to impose a specific per-layer branching, eg a balanced tree with
  // layer multiplicities [4, 5, 3] -> [[3,3,3,3,3]*4] (60 leaves).
  fanout_tree(name, fantype, tree, data={})::
    build_out(tree, fantype, name, '', data),

  // --- fanin: N input ports -> 1 output port ---------------------------------
  local build_in(tree, fantype, basename, path, data, offset=0) =
    if std.type(tree) == 'number'
    then fan_node(fantype, basename + path, tree, offset, tree, 1, data)
    else
      local b = std.length(tree);
      local spans = port_spans(tree, offset);
      local root = fan_node(fantype, basename + path, tree, offset, b, 1, data);
      local kids = [ build_in(tree[i], fantype, basename, path + '_%d' % i, data,
                              spans[i][0])
                     for i in std.range(0, b - 1) ];
      g.intern(innodes=kids, outnodes=[root],
               edges=[ g.edge(kids[i], root, 0, i) for i in std.range(0, b - 1) ],
               name=basename + path),

  fanin(name, fantype, multiplicity, max=10, data={})::
    build_in(self.factor(multiplicity, max), fantype, name, '', data),

  // As fanin() but driven by an EXPLICIT factor tree instead of a count.
  fanin_tree(name, fantype, tree, data={})::
    build_in(tree, fantype, name, '', data),
}
