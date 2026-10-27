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
// multiplicity <= max.  `data` is merged into each fan node's `data` (on top of
// the computed `multiplicity`), e.g. {tag_rules: [...]} for frame fans.

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
  local fan_node(fantype, name, mult, nin, nout, data) =
    g.pnode({ type: fantype, name: name, data: { multiplicity: mult } + data },
            nin=nin, nout=nout),

  // --- fanout: 1 input port -> N output ports --------------------------------
  local build_out(tree, fantype, basename, path, data) =
    if std.type(tree) == 'number'
    then fan_node(fantype, basename + path, tree, 1, tree, data)
    else
      local b = std.length(tree);
      local root = fan_node(fantype, basename + path, b, 1, b, data);
      local kids = [ build_out(tree[i], fantype, basename, path + '_%d' % i, data)
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
  local build_in(tree, fantype, basename, path, data) =
    if std.type(tree) == 'number'
    then fan_node(fantype, basename + path, tree, tree, 1, data)
    else
      local b = std.length(tree);
      local root = fan_node(fantype, basename + path, b, b, 1, data);
      local kids = [ build_in(tree[i], fantype, basename, path + '_%d' % i, data)
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
