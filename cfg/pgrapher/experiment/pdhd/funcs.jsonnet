// This provides some util functions.

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

{
  // Build a fanout-[pipelines]-fanin graph.  pipelines is a list of
  // pnode objects, one for each spine of the fan.
  fanpipe:: function(fout, pipelines, fin, name='fanpipe', outtags=[]) {

    local fanmult = std.length(pipelines),
    local fannums = std.range(0, fanmult - 1),

    local fanout = g.pnode({
      type: fout,
      name: name,
      data: {
        multiplicity: fanmult,
        tag_rules: [  // example in gen/test/test_fans.jsonnet
          {
            frame: {
              //'.*': 'number%d' % n,
              //'.*': 'gauss%d' % n,
              //'.*': 'framefanout%d ' % n,
              // '.*': 'orig%d' % n,
              '.*': 'orig',
            },
            trace: {
              // fake doing Nmult SP pipelines
              //orig: ['wiener', 'gauss'],
              //'.*': 'orig',
            },
          }
          for n in fannums
        ],
      },
    }, nin=1, nout=fanmult),


    local fanin = g.pnode({
      type: fin,
      name: name,
      data: {
        multiplicity: fanmult,
        tag_rules: [
          {
            frame: {
              //['number%d' % n]: ['output%d' % n, 'output'],
              '.*': 'framefanin',
            },
            trace: {
              //gauss: 'gauss%d' % n,
              //wiener: 'wiener%d' % n,
              ['raw%d' % n]: ['raw%d' % n],
              ['gauss%d' % n]: ['gauss%d' % n],
              ['wiener%d' % n]: ['wiener%d' % n],
              ['threshold%d' % n]: ['threshold%d' % n],
              ['dnnsp%d' % n]: ['dnnsp%d' % n],
            },

          }
          for n in fannums
        ],

        //tags: if outtags == [] then ['from-pipeline-%d' % n for n in fannums] else outtags,
        //tags: ['from-pipeline-%d' % n for n in fannums],
        tags: outtags,
      },
    }, nin=fanmult, nout=1),

    ret: g.intern(innodes=[fanout],
                  outnodes=[fanin],
                  centernodes=pipelines,
                  edges=
                  [g.edge(fanout, pipelines[n], n, 0) for n in std.range(0, fanmult - 1)] +
                  [g.edge(pipelines[n], fanin, 0, n) for n in std.range(0, fanmult - 1)],
                  name=name),
  }.ret,

  // The approximated sim+sigproc
  splat:: function(params, tools, anode, name=null) {
    local apaid = anode.data.ident,
    local sufix = if std.type(name) == "null" then apaid else name,
    local bg = g.pnode({
        type:'DepoBagger',
        name: sufix,
        data: {
            gate: [params.sim.ductor.start_time,
                   params.sim.ductor.start_time+params.sim.ductor.readout_time],
        },
    }, nin=1, nout=1),
    local sp = g.pnode({
        type: 'DepoFluxSplat',
        name: sufix,
        data: {
            anode: wc.tn(anode),
            field_response: wc.tn(tools.field), // for speed and origin
            sparse: true,
            tick: params.daq.tick,
            window_start: params.sim.ductor.start_time,
            window_duration: params.sim.ductor.readout_time,
            reference_time: 0.0,
            // Run wirecell-gen morse-* to find these numbers that match the extra
            // spread the sigproc induces.
            "smear_long": [
                2.691862363980221,
                2.6750200122535057,
                2.7137567141154055
            ],
            "smear_tran": [
                0.7377218875719689,
                0.7157764520393882,
                0.13980698710556544
            ]
        },
    }, nin=1, nout=1, uses=[anode, tools.field]),
    local rt = g.pnode({
        type: 'Retagger',
        name: sufix,
        data: {
            // Note: retagger keeps tag_rules an array to be like frame fanin/fanout.
            tag_rules: [{
                // Retagger also handles "frame" and "trace" like fanin/fanout
                // merge separately all traces like gaussN to gauss.
                frame: {
                ".*": "deposplat%d" % apaid
                },
                merge: {
                ".*": "deposplat%d" % apaid
                },
            }],
        },
    }, nin=1, nout=1),
    ret: g.pipeline([bg, sp, rt],"%s-%s" % [bg.name, sp.name]),
}.ret,


}
