// This provides multiple MagnifySink for e.g. protoDUNE

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

// multiple MagnifySink
// tagn (n = 0, 1, ... 5) for anode[n]
// FrameFanin tags configured in sim.jsonnet
function(tag, tools, outputfile) {

  local nanodes = std.length(tools.anodes),

  //local magnify = function(tag, index, tools) g.pnode({
  //  type: 'MagnifySink',
  //  name: 'mag%s%d' % [tag, index],
  //  data: {
  //    output_filename: outputfile,
  //    // root_file_mode: if index == 0 then "RECREATE" else "UPDATE",
  //    root_file_mode: 'RECREATE',
  //    frames: ['%s%d' % [tag, index]],  // note that if tag set, each apa should have a tag set for FrameFanin
  //    anode: wc.tn(tools.anodes[index]),
  //  },
  //}, nin=0, nout=0),

  //local multimagnify = [magnify(tag, n, tools) for n in std.range(0, nanodes-1)],

  local multimagnify = [
    g.pnode({
      type: 'MagnifySink',
      name: 'mag%s%d' % [tag, n],  // traces from source have NO tag
      data: {
        output_filename: outputfile,
        // root_file_mode: if (n == 0 && std.startsWith(tag,"orig")) then "RECREATE" else "UPDATE",
        root_file_mode: 'UPDATE',
        frames: ['%s%d' % [tag, n]],
        trace_has_tag: if tag == 'orig' then false else true,
        anode: wc.tn(tools.anodes[n]),
      },
    }, nin=1, nout=1)
    for n in std.range(0, nanodes - 1)
  ],

  local multimagnifysummaries = [
    g.pnode({
      type: 'MagnifySink',
      name: 'mag%s%d' % [tag, n],
      data: {
        output_filename: outputfile,
        //root_file_mode: if (n == 0 && std.startsWith(tag,"orig")) then "RECREATE" else "UPDATE",
        root_file_mode: 'UPDATE',
        summaries: ['%s%d' % [tag, n]],  // note that if tag set, each apa should have a tag set for FrameFanin
        summary_operator: { ['threshold%d' % n]: 'set' },  // []: obj comprehension
        anode: wc.tn(tools.anodes[n]),
      },
    }, nin=1, nout=1)
    for n in std.range(0, nanodes - 1)
  ],


  //return: g.pipeline([multimagnify[0], multimagnify[1], multimagnify[2], multimagnify[3], multimagnify[4], multimagnify[5]]),
  return: {
    magnify_pipelines: [g.pipeline([multimagnify[n]], name='magnifypipes%d' % n) for n in std.range(0, nanodes - 1)],
    magnifysummaries_pipelines: [g.pipeline([multimagnifysummaries[n]], name='magnifysumpipes%d' % n) for n in std.range(0, nanodes - 1)],
    magnify_pipeline: g.pipeline([multimagnify[n] for n in std.range(0, nanodes - 1)]),
  },


}.return
