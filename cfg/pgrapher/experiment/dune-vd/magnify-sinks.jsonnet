// This provides multiple MagnifySink for e.g. protoDUNE

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

// multiple MagnifySink
// tagn (n = 0, 1, ... 5) for anode[n]
// FrameFanin tags configured in sim.jsonnet
function(tools, outputfile) {

  local nanodes = std.length(tools.anodes),

  local magorig = [
    g.pnode({
      type: 'MagnifySink',
      name: 'magorig%d' % anode.data.ident,
      data: {
        output_filename: outputfile,
        root_file_mode: 'UPDATE',
        frames: ['orig%d' % anode.data.ident],
        trace_has_tag: false,   // traces from source have NO tag
        anode: wc.tn(anode),
      },
    }, nin=1, nout=1)
    for anode in tools.anodes
  ],

  local magraw = [
    g.pnode({
      type: 'MagnifySink',
      name: 'magraw%d' % anode.data.ident,
      data: {
        output_filename: outputfile,
        root_file_mode: 'UPDATE',
        frames: ['raw%d' % anode.data.ident],
        trace_has_tag: true,
        cmmtree: [["noisy", "T_noisy%d"%anode.data.ident],
                  ["sticky", "T_stky%d"%anode.data.ident],
                  ["ledge", "T_ldg%d"%anode.data.ident],
                  ["harmonic", "T_hm%d"%anode.data.ident] ], // maskmap in nf.jsonnet 
        anode: wc.tn(anode),
      },
    }, nin=1, nout=1)
    for anode in tools.anodes
  ],

  local magdecon = [
    g.pnode({
      type: 'MagnifySink',
      name: 'magdecon%d' % anode.data.ident,
      data: {
        output_filename: outputfile,
        root_file_mode: 'UPDATE',
        frames: ['gauss%d' % anode.data.ident, 'wiener%d' % anode.data.ident],
        trace_has_tag: true,
        anode: wc.tn(anode),
      },
    }, nin=1, nout=1)
    for anode in tools.anodes
  ],

  local magdebug = [
    g.pnode({
      type: 'MagnifySink',
      name: 'magdebug%d' % anode.data.ident,
      data: {
        output_filename: outputfile,
        root_file_mode: 'UPDATE',
        frames: ['tight_lf%d' %anode.data.ident, 'loose_lf%d' %anode.data.ident, 'cleanup_roi%d' %anode.data.ident,
                 'break_roi_1st%d' %anode.data.ident, 'break_roi_2nd%d' %anode.data.ident, 'decon_charge%d' %anode.data.ident,
                 'shrink_roi%d' %anode.data.ident, 'extend_roi%d' %anode.data.ident, 'mp2_roi%d' %anode.data.ident, 'mp3_roi%d' %anode.data.ident],
        trace_has_tag: true,
        anode: wc.tn(anode),
      },
    }, nin=1, nout=1)
    for anode in tools.anodes
  ],

  local magthr = [
    g.pnode({
      type: 'MagnifySink',
      name: 'magthr%d' % anode.data.ident,
      data: {
        output_filename: outputfile,
        root_file_mode: 'UPDATE',
        summaries: ['wiener%d' % anode.data.ident],  // note that if tag set, each apa should have a tag set for FrameFanin
        summary_operator: { ['wiener%d' % anode.data.ident]: 'set' },  // []: obj comprehension
        anode: wc.tn(anode),
      },
    }, nin=1, nout=1)
    for anode in tools.anodes
  ],

  local magdnnroi = [
    g.pnode({
      type: 'MagnifySink',
      name: 'magdnnroi%d' % anode.data.ident,
      data: {
        output_filename: outputfile,
        root_file_mode: 'UPDATE',
        frames: ['dnnsp%d' %anode.data.ident],
        trace_has_tag: true,
        anode: wc.tn(anode),
      },
    }, nin=1, nout=1)
    for anode in tools.anodes
  ],


  return: {
    orig_pipe: [g.pipeline([magorig[n]], name='magorigpipe%d' % n) for n in std.range(0, nanodes - 1)],
    raw_pipe: [g.pipeline([magraw[n]], name='magrawpipe%d' % n) for n in std.range(0, nanodes - 1)],
    decon_pipe: [g.pipeline([magdecon[n]], name='magdeconpipe%d' % n) for n in std.range(0, nanodes - 1)],
    debug_pipe: [g.pipeline([magdebug[n]], name='magdebugpipe%d' % n) for n in std.range(0, nanodes - 1)],
    threshold_pipe: [g.pipeline([magthr[n]], name='magthrpipe%d' % n) for n in std.range(0, nanodes - 1)],
    dnnroi_pipe: [g.pipeline([magdnnroi[n]], name='magdnnroipipe%d' % n) for n in std.range(0, nanodes - 1)],
  },


}.return
