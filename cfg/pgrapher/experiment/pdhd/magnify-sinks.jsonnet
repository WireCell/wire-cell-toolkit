// This provides multiple MagnifySink for e.g. protoDUNE

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

// multiple MagnifySink
// tagn (n = 0, 1, ... 5) for anode[n]
// FrameFanin tags configured in sim.jsonnet
// runinfo: optional {runNo, subRunNo, eventNo, total_time_bin} injected into
//          the first decon_trun_pipe sink's Trun tree; null skips Trun writing.
//          anodeNo is added automatically from the anode's ident.
function(tools, outputfile, runinfo=null) {

  local nanodes = std.length(tools.anodes),

  local magorig = [
    g.pnode({
      type: 'MagnifySink',
      name: 'magorig%d' % n,
      data: {
        output_filename: outputfile,
        root_file_mode: 'UPDATE',
        frames: ['orig%d' % n],
        trace_has_tag: false,   // traces from source have NO tag
        anode: wc.tn(tools.anodes[n]),
      },
    }, nin=1, nout=1)
    for n in std.range(0, nanodes - 1)
  ],

  local magraw = [
    g.pnode({
      type: 'MagnifySink',
      name: 'magraw%d' % n,
      data: {
        output_filename: outputfile,
        root_file_mode: 'UPDATE',
        frames: ['raw%d' % n],
        trace_has_tag: true,
        cmmtree: [["noisy", "T_noisy%d"%n],
                  ["sticky", "T_stky%d"%n],
                  ["ledge", "T_ldg%d"%n],
                  ["harmonic", "T_hm%d"%n] ], // maskmap in nf.jsonnet 
        anode: wc.tn(tools.anodes[n]),
      },
    }, nin=1, nout=1)
    for n in std.range(0, nanodes - 1)
  ],

  local magdecon = [
    g.pnode({
      type: 'MagnifySink',
      name: 'magdecon%d' % n,
      data: {
        output_filename: outputfile,
        root_file_mode: 'UPDATE',
        // 'rawdecon%d' is a special-mode tag emitted by OmnibusSigProc only
        // when its rawdecon_tag config is non-empty.  Absent in production runs;
        // MagnifySink silently skips missing tags, so leaving it here is safe.
        frames: ['gauss%d' % n, 'wiener%d' % n, 'rawdecon%d' % n],
        trace_has_tag: true,
        anode: wc.tn(tools.anodes[n]),
      },
    }, nin=1, nout=1)
    for n in std.range(0, nanodes - 1)
  ],

  local magdebug = [
    g.pnode({
      type: 'MagnifySink',
      name: 'magdebug%d' % n,
      data: {
        output_filename: outputfile,
        root_file_mode: 'UPDATE',
        frames: ['tight_lf%d' %n, 'loose_lf%d' %n, 'cleanup_roi%d' %n,
                 'break_roi_1st%d' %n, 'break_roi_2nd%d' %n,
                 'shrink_roi%d' %n, 'extend_roi%d' %n, 'mp2_roi%d' %n, 'mp3_roi%d' %n],
        trace_has_tag: true,
        anode: wc.tn(tools.anodes[n]),
      },
    }, nin=1, nout=1)
    for n in std.range(0, nanodes - 1)
  ],

  local magtruth = [
    g.pnode({
      type: 'MagnifySink',
      name: 'magtruth%d' % n,
      data: {
        output_filename: outputfile,
        root_file_mode: 'UPDATE',
        frames: ['deposplat%d' % n],
        trace_has_tag: true,
        anode: wc.tn(tools.anodes[n]),
      },
    }, nin=1, nout=1)
    for n in std.range(0, nanodes - 1)
  ],


  local magthr = [
    g.pnode({
      type: 'MagnifySink',
      name: 'magthr%d' % n,
      data: {
        output_filename: outputfile,
        root_file_mode: 'UPDATE',
        summaries: ['threshold%d' % n],  // note that if tag set, each apa should have a tag set for FrameFanin
        summary_operator: { ['threshold%d' % n]: 'set' },  // []: obj comprehension
        anode: wc.tn(tools.anodes[n]),
      },
    }, nin=1, nout=1)
    for n in std.range(0, nanodes - 1)
  ],

  // Trun-aware decon sink for the standalone SP->Magnify conversion
  // (pdhd/wct-sp-to-magnify.jsonnet).  First sink RECREATEs the file;
  // subsequent sinks UPDATE it.  Only the first sink carries runinfo to avoid
  // duplicate Trun cycles.  NB: sinks are named 'magdecon%d' by anode IDENT
  // while decon_pipe's are 'magdecon%d' by INDEX n -- do not wire decon_pipe
  // and decon_trun_pipe into the same job or the names may collide.
  local mksink_trun(n, anode) = g.pnode({
    type: 'MagnifySink',
    name: 'magdecon%d' % anode.data.ident,
    data: {
      output_filename: outputfile,
      root_file_mode: if n == 0 then 'RECREATE' else 'UPDATE',
      // 'rawdecon%d' is a special-mode tag (off in production);
      // MagnifySink silently skips tags absent from the input frame.
      frames: ['gauss%d' % anode.data.ident, 'wiener%d' % anode.data.ident,
               'rawdecon%d' % anode.data.ident],
      // Retagger (inserted upstream) copies wiener<N> → threshold<N>, so
      // summaries get written as h[uvw]_threshold<N> as Magnify expects.
      summaries: ['threshold%d' % anode.data.ident],
      summary_operator: { ['threshold%d' % anode.data.ident]: 'set' },
      cmmtree: [['bad', 'T_bad%d' % anode.data.ident]],
      trace_has_tag: true,
      anode: wc.tn(anode),
    } + (if n == 0 && runinfo != null
         then { runinfo: runinfo { anodeNo: anode.data.ident },
                geo_tree: 'T_geo%d' % anode.data.ident }
         else {}),
  }, nin=1, nout=1, uses=[anode]),

  local magdnnsp = [
    g.pnode({
      type: 'MagnifySink',
      name: 'magdnnsp%d' % n,
      data: {
        output_filename: outputfile,
        root_file_mode: 'UPDATE',
        frames: ['dnnsp%d' %n],
        trace_has_tag: true,
        anode: wc.tn(tools.anodes[n]),
      },
    }, nin=1, nout=1)
    for n in std.range(0, nanodes - 1)
  ],


  return: {
    truth_pipe: [g.pipeline([magtruth[n]], name='magtruthpipe%d' % n) for n in std.range(0, nanodes - 1)],
    orig_pipe: [g.pipeline([magorig[n]], name='magorigpipe%d' % n) for n in std.range(0, nanodes - 1)],
    raw_pipe: [g.pipeline([magraw[n]], name='magrawpipe%d' % n) for n in std.range(0, nanodes - 1)],
    decon_pipe: [g.pipeline([magdecon[n]], name='magdeconpipe%d' % n) for n in std.range(0, nanodes - 1)],
    decon_trun_pipe: [g.pipeline([mksink_trun(n, tools.anodes[n])], name='magpipe%d' % n) for n in std.range(0, nanodes - 1)],
    debug_pipe: [g.pipeline([magdebug[n]], name='magdebugpipe%d' % n) for n in std.range(0, nanodes - 1)],
    threshold_pipe: [g.pipeline([magthr[n]], name='magthrpipe%d' % n) for n in std.range(0, nanodes - 1)],
    dnnsp_pipe: [g.pipeline([magdnnsp[n]], name='magdnnsppipe%d' % n) for n in std.range(0, nanodes - 1)],
  },


}.return
