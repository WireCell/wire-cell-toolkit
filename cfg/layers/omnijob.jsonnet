// This is a top level config file given to wire-cell.
//
// It can configure a wide variety of types of WCT job graphs for a variety
// of supported detectors.  A detector is supported if it has configuration
// under the so called "layers" schema (wire-cell-toolkit/cfg/layers/).
// 
// CAVEAT and FIXME: the configuration currently configures a single-APA detector.
//
// This configuration produces a "top level function" with its arguments
// ("TLA"s) exposed to the "wire-cell" command line.  The user must provide one
// input file.  Each stage (or "task") in the pipeline may be a source of output
// files that the user may supply.
//
// Example usage:
// 
// wire-cell -c omnijob.jsonnet \
//           -A detector=pdsp \
//           -A variant=nominal \
//           -A input=depos.npz \
//           -A output=sp.npz \
//           -A tasks=sim,nf,sp
//
// - detector :: canonical name of supported detector (pdsp, uboone, etc).
// - variant :: the layers "mids" detector variant name .
// - input :: name of file provding data to input to first task.
// - ouput :: describe the output file(s) (see below).
// - tasks :: array or comma separated list of task names (see below).
//
// The "output" can be a string giving a single file name which will receive
// output from the end of the pipeline or it may be a map from "task" (see
// below) to an output file name.
//
// For example:
//
// wire-cell [...] \
//           --tla-code output="{sim:\"digits.npz\",sp:\"signals.npz\"}" \
//           -A tasks=sim,nf,sp \
//           [...]
//
// will produce "digits.npz" output from the simulation "sim" task and
// "signals.npz" from the signal processing "sp" task.  
//
// The available "task" pipeline items are:
//
// - drift :: drift depos through the bulk TPC volume
// - splat :: apply the "splat" sim+sigproc model to depos to produce frame
// - sim :: apply full simulation to depos to produce frame, include noise and digitization
// - sig :: apply signal simulation alone to convert depo to frame
// - noi :: apply noise simulation alone to frame
// - dig :: apply digitization alone to frame
// - nf :: apply noise filtering to frame
// - sp :: apply signal processing to frame
//
// Tasks are aplied in the order given by the "tasks" list.  

local high = import "layers/high.jsonnet";
local wc = high.wc;
local pg = high.pg;

// Main pipeline elements.
local builder(mid, anode, stages, outputs, dense=true) = {
    local last_stage = stages[std.length(stages)-1],
    main : {
        drift: [
            // explicitly filter depos in anode.
            pg.pnode({          
                type:"DepoSetFilter",
                name: "predrift",
                data:{
                    anode: wc.tn(anode),
                },
            }, nin=1, nout=1, uses=[anode]),
            mid.drifter(),
        ],
        splat: [
            // nothing inline, the pipline with the splat itself is part of the
            // sink branch of the tap
        ],
        sim: [
            mid.signal(anode),
            mid.noise(anode),
            mid.digitizer(anode),
        ],
        sig: [
            mid.signal(anode),
        ],
        noi: [
            mid.noise(anode),
        ],
        dig: [
            mid.digitizer(anode),
        ],
        nf: [
            mid.nf(anode),
        ],
        sp: [
            mid.sp(anode),
        ]
    },
    pre_sink(stage) :: 
        if stage == "splat"
        then [ mid.splat(anode) ]
        else [],

    reframer(stage) ::
        local reframers = {
            splat: [mid.reframer(anode, name=outputs[stage])],
            sp: [mid.reframer(anode, name=outputs[stage], tags=["gauss"])],
        };
        if dense
        then std.get(reframers, stage, [])
        else [],

    file_sink(stage) :: [
        if stage == "drift"
        then high.fio.depo_file_sink(outputs.drift)
        else high.fio.frame_file_sink(outputs[stage])
    ],

    sink(stage) ::
        pg.pipeline(self.pre_sink(stage) + self.reframer(stage) + self.file_sink(stage)),

    tap_or_sink(stage):
        if stage == last_stage then
            self.sink(stage)
        else
            high.fio.tap(if stage == "drift" || stage == "splat"
                         then "DepoSetFanout"
                         else "FrameFanout",
                         self.sink(stage), name=outputs[stage]),
        
    get_stage(stage):
        if std.objectHas(outputs, stage) then
            self.main[stage] + [self.tap_or_sink(stage)]
        else
            self.main[stage],
};
    

local source(stage, input) = 
    if std.member(["drift","splat","sim"], stage)
    then high.fio.depo_file_source(input)
    else high.fio.frame_file_source(input);


local output_objectify(stages, output) =
    local last_stage = stages[std.length(stages)-1];
    if std.type(output) == "string" then
        {[last_stage]: output}
    else output;

// Return configuration for single-APA job.
// - detector :: canonical name of supported detector (pdsp, uboone, etc)
// - input :: name of file provding data to input to first task
// - ouput :: name of file to receive output of last task or object mapping task to output file 
// - tasks :: array or comma separated list of task names
// - dense :: if false, save frames sparsely, else add reframers to make dense 
// - variant :: the layers mids detector variant name 
function (detector, input, output, tasks="drift,splat,sim,nf,sp", dense=true, variant="nominal")
    local params = high.params(detector, variant);
    local mid = high.api(detector, params);
    local stages = wc.listify(tasks);
    local outfiles = output_objectify(stages, output); // stage->filename

    local anodes = mid.anodes();
    local anode = anodes[0];

    local b = builder(mid, anode, stages, outfiles, dense);

    local head = [source(stages[0], input)];
    local guts = [b.get_stage(stage) for stage in stages];
    local body = std.flattenArrays(guts);

    local graph = pg.pipeline(head + body);
    high.main(graph, "Pgrapher")
