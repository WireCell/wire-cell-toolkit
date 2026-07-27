// This provides helpers to give TLA support to a Jsonnet file.

local detconf = import "detconf.jsonnet";
local pg = import "pgraph.jsonnet";

// Define default "control" object.  This defines things not relevant to
// detectors nor to high level construction.  Rather, relevant to how the user
// wants to run a job.
//
// Fixme: as-is, this is too rigid.
local control = {

    dft: {
        type: "FftwDFT",
        name: "",
        data: {}
    },

    rng: {
        type: 'Random',
        name: "",
        data: {
            generator: "default",
            seeds: [0,1,2,3,4],
        }
    },

};

// A Jsonnet can provide a flexible CLI + library interface for a detector
// of one TPC when it does like this:
//
// @code{jsonnet}
// local tlas = import "tlas.jsonnet";
// local one_tpc(tpc, control) { /* ... define "stage" attributes */ };
// tlas.single_tpc(one_tpc)
// @endcode
//
// User may then run the Jsonnet like:
//
// @code{sh}
// wcsonnet -A detname=pdhd -A tpcid=tpc0 -A stage="all" -A finalize=true the.jsonnet
// @endcode
//
// Or it may import the.jsonnet and call with finalize=false to get a subgraph.
local single_tpc(stages_function) = 
    function(stage="", detname="pdhd", tpcid="tpc0", finalize=false)
        local det = detconf[detname];
        local tpc = det.tpc[tpcid];
        local stages = stages_function(tpc, control);
        local graph = if stage == ""
                      then stages
                      else stages[stage];
        if finalize == true || finalize == "true"
        then pg.main(graph, 'Pgrapher', plugins=["WireCellSpng"])
        else graph;



local whole_detector(stages_function) = 
    function(stage="", detname="pdhd", finalize=false)
        local det = detconf[detname];
        local stages = stages_function(det, control);
        local graph = if stage == ""
                      then stages
                      else stages[stage];
        if finalize == true || finalize == "true"
        then pg.main(graph, 'Pgrapher', plugins=["WireCellSpng"])
        else graph;

local focus_tpc(stages_function) = 
    function(stage="", detname="pdhd", tpcid="tpc0", finalize=false)
        local det = detconf[detname];
        local tpc = det.tpc[tpcid];
        local stages = stages_function(det, tpc, control);
        local graph = if stage == ""
                      then stages
                      else stages[stage];
        if finalize == true || finalize == "true"
        then pg.main(graph, 'Pgrapher', plugins=["WireCellSpng"])
        else graph;

{
    control: control,
    single_tpc: single_tpc,
    whole_detector: whole_detector,
    focus_tpc: focus_tpc,
}
