local sys = import "system.jsonnet";
local jobs = import "jobs.jsonnet";
local pg = import "pgraph.jsonnet";

// tlas
function(input_filename, output_filename_pattern, detname='pdhd', engine='Pgrapher')
    local sysobj = sys.nominal(detname);
    local graph = jobs.depos_to_adc(sysobj, input_filename, output_filename_pattern);
    pg.main(graph, engine)



    
