local high = import "layers/high.jsonnet";
local wc = high.wc;
local pg = high.pg;

function(detector, variant="nominal",
         indepos=null,
         outdepos="depos-drifted.npz",
         frames="frames-adc-%(anode)s.npz")

    local mid = high.mid(detector, variant, options={sparse:false});

    local source = if std.type(indepos) == "null"
                  then mid.sim.track_depos()
                  else high.fio.depo_file_source(indepos);

    local drifter = mid.drifter();

    local anodes = mid.anodes();
    local nanodes = std.length(anodes);
    local anode_iota = std.range(0, nanodes-1);

    local apipes = [
        local anode = anodes[aid];
        pg.pipeline([
            mid.sim.signal(anode),
            // mid.sim.noise(anode),
            mid.sim.digitizer(anode),

            high.fio.frame_file_sink(
                std.format(frames, { anode: anode.data.ident}),
                digitize=true,
                tags=["orig%d" % anode.data.ident]),

        ]) for aid in anode_iota];

    local body = pg.fan.fanout('DepoSetFanout', apipes, "work");

    local graph = pg.pipeline([source, drifter, body], "main");
    local executor = "TbbFlow";
    // local executor = "Pgrapher";
    high.main(graph, executor)

