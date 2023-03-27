// This generates noise, adds it to ADC voltage level input from file
// and runs NF+SP.
// 
// See "wirecell-util npz-to-wct" as one way to produce ADC voltage
// level from ADC array.


local high = import "layers/high.jsonnet";
local wc = high.wc;
local pg = high.pg;

function(adcvoltfile, spfrfile,
         outprefix="", aid="0",
         detector="pdsp", variant="nominal",
         adcname="*", signame="gauss0")

    local options = {params: {sp:{field_file:spfrfile},
                              nf:{field_file:spfrfile} }};
    local mid = high.mid(detector, variant, options=options);
    local anodes = mid.anodes();
    local anode = anodes[std.parseInt(aid)];

    local source = high.fio.frame_file_source(adcvoltfile, tags=[adcname]);
    local sink = high.fio.frame_file_sink(outprefix + "sig.npz", tags=[signame],
                                          dense={chbeg:0,chend:800,tbbeg:0,tbend:6000});



    local graph = pg.pipeline([
        source,
        high.fio.tap('FrameFanout',
                     high.fio.frame_file_sink(outprefix + "src.npz", tags=[])),

        mid.sim.reframer(anode),
        high.fio.tap('FrameFanout',
                     high.fio.frame_file_sink(outprefix + "ref.npz", tags=[])),

        mid.sim.noise(anode),
        high.fio.tap('FrameFanout',
                     high.fio.frame_file_sink(outprefix + "noi.npz", tags=[])),
        
        mid.sim.digitizer(anode),
        high.fio.tap('FrameFanout',
                     high.fio.frame_file_sink(outprefix + "dig.npz", tags=[])),
        
        mid.sigproc.nf(anode),
        mid.sigproc.sp(anode),
        sink
    ], "main");

    local executor = "TbbFlow";
    high.main(graph, executor)

