// This is a top level config file given to wire-cell.
//
// It configures a very specific WCT job that applies the Resampler.
//
// It yields a top level function taking arguments, as example:
//
// wire-cell -c resample.jsonnet
//           -A input=sampled.npz \
//           -A output=resampled.npz \
//           -A period_ns=500
//

local high = import "layers/high.jsonnet";
local wc = high.wc;
local pg = high.pg;

function (input, output, period_ns="500", time_pad='linear', dft="cpu")
    local src = high.fio.frame_file_source(input);
    local sink = high.fio.frame_file_sink(output);
    local idft = high.services(dft).dft;
    local rs = pg.pnode({
        type: 'Resampler',
        name: "",
        data: {
            period: if std.type(period_ns) == "number"
                    then period_ns
                    else std.parseJson(period_ns),
            time_pad: time_pad,
            dft: wc.tn(idft)
        }
    }, nin=1, nout=1, uses=[idft]);
    local graph = pg.pipeline([src, rs, sink]);
    high.main(graph, "Pgrapher")
                        

