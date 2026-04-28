// local g = import 'pgraph.jsonnet';
// local wc = import 'wirecell.jsonnet';

function(g, wc, tools) {
    resamplers : [
        g.pnode({
            type: 'Resampler',
            name: 'resmp%d' %n,
            data: {
                period: 500*wc.ns,
                time_pad: "linear",
            }
        }, nin=1, nout=1,
        uses=[tools.dft, tools.anodes[n]])
        for n in std.range(0, std.length(tools.anodes) - 1)
    ]
}