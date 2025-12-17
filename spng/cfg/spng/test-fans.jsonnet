local fans_mod = import "spng/fans.jsonnet";
local fans = fans_mod();
{
    forfilter: fans.fanout_select("test", 3, [["w","g","d"], ["w","g","d"], ["w","g"]]),

    local names = ["a","b","c","d"],
    local fangen = fans.fanout_cross_gen(3, 2, function(num, M) {
        type: "FrameFanout",
        name: 'fan_' + names[num],
        data: {
            multiplicity: M
        }
    }),
    fangen_source0: fangen[1][0]
}
