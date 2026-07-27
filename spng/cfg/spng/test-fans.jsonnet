local fans_js = import "spng/fans.jsonnet";
local fans = fans_js();
{
    local crossed_views = [1,1,0],
    crossed: fans.fanout_select("crossed",
                                std.length(crossed_views), targets_list=[
                                    if is_crossed == 1
                                    then ["wiener", "dense"]
                                    else ["wiener"]
                                    for is_crossed in crossed_views]),


    forfilter: fans.fanout_select("test", 3, [["w","g","d"], ["w","g","d"], ["w","g"]]),

    local names = ["a","b","c","d"],
    local fangen = fans.fanout_cross_gen(3, 2, function(num, M) {
        type: "FrameFanout",
        name: 'fan_' + names[num],
        data: {
            multiplicity: M
        }
    }),
    fangen_source0: fangen[1][0],

    fosl: fans.fanout_shuntline(3),

}.fosl
