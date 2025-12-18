local wc = import "wirecell.jsonnet";

local ordered_uniq(seq) =
    local index = std.mapWithIndex(function (ind, obj) {ind:ind, obj:obj}, seq);
    local ulist = std.set(index, function (ent) std.toString(ent.obj));
    local ordered = std.sort(ulist, function(ent) ent.ind);
    [ent.obj for ent in ordered];


local fodder = [{c:"c"}, {a: "a"}, {a: "a"}, {a: "b"}, {a:"a"}, {b: "b"}, {c:"c"}];
{
    ul: wc.unique_list(fodder),
    ou: ordered_uniq(fodder),
}
