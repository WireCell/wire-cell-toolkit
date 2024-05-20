local foo = {
    fr: "bar",
    nf: {
        fr: $.fr,
    },
    sp: {
        fr: $.fr,
    }
};

local bar = std.mergePatch(foo, {fr:"baz"});
bar
