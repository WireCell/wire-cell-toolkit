// Exhaustive array of channel IDs grouped into coherent groups.
// This must correspond to the pattern in the wires-file!
function(n)
    [std.range(n * 2560 + u * 40, n * 2560 + (u + 1) * 40 - 1) for u in std.range(0, 19)] +
    [std.range(n * 2560 + 800 + v * 40, n * 2560 + 800 + (v + 1) * 40 - 1) for v in std.range(0, 19)] +
    [std.range(n * 2560 + 1600 + w * 48, n * 2560 + 1600 + (w + 1) * 48 - 1) for w in std.range(0, 19)]


