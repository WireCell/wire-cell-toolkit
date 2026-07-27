// Enforce schema via function calls.

// generic performance test schema
local perf(tech, kind, repeat, device, shape) = {
    tech: tech,
    kind: kind,
    repeat: repeat,
    device: device,
    shape: shape,
};

local repeats = {
    convo: { gpu:5000, cpu:50 },
    median: { gpu:5000, cpu:10 },
    sort: { gpu:10000, cpu:10 },
    arith: { gpu:10000, cpu:100 },
};

local shapes = [ [1024,8192], ];

function(techs="af,lt,ei", devices="gpu,cpu", tests="convo,median,sort,arith")
[
    perf(tech, tname, std.get(repeats[tname], dev, 10), dev, shape)
    for dev in std.split(devices,",")
    for tech in std.split(techs,",")
    for shape in shapes
    for tname in std.split(tests,",")
]
