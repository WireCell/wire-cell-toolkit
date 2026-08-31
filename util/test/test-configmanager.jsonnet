local rnd1 = {
    type: "Random",
    name: "dup",                // identical name
    data: {
        seeds: [1],
    }
};
local rnd2 = {
    type: "Random",
    name: "dup",                // identical name
    data: {
        seeds: [2],
    }
};
local cd = {
    type: "ConfigDumper",
    data: {
        components: ["Random:dup"]
    }
};
local cli = {
    type: "wire-cell",
    data: {
        plugins: ["WireCellGen", "WireCellApps"],
        apps: ["ConfigDumper"],
    }
};

[cli, rnd1, rnd2, cd]

    
