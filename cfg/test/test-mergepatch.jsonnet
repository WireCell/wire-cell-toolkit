local wc = import "wirecell.jsonnet";

local g1 = {
    wires_file: 'microboone-celltree-wires-v2.1.json.bz2',

    // distance between collection wire plane and a plane.
    xplanes: {
        danode: 10.0 * wc.mm,
        dresponse: 100.0 * wc.mm,
        dcathode: 1000.0 * wc.mm,
    },
    local xplanes = self.xplanes,  // to make available below

    volumes: [

        {
            wires: 0,
            xcenter: -6.0 * wc.mm,  // absolute center of APA
            local xcenter = self.xcenter,  // to make available below.
            faces: [

                {
                    local dcollection = 6.0 * wc.mm,
                    anode: xcenter + (xplanes.danode + dcollection),
                    response: xcenter + (xplanes.dresponse + dcollection),
                    cathode: xcenter + (xplanes.dcathode + dcollection),
                },
                null,
            ],
        },

    ],
};

local g2 =  {
    xplanes: {
        // Distances between collection plane and a logical plane.
        danode: 0.0,
        dresponse: 10*wc.cm,
        // Location of cathode measured from collection is based on a dump
        // of ubcore/v06_83_00/gdml/microboonev11.gdml by Matt Toups
        dcathode: 2.5480*wc.m,
    }
};

{
    "g1+g2": g1+g2,
    "g2+g1": g2+g1,
    "merge(g1,g2)": std.mergePatch(g1,g2),
    "merge(g2,g1)": std.mergePatch(g2,g1),
}
