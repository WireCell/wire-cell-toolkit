    local wc = import "wirecell.jsonnet";
    {
    
    wires_file: "protodunehd-wires-larsoft-v1.json.bz2",
    
    // distance between collection wire plane and a plane.
    xplanes: {
        danode: 10.0*wc.mm,
        dresponse: 100.0*wc.mm,
        dcathode: 1000.0*wc.mm
    },
    local xplanes = self.xplanes, // to make available below

    volumes: [
    
    {
        wires: 0,
        xcenter: -3575.0*wc.mm, // absolute center of APA
        local xcenter = self.xcenter, // to make available below.
        faces: [

        {
            local dcollection = 52.8*wc.mm,
            anode: xcenter + (xplanes.danode + dcollection),
            response: xcenter + (xplanes.dresponse + dcollection),
            cathode: xcenter + (xplanes.dcathode + dcollection),
        },

        {
            local dcollection = 52.8*wc.mm,
            anode: xcenter - (xplanes.danode + dcollection),
            response: xcenter - (xplanes.dresponse + dcollection),
            cathode: xcenter - (xplanes.dcathode + dcollection),
        },

        ],
    },


    {
        wires: 1,
        xcenter: 3573.0*wc.mm, // absolute center of APA
        local xcenter = self.xcenter, // to make available below.
        faces: [

        {
            local dcollection = 52.8*wc.mm,
            anode: xcenter + (xplanes.danode + dcollection),
            response: xcenter + (xplanes.dresponse + dcollection),
            cathode: xcenter + (xplanes.dcathode + dcollection),
        },

        {
            local dcollection = 52.8*wc.mm,
            anode: xcenter - (xplanes.danode + dcollection),
            response: xcenter - (xplanes.dresponse + dcollection),
            cathode: xcenter - (xplanes.dcathode + dcollection),
        },

        ],
    },


    {
        wires: 2,
        xcenter: -3575.0*wc.mm, // absolute center of APA
        local xcenter = self.xcenter, // to make available below.
        faces: [

        {
            local dcollection = 52.8*wc.mm,
            anode: xcenter + (xplanes.danode + dcollection),
            response: xcenter + (xplanes.dresponse + dcollection),
            cathode: xcenter + (xplanes.dcathode + dcollection),
        },

        {
            local dcollection = 52.8*wc.mm,
            anode: xcenter - (xplanes.danode + dcollection),
            response: xcenter - (xplanes.dresponse + dcollection),
            cathode: xcenter - (xplanes.dcathode + dcollection),
        },

        ],
    },


    {
        wires: 3,
        xcenter: 3573.0*wc.mm, // absolute center of APA
        local xcenter = self.xcenter, // to make available below.
        faces: [

        {
            local dcollection = 52.8*wc.mm,
            anode: xcenter + (xplanes.danode + dcollection),
            response: xcenter + (xplanes.dresponse + dcollection),
            cathode: xcenter + (xplanes.dcathode + dcollection),
        },

        {
            local dcollection = 52.8*wc.mm,
            anode: xcenter - (xplanes.danode + dcollection),
            response: xcenter - (xplanes.dresponse + dcollection),
            cathode: xcenter - (xplanes.dcathode + dcollection),
        },

        ],
    },

    ]
    }
    
