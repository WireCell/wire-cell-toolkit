    local wc = import "wirecell.jsonnet";
    {
    
    wires_file: "protodunevd-wires-larsoft-v3.json.bz2",
    
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
        xcenter: -3415.5*wc.mm, // absolute center of APA
        local xcenter = self.xcenter, // to make available below.
        faces: [

        {
            local dcollection = 0.4*wc.mm,
            anode: xcenter + (xplanes.danode + dcollection),
            response: xcenter + (xplanes.dresponse + dcollection),
            cathode: xcenter + (xplanes.dcathode + dcollection),
        },
null,
        ],
    },


    {
        wires: 1,
        xcenter: -3415.5*wc.mm, // absolute center of APA
        local xcenter = self.xcenter, // to make available below.
        faces: [

        {
            local dcollection = 0.4*wc.mm,
            anode: xcenter + (xplanes.danode + dcollection),
            response: xcenter + (xplanes.dresponse + dcollection),
            cathode: xcenter + (xplanes.dcathode + dcollection),
        },
null,
        ],
    },


    {
        wires: 2,
        xcenter: -3415.5*wc.mm, // absolute center of APA
        local xcenter = self.xcenter, // to make available below.
        faces: [

        {
            local dcollection = 0.4*wc.mm,
            anode: xcenter + (xplanes.danode + dcollection),
            response: xcenter + (xplanes.dresponse + dcollection),
            cathode: xcenter + (xplanes.dcathode + dcollection),
        },
null,
        ],
    },


    {
        wires: 3,
        xcenter: -3415.5*wc.mm, // absolute center of APA
        local xcenter = self.xcenter, // to make available below.
        faces: [

        {
            local dcollection = 0.4*wc.mm,
            anode: xcenter + (xplanes.danode + dcollection),
            response: xcenter + (xplanes.dresponse + dcollection),
            cathode: xcenter + (xplanes.dcathode + dcollection),
        },
null,
        ],
    },


    {
        wires: 4,
        xcenter: 3415.5*wc.mm, // absolute center of APA
        local xcenter = self.xcenter, // to make available below.
        faces: [
null,
        {
            local dcollection = 0.4*wc.mm,
            anode: xcenter - (xplanes.danode + dcollection),
            response: xcenter - (xplanes.dresponse + dcollection),
            cathode: xcenter - (xplanes.dcathode + dcollection),
        },

        ],
    },


    {
        wires: 5,
        xcenter: 3415.5*wc.mm, // absolute center of APA
        local xcenter = self.xcenter, // to make available below.
        faces: [
null,
        {
            local dcollection = 0.4*wc.mm,
            anode: xcenter - (xplanes.danode + dcollection),
            response: xcenter - (xplanes.dresponse + dcollection),
            cathode: xcenter - (xplanes.dcathode + dcollection),
        },

        ],
    },


    {
        wires: 6,
        xcenter: 3415.5*wc.mm, // absolute center of APA
        local xcenter = self.xcenter, // to make available below.
        faces: [
null,
        {
            local dcollection = 0.4*wc.mm,
            anode: xcenter - (xplanes.danode + dcollection),
            response: xcenter - (xplanes.dresponse + dcollection),
            cathode: xcenter - (xplanes.dcathode + dcollection),
        },

        ],
    },


    {
        wires: 7,
        xcenter: 3415.5*wc.mm, // absolute center of APA
        local xcenter = self.xcenter, // to make available below.
        faces: [
null,
        {
            local dcollection = 0.4*wc.mm,
            anode: xcenter - (xplanes.danode + dcollection),
            response: xcenter - (xplanes.dresponse + dcollection),
            cathode: xcenter - (xplanes.dcathode + dcollection),
        },

        ],
    },

    ]
    }
    
