// Add schema for configuration object attributes that the C++ code fail to
// create in the object returned by default_configuration().
//
// See test-config-export.bats
function( pkg, base )
    std.mergePatch(base, std.get({
        WireCellGen: {
            MultiDuctor: {
                config: {
                    fields: base.MultiDuctor.config.fields + [
                        {           // a complex sequence of objects
                            name: "chains", 
                            item: {
                                schema: "sequence"
                            },
                        },
                    ]
                }
            },
            WireBoundedDepos: {
                config: {
                    fields: base.WireBoundedDepos.config.fields + [
                        {           // a complex sequence of objects
                            name: "regions",
                            item: {
                                schema: "sequence",
                            },
                        },
                    ],
                },
            },
        },
    }, pkg, {}))
            
