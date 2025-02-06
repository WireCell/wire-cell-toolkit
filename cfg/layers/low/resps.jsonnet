function(params, options = {}) {
    sim: {
        er: {
            type: "ColdElecResponse",
            name: "sim",
            data: params.elec + params.ductor.binning,
        },
        fr: {
            type: "FieldResponse",
            name: "sim",
            data: { filename: std.get(options, "fields", params.ductor.field_file) }
        },
        rc: {
            type: 'RCResponse',
            data: params.ductor.binning { width: params.rc.width }
        },
    },
    nf: {
        fr: {
            type: "FieldResponse",
            name: "nf",
            data: { filename: std.get(options, "fields", params.nf.field_file) }
        },
    },
    sp: {
        er: {
            type: "ColdElecResponse",
            name: "sp",
            data: params.elec + params.sp.binning,
        },
        fr: {
            type: "FieldResponse",
            name: "sp",
            data: { filename: std.get(options, "fields", params.sp.field_file) }
        },
    },
        
}

