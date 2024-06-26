function(params) {
    sp: {
        type: "FieldResponse",
        name: "sp",
        data: { filename: params.sp.field_file } // singular
    },
    nf: {
        type: "FieldResponse",
        name: "nf",
        data: { filename: params.nf.field_file } // singular
    },
    sim : {
        type: "FieldResponse",
        name: params.ductor.field_file,
        data: { filename: params.ductor.field_file }
    },
    sims : [{
        type: "FieldResponse",
        name: fn,
        data: { filename: fn }
    } for fn in params.ductor.field_files], // plural
}

