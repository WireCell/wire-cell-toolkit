//  This defines basic parameters for known detectors.
// See also detectors.jsonnet from wire-cell-data for known files.
{
    hd: {
        connections: [2,2,0]
        faces: [0,1],
    },
    pdhd: hd,
    fdhd: hd,

    vd: {
        connections: [1,1,0],
        faces: [0,1],
    },
    pdvd: vd,
    fdvd: vd,

    uboone: {
        connections: [0,0,0],
        faces: [0],
    }
}

