local detector = import "spng/detector.jsonnet";
local pdvd = import "spng/detconfs/pdvd.jsonnet";
local pdhd = import "spng/detconfs/pdhd.jsonnet";

{
    pdvd: pdvd { name: "pdvd"},
    pdhd: pdhd { name: "pdhd"},

    /// Return the detector configuration object, reduced by listed tpcids.  If
    /// tpcids is not empty it selects a subset of TPCs to include, o.w. all are
    /// included.
    get(detname, tpcids=[])::
        detector.subset(self[detname], tpcids),

}
