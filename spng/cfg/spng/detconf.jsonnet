
local pdvd = import "detconfs/pdvd.jsonnet";
local pdhd = import "detconfs/pdhd.jsonnet";

{
    pdvd: pdvd { name: "pdvd"},
    pdhd: pdhd { name: "pdhd"},
}
