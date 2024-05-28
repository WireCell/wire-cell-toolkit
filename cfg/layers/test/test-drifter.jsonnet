
local high = import "layers/high.jsonnet";
local low = import "layers/low.jsonnet";
local dr = import "layers/low/drifter.jsonnet";

{
    // params :: high.params("pdsp", "nominal"),

    // api :: high.api("pdsp", self.params),

    ran : high.services.random,
    // xreg : low.util.driftsToXregions(self.params.geometry.drifts),
    // lar: self.params.lar,

    // drifter: self.api.drifter()
    // drifter : dr(high.services.random,
    //                 low.util.driftsToXregions(self.params.geometry.drifts),
    //                 self.params.lar),
}
