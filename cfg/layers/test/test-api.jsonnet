local high = import "layers/high.jsonnet";
local detectors = ["pdsp", "uboone"];
local variants = ["nominal", "spdir_hifr"];

{
    [det] : {
        params :: high.params(det, "nominal"),
        api :: high.api(det, self.params),
        methods : std.objectFieldsAll(self.api),
        nmethods: std.length(self.methods),
        nvariants: std.length([high.params(det, var) for var in variants]),
    } for det in detectors    
}
