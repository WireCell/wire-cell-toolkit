// Import all variant parameter packs to allow for dict like lookup.

local low = import "../../low.jsonnet";
//local ssss = import "variants/ssss.jsonnet";

local nominal = import "variants/nominal.jsonnet";
local actual = import "variants/actual.jsonnet";

{
    nominal: nominal,
    actual: actual,
} + low.ssss(nominal)
