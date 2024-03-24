// Import all variant parameter packs to allow for dict like lookup.

local ssss = import "variants/ssss.jsonnet";
local resample = import "variants/resample.jsonnet" ;
local nominal = import "variants/nominal.jsonnet";

{nominal: nominal} + ssss + resample
