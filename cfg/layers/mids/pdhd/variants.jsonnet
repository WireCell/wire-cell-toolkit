// Import all pdhd variant parameter packs to allow dict-like lookup.

local nominal = import "variants/nominal.jsonnet";
local ssss = import "variants/ssss.jsonnet";

{
    nominal: nominal,
    // No distinct "actual" yet; alias to nominal so lookups succeed.
    actual: nominal,
} + ssss
