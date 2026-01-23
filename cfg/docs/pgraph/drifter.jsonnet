// drift.jsonnet
local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";
{  
  drifter(electron_lifetime=8*units::ms, drift_speed=1.6*wc.mm/wc.us, name=""):: pg.pnode({
    type: "Drifter",
    name: name,
    data: {
      lifetime: electron_lifetime,
      drift_speed: drift_speed,
    }
  }, nin=1, nout=1),
}

