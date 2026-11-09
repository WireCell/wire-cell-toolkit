// Radioactive decays in the liquid argon, as LowEDepoSource "decay" sources.
//
// Activities are per kg of argon.  39Ar is well measured.  The others
// follow the order of magnitude of the DUNE FD radiological model and must
// be checked against the current DUNE background model before results are
// quoted.  The 222Rn daughters are put in equilibrium with 222Rn and uniform
// in the argon, which ignores that the ions drift to the cathode, except
// that 214Po decays only as the daughter of 214Bi, at its place (BiPo).  Decay
// schemes are simplified to their main branches (NNDC).
//
// Not included yet: neutrons (captures give a 6.1 MeV gamma cascade),
// external gammas from the cavern and cryostat, surface activity on the
// APAs and CPAs.

local wc = import "wirecell.jsonnet";

local beta(q, z, shape="allowed") = {q: q*wc.MeV, z: z, shape: shape};

{
    Ar39: {
        name: "Ar39", kind: "decay", activity_bq_per_kg: 1.01,
        branches: [{beta: beta(0.565, 19, "unique_first")}],
    },
    Ar42: {
        name: "Ar42", kind: "decay", activity_bq_per_kg: 1.28e-4,
        branches: [{beta: beta(0.599, 19, "unique_first")}],
    },
    K42: {                      // daughter of 42Ar
        name: "K42", kind: "decay", activity_bq_per_kg: 1.28e-4,
        branches: [
            {probability: 0.82, beta: beta(3.525, 20, "unique_first")},
            {probability: 0.18, beta: beta(2.000, 20), gammas: [1.5247*wc.MeV]},
        ],
    },
    Kr85: {
        name: "Kr85", kind: "decay", activity_bq_per_kg: 0.115,
        branches: [
            {probability: 0.9957, beta: beta(0.687, 37, "unique_first")},
            {probability: 0.0043, beta: beta(0.173, 37), gammas: [0.514*wc.MeV]},
        ],
    },
    Rn222: {
        name: "Rn222", kind: "decay", activity_bq_per_kg: 1.0e-3,
        branches: [{alpha: 5.4895*wc.MeV}],
    },
    Po218: {
        name: "Po218", kind: "decay", activity_bq_per_kg: 1.0e-3,
        branches: [{alpha: 6.0024*wc.MeV}],
    },
    Pb214: {
        name: "Pb214", kind: "decay", activity_bq_per_kg: 1.0e-3,
        branches: [
            {probability: 0.46, beta: beta(0.672, 83), gammas: [0.3519*wc.MeV]},
            {probability: 0.41, beta: beta(0.729, 83), gammas: [0.2952*wc.MeV]},
            {probability: 0.13, beta: beta(1.018, 83)},
        ],
    },
    Bi214: {
        name: "Bi214", kind: "decay", activity_bq_per_kg: 1.0e-3,
        // 214Po follows at the same place, half-life 164.3 us
        daughter: {name: "Po214", mean_life: 164.3*wc.us/std.log(2)},
        branches: [
            {probability: 0.19, beta: beta(3.270, 84)},
            {probability: 0.45, beta: beta(2.661, 84), gammas: [0.6093*wc.MeV]},
            {probability: 0.17, beta: beta(1.505, 84), gammas: [1.7645*wc.MeV]},
            {probability: 0.19, beta: beta(1.540, 84), gammas: [1.1203*wc.MeV, 0.6093*wc.MeV]},
        ],
    },
    Po214: {                    // only as the daughter of 214Bi
        name: "Po214", kind: "decay", rate_hz: 0,
        branches: [{alpha: 7.6869*wc.MeV}],
    },

    // A stress test, not a background estimate: thermal neutron captures on
    // 40Ar, as gamma cascades of 6.099 MeV (illustrative lines).  The
    // activity is arbitrary; use the per-capture pass probability.
    nCapture: {
        name: "nCapture", kind: "decay", activity_bq_per_kg: 5e-4,
        branches: [
            {probability: 0.50, gammas: [4.7453*wc.MeV, 1.1868*wc.MeV, 0.1673*wc.MeV]},
            {probability: 0.25, gammas: [3.7000*wc.MeV, 2.3990*wc.MeV]},
            {probability: 0.25, gammas: [2.7710*wc.MeV, 1.8810*wc.MeV, 1.1868*wc.MeV, 0.2602*wc.MeV]},
        ],
    },

    // The default list.
    all: [self.Ar39, self.Ar42, self.K42, self.Kr85,
          self.Rn222, self.Po218, self.Pb214, self.Bi214, self.Po214],
    // Sources for stress tests only.
    stress: [self.nCapture],
}
