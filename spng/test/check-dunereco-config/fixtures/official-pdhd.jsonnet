// Official-side fixture TEMPLATE for the check-dunereco-config checker.
//
// This is a TEMPLATE, not directly evaluatable: the check-dunereco-config
// script substitutes the placeholders
//
//     @REFERENCE@   absolute path to the dunereco reference DUNEWireCell dir
//     @ANODE@       anode/APA index to instantiate (integer)
//
// and evaluates the resulting file with wcsonnet.  Substitution (rather than a
// Jsonnet TLA/ext-var) is used for two reasons: Jsonnet `import` requires a
// string literal (so the reference path cannot be an ext-var), and the standard
// toolkit/cfg tree must remain ahead of the reference tree on WIRECELL_PATH to
// avoid shadowing shared files (wirecell.jsonnet, pgraph.jsonnet, pgrapher/...).
//
// It instantiates the "official" dunereco OSP signal-processing node for one
// PDHD anode so it can be evaluated to a flat JSON config list and compared,
// like-for-like, against the SPNG OSP configuration.
//
// params + tools are built from the standard pgrapher tree exactly as the
// official jobs (wct-sim-check.jsonnet, wcls-sp.jsonnet) do.  We emit the bare
// OmnibusSigProc (l1sp_pd_mode='') so the comparison targets the core OSP
// configuration that SPNG reimplements, not the PDHD-specific L1SP wrapper.
//
// Required ext-var (see pgrapher/experiment/pdhd/params.jsonnet):
//   elecGain   FE amplifier gain in mV/fC, e.g. -V elecGain=14.0

local g = import 'pgraph.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local params = import 'pgrapher/experiment/pdhd/params.jsonnet';
local sp_maker = import '@REFERENCE@/pdhd/sp.jsonnet';  // dunereco reference tree

local tools = tools_maker(params);
local sp = sp_maker(params, tools);
// l1sp_pd_mode='' yields the bare OmnibusSigProc (no L1SPFilterPD wrapper).
local node = sp.make_sigproc(tools.anodes[@ANODE@], l1sp_pd_mode='');
g.main(node, 'Pgrapher')
