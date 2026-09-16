// Official-side fixture TEMPLATE for check-dunereco-config -- protodunevd (PDVD).
//
// See fixtures/official-pdhd.jsonnet for the full explanation of the template
// mechanism (@REFERENCE@ / @ANODE@ substitution, absolute reference import,
// WIRECELL_PATH ordering, required elecGain ext-var).
//
// protodunevd has 8 anodes: bottom = ident 0..3, top = ident 4..7, each with
// its own filter suffix (_b / _t).  @ANODE@ selects which; index 0 is a bottom
// anode (match the SPNG pdvd tpcid used in the registry).
//
// The reference sp.jsonnet exposes make_sigproc(anode, l1sp_pd_mode='dump',...);
// we pass l1sp_pd_mode='' to emit the bare OmnibusSigProc that SPNG reimplements.

local g = import 'pgraph.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local params = import 'pgrapher/experiment/protodunevd/params.jsonnet';
local sp_maker = import '@REFERENCE@/protodunevd/sp.jsonnet';  // dunereco reference

local tools = tools_maker(params);
local sp = sp_maker(params, tools);
local node = sp.make_sigproc(tools.anodes[@ANODE@], l1sp_pd_mode='');
g.main(node, 'Pgrapher')
