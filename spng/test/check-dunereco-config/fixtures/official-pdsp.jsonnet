// Official-side fixture TEMPLATE for check-dunereco-config -- PDSP.
//
// See fixtures/official-pdhd.jsonnet for the full explanation of the template
// mechanism (@REFERENCE@ / @ANODE@ substitution, absolute reference import,
// WIRECELL_PATH ordering, required elecGain ext-var).
//
// PDSP is the *nominal* ProtoDUNE single-phase design: same plane layout and
// response philosophy as PDHD but with 6 APAs and, crucially, NONE of the PDHD
// "bad APA1" tunings (plane swaps, per-APA Wiener filters, etc).  It is
// therefore the better faithfulness reference for an idealised SPNG PDHD that
// assumes all APAs are good -- use the "pdhd-x-pdsp" comparison for that.
//
// The PDSP reference sp.jsonnet exposes make_sigproc(anode, name=null) (no L1SP
// wrapper), so we emit its OmnibusSigProc directly.

local g = import 'pgraph.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local params = import 'pgrapher/experiment/pdsp/params.jsonnet';
local sp_maker = import '@REFERENCE@/pdsp/sp.jsonnet';  // dunereco reference tree

local tools = tools_maker(params);
local sp = sp_maker(params, tools);
local node = sp.make_sigproc(tools.anodes[@ANODE@]);
g.main(node, 'Pgrapher')
