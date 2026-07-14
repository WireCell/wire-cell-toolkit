# CMMModifier

Modifies the Channel Mask Map (CMM) of an input frame to add or extend bad-channel entries based on configurable continuity-shorted, veto, dead-channel, and boundary-alignment rules, then passes all traces and tags through unchanged.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `CMMModifier` |
| Concrete class | `WireCell::Img::CMMModifier` |
| Node category | function |
| Primary interface | `IFrameFilter` |
| Input type(s) | `IFrame` |
| Output type(s) | `IFrame` |
| Configurable | yes |

## Configuration Parameters

| Parameter | Description |
| --- | --- |
| `anode` | type/name of the IAnodePlane component (default: "AnodePlane") |
| `cm_tag` | name of the ChannelMask entry in the CMM to process (default: "bad") |
| `trace_tag` | tag selecting traces used when evaluating dead-channel criteria (default: "gauss") |
| `m_start` | tick range start used when marking veto/dead channels as bad |
| `m_end` | tick range end used when marking veto/dead channels as bad |
| `ncount_cont_ch` | number of continuity (shorted-wire) channel ranges to process |
| `cont_ch_llimit` | array of lower channel bounds for continuity veto ranges |
| `cont_ch_hlimit` | array of upper channel bounds for continuity veto ranges |
| `ncount_veto_ch` | number of unconditional veto channel ranges |
| `veto_ch_llimit` | array of lower channel bounds for direct veto ranges |
| `veto_ch_hlimit` | array of upper channel bounds for direct veto ranges |
| `dead_ch_ncount` | minimum number of above-threshold ticks for a channel group to be considered alive |
| `dead_ch_charge` | charge threshold used when counting active ticks for dead-channel detection |
| `ncount_dead_ch` | number of dead-channel groups to check |
| `dead_ch_llimit` | array of lower channel bounds for dead-channel groups |
| `dead_ch_hlimit` | array of upper channel bounds for dead-channel groups |
| `ncount_org` | number of boundary-alignment intervals for dead-channel time ranges |
| `org_llimit` | array of lower tick boundaries for dead-range alignment |
| `org_hlimit` | array of upper tick boundaries for dead-range alignment |
