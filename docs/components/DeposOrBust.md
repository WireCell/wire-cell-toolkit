# DeposOrBust

Routes non-empty depo sets to output port 0 (IDepoSet) and replaces empty depo sets with a minimal empty IFrame on output port 1, providing a branch pattern to bypass downstream pipelines when there is nothing to process.

## Node Properties

| Property | Value |
| --- | --- |
| Factory name | `DeposOrBust` |
| Concrete class | `WireCell::Gen::DeposOrBust` |
| Node category | hydra |
| Primary interface | `IDepos2DeposOrFrame` |
| Input type(s) | `IDepoSet` |
| Output type(s) | `IDepoSet or IFrame` |
| Configurable | no |
